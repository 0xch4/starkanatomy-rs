//! Basic Fri proofs generation and verification.

use blake2::{Blake2b512, Digest};

use crate::{
    algebra::FieldElement,
    macros::{codeword_to_merkle, exact, try_match},
    merkle::{Hash, MerkleHash, MerkleTree},
    proof_stream::{ProofObject, ProofStream},
    univariate::Polynomial,
};

pub type Codeword = Vec<FieldElement>;

/// A parameterized Fri protocol.
/// `offset` here is refered to as `g` in the tutorial.
pub struct Fri {
    offset: FieldElement,
    omega: FieldElement,
    domain_length: usize,
    expansion_factor: usize,
    num_colinearity_tests: usize,
}

impl Fri {
    pub fn new(
        offset: FieldElement,
        omega: FieldElement,
        domain_length: usize,
        expansion_factor: usize,
        num_colinearity_tests: usize,
    ) -> Self {
        Self {
            offset,
            omega,
            domain_length,
            expansion_factor,
            num_colinearity_tests,
        }
    }

    pub fn get_domain_length(&self) -> usize {
        self.domain_length
    }

    fn num_rounds(&self) -> usize {
        let mut codeword_length = self.domain_length;
        let mut num_rounds = 0;

        while (codeword_length > self.expansion_factor)
            && (4 * self.num_colinearity_tests < codeword_length)
        {
            codeword_length = codeword_length / 2;
            num_rounds += 1;
        }

        num_rounds
    }

    pub fn evaluation_domain(&self) -> Codeword {
        let mut domain = Vec::new();

        for i in 0..self.domain_length {
            domain.push(self.offset * (self.omega.pow(i as _)));
        }

        domain
    }

    fn sample_indices(
        &self,
        seed: &[u8],
        size: usize,
        reduced_size: usize,
        n: usize,
    ) -> Vec<usize> {
        assert!(n <= reduced_size as usize);
        assert!(n <= 2 * reduced_size as usize);

        let mut indices = Vec::new();
        let mut reduced_indices = Vec::new();
        let mut counter: usize = 0;

        while indices.len() < n {
            let mut bytes = [0u8; 8];
            bytes
                .copy_from_slice(&Blake2b512::digest([seed, &counter.to_le_bytes()].concat())[..8]);
            let index = usize::from_be_bytes(bytes) % size;
            let reduced_index = index % reduced_size;
            counter += 1;
            if !reduced_indices.contains(&reduced_index) {
                indices.push(index);
                reduced_indices.push(reduced_index);
            }
        }

        indices
    }

    fn commit(&self, codeword: &Codeword, proof_stream: &mut ProofStream) -> Vec<Codeword> {
        const ONE: FieldElement = FieldElement::ONE;
        const TWO: FieldElement = FieldElement::new(2);
        let mut codewords = Vec::new();
        let mut codeword = codeword.clone();
        let mut omega = self.omega;
        let mut offset = self.offset;

        // Iterative randomized folding of the codeword onto itself
        for round in 0..self.num_rounds() {
            let leaves = codeword_to_merkle!(codeword);
            let tree = MerkleTree::commit(leaves);
            proof_stream.push(ProofObject::MerkleRoot(tree.get_root()));

            if round + 1 == self.num_rounds() {
                break;
            }

            let alpha =
                FieldElement::sample(proof_stream.prover_fiat_shamir()[..16].try_into().unwrap());

            codewords.push(codeword.clone());

            let mut new_codeword = Vec::new();
            for i in 0..codeword.len() / 2 {
                let left = codeword[i];
                let right = codeword[codeword.len() / 2 + i];
                let new_fe = TWO.inv()
                    * ((ONE + alpha / (offset * (omega.pow(i as _)))) * left
                        + (ONE - alpha / (offset * (omega.pow(i as _)))) * right);
                new_codeword.push(new_fe);
            }

            codeword = new_codeword.clone();
            omega = omega.pow(2);
            offset = offset.pow(2);
        }

        proof_stream.push(ProofObject::LastCodeword(codeword.clone()));
        codewords.push(codeword.clone());

        codewords
    }

    fn query(
        &self,
        cur_codeword: &Codeword,
        next_codeword: &Codeword,
        c_indices: &Vec<usize>,
        proof_stream: &mut ProofStream,
    ) -> Vec<usize> {
        let a_indices = c_indices.clone();
        let b_indices: Vec<usize> = c_indices
            .iter()
            .map(|i| i + cur_codeword.len() / 2)
            .collect();

        for s in 0..self.num_colinearity_tests {
            proof_stream.push(ProofObject::ColinearityTest((
                cur_codeword[a_indices[s]],
                cur_codeword[b_indices[s]],
                next_codeword[c_indices[s]],
            )));
        }

        let cur_tree = MerkleTree::commit(codeword_to_merkle!(cur_codeword));

        let next_tree = MerkleTree::commit(codeword_to_merkle!(next_codeword));

        for s in 0..self.num_colinearity_tests {
            proof_stream.push(ProofObject::MerklePath(cur_tree.open(a_indices[s])));
            proof_stream.push(ProofObject::MerklePath(cur_tree.open(b_indices[s])));
            proof_stream.push(ProofObject::MerklePath(next_tree.open(c_indices[s])));
        }

        [a_indices, b_indices].concat()
    }

    pub fn prove(&self, codeword: &Codeword, proof_stream: &mut ProofStream) -> Vec<usize> {
        assert!(self.domain_length == codeword.len());

        let codewords = self.commit(codeword, proof_stream);

        let top_level_indices = self.sample_indices(
            &proof_stream.prover_fiat_shamir(),
            codewords[1].len(),
            codewords[codewords.len() - 1].len(),
            self.num_colinearity_tests,
        );
        let mut indices = top_level_indices.clone();

        for i in 0..codewords.len() - 1 {
            indices = indices
                .iter()
                .map(|index| index % (codewords[i].len() / 2))
                .collect();

            self.query(&codewords[i], &codewords[i + 1], &indices, proof_stream);
        }

        top_level_indices
    }

    pub fn verify(&self, proof_stream: &mut ProofStream) -> Option<Vec<(usize, FieldElement)>> {
        let mut res = Vec::new();
        let mut omega = self.omega;
        let mut offset = self.offset;

        let mut roots = Vec::new();
        let mut alphas = Vec::new();
        for _ in 0..self.num_rounds() {
            let root = try_match!(proof_stream.pull().unwrap(), ProofObject::MerkleRoot)?;
            roots.push(root);

            alphas.push(FieldElement::sample(
                proof_stream.verifier_fiat_shamir()[..16]
                    .try_into()
                    .unwrap(),
            ));
        }

        let last_codeword = try_match!(proof_stream.pull().unwrap(), ProofObject::LastCodeword)?;

        exact!(
            roots[roots.len() - 1]
                == MerkleTree::commit(codeword_to_merkle!(last_codeword)).get_root()
        )
        .ok()?;

        let mut last_omega = omega;
        let mut last_offset = offset;
        for _ in 0..self.num_rounds() - 1 {
            last_omega = last_omega.pow(2);
            last_offset = last_offset.pow(2);
        }
        assert_eq!(
            last_omega.inv(),
            last_omega.pow(last_codeword.len() as u128 - 1)
        );

        let degree = (last_codeword.len() / self.expansion_factor) - 1;
        let mut last_domain = Vec::new();
        for i in 0..last_codeword.len() {
            last_domain.push((last_offset * last_omega.pow(i as _), last_codeword[i]));
        }
        let poly = Polynomial::interpolate(&last_domain);
        assert_eq!(
            poly.evaluate_on_domain(&last_domain.iter().map(|&(d, _)| d).collect()),
            last_codeword
        );

        exact!(poly.degree() <= degree).ok()?;

        let top_level_indices = self.sample_indices(
            &proof_stream.verifier_fiat_shamir(),
            self.domain_length >> 1,
            self.domain_length >> (self.num_rounds() - 1),
            self.num_colinearity_tests,
        );

        for r in 0..self.num_rounds() - 1 {
            let c_indices: Vec<usize> = top_level_indices
                .iter()
                .map(|index| index % (self.domain_length >> (r + 1)))
                .collect();
            let a_indices = c_indices.clone();
            let b_indices: Vec<usize> = a_indices
                .iter()
                .map(|index| index + (self.domain_length >> (r + 1)))
                .collect();

            // Verify colinearity criteria
            let mut aa = Vec::new();
            let mut bb = Vec::new();
            let mut cc = Vec::new();
            for s in 0..self.num_colinearity_tests {
                let (ay, by, cy) =
                    try_match!(proof_stream.pull().unwrap(), ProofObject::ColinearityTest)?;
                aa.push(ay);
                bb.push(by);
                cc.push(cy);

                if r == 0 {
                    res.push((a_indices[s], ay));
                    res.push((b_indices[s], by));
                }

                let ax = offset * omega.pow(a_indices[s] as _);
                let bx = offset * omega.pow(b_indices[s] as _);
                let cx = alphas[r];

                exact!(Polynomial::test_colinearity(vec![
                    (ax, ay),
                    (bx, by),
                    (cx, cy)
                ]))
                .ok()?;
            }

            for i in 0..self.num_colinearity_tests {
                let a_path = try_match!(proof_stream.pull().unwrap(), ProofObject::MerklePath)?;
                let a_hash: Hash = Blake2b512::digest(aa[i].to_le_bytes()).into();
                exact!(MerkleTree::verify_for_root(
                    roots[r],
                    MerkleHash::from(a_hash),
                    &a_path,
                    a_indices[i],
                ))
                .ok()?;

                let b_path = try_match!(proof_stream.pull().unwrap(), ProofObject::MerklePath)?;
                let b_hash: Hash = Blake2b512::digest(bb[i].to_le_bytes()).into();
                exact!(MerkleTree::verify_for_root(
                    roots[r],
                    MerkleHash::from(b_hash),
                    &b_path,
                    b_indices[i],
                ))
                .ok()?;

                let c_path = try_match!(proof_stream.pull().unwrap(), ProofObject::MerklePath)?;
                let c_hash: Hash = Blake2b512::digest(cc[i].to_le_bytes()).into();
                exact!(MerkleTree::verify_for_root(
                    roots[r + 1],
                    MerkleHash::from(c_hash),
                    &c_path,
                    c_indices[i],
                ))
                .ok()?;
            }

            omega = omega.pow(2);
            offset = offset.pow(2);
        }

        Some(res)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn fri() {
        let degree = 63;
        let expansion_factor = 4;
        let num_colinearity_tests = 17;

        let initial_codeword_length = (degree + 1) * expansion_factor;
        let mut log_codeword_length = 0;
        let mut codeword_length = initial_codeword_length;
        while codeword_length > 1 {
            codeword_length /= 2;
            log_codeword_length += 1;
        }
        assert_eq!(1 << log_codeword_length, initial_codeword_length);

        let omega = FieldElement::primitive_nth_root(initial_codeword_length as _);
        let generator = FieldElement::GENERATOR;
        assert_eq!(omega.pow(1 << log_codeword_length), FieldElement::ONE);
        assert_ne!(omega.pow(1 << (log_codeword_length - 1)), FieldElement::ONE);

        let fri = Fri::new(
            generator,
            omega,
            initial_codeword_length,
            expansion_factor,
            num_colinearity_tests,
        );

        let poly = Polynomial::new_non_zero(
            (0..degree + 1)
                .map(|i| FieldElement::new(i as u128))
                .collect(),
        );
        let domain = (0..initial_codeword_length)
            .map(|i| omega.pow(i as _))
            .collect();
        let mut codeword = poly.evaluate_on_domain(&domain);
        let mut proof_stream = ProofStream::new();

        fri.prove(&codeword, &mut proof_stream);
        let verdict = fri.verify(&mut proof_stream).unwrap();
        for (x, y) in verdict {
            assert_eq!(poly.evaluate(omega.pow(x as _)), y);
        }

        let mut proof_stream = ProofStream::new();
        for i in 0..degree / 3 {
            codeword[i] = FieldElement::ZERO;
        }
        fri.prove(&codeword, &mut proof_stream);
        assert!(fri.verify(&mut proof_stream).is_none());
    }
}
