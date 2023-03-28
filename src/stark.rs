//! STARK protocol implementation.

use std::collections::BTreeMap;

use blake2::{Blake2b512, Digest};
use rand::{rngs::OsRng, RngCore};

use crate::{
    algebra::FieldElement,
    fri::Fri,
    macros::{codeword_to_merkle, try_match},
    merkle::{Hash, MerkleHash, MerkleTree},
    multivariate::MultiPolynomial,
    proof_stream::{ProofObject, ProofStream},
    univariate::Polynomial,
};

/// Parameterized STARK.
pub struct Stark {
    expansion_factor: usize,
    num_randomizers: usize,
    num_registers: usize,
    original_trace_length: usize,
    generator: FieldElement,
    omega: FieldElement,
    omicron: FieldElement,
    omicron_domain: Vec<FieldElement>,
    fri: Fri,
}

impl Stark {
    pub fn new(
        expansion_factor: usize,
        num_colinearity_checks: usize,
        security_level: usize,
        num_registers: usize,
        num_cycles: usize,
        transition_constraints_degree: usize,
    ) -> Self {
        assert!(expansion_factor.is_power_of_two());
        assert!(expansion_factor >= 4);
        assert!(num_colinearity_checks * 2 >= security_level);

        let num_randomizers = 4 * num_colinearity_checks;
        let randomized_trace_length = num_cycles + num_randomizers;
        let omicron_domain_length = 1
            << (usize::BITS
                - (randomized_trace_length * transition_constraints_degree).leading_zeros());
        let fri_domain_length = omicron_domain_length * expansion_factor;
        let omega = FieldElement::primitive_nth_root(fri_domain_length as _);
        let omicron = FieldElement::primitive_nth_root(omicron_domain_length as _);
        let mut omicron_domain = Vec::with_capacity(omicron_domain_length);
        for i in 0..omicron_domain_length {
            omicron_domain.push(omicron.pow(i as _));
        }

        Self {
            expansion_factor,
            num_randomizers,
            num_registers,
            original_trace_length: num_cycles,
            generator: FieldElement::GENERATOR,
            omega,
            omicron,
            omicron_domain,
            fri: Fri::new(
                FieldElement::GENERATOR,
                omega,
                fri_domain_length,
                expansion_factor,
                num_colinearity_checks,
            ),
        }
    }

    pub fn omicron(&self) -> FieldElement {
        self.omicron
    }

    fn transition_degree_bounds(
        &self,
        transition_constraints: &Vec<MultiPolynomial>,
    ) -> Vec<usize> {
        let point_degrees = [
            vec![1],
            vec![self.original_trace_length + self.num_randomizers - 1; 2 * self.num_registers],
        ]
        .concat();

        transition_constraints
            .iter()
            .map(|a| {
                a.clone()
                    .into_inner()
                    .iter()
                    .map(|(k, _)| k.iter().zip(point_degrees.iter()).map(|(r, l)| r * l).sum())
                    .max()
                    .unwrap()
            })
            .collect()
    }

    fn transition_quotient_degree_bounds(
        &self,
        transition_constraints: &Vec<MultiPolynomial>,
    ) -> Vec<usize> {
        self.transition_degree_bounds(transition_constraints)
            .iter()
            .map(|d| d - (self.original_trace_length - 1))
            .collect()
    }

    fn max_degree(&self, transition_constraints: &Vec<MultiPolynomial>) -> usize {
        let md = *self
            .transition_quotient_degree_bounds(transition_constraints)
            .iter()
            .max()
            .unwrap();
        (1 << (usize::BITS - md.leading_zeros())) - 1
    }

    fn transition_zerofier(&self) -> Polynomial {
        let domain = self.omicron_domain[0..(self.original_trace_length - 1)].to_vec();
        Polynomial::zerofier(&domain)
    }

    fn boundary_zerofiers(&self, boundary: &Vec<(usize, usize, FieldElement)>) -> Vec<Polynomial> {
        let mut zerofiers = Vec::new();
        for s in 0..self.num_registers {
            let points = boundary
                .iter()
                .filter(|(_, r, _)| s.eq(r))
                .map(|&(c, _, _)| self.omicron.pow(c as _))
                .collect();
            zerofiers.push(Polynomial::zerofier(&points));
        }
        zerofiers
    }

    fn boundary_interpolants(
        &self,
        boundary: &Vec<(usize, usize, FieldElement)>,
    ) -> Vec<Polynomial> {
        let mut interpolants = Vec::new();
        for s in 0..self.num_registers {
            let domain = boundary
                .iter()
                .filter(|(_, r, _)| s.eq(r))
                .map(|&(c, _, v)| (self.omicron.pow(c as _), v))
                .collect();
            interpolants.push(Polynomial::interpolate(&domain));
        }
        interpolants
    }

    fn boundary_quotient_degree_bounds(
        &self,
        boundary: &Vec<(usize, usize, FieldElement)>,
    ) -> Vec<usize> {
        let randomized_trace_degree = self.original_trace_length + self.num_randomizers - 1;
        self.boundary_zerofiers(boundary)
            .iter()
            .map(|z| randomized_trace_degree - z.degree())
            .collect()
    }

    fn sample_weights(&self, n: usize, seed: &[u8]) -> Vec<FieldElement> {
        (0..n)
            .map(|i| {
                let mut arr = [0u8; 16];
                arr.copy_from_slice(&Blake2b512::digest([seed, &i.to_le_bytes()].concat())[..16]);
                FieldElement::sample(arr)
            })
            .collect()
    }

    pub fn prove(
        &self,
        trace: &Vec<Vec<FieldElement>>,
        transition_constraints: &Vec<MultiPolynomial>,
        boundary: &Vec<(usize, usize, FieldElement)>,
        proof_stream: &mut ProofStream,
    ) -> Vec<u8> {
        let mut trace = trace.clone();
        for _ in 0..self.num_randomizers {
            trace.push(
                (0..self.num_registers)
                    .map(|_| {
                        let mut random = [0u8; 16];
                        OsRng.fill_bytes(random.as_mut());
                        FieldElement::sample(random)
                    })
                    .collect(),
            );
        }

        let trace_domain: Vec<FieldElement> =
            (0..trace.len()).map(|i| self.omicron.pow(i as _)).collect();
        let mut trace_polynomials = Vec::new();
        for s in 0..self.num_registers {
            let single_trace: Vec<FieldElement> = trace.iter().map(|row| row[s]).collect();
            trace_polynomials.push(Polynomial::interpolate(
                &trace_domain
                    .iter()
                    .zip(single_trace.iter())
                    .map(|(&d, &v)| (d, v))
                    .collect(),
            ));
        }

        let mut boundary_quotients = Vec::new();
        for s in 0..self.num_registers {
            let interpolants = self.boundary_interpolants(boundary);
            let zerofiers = self.boundary_zerofiers(boundary);
            let quotient = (trace_polynomials[s].sub(interpolants.get(s).unwrap()))
                .exact_div(zerofiers.get(s).unwrap())
                .unwrap();
            boundary_quotients.push(quotient);
        }

        let fri_domain = self.fri.evaluation_domain();
        let mut boundary_quotient_codewords = Vec::new();
        for s in 0..self.num_registers {
            let codeword = boundary_quotients[s].evaluate_on_domain(&fri_domain);
            let leaves = codeword_to_merkle!(codeword);
            let tree = MerkleTree::commit(leaves);
            proof_stream.push(ProofObject::MerkleRoot(tree.get_root()));
            boundary_quotient_codewords.push(codeword);
        }

        let polys: Vec<Polynomial> = vec![Polynomial::new_non_zero(vec![
            FieldElement::ZERO,
            FieldElement::ONE,
        ])]
        .into_iter()
        .chain(trace_polynomials.iter().cloned())
        .chain(trace_polynomials.iter().map(|tp| tp.scale_by(self.omicron)))
        .collect();
        let transition_polynomials: Vec<Polynomial> = transition_constraints
            .iter()
            .map(|a| a.evaluate_symbolic(&polys))
            .collect();

        let transition_quotients: Vec<Polynomial> = transition_polynomials
            .iter()
            .map(|tp| tp.exact_div(&self.transition_zerofier()).unwrap())
            .collect();

        let randomizer_polynomial = Polynomial::new_non_zero(
            (0..self.max_degree(transition_constraints) + 1)
                .map(|_| {
                    let mut random = [0u8; 16];
                    OsRng.fill_bytes(random.as_mut());
                    FieldElement::sample(random)
                })
                .collect(),
        );
        let randomizer_codeword = randomizer_polynomial.evaluate_on_domain(&fri_domain);
        let leaves = codeword_to_merkle!(randomizer_codeword);
        let tree = MerkleTree::commit(leaves);
        proof_stream.push(ProofObject::MerkleRoot(tree.get_root()));

        assert!(transition_quotients
            .iter()
            .zip(
                self.transition_quotient_degree_bounds(transition_constraints)
                    .iter()
            )
            .all(|(tq, &db)| tq.degree() == db));

        let x = Polynomial::new_non_zero(vec![FieldElement::ZERO, FieldElement::ONE]);
        let mut terms = vec![randomizer_polynomial];
        for i in 0..transition_quotients.len() {
            terms.push(transition_quotients[i].clone());
            let shift = self.max_degree(transition_constraints)
                - self.transition_quotient_degree_bounds(transition_constraints)[i];
            terms.push((x.clone().pow(shift as _)).mul(transition_quotients.get(i).unwrap()));
        }
        for i in 0..self.num_registers {
            terms.push(boundary_quotients[i].clone());
            let shift = self.max_degree(transition_constraints)
                - self.boundary_quotient_degree_bounds(boundary)[i];
            terms.push((x.clone().pow(shift as _)).mul(boundary_quotients.get(i).unwrap()));
        }

        let weights = self.sample_weights(
            1 + 2 * transition_quotients.len() + 2 * boundary_quotients.len(),
            &proof_stream.prover_fiat_shamir(),
        );

        let combination = terms
            .iter()
            .zip(weights.iter())
            .fold(Polynomial::zero(), |acc, (term, &weight)| {
                acc.add(&Polynomial::new_non_zero(vec![weight]).mul(term))
            });
        let combined_codeword = combination.evaluate_on_domain(&fri_domain);

        let indices = self.fri.prove(&combined_codeword, proof_stream);
        let duplicated_indices: Vec<usize> = indices
            .iter()
            .map(|i| *i)
            .chain(
                indices
                    .iter()
                    .map(|i| (i + self.expansion_factor) % self.fri.get_domain_length()),
            )
            .collect();
        let mut quadrupled_indices: Vec<usize> = duplicated_indices
            .iter()
            .map(|i| *i)
            .chain(
                duplicated_indices
                    .iter()
                    .map(|i| (i + self.fri.get_domain_length() / 2) % self.fri.get_domain_length()),
            )
            .collect();
        quadrupled_indices.sort();

        for bqc in boundary_quotient_codewords {
            for &i in quadrupled_indices.iter() {
                proof_stream.push(ProofObject::Point(bqc[i]));
                let leaves = codeword_to_merkle!(bqc);
                let tree = MerkleTree::commit(leaves);
                let path = tree.open(i);
                proof_stream.push(ProofObject::MerklePath(path));
            }
        }

        for &i in quadrupled_indices.iter() {
            proof_stream.push(ProofObject::Point(randomizer_codeword[i]));
            let leaves = codeword_to_merkle!(randomizer_codeword);
            let tree = MerkleTree::commit(leaves);
            let path = tree.open(i);
            proof_stream.push(ProofObject::MerklePath(path));
        }

        proof_stream.serialize()
    }

    pub fn verify(
        &self,
        proof: &Vec<u8>,
        transition_constraints: &Vec<MultiPolynomial>,
        boundary: &Vec<(usize, usize, FieldElement)>,
        prefix: Option<Vec<u8>>,
    ) -> bool {
        if self.original_trace_length != 1 + boundary.iter().map(|(c, _, _)| c).max().unwrap() {
            return false;
        }
        let mut proof_stream = ProofStream::deserialize(prefix, &proof);

        let mut boundary_quotient_roots = Vec::new();
        for _ in 0..self.num_registers {
            boundary_quotient_roots
                .push(try_match!(proof_stream.pull().unwrap(), ProofObject::MerkleRoot).unwrap());
        }

        let randomizer_root =
            try_match!(proof_stream.pull().unwrap(), ProofObject::MerkleRoot).unwrap();
        let weights = self.sample_weights(
            1 + 2 * transition_constraints.len() + 2 * self.boundary_interpolants(boundary).len(),
            &proof_stream.verifier_fiat_shamir(),
        );

        let opt_pv = self.fri.verify(&mut proof_stream);
        if opt_pv.is_none() {
            return false;
        }
        let mut polynomial_values = opt_pv.unwrap();
        polynomial_values.sort_by_key(|iv| iv.0);

        let indices: Vec<usize> = polynomial_values.iter().map(|&(i, _)| i).collect();
        let values: Vec<FieldElement> = polynomial_values.iter().map(|&(_, v)| v).collect();
        let mut duplicated_indices: Vec<usize> = indices
            .iter()
            .map(|i| *i)
            .chain(
                indices
                    .iter()
                    .map(|i| (i + self.expansion_factor) % self.fri.get_domain_length()),
            )
            .collect();
        duplicated_indices.sort();

        let mut leaves_vector = Vec::new();
        for &r in boundary_quotient_roots.iter() {
            let mut leaves = BTreeMap::new();
            for &i in duplicated_indices.iter() {
                let fe = try_match!(proof_stream.pull().unwrap(), ProofObject::Point).unwrap();
                let hash: Hash = Blake2b512::digest(fe.to_le_bytes()).into();
                let leaf = MerkleHash::from(hash);
                leaves.insert(i, fe);
                let path =
                    try_match!(proof_stream.pull().unwrap(), ProofObject::MerklePath).unwrap();
                if !MerkleTree::verify_for_root(r, leaf, &path, i) {
                    return false;
                }
            }
            leaves_vector.push(leaves);
        }

        let mut randomizer = BTreeMap::new();
        for &i in duplicated_indices.iter() {
            let fe = try_match!(proof_stream.pull().unwrap(), ProofObject::Point).unwrap();
            let hash: Hash = Blake2b512::digest(fe.to_le_bytes()).into();
            let leaf = MerkleHash::from(hash);
            randomizer.insert(i, fe);
            let path = try_match!(proof_stream.pull().unwrap(), ProofObject::MerklePath).unwrap();
            if !MerkleTree::verify_for_root(randomizer_root, leaf, &path, i) {
                return false;
            }
        }

        let zerofiers = self.boundary_zerofiers(boundary);
        let interpolants = self.boundary_interpolants(boundary);

        for (i, &cur_index) in indices.iter().enumerate() {
            let domain_current_index = self.generator * (self.omega.pow(cur_index as _));
            let next_index = (cur_index + self.expansion_factor) % self.fri.get_domain_length();
            let domain_next_index = self.generator * (self.omega.pow(next_index as _));
            let mut current_trace = vec![FieldElement::ZERO; self.num_registers];
            let mut next_trace = current_trace.clone();

            for s in 0..self.num_registers {
                current_trace[s] = leaves_vector[s][&cur_index]
                    * zerofiers[s].evaluate(domain_current_index)
                    + interpolants[s].evaluate(domain_current_index);
                next_trace[s] = leaves_vector[s][&next_index]
                    * zerofiers[s].evaluate(domain_next_index)
                    + interpolants[s].evaluate(domain_next_index);
            }

            let point = [vec![domain_current_index], current_trace, next_trace].concat();
            let transition_constraints_values: Vec<FieldElement> = transition_constraints
                .iter()
                .map(|p| p.evaluate(&point))
                .collect();

            let mut terms = vec![randomizer[&cur_index]];
            let transition_quotient_degree_bounds =
                self.transition_quotient_degree_bounds(&transition_constraints);
            let boundary_quotient_degree_bounds = self.boundary_quotient_degree_bounds(boundary);

            for (s, tcv) in transition_constraints_values.iter().enumerate() {
                let quotient = *tcv / self.transition_zerofier().evaluate(domain_current_index);
                terms.push(quotient);
                let shift =
                    self.max_degree(&transition_constraints) - transition_quotient_degree_bounds[s];
                terms.push(quotient * (domain_current_index.pow(shift as _)));
            }
            for (s, bqv) in leaves_vector.iter().map(|l| l[&cur_index]).enumerate() {
                terms.push(bqv);
                let shift =
                    self.max_degree(&transition_constraints) - boundary_quotient_degree_bounds[s];
                terms.push(bqv * (domain_current_index.pow(shift as _)));
            }

            let combination = terms
                .iter()
                .zip(weights.iter())
                .map(|(&t, &w)| t * w)
                .fold(FieldElement::ZERO, |a, b| a + b);

            if combination != *values.get(i).unwrap_or(&FieldElement::ZERO) {
                return false;
            }
        }

        true
    }
}

#[cfg(test)]
mod test {
    use crate::{
        proof_stream,
        rescue_prime::{RescuePrime, M, N},
    };

    use super::*;

    #[test]
    fn rescue_prime() {
        let expansion_factor = 4;
        let num_colinearity_checks = 2;
        let security_level = 2;
        let rp = RescuePrime::new();
        let mut output = FieldElement::new(228894434762048332457318);

        for _ in 0..20 {
            let input = output;
            output = rp.hash(input);
            let num_cycles = N + 1;
            let state_width = M;

            let stark = Stark::new(
                expansion_factor,
                num_colinearity_checks,
                security_level,
                state_width,
                num_cycles,
                2,
            );

            let trace = rp.trace(input);
            let air = rp.transition_constraints(stark.omicron());
            let boundary = rp.boundary_constraints(output);
            let mut proof_stream = proof_stream::ProofStream::new();
            let proof = stark.prove(&trace, &air, &boundary, &mut proof_stream);
            assert!(stark.verify(&proof, &air, &boundary, None));

            let false_output = output + FieldElement::ONE;
            let false_boundary = rp.boundary_constraints(false_output);
            assert!(!stark.verify(&proof, &air, &false_boundary, None));
        }
    }
}
