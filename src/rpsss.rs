//! Signature scheme implementation based on hashing and STARKs.

use rand::{rngs::OsRng, Rng};

use crate::{
    algebra::FieldElement,
    proof_stream::{self, ProofStream},
    rescue_prime::RescuePrime,
    stark::Stark,
};

/// Rescue Prime STARK Signature Scheme.
pub struct RPSSS {
    rp: RescuePrime,
    stark: Stark,
}

impl RPSSS {
    pub fn new() -> Self {
        let rp = RescuePrime::new();
        let stark = Stark::new(
            4,
            64,
            2 * 64,
            crate::rescue_prime::M,
            crate::rescue_prime::N + 1,
            3,
        );
        Self { rp, stark }
    }

    pub fn stark_prove(
        &self,
        input: FieldElement,
        proof_stream: &mut proof_stream::ProofStream,
    ) -> Vec<u8> {
        let output = self.rp.hash(input);
        let trace = self.rp.trace(input);
        let transition_constraints = self.rp.transition_constraints(self.stark.omicron());
        let boundary_constraints = self.rp.boundary_constraints(output);

        self.stark.prove(
            &trace,
            &transition_constraints,
            &boundary_constraints,
            proof_stream,
        )
    }

    pub fn stark_verify(
        &self,
        output: FieldElement,
        stark_proof: Vec<u8>,
        prefix: Vec<u8>,
    ) -> bool {
        let boundary_constraints = self.rp.boundary_constraints(output);
        let transition_constraints = self.rp.transition_constraints(self.stark.omicron());

        self.stark.verify(
            &stark_proof,
            &transition_constraints,
            &boundary_constraints,
            Some(prefix),
        )
    }

    pub fn keygen(&self) -> (FieldElement, FieldElement) {
        let mut arr = [0u8; 16];
        OsRng.fill(&mut arr);
        let sk = FieldElement::sample(arr);
        let pk = self.rp.hash(sk);
        (sk, pk)
    }

    pub fn sign(&self, sk: FieldElement, prefix: Vec<u8>) -> Vec<u8> {
        let mut sps = ProofStream::with_prefix(prefix);
        self.stark_prove(sk, &mut sps)
    }

    pub fn verify(&self, pk: FieldElement, signature: Vec<u8>, prefix: Vec<u8>) -> bool {
        self.stark_verify(pk, signature, prefix)
    }
}

#[cfg(test)]
mod test {
    use blake2::{Blake2s256, Digest};

    use super::*;

    #[test]
    fn rpsss() {
        let rpsss = RPSSS::new();
        let (sk, pk) = rpsss.keygen();
        let prefix = Blake2s256::digest("Hello world!").to_vec();
        let fake_prefix = Blake2s256::digest("Bad message!").to_vec();
        let signature = rpsss.sign(sk, prefix.clone());
        assert!(rpsss.verify(pk, signature.clone(), prefix));
        assert!(!rpsss.verify(pk, signature, fake_prefix));
    }
}
