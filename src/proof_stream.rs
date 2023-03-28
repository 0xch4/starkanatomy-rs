//! Proof stream implementation for proofs made non-interactive via the Fiat-Shamir transform.

#[cfg(test)]
use std::collections::BTreeMap;

use serde_derive::{Deserialize, Serialize};
use sha3::{Digest, Sha3_256};

use crate::{algebra::FieldElement, merkle::MerkleHash};

/// A proof object, pushed and pulled from a `ProofStream`.
/// The test version is only there to reproduce the tutorials' test cases.
#[cfg(test)]
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub enum ProofObject {
    Bytes(Vec<u8>),
    Map(BTreeMap<u8, u8>),
    BytesVec(Vec<Vec<u8>>),
    Point(FieldElement),
    MerkleRoot(MerkleHash),
    LastCodeword(Vec<FieldElement>),
    ColinearityTest((FieldElement, FieldElement, FieldElement)),
    MerklePath(Vec<MerkleHash>),
}

#[cfg(not(test))]
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub enum ProofObject {
    Point(FieldElement),
    MerkleRoot(MerkleHash),
    LastCodeword(Vec<FieldElement>),
    ColinearityTest((FieldElement, FieldElement, FieldElement)),
    MerklePath(Vec<MerkleHash>),
}

/// A proof stream, supporting pushing and pulling of `ProofObject`s, and generation of Fiat-Shamir challenges.
/// The optional `prefix` is used to tie the Fiat-Shamir challenges to a broader context, e.g. a specific document in the case of `RPSSS`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProofStream {
    prefix: Option<Vec<u8>>,
    objects: Vec<ProofObject>,
    read_index: usize,
}

impl ProofStream {
    pub fn new() -> Self {
        Self {
            prefix: None,
            objects: Vec::new(),
            read_index: 0,
        }
    }

    pub fn with_prefix(prefix: Vec<u8>) -> Self {
        Self {
            prefix: Some(prefix),
            objects: Vec::new(),
            read_index: 0,
        }
    }

    pub fn push(&mut self, object: ProofObject) {
        self.objects.push(object);
    }

    pub fn pull(&mut self) -> Option<ProofObject> {
        self.read_index += 1;
        self.objects.get(self.read_index - 1).cloned()
    }

    pub fn serialize(&self) -> Vec<u8> {
        bincode::serialize(&self.objects).unwrap()
    }

    pub fn deserialize(prefix: Option<Vec<u8>>, bytes: &Vec<u8>) -> Self {
        let objects: Vec<ProofObject> = bincode::deserialize(bytes).unwrap();
        Self {
            prefix,
            objects,
            read_index: 0,
        }
    }

    pub fn prover_fiat_shamir(&self) -> [u8; 32] {
        if let Some(prefix) = &self.prefix {
            Sha3_256::digest(&[prefix.clone(), self.serialize()].concat()).into()
        } else {
            Sha3_256::digest(&self.serialize()).into()
        }
    }

    pub fn verifier_fiat_shamir(&self) -> [u8; 32] {
        if let Some(prefix) = &self.prefix {
            Sha3_256::digest(
                &[
                    prefix.clone(),
                    bincode::serialize(&self.objects[..self.read_index].to_vec()).unwrap(),
                ]
                .concat(),
            )
            .into()
        } else {
            Sha3_256::digest(
                &bincode::serialize(&self.objects[..self.read_index].to_vec()).unwrap(),
            )
            .into()
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn serialize() {
        let mut proof1 = ProofStream::new();
        proof1.push(ProofObject::Bytes(vec![1]));
        proof1.push(ProofObject::Map(BTreeMap::from([(1, 1)])));
        proof1.push(ProofObject::BytesVec(vec![vec![1]]));
        proof1.push(ProofObject::Bytes(vec![2]));

        let serialized = proof1.serialize();
        let mut proof2 = ProofStream::deserialize(None, &serialized);

        assert_eq!(proof1.pull(), proof2.pull());
        assert_eq!(proof1.pull(), proof2.pull());
        assert_eq!(proof1.pull(), proof2.pull());
        assert_eq!(proof1.pull(), Some(ProofObject::Bytes(vec![2])));
        assert_eq!(proof2.pull(), Some(ProofObject::Bytes(vec![2])));
        assert_eq!(proof1.prover_fiat_shamir(), proof2.prover_fiat_shamir());
    }
}
