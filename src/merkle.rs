//! Unoptimized, recursive Merkle tree implementation.

use std::ops::Deref;

use blake2::{Blake2b512, Digest};
use serde_big_array::BigArray;
use serde_derive::{Deserialize, Serialize};

pub type Hash = [u8; 64];

/// A 512-bit hash, used for Merkle tree leaves and internal nodes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct MerkleHash {
    #[serde(with = "BigArray")]
    inner: Hash,
}

impl From<[u8; 64]> for MerkleHash {
    fn from(inner: [u8; 64]) -> Self {
        Self { inner }
    }
}

impl Deref for MerkleHash {
    type Target = [u8; 64];

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

/// A Merkle tree, which owns a vector of (already hashed) leaves, from which a root hash is derived.
pub struct MerkleTree {
    leaves: Vec<MerkleHash>,
    root: MerkleHash,
}

impl MerkleTree {
    pub fn get_root(&self) -> MerkleHash {
        self.root
    }

    pub fn commit(leaves: Vec<MerkleHash>) -> Self {
        assert!(leaves.len().is_power_of_two());
        if leaves.len() == 1 {
            let root = leaves[0];
            Self { leaves, root: root }
        } else {
            let first_half = leaves[..leaves.len() / 2].to_vec();
            let second_half = leaves[leaves.len() / 2..].to_vec();
            let root: Hash = Blake2b512::digest(
                &[
                    &MerkleTree::commit(first_half).root[..],
                    &MerkleTree::commit(second_half).root[..],
                ]
                .concat(),
            )
            .into();

            Self {
                leaves,
                root: MerkleHash::from(root),
            }
        }
    }

    pub fn open(&self, index: usize) -> Vec<MerkleHash> {
        fn open_rec(index: usize, leaves: &Vec<MerkleHash>) -> Vec<MerkleHash> {
            if leaves.len() == 2 {
                vec![leaves[1 - index]]
            } else if index < leaves.len() / 2 {
                let first_half = leaves[..leaves.len() / 2].to_vec();
                let second_half = leaves[leaves.len() / 2..].to_vec();
                let mut res = open_rec(index, &first_half);
                res.push(MerkleTree::commit(second_half).root);
                res
            } else {
                let first_half = leaves[..leaves.len() / 2].to_vec();
                let second_half = leaves[leaves.len() / 2..].to_vec();
                let mut res = open_rec(index - leaves.len() / 2, &second_half);
                res.push(MerkleTree::commit(first_half).root);
                res
            }
        }

        assert_eq!(self.leaves.len() & (self.leaves.len() - 1), 0);
        assert!(index < self.leaves.len());
        open_rec(index, &self.leaves)
    }

    pub fn verify(&self, leaf: MerkleHash, path: &[MerkleHash], index: usize) -> bool {
        Self::verify_for_root(self.root, leaf, path, index)
    }

    pub fn verify_for_root(
        root: MerkleHash,
        leaf: MerkleHash,
        path: &[MerkleHash],
        index: usize,
    ) -> bool {
        assert!(index < (1 << path.len()));
        assert!(path.len() > 0);
        if path.len() == 1 {
            let is_root: Hash = if index == 0 {
                Blake2b512::digest(&[&leaf[..], &path[0][..]].concat()).into()
            } else {
                Blake2b512::digest(&[&path[0][..], &leaf[..]].concat()).into()
            };
            root.as_ref() == is_root
        } else {
            if index % 2 == 0 {
                let new_root: Hash = Blake2b512::digest(&[&leaf[..], &path[0][..]].concat()).into();
                Self::verify_for_root(root, MerkleHash::from(new_root), &path[1..], index >> 1)
            } else {
                let new_root: Hash = Blake2b512::digest(&[&path[0][..], &leaf[..]].concat()).into();
                Self::verify_for_root(root, MerkleHash::from(new_root), &path[1..], index >> 1)
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::{rngs::OsRng, Rng};

    fn random_leaf() -> MerkleHash {
        let mut leaf = [0u8; 64];
        OsRng.fill(&mut leaf);
        MerkleHash::from(leaf)
    }

    #[test]
    fn merkle() {
        const N: usize = 64;

        let mut leaves = Vec::new();
        for _ in 0..N {
            leaves.push(random_leaf());
        }
        let mut tree = MerkleTree::commit(leaves.clone());

        for i in 0..N {
            let path = tree.open(i);
            assert!(tree.verify(leaves[i], &path, i));
        }

        for i in 0..N {
            let path = tree.open(i);
            assert!(!tree.verify(random_leaf(), &path, i));
        }

        for i in 0..N {
            let path = tree.open(i);
            let mut j = i;
            while j == i {
                j = (i + 1 + (OsRng.gen_range(0..N))) % N;
            }
            assert!(!tree.verify(leaves[j], &path, i));
        }

        for i in 0..N {
            let path = tree.open(i);
            let mut j = i;
            while j == i {
                j = (i + 1 + (OsRng.gen_range(0..N))) % N;
            }
            assert!(!tree.verify(leaves[i], &path, j));
        }

        let root = tree.root;
        for i in 0..N {
            let path = tree.open(i);
            tree.root = random_leaf();
            assert!(!tree.verify(leaves[i], &path, i));
        }
        tree.root = root;

        for i in 0..N {
            let path = tree.open(i);
            for j in 0..path.len() {
                let mut fake_path = path.clone();
                fake_path[j] = random_leaf();
                assert!(!tree.verify(leaves[i], &fake_path, i));
            }
        }

        let fake_tree = MerkleTree::commit(vec![random_leaf(); N]);
        for i in 0..N {
            let path = tree.open(i);
            assert!(!fake_tree.verify(leaves[i], &path, i));
        }
    }
}
