/// This macro is used to try and deconstruct an enum variant.
macro_rules! try_match {
    ($var:expr, $enum_variant:path) => {
        match $var {
            $enum_variant(v) => Some(v),
            _ => None,
        }
    };
}

/// This macro is used to check that a boolean expression is true.
/// It returns a 'Result' rather than an 'Option' so that it is '#[must_use]'.
macro_rules! exact {
    ($b:expr) => {
        $b.then_some(()).ok_or(Box::<dyn std::error::Error>::from(
            "statement should be exact, but is false",
        ))
    };
}

/// This macro is used to convert a codeword into a merkle tree leaves (hashes).
macro_rules! codeword_to_merkle {
    ($codeword:expr) => {
        $codeword
            .iter()
            .map(|fe| {
                let hash: Hash = Blake2b512::digest(fe.to_le_bytes()).into();
                MerkleHash::from(hash)
            })
            .collect()
    };
}

pub(crate) use codeword_to_merkle;
pub(crate) use exact;
pub(crate) use try_match;
