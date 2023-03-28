//! Elementary algebra over the fixed finite field (F_p) with p = 407 * 2^119 + 1.

use std::{
    iter::Sum,
    ops::{Add, Div, Mul, Neg, Sub},
};

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::cast::ToPrimitive;
use serde_derive::{Deserialize, Serialize};

const P: u128 = 407 * (1 << 119) + 1;

fn xgcd(a: u128, b: u128) -> (u128, u128, u128) {
    let egcd = BigInt::from(a).extended_gcd(&BigInt::from(b));
    (
        egcd.gcd.mod_floor(&BigInt::from(P)).to_u128().unwrap(),
        egcd.x.mod_floor(&BigInt::from(P)).to_u128().unwrap(),
        egcd.y.mod_floor(&BigInt::from(P)).to_u128().unwrap(),
    )
}

/// A field element, represented as a 128-bit unsigned integer in [0, p).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct FieldElement(u128);

impl FieldElement {
    pub const fn new(n: u128) -> Self {
        Self(n % P)
    }

    pub const ZERO: Self = Self(0);
    pub const ONE: Self = Self(1);
    pub const GENERATOR: Self = Self(85408008396924667383611388730472331217);

    pub fn is_zero(&self) -> bool {
        self.0 == 0
    }

    pub fn to_le_bytes(self) -> [u8; 16] {
        self.0.to_le_bytes()
    }

    pub fn to_u128(self) -> u128 {
        self.0
    }

    pub fn inv(self) -> Self {
        let (_, a, _) = xgcd(self.0, P);
        assert_eq!(self.mul(Self(a % P)), Self::ONE);
        Self(a % P)
    }

    pub fn primitive_nth_root(n: u128) -> Self {
        assert!(n <= (1 << 119));
        assert!(n.is_power_of_two());
        let mut root = Self::GENERATOR;
        let mut order = 1 << 119;

        while order != n {
            root = root.mul(root);
            order = order / 2;
        }

        root
    }

    pub fn sample(random_bytes: [u8; 16]) -> Self {
        FieldElement::new(u128::from_be_bytes(random_bytes))
    }

    pub fn pow(&self, n: u128) -> Self {
        let mut acc = Self::ONE;

        for b in (0..u128::BITS).rev().map(|i| (n >> i) & 1) {
            acc = acc * acc;
            if b != 0 {
                acc = acc * (*self);
            }
        }

        acc
    }
}

impl Mul for FieldElement {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(
            BigInt::from(self.0)
                .mul(BigInt::from(rhs.0))
                .mod_floor(&BigInt::from(P))
                .to_u128()
                .unwrap(),
        )
    }
}

impl Add for FieldElement {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(
            BigInt::from(self.0)
                .add(BigInt::from(rhs.0))
                .mod_floor(&BigInt::from(P))
                .to_u128()
                .unwrap(),
        )
    }
}

impl Sub for FieldElement {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(
            BigInt::from(self.0)
                .sub(BigInt::from(rhs.0))
                .mod_floor(&BigInt::from(P))
                .to_u128()
                .unwrap(),
        )
    }
}

impl Sum for FieldElement {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ZERO, |acc, x| acc + x)
    }
}

impl Neg for FieldElement {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(
            BigInt::from(self.0)
                .neg()
                .mod_floor(&BigInt::from(P))
                .to_u128()
                .unwrap(),
        )
    }
}

impl Div for FieldElement {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self.mul(rhs.inv())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn field_ops() {
        assert_eq!(
            FieldElement(2395879873598273984792873)
                .sub(FieldElement(39959992))
                .0,
            2395879873598273944832881
        );
        assert_eq!(
            FieldElement(39959992)
                .sub(FieldElement(2395879873598273984792873))
                .0,
            270497897142227984256051138493105288336
        );
        assert_eq!(
            FieldElement(2395879873598273984792873).inv().0,
            124611250983334475158786519933932092045
        );
        assert_eq!(
            FieldElement(39959992).inv().0,
            154697311270125722245319655849728864008
        );
    }
}
