//! Basic univariate polynomial algebra.

use std::{collections::BTreeSet, ops::Neg};

use crate::algebra::*;

/// A univariate polynomial, internally represented as a vector of coefficients.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polynomial(Vec<FieldElement>);

impl Polynomial {
    pub fn zero() -> Self {
        Self(vec![FieldElement::ZERO])
    }

    pub fn one() -> Self {
        Self(vec![FieldElement::ONE])
    }

    pub fn into_inner(self) -> Vec<FieldElement> {
        self.0
    }

    fn canonize(&mut self) {
        while self.0.len() > 1 && self.0.last().is_some_and(|l| l.is_zero()) {
            self.0.pop();
        }
    }

    pub fn is_zero(&self) -> bool {
        self.0.len() == 1 && self.0.first().is_some_and(|fe| fe.is_zero())
    }

    pub fn new_non_zero(coeffs: Vec<FieldElement>) -> Self {
        assert!(*coeffs.last().unwrap() != FieldElement::ZERO);
        let mut res = Self(coeffs);
        res.canonize();
        res
    }

    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }

    pub fn leading_coeff(&self) -> FieldElement {
        *self.0.last().unwrap()
    }

    pub fn evaluate(&self, x: FieldElement) -> FieldElement {
        let mut monomial = FieldElement::ONE;
        let mut res = FieldElement::ZERO;

        for &c in self.0.iter() {
            res = res + c * monomial;
            monomial = monomial * x
        }

        res
    }

    pub fn evaluate_on_domain(&self, domain: &Vec<FieldElement>) -> Vec<FieldElement> {
        domain.into_iter().map(|&fe| self.evaluate(fe)).collect()
    }

    pub fn interpolate(domain: &Vec<(FieldElement, FieldElement)>) -> Self {
        assert!(!domain.is_empty());
        assert_eq!(
            BTreeSet::from_iter(domain.iter().map(|(d, _)| d)).len(),
            domain.len()
        );
        let x = Self(vec![FieldElement::ZERO, FieldElement::ONE]);
        let mut acc = Self::zero();

        for i in 0..domain.len() {
            let mut prod = Self(vec![domain[i].1]);
            for j in 0..domain.len() {
                if i != j {
                    prod = prod.mul(
                        &x.sub(&Self(vec![domain[j].0]))
                            .mul(&Self(vec![(domain[i].0 - domain[j].0).inv()])),
                    );
                }
            }
            acc = acc.add(&prod)
        }

        acc.canonize();
        acc
    }

    pub fn zerofier(domain: &Vec<FieldElement>) -> Self {
        let x = Self(vec![FieldElement::ZERO, FieldElement::ONE]);
        let mut acc = Self(vec![FieldElement::ONE]);

        for &d in domain.iter() {
            acc = acc.mul(&x.sub(&Self(vec![d])));
        }

        acc.canonize();
        acc
    }

    pub fn scale_by(&self, factor: FieldElement) -> Self {
        let mut res = self.0.clone();
        let mut f = FieldElement::ONE;

        for i in 0..res.len() {
            res[i] = res[i] * f;
            f = f * factor;
        }

        let mut res = Self(res);
        res.canonize();
        res
    }

    pub fn test_colinearity(domain: Vec<(FieldElement, FieldElement)>) -> bool {
        Self::interpolate(&domain).degree() <= 1
    }

    pub fn divide_with_remainder(&self, d: &Self) -> (Self, Self) {
        assert!(!d.is_zero());
        if self.degree() < d.degree() {
            return (Self::zero(), self.clone());
        }

        let mut r = self.clone();
        let mut q_coeffs = vec![FieldElement::ZERO; self.degree() - d.degree() + 1];
        for _ in 0..q_coeffs.len() {
            if r.degree() < d.degree() {
                break;
            }
            let c = r.leading_coeff() / d.leading_coeff();
            let mut to_sub = vec![FieldElement::ZERO; r.degree() - d.degree() + 1];
            to_sub[r.degree() - d.degree()] = c;
            let to_sub = Self(to_sub).mul(d);
            q_coeffs[r.degree() - d.degree()] = c;
            r = r.sub(&to_sub);
            r.canonize();
        }

        let mut q = Self(q_coeffs);
        q.canonize();
        r.canonize();
        (q, r)
    }

    pub fn modulo(&self, rhs: &Self) -> Self {
        let (_, r) = self.divide_with_remainder(rhs);
        r
    }

    pub fn neg(&self) -> Self {
        Self(self.0.clone().iter_mut().map(|c| c.neg()).collect())
    }

    pub fn add(&self, rhs: &Self) -> Self {
        if self.0.len() >= rhs.0.len() {
            let mut coeffs = self.0.clone();
            for i in 0..rhs.0.len() {
                coeffs[i] = coeffs[i] + rhs.0[i];
            }

            let mut res = Self(coeffs);
            res.canonize();
            res
        } else {
            let mut res = rhs.add(self);
            res.canonize();
            res
        }
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        let mut coeffs = vec![FieldElement::ZERO; self.0.len() + rhs.0.len() - 1];
        for i in 0..self.0.len() {
            if !self.0[i].is_zero() {
                for j in 0..rhs.0.len() {
                    coeffs[i + j] = coeffs[i + j] + (self.0[i] * rhs.0[j]);
                }
            }
        }

        let mut res = Self(coeffs);
        res.canonize();
        res
    }

    pub fn exact_div(&self, rhs: &Self) -> Option<Self> {
        let (q, r) = self.divide_with_remainder(rhs);
        r.is_zero().then(|| q)
    }

    pub fn pow(&self, n: usize) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        if n == 0 {
            return Self::one();
        }

        let mut acc = Self::one();
        for b in (0..usize::BITS).rev().map(|i| (n >> i) & 1) {
            acc = acc.mul(&acc);
            if b != 0 {
                acc = acc.mul(self);
            }
        }

        acc.canonize();
        acc
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::OsRng, Rng, RngCore};

    const ZERO: FieldElement = FieldElement::ZERO;
    const ONE: FieldElement = FieldElement::ONE;
    const TWO: FieldElement = FieldElement::new(2);
    const THREE: FieldElement = FieldElement::new(3);
    const FOUR: FieldElement = FieldElement::new(4);
    const FIVE: FieldElement = FieldElement::new(5);

    #[test]
    fn distributivity() {
        let a = Polynomial::new_non_zero(vec![ONE, ZERO, FIVE, TWO]);
        let b = Polynomial::new_non_zero(vec![TWO, TWO, ONE]);
        let c = Polynomial::new_non_zero(vec![ZERO, FIVE, TWO, FIVE, FIVE, ONE]);
        assert_eq!(a.mul(&b.add(&c)), a.mul(&b).add(&a.mul(&c)));
    }

    #[test]
    fn division() {
        let a = Polynomial::new_non_zero(vec![ONE, ZERO, FIVE, TWO]);
        let b = Polynomial::new_non_zero(vec![TWO, TWO, ONE]);
        let c = Polynomial::new_non_zero(vec![ZERO, FIVE, TWO, FIVE, FIVE, ONE]);

        let (q, r) = a.mul(&b).divide_with_remainder(&a);
        assert!(r.is_zero());
        assert_eq!(q, b);

        let (q, r) = a.mul(&b).divide_with_remainder(&b);
        assert!(r.is_zero());
        assert_eq!(q, a);

        let (q, r) = a.mul(&b).divide_with_remainder(&c);
        assert!(!r.is_zero());
        assert_eq!(q.mul(&c).add(&r), a.mul(&b));
    }

    #[test]
    fn interpolate() {
        let mut domain = vec![
            (ZERO, FIVE),
            (ONE, TWO),
            (TWO, TWO),
            (THREE, ONE),
            (FOUR, FIVE),
            (FIVE, ZERO),
        ];
        let lagrange = Polynomial::interpolate(&domain);

        for &(x, y) in domain.iter() {
            assert_eq!(lagrange.evaluate(x), y);
        }

        assert!(lagrange.evaluate(FieldElement::new(363)) != ZERO);
        assert_eq!(lagrange.degree(), domain.len() - 1);

        domain.push((
            FieldElement::new(363),
            lagrange.evaluate(FieldElement::new(363)),
        ));
        let lagrange = Polynomial::interpolate(&domain);
        assert_eq!(lagrange.degree(), domain.len() - 2);
    }

    #[test]
    fn random_fe() {
        let mut random_bytes = [0u8; 16];
        OsRng.fill(&mut random_bytes);
    }

    #[test]
    fn zerofier() {
        fn random_fe() -> FieldElement {
            let mut random_bytes = [0u8; 16];
            OsRng.fill(&mut random_bytes);
            FieldElement::sample(random_bytes)
        }

        for _ in 0..10 {
            let degree = OsRng.next_u64() as u8;
            let mut domain = Vec::new();

            while domain.len() != degree as usize {
                let fe = random_fe();
                if !domain.contains(&fe) {
                    domain.push(fe);
                }
            }

            let zerofier = Polynomial::zerofier(&domain);
            assert_eq!(zerofier.degree(), degree as usize);

            for &d in domain.iter() {
                assert_eq!(zerofier.evaluate(d), FieldElement::ZERO);
            }

            let mut randfe = random_fe();
            while domain.contains(&randfe) {
                randfe = random_fe();
            }
            assert_ne!(zerofier.evaluate(randfe), FieldElement::ZERO);
        }
    }
}
