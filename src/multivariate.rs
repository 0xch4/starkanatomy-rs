//! Basic multivariate polynomial algebra.

use std::{cmp::max, collections::BTreeMap, ops::Neg, vec};

use crate::{algebra::*, univariate::Polynomial};

/// A multivariate polynomial, represented internally as a map from vectors of exponents
/// (each of length the number of variables) to the corresponding scalar coefficients.
/// Keys (exponents) that would otherwise be mapped to a ZERO coefficient are omitted.
#[derive(Clone, Debug)]
pub struct MultiPolynomial(BTreeMap<Vec<usize>, FieldElement>);

impl MultiPolynomial {
    pub const ZERO: Self = Self(BTreeMap::new());

    pub fn constant(fe: FieldElement, nvar: usize) -> Self {
        Self(BTreeMap::from([(vec![0; nvar], fe)]))
    }

    pub fn into_inner(self) -> BTreeMap<Vec<usize>, FieldElement> {
        self.0
    }

    fn num_variables(&self) -> usize {
        self.0.first_key_value().map_or(0, |(k, _)| k.len())
    }

    fn canonize(&mut self) {
        self.0.retain(|_, &mut fe| fe != FieldElement::ZERO);
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    pub fn univariate_variables(nvar: usize) -> Vec<Self> {
        assert!(nvar > 0);
        let mut res = vec![Self::ZERO; nvar];

        for i in 0..nvar {
            let mut ith_var = vec![0; nvar];
            ith_var[i] = 1;
            res[i].0.insert(ith_var, FieldElement::ONE);
        }

        res
    }

    pub fn lift_univariate(poly: &Polynomial, nvar: usize, index: usize) -> Self {
        assert!(nvar > 0);
        if poly.is_zero() {
            return Self::ZERO;
        }
        let univariates = Self::univariate_variables(nvar);
        let x = univariates[index].clone();
        let mut acc = Self::ZERO;
        let coeffs = poly.clone().into_inner();

        for i in 0..coeffs.len() {
            acc = acc.add(&MultiPolynomial::constant(coeffs[i], nvar).mul(&x.pow(i as _)));
        }

        acc
    }

    pub fn evaluate(&self, fes: &Vec<FieldElement>) -> FieldElement {
        let mut acc = FieldElement::ZERO;

        for (k, &v) in self.0.iter() {
            let mut prod = v;
            for i in 0..k.len() {
                prod = prod * fes[i].pow(k[i] as _);
            }
            acc = acc + prod;
        }

        acc
    }

    pub fn evaluate_symbolic(&self, pols: &Vec<Polynomial>) -> Polynomial {
        let mut acc = Polynomial::zero();

        for (k, &v) in self.0.iter() {
            let mut prod = Polynomial::new_non_zero(vec![v]);
            for i in 0..k.len() {
                prod = prod.mul(&pols[i].pow(k[i]));
            }
            acc = acc.add(&prod);
        }

        acc
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        let mut btree = BTreeMap::new();
        let nvar = max(self.num_variables(), rhs.num_variables());

        for (lk, &lv) in self.0.iter() {
            for (rk, &rv) in rhs.0.iter() {
                let mut vars = vec![0; nvar];
                for i in 0..lk.len() {
                    vars[i] += lk[i];
                }
                for i in 0..rk.len() {
                    vars[i] += rk[i];
                }
                if let Some(&prev) = btree.get(&vars) {
                    btree.insert(vars, prev + lv * rv);
                } else {
                    btree.insert(vars, lv * rv);
                }
            }
        }

        let mut res = Self(btree);
        res.canonize();
        res
    }

    pub fn add(&self, rhs: &Self) -> Self {
        let mut btree = BTreeMap::new();
        let nvar = max(self.num_variables(), rhs.num_variables());

        for (lk, &lv) in self.0.iter() {
            let mut newk = lk.clone();
            newk.append(&mut vec![0; nvar - lk.len()]);
            btree.insert(newk, lv);
        }
        for (rk, &rv) in rhs.0.iter() {
            let mut newk = rk.clone();
            newk.append(&mut vec![0; nvar - rk.len()]);
            if let Some(&prev) = btree.get(&newk) {
                btree.insert(newk, prev + rv);
            } else {
                btree.insert(newk, rv);
            }
        }

        let mut res = Self(btree);
        res.canonize();
        res
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }

    pub fn neg(&self) -> Self {
        let mut btree = BTreeMap::new();
        for (k, v) in self.0.iter() {
            btree.insert(k.clone(), v.neg());
        }
        Self(btree)
    }

    pub fn pow(&self, n: u128) -> Self {
        if self.is_zero() {
            return Self::ZERO;
        }
        let nvar = self.num_variables();
        let vars = vec![0; nvar];
        let mut acc = Self::ZERO;
        acc.0.insert(vars, FieldElement::ONE);

        for b in (0..u128::BITS).rev().map(|i| (n >> i) & 1) {
            acc = acc.mul(&acc);
            if b != 0 {
                acc = acc.mul(self);
            }
        }

        acc
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const ZERO: FieldElement = FieldElement::ZERO;
    const ONE: FieldElement = FieldElement::ONE;
    const TWO: FieldElement = FieldElement::new(2);
    const FIVE: FieldElement = FieldElement::new(5);

    #[test]
    fn evaluate() {
        let univariate_vars = MultiPolynomial::univariate_variables(4);

        let mpoly_a = MultiPolynomial::constant(ONE, 4)
            .mul(&univariate_vars[0].clone())
            .add(&MultiPolynomial::constant(TWO, 4).mul(&univariate_vars[1].clone()))
            .add(&MultiPolynomial::constant(FIVE, 4).mul(&univariate_vars[2].clone().pow(3)));

        let mpoly_b = MultiPolynomial::constant(ONE, 4)
            .mul(&univariate_vars[0].clone())
            .mul(&univariate_vars[3].clone())
            .add(&MultiPolynomial::constant(FIVE, 4).mul(&univariate_vars[3].clone().pow(3)))
            .add(&MultiPolynomial::constant(FIVE, 4));

        let mpoly_c = mpoly_a.mul(&mpoly_b);
        let fes = vec![ZERO, FIVE, FIVE, ZERO];

        let eval_a = mpoly_a.evaluate(&fes);
        let eval_b = mpoly_b.evaluate(&fes);
        let eval_c = mpoly_c.evaluate(&fes);

        assert_eq!(eval_a * eval_b, eval_c);
        assert_eq!(eval_a + eval_b, mpoly_a.add(&mpoly_b).evaluate(&fes));
    }

    #[test]
    fn lift() {
        let upoly = Polynomial::interpolate(&vec![(ZERO, TWO), (ONE, FIVE), (TWO, FIVE)]);
        let mpoly = MultiPolynomial::lift_univariate(&upoly, 4, 3);

        assert_eq!(
            upoly.evaluate(FIVE),
            mpoly.evaluate(&vec![ZERO, ZERO, ZERO, FIVE])
        );
    }
}
