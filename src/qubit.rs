use approx::abs_diff_eq;
use num::complex::Complex64;
use num::Complex;
use std::fmt::{Display, Formatter};

#[derive(Debug, PartialEq)]
pub enum Binary {
    Zero,
    One,
}

impl Display for Binary {
    fn fmt(&self, f: &mut Formatter) -> Result<(), std::fmt::Error> {
        match self {
            Binary::Zero => write!(f, "0"),
            Binary::One => write!(f, "1"),
        }
    }
}

/// Representation of a qubit.
#[derive(Debug)]
pub struct Qubit(pub Complex64, pub Complex64);

impl Qubit {
    pub fn zero() -> Self {
        Qubit(Complex::new(1.0, 0.0), Complex::new(0.0, 0.0))
    }

    pub fn one() -> Self {
        Qubit(Complex::new(0.0, 0.0), Complex::new(1.0, 0.0))
    }

    pub fn read(&mut self, measurement: f64) -> Binary {
        if measurement < self.0.norm_sqr() {
            self.0 = Complex::new(1.0, 0.0);
            self.1 = Complex::new(0.0, 0.0);

            Binary::Zero
        } else {
            self.0 = Complex::new(0.0, 0.0);
            self.1 = Complex::new(1.0, 0.0);

            Binary::One
        }
    }
}

impl PartialEq for Qubit {
    fn eq(&self, other: &Self) -> bool {
        abs_diff_eq!(self.0.re, other.0.re)
            && abs_diff_eq!(self.0.im, other.0.im)
            && abs_diff_eq!(self.1.re, other.1.re)
            && abs_diff_eq!(self.1.im, other.1.im)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_1_SQRT_2;

    #[test]
    fn init_qubits() {
        assert_eq!(Qubit::zero().0, Complex::new(1.0, 0.0));
        assert_eq!(Qubit::zero().1, Complex::new(0.0, 0.0));
        assert_eq!(Qubit::one().0, Complex::new(0.0, 0.0));
        assert_eq!(Qubit::one().1, Complex::new(1.0, 0.0));
    }

    #[test]
    fn read() {
        let mut q = Qubit(
            Complex::new(FRAC_1_SQRT_2, 0.0),
            Complex::new(FRAC_1_SQRT_2, 0.0),
        );

        assert_eq!(q.read(0.49), Binary::Zero);

        let mut q = Qubit(
            Complex::new(FRAC_1_SQRT_2, 0.0),
            Complex::new(FRAC_1_SQRT_2, 0.0),
        );

        assert_eq!(q.read(0.51), Binary::One);
    }
}
