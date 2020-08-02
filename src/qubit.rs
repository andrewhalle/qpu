use approx::abs_diff_eq;
use num::complex::Complex64;
use num::FromPrimitive;

/// Representation of a qubit.
#[derive(Debug)]
pub struct Qubit(pub Complex64, pub Complex64);

impl Qubit {
    pub fn zero() -> Self {
        Qubit(
            FromPrimitive::from_f64(1.0).unwrap(),
            FromPrimitive::from_f64(0.0).unwrap(),
        )
    }

    pub fn one() -> Self {
        Qubit(
            FromPrimitive::from_f64(0.0).unwrap(),
            FromPrimitive::from_f64(1.0).unwrap(),
        )
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

    #[test]
    fn init_qubits() {
        assert_eq!(Qubit::zero().0, FromPrimitive::from_f64(1.0).unwrap());
        assert_eq!(Qubit::zero().1, FromPrimitive::from_f64(0.0).unwrap());
        assert_eq!(Qubit::one().0, FromPrimitive::from_f64(0.0).unwrap());
        assert_eq!(Qubit::one().1, FromPrimitive::from_f64(1.0).unwrap());
    }
}
