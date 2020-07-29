use num::complex::Complex64;
use num::FromPrimitive;

/// Representation of a qubit.
#[derive(Debug, PartialEq)]
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
