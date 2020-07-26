/// Representation of a qubit.
#[derive(Debug, PartialEq)]
pub struct Qubit(pub f64, pub f64);

impl Qubit {
    pub fn zero() -> Self {
        Qubit(1.0, 0.0)
    }

    pub fn one() -> Self {
        Qubit(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_qubits() {
        assert_eq!(Qubit::zero().0, 1.0);
        assert_eq!(Qubit::zero().1, 0.0);
        assert_eq!(Qubit::one().0, 0.0);
        assert_eq!(Qubit::one().1, 1.0);
    }
}
