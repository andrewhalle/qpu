use rand::prelude::*;

#[derive(Debug, PartialEq)]
pub enum Binary {
    Zero,
    One,
}

/// Representation of a qubit.
#[derive(Debug, PartialEq)]
pub struct Qubit(f64, f64);

impl Qubit {
    pub fn zero() -> Self {
        Qubit(1.0, 0.0)
    }

    pub fn one() -> Self {
        Qubit(0.0, 1.0)
    }
}

/// Operations supported on qubits (both single-qubit operations and multi-qubit operations).
pub mod ops {
    use super::*;
    use std::f64::consts::*;
    use std::mem::swap;

    /// Quantum NOT operator.
    pub fn not(q: &mut Qubit) {
        swap(&mut q.0, &mut q.1);
    }

    /// Quantum HAD operator (Hadamard gate).
    pub fn had(q: &mut Qubit) {
        let tmp = q.0;
        q.0 = FRAC_1_SQRT_2 * (q.0 + q.1);
        q.1 = FRAC_1_SQRT_2 * (tmp - q.1);
    }

    /// Quantum READ operator. Returns a random result with probability based on the
    /// magnitudes.
    pub fn read(q: &mut Qubit) -> Binary {
        let measurement: f64 = thread_rng().gen();
        read_deterministic(q, measurement)
    }

    fn read_deterministic(q: &mut Qubit, measurement: f64) -> Binary {
        if measurement < q.0.powi(2) {
            q.0 = 1.0;
            q.1 = 0.0;

            Binary::Zero
        } else {
            q.0 = 0.0;
            q.1 = 1.0;

            Binary::One
        }
    }

    /// Quantum WRITE operator. Deterministically sets the values of a qubit.
    pub fn write(target: &mut Qubit, value: Binary) {
        let mut q = match value {
            Binary::Zero => Qubit::zero(),
            Binary::One => Qubit::one(),
        };

        swap(target, &mut q);
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn not() {
            let mut q = Qubit::zero();
            ops::not(&mut q);
            assert_eq!(q, Qubit::one());
        }

        #[test]
        fn had() {
            let mut q = Qubit::zero();
            ops::had(&mut q);
            assert_eq!(q, Qubit(FRAC_1_SQRT_2, FRAC_1_SQRT_2));

            let mut q = Qubit::one();
            ops::had(&mut q);
            assert_eq!(q, Qubit(FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
        }

        #[test]
        fn read_deterministic() {
            let mut q = Qubit::zero();
            ops::had(&mut q);
            let res = ops::read_deterministic(&mut q, 0.49);
            assert_eq!(res, Binary::Zero);

            ops::had(&mut q);
            let res = ops::read_deterministic(&mut q, 0.51);
            assert_eq!(res, Binary::One);
        }

        #[test]
        fn write() {
            let mut q = Qubit::one();
            ops::write(&mut q, Binary::Zero);
            assert_eq!(q, Qubit::zero());
        }
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
