use rand::prelude::*;
use std::f64::consts::*;
use std::mem::swap;

#[derive(Debug, PartialEq)]
pub enum Binary {
    Zero,
    One,
}

/// Representation of a qubit.
#[derive(Debug, PartialEq)]
struct Qubit(f64, f64);

impl Qubit {
    fn zero() -> Self {
        Qubit(1.0, 0.0)
    }

    fn one() -> Self {
        Qubit(0.0, 1.0)
    }
}

/// An quantum computer prepared with a single qubit.
pub struct QuantumComputer {
    q: Qubit,
}

impl QuantumComputer {
    /// Initialize the quantum computer with a zero qubit.
    pub fn reset(_n: u8) -> Self {
        QuantumComputer { q: Qubit::zero() }
    }

    /// Quantum NOT operator.
    pub fn not(&mut self) {
        swap(&mut self.q.0, &mut self.q.1);
    }

    /// Quantum HAD operator (Hadamard gate).
    pub fn had(&mut self) {
        let tmp = self.q.0;
        self.q.0 = FRAC_1_SQRT_2 * (self.q.0 + self.q.1);
        self.q.1 = FRAC_1_SQRT_2 * (tmp - self.q.1);
    }

    /// Quantum READ operator. Returns a random result with probability based on the
    /// magnitudes.
    pub fn read(&mut self) -> Binary {
        let measurement: f64 = thread_rng().gen();
        self.read_deterministic(measurement)
    }

    fn read_deterministic(&mut self, measurement: f64) -> Binary {
        if measurement < self.q.0.powi(2) {
            self.q.0 = 1.0;
            self.q.1 = 0.0;

            Binary::Zero
        } else {
            self.q.0 = 0.0;
            self.q.1 = 1.0;

            Binary::One
        }
    }

    /// Quantum WRITE operator. Deterministically sets the values of a qubit.
    pub fn write(&mut self, value: Binary) {
        let mut q = match value {
            Binary::Zero => Qubit::zero(),
            Binary::One => Qubit::one(),
        };

        swap(&mut self.q, &mut q);
    }
}

#[cfg(test)]
mod qc_tests {
    use super::*;

    #[test]
    fn not() {
        let mut qc = QuantumComputer::reset(1);
        qc.not();
        assert_eq!(qc.q, Qubit::one());
    }

    #[test]
    fn had() {
        let mut qc = QuantumComputer::reset(1);
        qc.had();
        assert_eq!(qc.q, Qubit(FRAC_1_SQRT_2, FRAC_1_SQRT_2));

        qc = QuantumComputer::reset(1);
        qc.write(Binary::One);
        qc.had();
        assert_eq!(qc.q, Qubit(FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
    }

    #[test]
    fn read_deterministic() {
        let mut qc = QuantumComputer::reset(1);
        qc.had();
        let res = qc.read_deterministic(0.49);
        assert_eq!(res, Binary::Zero);

        qc.had();
        let res = qc.read_deterministic(0.51);
        assert_eq!(res, Binary::One);
    }

    #[test]
    fn write() {
        let mut qc = QuantumComputer::reset(0);
        qc.write(Binary::One);
        assert_eq!(qc.q, Qubit::one());
    }
}

#[cfg(test)]
mod qubit_tests {
    use super::*;

    #[test]
    fn init_qubits() {
        assert_eq!(Qubit::zero().0, 1.0);
        assert_eq!(Qubit::zero().1, 0.0);
        assert_eq!(Qubit::one().0, 0.0);
        assert_eq!(Qubit::one().1, 1.0);
    }
}
