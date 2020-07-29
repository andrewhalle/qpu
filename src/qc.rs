use super::qubit::Qubit;
use super::*;
use num::FromPrimitive;
use rand::prelude::*;
use std::f64::consts::*;
use std::mem::swap;

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
        if measurement < self.q.0.norm_sqr() {
            self.q.0 = FromPrimitive::from_f64(1.0).unwrap();
            self.q.1 = FromPrimitive::from_f64(0.0).unwrap();

            Binary::Zero
        } else {
            self.q.0 = FromPrimitive::from_f64(0.0).unwrap();
            self.q.1 = FromPrimitive::from_f64(1.0).unwrap();

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
    use num::FromPrimitive;

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
        assert_eq!(
            qc.q,
            Qubit(
                FromPrimitive::from_f64(FRAC_1_SQRT_2).unwrap(),
                FromPrimitive::from_f64(FRAC_1_SQRT_2).unwrap()
            )
        );

        qc = QuantumComputer::reset(1);
        qc.write(Binary::One);
        qc.had();
        assert_eq!(
            qc.q,
            Qubit(
                FromPrimitive::from_f64(FRAC_1_SQRT_2).unwrap(),
                FromPrimitive::from_f64(-FRAC_1_SQRT_2).unwrap()
            )
        );
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
