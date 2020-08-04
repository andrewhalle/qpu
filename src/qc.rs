use super::qubit::{Qubit, Binary};
use rand::prelude::*;
use std::f64::consts::*;
use std::mem::swap;
use num::Complex;

mod helpers {
    use super::*;
    use std::slice::IterMut;

    pub fn half_angle_factors(angle: f64) -> (f64, f64) {
        let half_angle = angle / 2.0;
        (half_angle.cos(), half_angle.sin())
    }

    pub fn bitmask_for_each<F: FnMut(&mut Qubit, usize)>(mask: usize, iterator: IterMut<Qubit>, mut f: F) {
        let mut curr = 1_usize;
        for q in iterator {
            if mask & curr != 0 {
                f(q, curr);
            }

            curr <<= 1;
        }
    }
}

/// An quantum computer prepared with a single qubit.
pub struct QuantumComputer {
    qs: Vec<Qubit>,
}

impl QuantumComputer {
    /// Initialize the quantum computer with n qubits, all zero.
    pub fn reset(n: u8) -> Self {
        QuantumComputer {
            qs: vec![Qubit::zero(); n.into()],
        }
    }

    /// Start a label block (for the circuit diagram).
    pub fn label(&mut self, _msg: &str) {}

    /// End a label block (for the circuit diagram).
    pub fn end_label(&mut self) {}

    /// Quantum NOT operator. Applies to all qubits identified by the mask `target`.
    pub fn not(&mut self, target: usize) {
        helpers::bitmask_for_each(target, self.qs.iter_mut(), |q, _| {
            swap(&mut q.0, &mut q.1);
        });
    }

    /// Quantum HAD operator (Hadamard gate). Applies to all qubits identified by the mask
    /// `target`.
    pub fn had(&mut self, target: usize) {
        helpers::bitmask_for_each(target, self.qs.iter_mut(), |q, _| {
            let tmp = q.0;
            q.0 = FRAC_1_SQRT_2 * (q.0 + q.1);
            q.1 = FRAC_1_SQRT_2 * (tmp - q.1);
        });
    }

    /// Quantum READ operator. Returns a random result with probability based on the
    /// magnitudes. Applies to all qubits identified by the mask `target`.
    pub fn read(&mut self, target: usize) -> usize {
        let mut rng = thread_rng();
        let mut res = 0;

        helpers::bitmask_for_each(target, self.qs.iter_mut(), |q, curr| {
            let measurement = rng.gen();
            if let Binary::One = q.read(measurement) {
                res |= curr;
            }
        });

        res
    }

    /// Quantum WRITE operator. Deterministically sets the values of a qubit based on `value`.
    /// Applies to the qubits specified by `target`.
    pub fn write(&mut self, value: usize, target: usize) {
        helpers::bitmask_for_each(target, self.qs.iter_mut(), |q, curr| {
            let mut new = if value & curr != 0 {
                Qubit::one()
            } else {
                Qubit::zero()
            };

            swap(q, &mut new);
        });
    }

    /// Quantum PHASE operator. Maps |1> to e^{i \phi} |1>. Applies to the qubits specified by
    /// `target`.
    pub fn phase(&mut self, angle: f64, target: usize) {
        helpers::bitmask_for_each(target, self.qs.iter_mut(), |q, _| {
            q.1 = (Complex::i() * angle).exp() * q.1;
        });
    }

    /// Quantum ROTX operator. Rotates in the X plane of the Bloch sphere. Applies to the qubits
    /// specified by `target`.
    pub fn rotx(&mut self, angle: f64, target: usize) {
        let (f1, f2) = helpers::half_angle_factors(angle);

        helpers::bitmask_for_each(target, self.qs.iter_mut(), |q, _| {
            let tmp = q.0;
            q.0 = (f1 * q.0) + (-Complex::i() * f2 * q.1);
            q.1 = (-Complex::i() * f2 * tmp) + (f1 * q.1);

            q.0 *= Complex::i();
            q.1 *= Complex::i();
        });
    }

    /// Quantum ROTY operator. Rotates in the Y plane of the Bloch sphere. Applies to the qubits
    /// specified by `target`.
    pub fn roty(&mut self, angle: f64, target: usize) {
        let (f1, f2) = helpers::half_angle_factors(angle);

        helpers::bitmask_for_each(target, self.qs.iter_mut(), |q, _| {
            let tmp = q.0;
            q.0 = (f1 * q.0) + (-f2 * q.1);
            q.1 = (f2 * tmp) + (f1 * q.1);

            q.0 *= Complex::i();
            q.1 *= Complex::i();
        });
    }

    /// Quantum ROOT-of-NOT operator. Two applications should equal a NOT. Applies to all qubits
    /// specified by `target`.
    pub fn root_of_not(&mut self, target: usize) {
        self.had(target);
        self.phase(-FRAC_PI_2, target);
        self.had(target);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn not() {
        let mut qc = QuantumComputer::reset(1);
        qc.not(0b1);
        assert_eq!(qc.qs[0], Qubit::one());
    }

    #[test]
    fn had() {
        let mut qc = QuantumComputer::reset(1);
        qc.had(0b1);
        assert_eq!(
            qc.qs[0],
            Qubit(
                Complex::new(FRAC_1_SQRT_2, 0.0),
                Complex::new(FRAC_1_SQRT_2, 0.0),
            )
        );

        qc = QuantumComputer::reset(1);
        qc.not(0b1);
        qc.had(0b1);
        assert_eq!(
            qc.qs[0],
            Qubit(
                Complex::new(FRAC_1_SQRT_2, 0.0),
                Complex::new(-FRAC_1_SQRT_2, 0.0),
            )
        );
    }

    #[test]
    fn read() {
        let mut qc = QuantumComputer::reset(4);
        qc.not(0b1010);
        assert_eq!(qc.read(0b1111), 10);
    }

    #[test]
    fn write() {
        let mut qc = QuantumComputer::reset(1);
        qc.write(0x1, 0x1);
        assert_eq!(qc.qs[0], Qubit::one());

        let mut qc = QuantumComputer::reset(4);
        qc.write(0b1010, 0b1111);
        assert_eq!(qc.read(0b1111), 10);

        let mut qc = QuantumComputer::reset(4);
        qc.write(0b1010, 0b11);
        assert_eq!(qc.read(0b1111), 2);

        let mut qc = QuantumComputer::reset(4);
        qc.write(0b1010, 0b1100);
        assert_eq!(qc.read(0b1111), 8);
    }

    #[test]
    fn phase() {
        let mut qc = QuantumComputer::reset(1);
        qc.had(0b1);
        qc.phase(FRAC_PI_4, 0b1);
        assert_eq!(qc.qs[0].read(0.49), Binary::Zero);

        let mut qc = QuantumComputer::reset(1);
        qc.had(0b1);
        qc.phase(FRAC_PI_4, 0b1);
        assert_eq!(qc.qs[0].read(0.51), Binary::One);

        let mut qc = QuantumComputer::reset(1);
        qc.had(0b1);
        qc.phase(FRAC_PI_4, 0b1);

        assert_eq!(
            qc.qs[0],
            Qubit(
                Complex::new(FRAC_1_SQRT_2, 0.0),
                (Complex::i() * FRAC_PI_4).exp() * Complex::new(FRAC_1_SQRT_2, 0.0)
            )
        );
    }

    #[test]
    fn rotx() {
        let mut qc = QuantumComputer::reset(1);
        qc.rotx(PI, 0b1);
        assert_eq!(qc.qs[0], Qubit::one());
    }

    #[test]
    fn roty() {
        let mut qc = QuantumComputer::reset(1);
        qc.roty(PI, 0b1);
        assert_eq!(qc.qs[0], Qubit(Complex::new(0.0, 0.0), Complex::new(0.0, 1.0)));
    }

    #[test]
    fn root_of_not() {
        let mut qc = QuantumComputer::reset(1);
        qc.root_of_not(0b1);
        qc.root_of_not(0b1);
        assert_eq!(qc.qs[0], Qubit::one());
    }
}
