use super::qubit::Qubit;
use super::*;
use num::{Complex, FromPrimitive};
use rand::prelude::*;
use std::f64::consts::*;
use std::mem::swap;

mod helpers {
    pub fn half_angle_factors(angle: f64) -> (f64, f64) {
        let half_angle = angle / 2.0;
        (half_angle.cos(), half_angle.sin())
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

    /// Quantum NOT operator. Applies to all qubits in the quantum computer.
    pub fn not(&mut self) {
        for q in self.qs.iter_mut() {
            swap(&mut q.0, &mut q.1);
        }
    }

    /// Quantum HAD operator (Hadamard gate). Applies to all qubits in the quantum computer.
    pub fn had(&mut self) {
        for q in self.qs.iter_mut() {
            let tmp = q.0;
            q.0 = FRAC_1_SQRT_2 * (q.0 + q.1);
            q.1 = FRAC_1_SQRT_2 * (tmp - q.1);
        }
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

    /// Quantum PHASE operator. Maps |1> to e^{i \phi} |1>.
    pub fn phase(&mut self, angle: f64) {
        self.q.1 = (Complex::i() * angle).exp() * self.q.1;
    }

    /// Quantum ROTX operator. Rotates in the X plane of the Bloch sphere.
    pub fn rotx(&mut self, angle: f64) {
        let (f1, f2) = helpers::half_angle_factors(angle);

        let tmp = self.q.0;
        self.q.0 = (f1 * self.q.0) + (-Complex::i() * f2 * self.q.1);
        self.q.1 = (-Complex::i() * f2 * tmp) + (f1 * self.q.1);

        self.q.0 *= Complex::i();
        self.q.1 *= Complex::i();
    }

    /// Quantum ROTY operator. Rotates in the Y plane of the Bloch sphere.
    pub fn roty(&mut self, angle: f64) {
        let (f1, f2) = helpers::half_angle_factors(angle);

        let tmp = self.q.0;
        self.q.0 = (f1 * self.q.0) + (-f2 * self.q.1);
        self.q.1 = (f2 * tmp) + (f1 * self.q.1);

        self.q.0 *= Complex::i();
        self.q.1 *= Complex::i();
    }

    /// Quantum ROOT-of-NOT operator. Two applications should equal a NOT.
    pub fn root_of_not(&mut self) {
        self.had();
        self.phase(-FRAC_PI_2);
        self.had();
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
        let mut qc = QuantumComputer::reset(1);
        qc.write(Binary::One);
        assert_eq!(qc.q, Qubit::one());
    }

    #[test]
    fn phase() {
        let mut qc = QuantumComputer::reset(1);
        qc.had();
        qc.phase(FRAC_PI_4);
        assert_eq!(qc.read_deterministic(0.49), Binary::Zero);

        let mut qc = QuantumComputer::reset(1);
        qc.had();
        qc.phase(FRAC_PI_4);
        assert_eq!(qc.read_deterministic(0.51), Binary::One);

        let mut qc = QuantumComputer::reset(1);
        qc.had();
        qc.phase(FRAC_PI_4);

        assert_eq!(
            qc.q,
            Qubit(
                FromPrimitive::from_f64(FRAC_1_SQRT_2).unwrap(),
                (Complex::i() * FRAC_PI_4).exp()
                    * <Complex<f64> as FromPrimitive>::from_f64(FRAC_1_SQRT_2).unwrap()
            )
        );
    }

    #[test]
    fn rotx() {
        let mut qc = QuantumComputer::reset(1);
        qc.rotx(PI);
        assert_eq!(qc.q, Qubit::one());
    }

    #[test]
    fn roty() {
        let mut qc = QuantumComputer::reset(1);
        qc.roty(PI);
        assert_eq!(qc.q, Qubit(Complex::new(0.0, 0.0), Complex::new(0.0, 1.0)));
    }

    #[test]
    fn root_of_not() {
        let mut qc = QuantumComputer::reset(1);
        qc.root_of_not();
        qc.root_of_not();
        assert_eq!(qc.q, Qubit::one());
    }
}
