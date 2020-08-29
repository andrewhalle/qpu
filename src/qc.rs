use approx::relative_ne;
use nalgebra::{Complex, DMatrix, DVector};
use rand::prelude::*;
use std::convert::TryInto;
use std::f64::consts::*;
use std::mem::swap;

mod helpers {
    use super::*;

    pub fn half_angle_factors(angle: f64) -> (C64, C64) {
        let half_angle = angle / 2.0;
        (
            C64::new(half_angle.cos(), 0.0),
            C64::new(half_angle.sin(), 0.0),
        )
    }

    pub fn build_operator_matrix(base: C64Matrix, n: u8, target: u8) -> C64Matrix {
        assert!(base.is_square());
        let n = Into::<usize>::into(n);
        let target = Into::<usize>::into(target);

        // the size of base determines the number of qubits being acted on
        let num_qubits: usize = base.nrows().trailing_zeros().try_into().unwrap();

        // set up list of matrices [I, I, .., I] we need n - num_qubits + 1
        // because our final matrix needs to be of size (2^n, 2^n)
        let mut matrices = vec![C64Matrix::identity(2, 2); n - num_qubits + 1];
        matrices[target] = base;

        let mut retval = matrices.pop().unwrap();
        while let Some(m) = matrices.pop() {
            retval = retval.kronecker(&m);
        }

        retval
    }
}

type C64 = Complex<f64>;
type C64Vector = DVector<C64>;
type C64Matrix = DMatrix<C64>;

const ZERO: C64 = C64::new(0.0, 0.0);
const ONE: C64 = C64::new(1.0, 0.0);
const FRAC2: C64 = C64::new(FRAC_1_SQRT_2, 0.0);

/// A quantum computer prepared with a number of qubits. The quantum computer does not store qubits
/// individually, but rather stores the complex amplitudes associated to each basis state (for the
/// number of qubits the QC is initialized with). There will be 2^n such amplitudes.
pub struct QuantumComputer {
    amplitudes: C64Vector,
    n: u8,
}

impl QuantumComputer {
    /// Initialize the quantum computer with n qubits, in state |00..0>
    pub fn reset(n: u8) -> Self {
        assert!(n <= 16);

        let mut amplitudes = C64Vector::from_element(2_usize.pow(n.into()), ZERO);

        amplitudes[0] = ONE;

        QuantumComputer { amplitudes, n }
    }

    /// Overall phase doesn't matter for quantum computing, so after applying an operator
    /// matrix, the amplitudes can be transformed into a standard form by multiplying by
    /// any complex number of unit magnitude.
    ///
    /// Here we choose to make the first non-zero element of the amplitudes vector
    /// real and positive.
    fn remove_phase(&mut self) {
        let first_nonzero_elem = self
            .amplitudes
            .iter()
            .find(|elem| relative_ne!(**elem, ZERO))
            .unwrap();

        let r = 1.0;
        let theta = -first_nonzero_elem.arg();
        let cancel_factor = C64::from_polar(&r, &theta);

        self.amplitudes *= cancel_factor;
    }

    fn apply_operator(&mut self, target: u8, op: C64Matrix) {
        assert!(target < self.n);

        let op = helpers::build_operator_matrix(op, self.n, target);

        let mut amps = C64Vector::zeros(0);
        swap(&mut amps, &mut self.amplitudes);
        self.amplitudes = op * amps;

        self.remove_phase();
    }

    /// Quantum NOT operator. Applies to the single qubit specified by `target`.
    pub fn not(&mut self, target: u8) {
        self.apply_operator(
            target,
            C64Matrix::from_row_slice(2, 2, &[ZERO, ONE, ONE, ZERO]),
        );
    }

    /// Quantum HAD operator (Hadamard gate). Applies to the single qubit specified by `target`.
    pub fn had(&mut self, target: u8) {
        self.apply_operator(
            target,
            C64Matrix::from_row_slice(2, 2, &[FRAC2, FRAC2, FRAC2, -FRAC2]),
        );
    }

    /// Quantum READ operator. Collapses the state of the qubit specified by `target` and returns
    /// the result. Will cause re-normalization.
    pub fn read(&mut self, target: u8) -> u8 {
        let mut rng = thread_rng();
        let measurement = rng.gen();

        self.read_deterministic(target, measurement)
    }

    fn read_deterministic(&mut self, target: u8, measurement: f64) -> u8 {
        // we add the norm squares of the amplitudes corresponding to states where the target
        // qubit is 0.
        let zero_probability = self
            .amplitudes
            .iter()
            .enumerate()
            .filter(|(i, _)| (i & (1 << target)) >> target == 0)
            .map(|(_, val)| val.norm_sqr())
            .sum();

        let result = if measurement < zero_probability { 0 } else { 1 };

        // eliminate states that don't agree with the measurement
        for (_, amp) in self.amplitudes.iter_mut().enumerate().filter(|(i, _)| {
            let bit = (i & (1 << target)) >> target;

            result != bit
        }) {
            *amp = ZERO;
        }

        self.renormalize();

        result.try_into().unwrap()
    }

    fn renormalize(&mut self) {
        let sum: f64 = self.amplitudes.iter().map(|v| v.norm_sqr()).sum();

        for a in self.amplitudes.iter_mut() {
            *a /= sum.sqrt();
        }
    }

    /// Quantum WRITE operator. Sets the value of the qubit specified by
    /// `target`. Equivlent to a read + a conditional NOT. Value must be 0 or 1.
    pub fn write(&mut self, target: u8, value: u8) {
        let mut rng = thread_rng();
        let measurement = rng.gen();

        self.write_deterministic(target, value, measurement);
    }

    fn write_deterministic(&mut self, target: u8, value: u8, measurement: f64) {
        assert!(value == 0 || value == 1);

        if self.read_deterministic(target, measurement) != value {
            self.not(target);
        }
    }

    /// Quantum PHASE operator. Maps |1> to e^{i \phi} |1>. Applies to `target`.
    pub fn phase(&mut self, angle: f64, target: u8) {
        self.apply_operator(
            target,
            C64Matrix::from_row_slice(2, 2, &[ONE, ZERO, ZERO, (C64::i() * angle).exp()]),
        );
    }

    /// Quantum ROTX operator. Rotates in the X plane of the Bloch sphere. Applies to `target`.
    pub fn rotx(&mut self, angle: f64, target: u8) {
        let (cos, sin) = helpers::half_angle_factors(angle);

        self.apply_operator(
            target,
            C64Matrix::from_row_slice(2, 2, &[cos, sin * -C64::i(), sin * -C64::i(), cos]),
        );
    }

    /// Quantum ROTY operator. Rotates in the Y plane of the Bloch sphere. Applies to `target`.
    pub fn roty(&mut self, angle: f64, target: u8) {
        let (cos, sin) = helpers::half_angle_factors(angle);

        self.apply_operator(
            target,
            C64Matrix::from_row_slice(2, 2, &[cos, -sin, sin, cos]),
        );
    }

    /// Quantum ROOT-of-NOT operator. Two applications should equal a NOT. Applies to all qubits
    /// specified by `target`.
    pub fn root_of_not(&mut self, target: u8) {
        self.had(target);
        self.phase(-FRAC_PI_2, target);
        self.had(target);
    }

    /// Swap 2 qubits. Qubits must be neighbors (only local interactions allowed).
    /// Applies to qubits `target` and `target + 1`.
    pub fn swap(&mut self, target: u8) {
        #[rustfmt::skip]
        self.apply_operator(
            target,
            C64Matrix::from_row_slice(4, 4, &[ ONE,  ZERO, ZERO, ZERO,
                                               ZERO, ZERO, ONE,  ZERO,
                                               ZERO, ONE,  ZERO, ZERO,
                                               ZERO, ZERO, ZERO, ONE  ]),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn remove_phase() {
        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes =
            C64Vector::from_column_slice(&[ZERO, C64::i() / 2.0, -ONE / 2.0, -C64::i() / 2.0]);
        qc.remove_phase();
        assert_relative_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[ZERO, ONE / 2.0, C64::i() / 2.0, -ONE / 2.0])
        );
    }

    #[test]
    fn not() {
        let mut qc = QuantumComputer::reset(1);
        qc.not(0);
        assert_relative_eq!(qc.amplitudes, C64Vector::from_column_slice(&[ZERO, ONE]));

        // the amplitudes I'm using here aren't realizable, but they do demonstrate the appropriate
        // swapping of states
        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = C64Vector::from_column_slice(&[ONE, ONE, ZERO, ONE]);
        qc.not(0);
        assert_relative_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[ONE, ONE, ONE, ZERO])
        );

        // the amplitudes I'm using here aren't realizable, but they do demonstrate the appropriate
        // swapping of states
        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = C64Vector::from_column_slice(&[ONE, ONE, ZERO, ONE]);
        qc.not(1);
        assert_relative_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[ZERO, ONE, ONE, ONE])
        );

        let mut qc = QuantumComputer::reset(3);
        qc.not(2);
        assert_relative_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO])
        );
    }

    #[test]
    fn had() {
        let mut qc = QuantumComputer::reset(1);
        qc.had(0);
        assert_relative_eq!(qc.amplitudes, C64Vector::from_column_slice(&[FRAC2, FRAC2]));

        qc = QuantumComputer::reset(1);
        qc.not(0);
        qc.had(0);
        assert_relative_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[FRAC2, -FRAC2])
        );
    }

    #[test]
    fn read_deterministic() {
        let mut qc = QuantumComputer::reset(1);
        qc.had(0);
        assert_eq!(qc.read_deterministic(0, 0.49), 0);

        let mut qc = QuantumComputer::reset(1);
        qc.had(0);
        assert_eq!(qc.read_deterministic(0, 0.51), 1);

        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = C64Vector::from_column_slice(&[FRAC2, ZERO, ZERO, FRAC2]);
        assert_eq!(qc.read_deterministic(0, 0.49), 0);
        assert_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[ONE, ZERO, ZERO, ZERO])
        );

        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = C64Vector::from_column_slice(&[ONE, ONE, ZERO, ONE]);
        qc.renormalize();
        qc.read_deterministic(0, 0.1);
        assert_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[ONE, ZERO, ZERO, ZERO])
        );

        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = C64Vector::from_column_slice(&[ONE, ONE, ZERO, ONE]);
        qc.renormalize();
        qc.read_deterministic(0, 0.9);
        assert_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[ZERO, FRAC2, ZERO, FRAC2])
        );
    }

    #[test]
    fn write_deterministic() {
        let mut qc = QuantumComputer::reset(1);
        qc.write_deterministic(0, 1, 0.0);
        assert_eq!(qc.amplitudes, C64Vector::from_column_slice(&[ZERO, ONE]));

        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = C64Vector::from_column_slice(&[FRAC2, ZERO, ZERO, FRAC2]);
        qc.write_deterministic(0, 1, 0.49);
        assert_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[ZERO, ONE, ZERO, ZERO])
        );
    }

    #[test]
    fn phase() {
        let mut qc = QuantumComputer::reset(1);
        qc.phase(FRAC_PI_4, 0);
        assert_relative_eq!(qc.amplitudes, C64Vector::from_column_slice(&[ONE, ZERO]));

        let e_pi_4 = (C64::i() * FRAC_PI_4).exp();
        let mut qc = QuantumComputer::reset(1);
        qc.had(0);
        qc.phase(FRAC_PI_4, 0);
        assert_relative_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[FRAC2, FRAC2 * e_pi_4])
        );
    }

    #[test]
    fn rotx() {
        let mut qc = QuantumComputer::reset(1);
        qc.rotx(PI, 0);
        assert_relative_eq!(qc.amplitudes, C64Vector::from_column_slice(&[ZERO, ONE]));
    }

    #[test]
    fn roty() {
        let mut qc = QuantumComputer::reset(1);
        qc.roty(PI, 0);
        assert_relative_eq!(qc.amplitudes, C64Vector::from_column_slice(&[ZERO, ONE]),);
    }

    #[test]
    fn root_of_not() {
        let mut qc = QuantumComputer::reset(1);
        qc.root_of_not(0);
        qc.root_of_not(0);
        assert_relative_eq!(qc.amplitudes, C64Vector::from_column_slice(&[ZERO, ONE]));
    }

    #[test]
    fn swap() {
        let mut qc = QuantumComputer::reset(4);
        qc.not(0);
        qc.swap(0);
        let mut expected = C64Vector::zeros(16);
        expected[2] = ONE;
        assert_relative_eq!(qc.amplitudes, expected);
    }
}
