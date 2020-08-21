use nalgebra::{Complex, DMatrix, DVector};
use rand::prelude::*;
use std::convert::TryInto;
use std::f64::consts::*;
use std::mem::swap;

mod helpers {
    use super::*;

    /*
    pub fn half_angle_factors(angle: f64) -> (f64, f64) {
        let half_angle = angle / 2.0;
        (half_angle.cos(), half_angle.sin())
    }
    */

    pub fn build_single_qubit_operator_matrix(
        base: C64Matrix,
        n: u8,
        target: PowerOfTwo,
    ) -> C64Matrix {
        let mut matrices = vec![C64Matrix::identity(2, 2); n.into()];
        matrices[target.index()] = base;

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

pub struct PowerOfTwo(usize);
impl From<usize> for PowerOfTwo {
    fn from(u: usize) -> Self {
        assert!(u.is_power_of_two());

        Self(u)
    }
}

impl PowerOfTwo {
    fn index(&self) -> usize {
        self.0.trailing_zeros().try_into().unwrap()
    }
}

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

    fn apply_single_qubit_operator(&mut self, target: PowerOfTwo, op: C64Matrix) {
        let op = helpers::build_single_qubit_operator_matrix(op, self.n, target);

        let mut amps = C64Vector::zeros(0);
        swap(&mut amps, &mut self.amplitudes);
        self.amplitudes = op * amps;
    }

    /// Quantum NOT operator. Applies to the single qubit specified by `target`.
    pub fn not<T: Into<PowerOfTwo>>(&mut self, target: T) {
        self.apply_single_qubit_operator(
            target.into(),
            C64Matrix::from_row_slice(2, 2, &[ZERO, ONE, ONE, ZERO]),
        );
    }

    /// Quantum HAD operator (Hadamard gate). Applies to the single qubit specified by `target`.
    pub fn had<T: Into<PowerOfTwo>>(&mut self, target: T) {
        self.apply_single_qubit_operator(
            target.into(),
            C64Matrix::from_row_slice(2, 2, &[FRAC2, FRAC2, FRAC2, -FRAC2]),
        );
    }

    /*
    /// Quantum READ operator. Collapses the state of the qubit specified by `target` and returns
    /// the result. Will cause re-normalization.
    pub fn read<T: Into<QubitAddress> + Copy>(&mut self, target: T) -> u8 {
        let mut rng = thread_rng();
        let measurement = rng.gen();

        self.read_deterministic(target, measurement)
    }

    fn read_deterministic<T>(&mut self, target: T, measurement: f64) -> u8
    where
        T: Into<QubitAddress> + Copy,
    {
        let mut zero_probability = 0.0;

        // the sum of norm squares of the left-hand side of operator pairs is equal to the
        // probability that this qubit measures 0.
        helpers::for_each_operator_pair(&mut self.amplitudes, target, |a1, _| {
            zero_probability += a1.norm_sqr();
        });

        let result = if measurement < zero_probability { 0 } else { 1 };

        // we then eliminate incompatible states (e.g. if we measure a 0, we eliminate all states
        // where this qubit would be a 1, by setting that state's amplitude to 0).
        helpers::for_each_operator_pair(&mut self.amplitudes, target, |a1, a2| {
            if result == 0 {
                *a2 = Complex64::new(0.0, 0.0);
            } else {
                *a1 = Complex64::new(0.0, 0.0);
            }
        });

        self.renormalize();

        result
    }

    // XXX check this implementation
    fn renormalize(&mut self) {
        let mut sum = 0.0;

        for a in self.amplitudes.iter() {
            sum += a.norm_sqr();
        }

        for a in self.amplitudes.iter_mut() {
            *a /= sum.sqrt();
        }
    }

    /// Quantum WRITE operator. Sets the value of the qubit specified by
    /// `target`. Equivlent to a read + a conditional NOT. Value must be 0 or 1.
    pub fn write<T: Into<QubitAddress> + Copy>(&mut self, target: T, value: u8) {
        let mut rng = thread_rng();
        let measurement = rng.gen();

        self.write_deterministic(target, value, measurement);
    }

    fn write_deterministic<T>(&mut self, target: T, value: u8, measurement: f64)
    where
        T: Into<QubitAddress> + Copy,
    {
        assert!(value == 0 || value == 1);

        if self.read_deterministic(target, measurement) != value {
            self.not(target);
        }
    }

    /// Quantum PHASE operator. Maps |1> to e^{i \phi} |1>. Applies to the qubits specified by
    /// `target`.
    pub fn phase<T: Into<QubitTarget>>(&mut self, angle: f64, target: T) {
        helpers::bitmask_for_each(target, self.qs.iter_mut(), |q, _| {
            q.1 = (Complex::i() * angle).exp() * q.1;
        });
    }

    /// Quantum ROTX operator. Rotates in the X plane of the Bloch sphere. Applies to the qubits
    /// specified by `target`.
    pub fn rotx<T: Into<QubitTarget>>(&mut self, angle: f64, target: T) {
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
    pub fn roty<T: Into<QubitTarget>>(&mut self, angle: f64, target: T) {
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
    pub fn root_of_not<T: Into<QubitTarget>>(&mut self, target: T) {
        let target = target.into();
        self.had(target);
        self.phase(-FRAC_PI_2, target);
        self.had(target);
    }

    /// Exchange 2 targets. The targets must specify the same number of qubits.
    pub fn exchange<T: Into<QubitTarget>>(&mut self, src: T, dest: T) {
        let src = src.into();
        let dest = dest.into();

        assert_eq!(src.count(), dest.count());

        let mut qubits = Vec::new();
        helpers::bitmask_for_each(src, self.qs.iter_mut(), |q, _| {
            qubits.push(Qubit::zero());
            swap(qubits.last_mut().unwrap(), q);
        });

        let mut curr = 0;
        helpers::bitmask_for_each(dest, self.qs.iter_mut(), |q, _| {
            swap(&mut qubits[curr], q);
            curr += 1;
        });

        let mut curr = 0;
        helpers::bitmask_for_each(src, self.qs.iter_mut(), |q, _| {
            swap(&mut qubits[curr], q);
            curr += 1;
        });
    }
    */
}

#[cfg(test)]
mod tests {
    use super::*;

    /*
    #[test]
    fn qint() {
        let mut qc = QuantumComputer::reset(5);
        assert_eq!(qc.qint(3, ""), QubitTarget::Register { skip: 0, take: 3 });
        assert_eq!(qc.qint(2, ""), QubitTarget::Register { skip: 3, take: 2 });
    }
    */

    #[test]
    fn not() {
        let mut qc = QuantumComputer::reset(1);
        qc.not(0b1);
        assert_eq!(qc.amplitudes, C64Vector::from_column_slice(&[ZERO, ONE]));

        // the amplitudes I'm using here aren't realizable, but they do demonstrate the appropriate
        // swapping of states
        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = C64Vector::from_column_slice(&[ONE, ONE, ZERO, ONE]);
        qc.not(0b1);
        assert_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[ONE, ONE, ONE, ZERO])
        );

        // the amplitudes I'm using here aren't realizable, but they do demonstrate the appropriate
        // swapping of states
        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = C64Vector::from_column_slice(&[ONE, ONE, ZERO, ONE]);
        qc.not(0b10);
        assert_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[ZERO, ONE, ONE, ONE])
        );

        let mut qc = QuantumComputer::reset(3);
        qc.not(0b100);
        assert_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO])
        );
    }

    #[test]
    fn had() {
        let mut qc = QuantumComputer::reset(1);
        qc.had(0b1);
        assert_eq!(qc.amplitudes, C64Vector::from_column_slice(&[FRAC2, FRAC2]));

        qc = QuantumComputer::reset(1);
        qc.not(0b1);
        qc.had(0b1);
        assert_eq!(
            qc.amplitudes,
            C64Vector::from_column_slice(&[FRAC2, -FRAC2])
        );
    }

    /*
    #[test]
    fn read_deterministic() {
        let mut qc = QuantumComputer::reset(1);
        qc.had(0b1);
        assert_eq!(qc.read_deterministic(0b1, 0.49), 0);

        let mut qc = QuantumComputer::reset(1);
        qc.had(0b1);
        assert_eq!(qc.read_deterministic(0b1, 0.51), 1);

        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = vec![frac2, zero, zero, frac2];
        assert_eq!(qc.read_deterministic(0b1, 0.49), 0);
        assert_eq!(qc.amplitudes, vec![one, zero, zero, zero]);

        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = vec![one, one, zero, one];
        qc.renormalize();
        qc.read_deterministic(0b1, 0.1);
        assert_eq!(qc.amplitudes, vec![one, zero, zero, zero]);

        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = vec![one, one, zero, one];
        qc.renormalize();
        qc.read_deterministic(0b1, 0.9);
        assert_eq!(qc.amplitudes, vec![zero, frac2, zero, frac2]);
    }

    #[test]
    fn write() {
        let mut qc = QuantumComputer::reset(1);
        qc.write(0b1, 1);
        assert_eq!(qc.amplitudes, vec![zero, one]);

        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = vec![frac2, zero, zero, frac2];
        qc.write_deterministic(0b1, 1, 0.49);
        assert_eq!(qc.amplitudes, vec![zero, one, zero, zero]);
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
        assert_eq!(
            qc.qs[0],
            Qubit(Complex::new(0.0, 0.0), Complex::new(0.0, 1.0))
        );
    }

    #[test]
    fn root_of_not() {
        let mut qc = QuantumComputer::reset(1);
        qc.root_of_not(0b1);
        qc.root_of_not(0b1);
        assert_eq!(qc.qs[0], Qubit::one());
    }

    #[test]
    fn exchange() {
        let mut qc = QuantumComputer::reset(4);
        let a = qc.qint(2, "");
        let b = qc.qint(2, "");
        qc.write(0b1001, 0b1111);
        qc.exchange(a, b);
        assert_eq!(qc.read(0b1111), 6);

        let mut qc = QuantumComputer::reset(6);
        let a = qc.qint(3, "");
        let b = qc.qint(3, "");
        qc.write(0b101000, 0b111111);
        qc.exchange(a, b);
        assert_eq!(qc.read(0b111111), 5);
    }
    */
}
