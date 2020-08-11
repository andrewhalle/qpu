use num::complex::Complex64;
use rand::prelude::*;
use std::collections::HashSet;
use std::convert::TryInto;
use std::f64::consts::*;
use std::mem::swap;

mod helpers {
    use super::*;
    /*
        use std::slice::IterMut;
    */

    pub fn half_angle_factors(angle: f64) -> (f64, f64) {
        let half_angle = angle / 2.0;
        (half_angle.cos(), half_angle.sin())
    }

    pub fn for_each_operator_pair<
        T: Into<QubitAddress>,
        F: FnMut(&mut Complex64, &mut Complex64),
    >(
        amplitudes: &mut Vec<Complex64>,
        target: T,
        mut f: F,
    ) {
        let mut processed = HashSet::new();
        let target = target.into().0;

        for i in 0..amplitudes.len() {
            if !processed.contains(&i) {
                let pair = i + target;

                let (first, last) = amplitudes.split_at_mut(pair);
                f(&mut first[i], &mut last[0]);
                processed.insert(i);
                processed.insert(pair);
            }
        }
    }

    /*
        pub fn bitmask_for_each<T: Into<QubitTarget>, F: FnMut(&mut Qubit, usize)>(
            target: T,
            iterator: IterMut<Qubit>,
            mut f: F,
        ) {
            let target = target.into();
            let (skip, mask) = match target {
                QubitTarget::Absolute(t) => (0, t),
                QubitTarget::Register { skip, take } => (skip, (1_usize << take) - 1),
            };

            let mut curr = 1_usize;
            for q in iterator.skip(skip) {
                if mask & curr != 0 {
                    f(q, curr);
                }

                curr <<= 1;
            }
        }
    */
}
/*

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum QubitTarget {
    Absolute(usize),
    Register { skip: usize, take: usize },
}

impl From<usize> for QubitTarget {
    fn from(u: usize) -> QubitTarget {
        QubitTarget::Absolute(u)
    }
}

impl QubitTarget {
    fn count(&self) -> usize {
        match self {
            QubitTarget::Absolute(t) => t.count_ones().try_into().unwrap(),
            QubitTarget::Register { take, .. } => *take,
        }
    }
}
*/
pub struct QubitAddress(usize);

impl From<usize> for QubitAddress {
    fn from(u: usize) -> Self {
        if !u.is_power_of_two() {
            panic!("QubitAddress specifies a single qubit.");
        }

        QubitAddress(u)
    }
}

/// A quantum computer prepared with a number of qubits. The quantum computer does not store qubits
/// individually, but rather stores the complex amplitudes associated to each basis state (for the
/// number of qubits the QC is initialized with). There will be 2^n such amplitudes.
pub struct QuantumComputer {
    amplitudes: Vec<Complex64>,
    allocated: usize,
}

impl QuantumComputer {
    /// Initialize the quantum computer with n qubits, in state |00..0>
    pub fn reset(n: u8) -> Self {
        let mut amplitudes = vec![Complex64::new(0.0, 0.0); 2_usize.pow(n.into())];

        amplitudes[0] = Complex64::new(1.0, 0.0);

        QuantumComputer {
            amplitudes,
            allocated: 0,
        }
    }

    /// Start a label block (for the circuit diagram).
    pub fn label(&mut self, _msg: &str) {}

    /// End a label block (for the circuit diagram).
    pub fn end_label(&mut self) {}

    /*
    /// Create a register using the bottom `n` unallocated bits.
    pub fn qint(&mut self, take: usize, _name: &str) -> QubitTarget {
        let skip = self.allocated;
        self.allocated += take;

        QubitTarget::Register { skip, take }
    }
    */

    /// Quantum NOT operator. Applies to the single qubit specified by `target`.
    pub fn not<T: Into<QubitAddress>>(&mut self, target: T) {
        helpers::for_each_operator_pair(&mut self.amplitudes, target, |a1, a2| {
            swap(a1, a2);
        });
    }

    /// Quantum HAD operator (Hadamard gate). Applies to the single qubit specified by `target`.
    pub fn had<T: Into<QubitAddress>>(&mut self, target: T) {
        helpers::for_each_operator_pair(&mut self.amplitudes, target, |a1, a2| {
            let tmp = *a1;
            *a1 = FRAC_1_SQRT_2 * (*a1 + *a2);
            *a2 = FRAC_1_SQRT_2 * (tmp - *a2);
        });
    }

    /*
    /// Quantum READ operator. Returns a random result with probability based on the
    /// magnitudes. Applies to all qubits identified by the mask `target`.
    pub fn read<T: Into<QubitTarget>>(&mut self, target: T) -> usize {
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
    pub fn write<T: Into<QubitTarget>>(&mut self, value: usize, target: T) {
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
        assert_eq!(
            qc.amplitudes,
            vec![Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)]
        );

        // the amplitudes I'm using here aren't realizable, but they do demonstrate the appropriate
        // swapping of states
        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        ];
        qc.not(0b1);
        assert_eq!(
            qc.amplitudes,
            vec![
                Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(0.0, 0.0),
            ],
        );

        // the amplitudes I'm using here aren't realizable, but they do demonstrate the appropriate
        // swapping of states
        let mut qc = QuantumComputer::reset(2);
        qc.amplitudes = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        ];
        qc.not(0b10);
        assert_eq!(
            qc.amplitudes,
            vec![
                Complex64::new(0.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0),
            ],
        );
    }

    #[test]
    fn had() {
        let mut qc = QuantumComputer::reset(1);
        qc.had(0b1);
        assert_eq!(
            qc.amplitudes,
            vec![
                Complex64::new(FRAC_1_SQRT_2, 0.0),
                Complex64::new(FRAC_1_SQRT_2, 0.0),
            ]
        );

        qc = QuantumComputer::reset(1);
        qc.not(0b1);
        qc.had(0b1);
        assert_eq!(
            qc.amplitudes,
            vec![
                Complex64::new(FRAC_1_SQRT_2, 0.0),
                Complex64::new(-FRAC_1_SQRT_2, 0.0),
            ]
        );
    }

    /*
    #[test]
    fn read() {
        let mut qc = QuantumComputer::reset(4);
        qc.not(0b1010);
        assert_eq!(qc.read(0b1111), 10);

        let mut qc = QuantumComputer::reset(4);
        qc.not(0b1010);
        assert_eq!(qc.read(QubitTarget::Register { skip: 2, take: 2 }), 2);
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

        let mut qc = QuantumComputer::reset(4);
        qc.write(0b10, QubitTarget::Register { skip: 2, take: 2 });
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
