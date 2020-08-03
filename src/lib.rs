use std::fmt::{Display, Formatter};

#[derive(Debug, PartialEq)]
pub enum Binary {
    Zero,
    One,
}

impl Display for Binary {
    fn fmt(&self, f: &mut Formatter) -> Result<(), std::fmt::Error> {
        match self {
            Binary::Zero => write!(f, "0"),
            Binary::One => write!(f, "1"),
        }
    }
}

impl Binary {
    pub fn is_one(&self) -> bool {
        match self {
            Binary::Zero => false,
            Binary::One => true,
        }
    }
}

mod qc;
mod qint;
mod qubit;

pub use qc::QuantumComputer;
