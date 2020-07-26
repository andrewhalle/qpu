#[derive(Debug, PartialEq)]
pub enum Binary {
    Zero,
    One,
}

mod qc;
mod qubit;

pub use qc::QuantumComputer;
