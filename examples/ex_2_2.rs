// Ex 2-2: Random Byte
use qpu::QuantumComputer;

fn main() {
    let mut qc = QuantumComputer::reset(8);
    qc.had(0xff);
    println!("{:#010b}", qc.read(0xff));
}
