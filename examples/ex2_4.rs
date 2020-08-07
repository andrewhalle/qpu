// Ex 2-4: Quantum Spy Hunter
use qpu::{QuantumComputer, QubitTarget};

fn get_random_bit<T: Into<QubitTarget>>(qc: &mut QuantumComputer, target: T) -> usize {
    let target = target.into();

    qc.write(0, target);
    qc.had(target);

    qc.read(target)
}

// Runs the spy hunter protocol. Returns true if a spy is caught.
fn run_spy_hunter_protocol() -> bool {
    let mut qc = QuantumComputer::reset(3);

    let a = qc.qint(1, "alice");
    let fiber = qc.qint(1, "fiber");
    let b = qc.qint(1, "bob");

    qc.label("Generate two random bits.");
    let send_had = get_random_bit(&mut qc, a);
    let send_value = get_random_bit(&mut qc, a);

    let spy_is_present = get_random_bit(&mut qc, fiber);

    let recv_had = get_random_bit(&mut qc, b);

    // Prepare Alice's qubit.
    qc.write(0, a);
    qc.label("set value");
    if send_value == 1 {
        qc.not(a);
    }
    qc.end_label();

    qc.label("apply had");
    if send_had == 1 {
        qc.had(a);
    }
    qc.end_label();

    // Send the qubit!
    qc.exchange(a, fiber);

    // Activate the spy
    if spy_is_present == 1 {
        qc.label("spy");
        qc.had(fiber);

        let stolen_data = qc.read(fiber);

        qc.write(0, fiber);
        if stolen_data == 1 {
            qc.not(fiber);
        }
        qc.had(fiber);
        qc.end_label();
    }

    // Receive the qubit!
    qc.exchange(fiber, b);

    qc.label("apply had");
    if recv_had == 1 {
        qc.had(b);
    }
    qc.end_label();

    qc.label("read value");
    let recv_val = qc.read(b);
    qc.end_label();

    // Now Alice emails Bob to tell
    // him her had setting and value.
    // If the had setting matches and the
    // value does not, there's a spy!
    if send_had == recv_had {
        if send_value != recv_val {
            return true;
        }
    }

    false
}

fn main() {
    let mut caught = 0;
    let mut total = 0;
    while total < 1000 {
        if run_spy_hunter_protocol() {
            caught += 1;
        }

        total += 1;
    }

    println!("Caught {} spies in {} trials.", caught, total);
}
