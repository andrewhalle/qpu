// Ex 2-4: Quantum Spy Hunter
use qpu::{Binary, QuantumComputer};

fn get_random_bit(q: &mut QInt) -> Binary {
    q.write(0);
    q.had();

    q.read()
}

// Runs the spy hunter protocol. Returns true if a spy is caught.
fn run_spy_hunter_protocol() -> bool {
    let mut qc = QuantumComputer::reset(3);

    let mut a = qc.qint(1, "alice");
    let mut fiber = qc.qint(1, "fiber");
    let mut b = qc.qint(1, "bob");

    qc.label("Generate two random bits.");
    let send_had = get_random_bit(&mut a);
    let send_value = get_random_bit(&mut a);

    let spy_is_present = get_random_bit(&mut fiber);

    let recv_had = get_random_bit(&mut b);

    // Prepare Alice's qubit.
    a.write(0);
    qc.label("set value");
    if send_value.is_one() {
        a.not();
    }
    qc.end_label();

    qc.label("apply had");
    if send_had.is_one() {
        a.had();
    }
    qc.end_label();

    // Send the qubit!
    fiber.exchange(a);

    // Activate the spy
    if spy_is_present.is_one() {
        qc.label("spy");
        fiber.had();

        let stolen_data = fiber.read();

        fiber.write(0);
        if stolen_data.is_one() {
            fiber.not();
        }
        fiber.had();
        qc.end_label();
    }

    // Receive the qubit!
    fiber.exchange(b);

    qc.label("apply had");
    if recv_had.is_one() {
        b.had();
    }
    qc.end_label();

    qc.label("read value");
    let recv_val = b.read();
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

fn main() {}
