// Ex 2-4: Quantum Spy Hunter
use qpu::QuantumComputer;
use rand::prelude::*;

enum SimulationResult {
    UnusableBit,
    UsableBit,
    SpyDetected,
}

// Runs the spy hunter protocol. Returns true if a spy is caught.
fn send_one_qubit() -> SimulationResult {
    let mut rng = thread_rng();

    // -- simulation parameters
    let alice_not = rng.gen_bool(0.5);
    let alice_val = if alice_not { 1 } else { 0 };
    let alice_had = rng.gen_bool(0.5);

    // being a spy-hunter qubit means we send the value over the authenticated channel
    // to try to detect a spy.
    let is_spy_hunter = rng.gen_bool(0.01);
    let spy_is_present = rng.gen_bool(0.5);
    let spy_had = rng.gen_bool(0.5);

    let bob_had = rng.gen_bool(0.5);

    let mut qc = QuantumComputer::reset(3);

    let a = qc.qint(1, "alice");
    let fiber = qc.qint(1, "fiber");
    let b = qc.qint(1, "bob");

    // Prepare Alice's qubit.
    qc.write(0, a);
    qc.label("set value");
    if alice_not {
        qc.not(a);
    }
    qc.end_label();

    qc.label("apply had");
    if alice_had {
        qc.had(a);
    }
    qc.end_label();

    // Send the qubit!
    qc.exchange(a, fiber);

    // Activate the spy
    if spy_is_present {
        qc.label("spy");

        let stolen_data = qc.read(fiber);

        qc.write(0, fiber);
        if stolen_data == 1 {
            qc.not(fiber);
        }

        if spy_had {
            qc.had(fiber);
        }
        qc.end_label();
    }

    // Receive the qubit!
    qc.exchange(fiber, b);

    qc.label("apply had");
    if bob_had {
        qc.had(b);
    }
    qc.end_label();

    qc.label("read value");
    let bob_val = qc.read(b);
    qc.end_label();

    // basis setting is always shared
    if alice_had == bob_had {
        // only check values if this is a spy-hunter qubit, otherwise we don't share
        // the value
        if is_spy_hunter && alice_val != bob_val {
            return SimulationResult::SpyDetected;
        } else if is_spy_hunter {
            // spy hunter qubits aren't usable because we share the value over a public channel.
            return SimulationResult::UnusableBit;
        } else {
            return SimulationResult::UsableBit;
        }
    }

    SimulationResult::UnusableBit
}

fn main() {
    let mut usable = 0;
    let mut caught = 0;
    let mut total = 0;

    while total < 1000 {
        match send_one_qubit() {
            SimulationResult::UnusableBit => {}
            SimulationResult::UsableBit => usable += 1,
            SimulationResult::SpyDetected => {
                caught += 1;
                break;
            }
        }

        total += 1;
    }

    if caught != 0 {
        println!("Caught a spy!");
    } else {
        println!("{} usable bits extracted out of {}", usable, total);
    }
}
