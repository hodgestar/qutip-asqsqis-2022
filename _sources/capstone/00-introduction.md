# Capstone project: Simulating a Neutral Atom Quantum Device

Welcome to the capstone project at the end of this tutorial series. In the
project you're going to use everything you've learned this week to write
your own simulated quantum computer that can run simple circuits.

The goal is to eventually be able to write code that looks like this:

```
circuit = NeutralAtomCircuit(atoms)

# Pretty sequence of rotations the form a closed loop (on qubit 0)
circuit.rx(0, np.pi / 4)
circuit.rz(0, np.pi / 4)
circuit.ry(0, np.pi / 4)
```

And perhaps to add support for two qubit gates at the end.

You'll be given an outline of the `NeutralAtomCircuit` class and you'll
need to implement the physics for the gates.

Don't worry about going slowly. In the longer term its much more beneficial to
do a little bit that you understand well, rather than a lot that you understand
poorly.

## Guideline

This is quite a complex project, so you'll want to break it down into steps.
Don't be afraid to update the list of steps as you go. You'll be learning
more about the task at every step. Here is an initial list of steps for you
to start with:

- Reading the provided lecture notes:

  - [Neutral Atom Devices](./01-neutral-atom-devices)
  - [A Single Atom](./02-a-single-atom)

- Read section 2, *Description of the physical system* from
  [Pulser: An open-source package for the design of pulse sequences in programmable neutral-atom arrays](https://arxiv.org/abs/2104.15044). Don't try read the whole paper!

- Open the [solution outline](./03-neutral-atom-simulator) and read
  it thoroughly. Don't write any code -- just look at it and think about
  each method you'll need to implement.

- In another notebook, write the control Hamiltonian for a *single* qubit.
  You can start by just writing one term for the Hamiltonian, and then add
  the other term later.

- Use `sesolve` to evolution the single qubit in time with different control
  pulses. Plot the evolution either on the Bloch sphere or by taking
  expectation values.

- Implement `RX(theta)`, `RY(theta)`, and `RZ(theta)` for your single qubit.
  Try them out and check that they work individually.

- Run the simple circuit from the example we gave right at the top. If you
  display the states on the Bloch sphere, you should see a pretty pattern.

- Triple check that everything is working perfectly with a single qubit!
  Put many different initial states through many different combinations of
  operations. Keep trying to find errors in your implementation until no
  matter how hard you try, you can't think of new things to check.

- Do the step above again.

- Now add a second qubit. Add the required terms to the Hamiltonian. Update
  `rx`, `ry`, and `rz` to operate on the specified qubit.

- Check everything again!

If you get this far, you've done great! Take a moment to celebrate.


## Worked solutions

Once you have made a solid attempt on your own, here are two worked solutions
that you can look at if you get stuck or if you would like to read a
different solution for comparison.

The first is an exploratory notebook that shows how you might go about trying
things out and testing and improving your understanding of the system.

The second is a more complete cleaned up solution.

- [Neutral Atom Simulator (Exploratorion)](./04-neutral-atom-simulator-exploration.md)
- [Neutral Atom Simulator (Worked Solution)](./05-neutral-atom-simulator-worked-solution.md)


## Further things to do

It's unlikely that you'll get to everything above in a single day. But
if you do get to them, whether it's during the project day or afterwards,
here are some things you could try next:

- Implement the `CZ` (controlled Z-gate) on two qubits. This itself is an
  involved task:

  - Turn your qubits into a *qutrits* -- i.e.
    make each a state formed from the three basis states, $|0>$, $|1>$, and
    $|r>$ where $|r>$ is the Ryberg state.

  - Update you Hamiltonian and gates to act correctly on the new *qutrits*.

  - Check everything again!

  - Add the Rydberg interaction terms to your Hamiltonian.

  - Add control terms to the Hamiltonian for driving the *qutrits* into their
    $|r>$ states.

  - Implement the CZ gate using the protocol suggested in the paper above.

  - Try out simple circuits using `CNOT`.

  - Check everything as thoroughly as you can.

  - Celebrate! This was a big task.

Even further extensions:

- Implement more gates.

- Add dissipaton and other environmental noise to the implementation using
  `mesolve` and collapse operators. Look at how this affects the performance
  of gates and circuits. Try design better control pulses for these noisy gates.

- The QASM (Quantum Assembly) is a very simple programming language for
  describing quantum circuits. Learn how it works and try a program that
  will take QASM as input and run the circuit on your simulator.

- Implement a visual simulator that shows what is happening as a circuit is
  run.

All of the tasks in this section are quite significant projects themselves.
You are most definitely not required to work on them. Only start one if
you're very excited to work on it and have the free time to spend on it.

And that's the end! Keep using QuTiP, learning physics and asking questions.
