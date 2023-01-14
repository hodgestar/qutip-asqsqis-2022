---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# 03: Neutral Atom Simulator (Outline)

+++

## Imports

```{code-cell} ipython3
import qutip
import numpy as np
```

```{code-cell} ipython3
%matplotlib inline
import matplotlib.pyplot as plt
```

## Outline of an implementation

```{code-cell} ipython3
class NeutralAtomCircuit:
    """ A neutral atom circuit laser pulse simulator. """
    def __init__(self, atoms):
        self.n = len(atoms)
        self.atoms = atoms
        self.pulses = []
        self.t = 0

    def state(self, s):
        """ Return a ket corresponding to the given string description. """

    def run(self, psi, tlist=None, step_rate=20):
        """ Run the circuit. """        

    # Gates implemented directly using control pulses:

    def rx(self, i, theta, dt=1):
        """ Add a laser pulse that rotates qubit i an angle theta around the X axis. """

    def ry(self, i, theta, dt=1):
        """ Add a laser pulse that rotates qubit i an angle theta around the Y axis. """

    def rz(self, i, theta, dt=1):
        """ Add a laser pulse that rotates qubit i an angle theta around the Z axis. """

    def cz(self, control, target, dt=1):
        """ Add a pulse that performs a controlled-Z phase flip operation on two qubits. """

    # Gates implemented on top of other gates:

    def h(self, i):
        """ Apply a hadamard gate to qubit i. """

    def cnot(self, i):
        """ Apply a CNOT gate to qubit i. """
```

# Implementation

```{code-cell} ipython3
# Your code goes here
```

# Example uses

```{code-cell} ipython3
atoms = [(1, 0), (0, 1)]
circuit = NeutralAtomCircuit(atoms)

# Example circuit 1: Pretty sequence of rotations the form a closed loop (qubit 0)
#circuit.rx(0, np.pi / 4)
#circuit.rz(0, np.pi / 4)
#circuit.ry(0, np.pi / 4)

# Example circuit 2: Pretty sequence of rotations the form a closed loop (qubit 1)
#circuit.rx(1, np.pi / 4)
#circuit.rz(1, np.pi / 4)
#circuit.ry(1, np.pi / 4)

# Example circuit 3:
# Combine the two pervious circuits into one

# Example circuit 4: Apply a CNOT
# CNOT = (I ⊗ H)CZ(I ⊗ H)
# H = X Ry(pi/2)
circuit.rx(1, np.pi / 2)
circuit.ry(1, np.pi / 4)
circuit.cz(0, 1)
circuit.rx(1, np.pi / 2)
circuit.ry(1, np.pi / 4)

psi = circuit.state("11")
result = circuit.run(psi, step_rate=20)
```
