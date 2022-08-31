---
jupytext:
  formats: md:myst,ipynb
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

# 00: Introduction to QuTiP

+++

## Tasks

- [Pre-flight check](#Pre-flight-check)
  - Open Jupyter
  - Run the smoke test notebook
- [Import qutip](#Imports)
  - Imports
- [States](#States)
  - Creating and examining quantum states
  - Arithmetic with states
  - What is a Qobj?
- [Bloch sphere](#Bloch-sphere)
  - Plotting 2-level states on the Bloch sphere
  - Plotting many states at once
  - Changing the colors
  - Changing the transparency
- [Operators](#Operators)
  - Creating and examining quantum operators
  - Calculating expectation values
  - Making projectors
- [Eigenvalues and eigenstates](#Eigenvalues-and-eigenstates)
  - Eigenvalues
  - Eigenstates
- [Density matrices](#Density-matrices)
  - Creating density matrices
  - Plotting 2-level density matrices on the Bloch sphere
  - Plot some pure states
  - Plot some mixed states
- [Time-dependent operators](#Time-dependent-operators)
  - Create a time-dependent operator
  - What is a QobjEvo?
  - Adding arguments
- [Tensor products and partial traces](#Tensor-products-and-partial-traces)
  - Construct a state with two or three qubits
  - Construct an operator on multiple qubits
  - Construct an entangled state
  - Take the partial trace of a product state
  - Take the partial trace of an entangled state
- [Hamiltonians](#Hamiltonians)
  - Building Hamiltonians
  - Determining eigenvalues and eigenstates
  - Plotting eigenvalues
- [Solving the Schrödinger Equation using sesolve](#Solving-the-Schrödinger-Equation-using-sesolve)
  - Plotting the states
  - Plotting the states on the Bloch sphere
  - Plotting expectation values
- [Using functions to make your notebooks neater](#Using-functions-to-make-your-notebooks-neater)
  - Write a function for displaying the eigenstates of an operator
  - Write a function for displaying states on the Bloch sphere

+++

## Pre-flight check

+++

If you haven't installed Python, QuTiP, Jupyter and the other requirements yet, do so now by following one of the [sets of instructions](https://github.com/hodgestar/qutip-asqsqis-2022/).

You can test that your setup is working by running the [smoke test notebook](https://github.com/hodgestar/qutip-asqsqis-2022/blob/main/smoke-test.ipynb).

+++

## Imports

```{code-cell}
%matplotlib inline
import matplotlib.pyplot as plt
```

```{code-cell}
import qutip
import numpy as np
```

+++ {"tags": ["hide-cell"]}

## Helper functions

```{code-cell}
:tags: [hide-cell]

def display_list(a, title=None):
    """ Display a list of items using Jupyter's display function. """
    if title:
        print(title)
        print("=" * len(title))
    for item in a:
        display(item)
```

```{code-cell}
:tags: [hide-cell]

def print_attrs(obj, name, attrs):
    """ Print the listed attributes of an object. """
    for attr in attrs:
        print(f"{name}.{attr}: {getattr(obj, attr)}")
```

```{code-cell}
:tags: [hide-cell]

def show_bloch(states, vector_color=None, vector_alpha=None):
    """ Show states or density matrices on the Bloch sphere. """
    bloch = qutip.Bloch()
    bloch.add_states(states)
    if vector_color is not None:
        if vector_color == "rb":
            vector_color = [
                (c, 0, 1 - c) for c in np.linspace(0.1, 1.0, len(states))
            ]
        bloch.vector_color = vector_color
    if vector_alpha is not None:
        bloch.vector_alpha = vector_alpha
    bloch.show()
```

```{code-cell}
:tags: [hide-cell]

def display_eigenstates(op):
    """ Display the eigenvalues and eigenstates of an operator. """
    evals, evecs = op.eigenstates()
    print("Eigenvalues:", evals)
    print()
    print("Eigenstates")
    print("===========")
    for v in evecs:
        display(v)
```

## States

- Creating and examining quantum states
- Arithmetic with states
- What is a Qobj?

```{code-cell}
:tags: [hide-cell]

# a qubit in the |0> state
ket = qutip.ket("0")
```

```{code-cell}
:tags: [hide-cell]

# show some of the properties relevant to kets
display(ket)
print_attrs(ket, "ket", ["type", "dims", "shape", "isket", "isbra"])
```

```{code-cell}
:tags: [hide-cell]

# show the adjoint of the ket -- a bra
ket.dag()
```

```{code-cell}
:tags: [hide-cell]

# construct the bra directly
bra = qutip.bra("0")
```

```{code-cell}
:tags: [hide-cell]

# a qubit in the |1> state
qutip.ket("1")
```

```{code-cell}
:tags: [hide-cell]

# construct |0> and |1> using qutip.basis
q0 = qutip.basis(2, 0)
q1 = qutip.basis(2, 1)
display_list([q0, q1])
```

```{code-cell}
:tags: [hide-cell]

# show that Qobj's with the same values are equal
q0 == q0
```

```{code-cell}
:tags: [hide-cell]

# show that Qobj's with different global phases are not equal
display(q0 * 1j)
print("Equal:", q0 == q0 * 1j)
print("Overlap:", q0.overlap(q0 * 1j))
print("Overlap is 1:", np.abs(q0.overlap(q0 * 1j)) == 1.0)
```

```{code-cell}
:tags: [hide-cell]

# orthogonal vectors have 0 overlap
q0.overlap(q1)
```

```{code-cell}
:tags: [hide-cell]

# for states, the overlap is the inner product
print("q0 · q1:", (q0.dag() * q1).full()[0, 0])
print("q0 · q1:", (q1.dag() * q0).full()[0, 0])
print("(1j * q0) · q0:", ((1j * q0).dag() * q0).full()[0, 0])
```

```{code-cell}
:tags: [hide-cell]

# adding states
q0 + q1
```

```{code-cell}
:tags: [hide-cell]

# calculating the norm
(q0 + q1).norm()
```

```{code-cell}
:tags: [hide-cell]

# manually normalizing a state to unit length
(q0 + q1) / (q0 + q1).norm()
```

```{code-cell}
:tags: [hide-cell]

# using .unit() to normalize the length of a state
(q0 + q1).unit()
```

## Bloch sphere

- Plotting 2-level states on the Bloch sphere
- Plotting many states at once
- Changing the colors
- Changing the transparency

```{code-cell}
:tags: [hide-cell]

# show a single state
show_bloch(q0)
```

```{code-cell}
:tags: [hide-cell]

# show both basis states
show_bloch([q0, q1])
```

```{code-cell}
:tags: [hide-cell]

# show many states
theta = np.linspace(0, np.pi, 20)
states = [q0 + np.exp(1j * th) * q1 for th in theta]
states = [s.unit() for s in states]
show_bloch(states)
```

```{code-cell}
:tags: [hide-cell]

# How to remove the colours:
show_bloch(
    states,
    vector_color=[(0, 0, 1.0)] * len(states),
)
```

```{code-cell}
:tags: [hide-cell]

# Use transparency to show where states are in the list
show_bloch(
    states,
    vector_color=[(0, 0, 1.0)] * len(states),
    vector_alpha=np.linspace(0.1, 1.0, len(states)),
)
```

```{code-cell}
:tags: [hide-cell]

# Use color map to show where states are in the list
show_bloch(
    states,
    vector_color=[(c, 0, 1 - c) for c in np.linspace(0.1, 1.0, len(states))],
)
```

## Operators

- Creating and examining quantum operators
- Calculating expectation values
- Making projectors

```{code-cell}
:tags: [hide-cell]

# the Pauli X operator
sx = qutip.sigmax()
```

```{code-cell}
:tags: [hide-cell]

# show some of the properties relevant to operators
display(sx)
print_attrs(sx, "sx", ["type", "isoper", "dims", "shape", "isherm"])
```

```{code-cell}
:tags: [hide-cell]

# show the adjoint
sx.dag()
```

```{code-cell}
:tags: [hide-cell]

# show the diagonal entries
sx.diag()
```

```{code-cell}
:tags: [hide-cell]

# show the trace
sx.tr()
```

```{code-cell}
:tags: [hide-cell]

# show the (trace) norm
sx.norm()
```

```{code-cell}
:tags: [hide-cell]

# apply sx to the basis states -- it swaps them
display_list([sx * q0, sx * q1])
```

```{code-cell}
:tags: [hide-cell]

# construct the other Pauli matrices and the identity matrix
sy = qutip.sigmay()
sz = qutip.sigmaz()
I = qutip.qeye(2)
display_list([I, sx, sy, sz])
```

```{code-cell}
:tags: [hide-cell]

# Hadamard gate
H = qutip.qip.operations.hadamard_transform()
H
```

```{code-cell}
:tags: [hide-cell]

# Making projectors
q0.proj()
```

```{code-cell}
:tags: [hide-cell]

q0 * q0.dag()
```

## Eigenvalues and eigenstates

- Eigenvalues
- Eigenstates

```{code-cell}
:tags: [hide-cell]

display_eigenstates(sz)
```

```{code-cell}
:tags: [hide-cell]

display_eigenstates(sx)
```

```{code-cell}
:tags: [hide-cell]

display_eigenstates(sy)
```

## Density matrices

- Creating density matrices
- Plotting 2-level density matrices on the Bloch sphere
- Plot some pure states
- Plot some mixed states

```{code-cell}
:tags: [hide-cell]

rho0 = qutip.ket2dm(q0)
rho0
```

```{code-cell}
:tags: [hide-cell]

qutip.ket2dm(q0 * np.exp(0.1j))
```

```{code-cell}
:tags: [hide-cell]

q0 * q0.dag()  # global phases cancel
```

```{code-cell}
:tags: [hide-cell]

rho = q0 * q0.dag() + q1 * q1.dag()
rho
```

```{code-cell}
:tags: [hide-cell]

rho.norm()
```

```{code-cell}
:tags: [hide-cell]

rho = rho.unit()
rho
```

```{code-cell}
:tags: [hide-cell]

bloch = qutip.Bloch()

bloch.add_states([qutip.ket2dm(q0), qutip.ket2dm(q1)])
bloch.add_states([rho])

bloch.show()
```

```{code-cell}
:tags: [hide-cell]

# p = 0.5 * (I + n . (sx, sy, sz))
# if p == 0.5 I, then n = (0, 0, 0)
rho == 0.5 * qutip.qeye(2)
```

```{code-cell}
:tags: [hide-cell]

basis = [qutip.qeye(2), sx, sy, sz]
for v in basis:
    display(v)
for v in basis:
    assert (v * v).tr() == 2.0
for v in basis:
    for w in basis:
        if v != w:
            assert (v * w).tr() == 0.0
```

## Time-dependent operators

- Create a time-dependent operator
- What is a QobjEvo?
- Adding arguments

```{code-cell}
:tags: [hide-cell]

# Create a simple quantum object that represents sx * t
sigmax_t = qutip.QobjEvo([[sx, "t"]])
```

```{code-cell}
:tags: [hide-cell]

display(sigmax_t)
display_list([
    sigmax_t(t) for t in [0.0, 0.25, 0.5, 1.0]
])
```

```{code-cell}
:tags: [hide-cell]

# Multiplication:
display_list([
    (0.1 * sigmax_t)(1.0),
    (0.5 * sigmax_t)(1.0),
])
```

```{code-cell}
:tags: [hide-cell]

# Addition:
sigmay_t = qutip.QobjEvo([[sy, "t"]])
display_list([
    (sigmax_t * sigmay_t)(1.0),
])
```

```{code-cell}
:tags: [hide-cell]

# Take the conjugate and the adjoint:
display_list([
    sigmay_t.conj()(0.5),
    sigmay_t.dag()(0.5),
])
```

```{code-cell}
:tags: [hide-cell]

# Plotting the evolution under sigmax_t
states = [(1j * sigmax_t(t)).expm() * q0 for t in np.linspace(0, np.pi / 2, 10)]
show_bloch(states, vector_color="rb")
```

```{code-cell}
:tags: [hide-cell]

# Including a constant term
sz_plus_sy_t = qutip.QobjEvo([sz, [sy, "t"]])
display_list([
    sz_plus_sy_t(0.1),
    sz_plus_sy_t(1.0),
])
```

```{code-cell}
:tags: [hide-cell]

# Using a function like sin()
sz_sin_t = qutip.QobjEvo([[sz, "sin(t)"]])
display_list([
    sz_sin_t(t) for t in np.linspace(0, np.pi / 2, 4)
])
```

```{code-cell}
:tags: [hide-cell]

# Supplying arguments other than time
sz_g_sin_t = qutip.QobjEvo([[sz, "g * sin(t)"]], args={"g": 2.0})
display_list([
    sz_g_sin_t(np.pi / 2, args={"g" : g}) for g in [0, 0.2, 0.4]
])
```

## Tensor products and partial traces

- Construct a state with two or three qubits
- Construct an operator on multiple qubits
- Construct an entangled state
- Take the partial trace of a product state
- Take the partial trace of an entangled state

```{code-cell}
:tags: [hide-cell]

q00 = qutip.tensor([q0, q0])
q00
```

```{code-cell}
:tags: [hide-cell]

pt_q0 = q00.ptrace(0)
pt_q0
```

```{code-cell}
:tags: [hide-cell]

q01 = qutip.tensor([q0, q1])
q01
```

```{code-cell}
:tags: [hide-cell]

display_list([
    q01.ptrace(0),
    q01.ptrace(1),
])
```

## Hamiltonians

- Building Hamiltonians
- Plotting eigenvalues

```{code-cell}
:tags: [hide-cell]

H = qutip.QobjEvo([[sx, "a * cos(t)"], [sy, "b * sin(t)"]], args={"a": 0.2, "b": 0.1})
```

```{code-cell}
:tags: [hide-cell]

tlist = np.linspace(0, np.pi, 40)
evals = np.array([H(t).eigenenergies() for t in tlist])

plt.plot(tlist, evals[:, 0])
plt.plot(tlist, evals[:, 1])
plt.ylabel("Energy level")
plt.xlabel("t");
```

## Solving the Schrödinger Equation using sesolve

- Plotting the states
- Plotting the states on the Bloch sphere
- Plotting expectation values

```{code-cell}
:tags: [hide-cell]

tlist = np.linspace(0, np.pi, 40)
result = qutip.sesolve(H, q0, tlist)
```

```{code-cell}
:tags: [hide-cell]

show_bloch(result.states, vector_color="rb")
```

## Using functions to make your notebooks neater

- Write a function for displaying the eigenstates of an operator
- Write a function for displaying states on the Bloch sphere

```{code-cell}
:tags: [hide-cell]

# See the helper functions at the top for examples.
```
