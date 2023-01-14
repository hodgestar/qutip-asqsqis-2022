---
jupytext:
  formats: md:myst,ipynb
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# 00: Introduction to QuTiP

+++

In this first tutorial, you'll play around with the core of QuTiP.

This includes constructing quantum states (e.g. qubits), quantum operators (e.g.
gates in a quantum circuit), Hamiltonians, and density matrices. Near the end
you'll simulate the evolution of a quantum system using QuTiP's numerical solver
for the Schrödinger equation and visualize the result. Along the way you'll find
eigenvalues and eigenvectors and plot states on the Bloch sphere.

Individually none of the tasks should require more than a few lines of code, but
there is a lot to explore. Take your time. Tackle individual steps in multiple
ways. Read the documentation on the functions you're using. Explore the objects
you're creating. Ask questions.

The topics covered in this tutorial will form the foundation we use during the
rest of the week -- and that you'll use when using QuTiP yourself once the
Summer school is over.

+++

## Tasks

- [Pre-flight check](#Pre-flight-check)
- [Import qutip](#Imports)
- [Accessing documentation](#Accessing-documentation)
- [States](#States)
- [Bloch sphere](#Bloch-sphere)
- [Operators](#Operators)
- [Eigenvalues and eigenstates](#Eigenvalues-and-eigenstates)
- [Density matrices](#Density-matrices)
- [Time-dependent operators](#Time-dependent-operators)
- [Tensor products and partial traces](#Tensor-products-and-partial-traces)
- [Hamiltonians](#Hamiltonians)
- [Solving the Schrödinger Equation using sesolve](#Solving-the-Schrödinger-Equation-using-sesolve)
- [Using functions to make your notebooks neater](#Using-functions-to-make-your-notebooks-neater)

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

## Accessing documentation

Once you're back home after the summer school, the QuTiP documentation will be your primary means of learning how to use QuTiP. In these tutorials, we'll encourage you to read the documentation first, so that by the end of the week, you'll be familiar with it and able to answer many of your own questions by reading it.

You'll read the documentation in two main ways:

- You can type a `?` after any Python object to read its documentation. You can try it right now by typing `qutip.ket?` into notebook cell and pressing enter.

- You can find the documentation online at [qutip.org](https://qutip.org/documentation.html). Open it now and have a look. You can also download it as a PDF in case you need to work offline.

+++

## States

- Creating and examining quantum states
- Arithmetic with states
- Checking whether states are equal and whether states with a different global phase are equal
- What is a Qobj?

Functions you'll want to try out in this section are:

- `qutip.ket`
- `qutip.bra`
- `qutip.basis`

These functions all return a "quantum object" or `Qobj` which is the data type QuTiP uses to store all kinds of quantum objects (states, operators, density matrices). A `Qobj` has some many methods and attributes that are useful.

Some attributes to look at include `type`, `dims`, `shape`, `isket`, and `isbra`.

A few important methods are `.norm()`, `.unit()`, `.dag()`, and `.overlap(other)`.

You can use the arithmetic operators `*` and `+` as you would expect. Try out some simple numerical expressions from this morning's theory lecture.

There are many more methods and attributes, but these are enough to start with.

```{code-cell}

```

## Bloch sphere

- Plotting 2-level states on the Bloch sphere
- Plotting many states at once
- Changing the colors
- Changing the transparency

QuTiP provides a class called `qutip.Bloch`. You can create an instance of it and add states to plot with `.add_states(...)`. Once you've added all the states you want to, you can all `.show()` to render it.

The Bloch sphere is very useful for visualizing the state of a two-level system, and you'll use it regularly.

You can change the colors and transparency of the arrow shown for each state by setting the `.vector_color` attribute to a list of colors, and the `.vector_alpha` attribute to a list of transparency values.

```{code-cell}

```

## Operators

- Creating and examining quantum operators
- Calculating expectation values
- Making projectors

Functions to read about and try out in this section are:

- `qutip.sigmax`, `qutip.sigmay`, `qutip.sigmaz`
- `qutip.qeye`
- `qutip.qip.operations.hadamard_transform`

Operators are also represented by `Qobj` instances and there are some attitional `Qobj` attributes that are relevant now: `type`, `isoper`, `dims`, `shape`, `isherm`. Look at those you've seen already and compare them with what you've already seen for states.

There are some additional `Qobj` methods to try out on operators -- `.diag()`, `.tr()` -- and some to revisit -- `.dag()`, `.norm()`, `.unit()`.

You can create a projection operator by calling `.proj()` on a state. Try this out. Also try creating the projector yourself without calling `.proj()`.

```{code-cell}

```

## Eigenvalues and eigenstates

- Eigenvalues
- Eigenstates

It's often important to know the eigenvalues and eigenstates of operators. QuTiP makes this easy. Read about and try out the `.eigenenergies()` and `.eigenstates()` methods of `Qobj`.

```{code-cell}

```

## Density matrices

- Creating density matrices
- Plotting 2-level density matrices on the Bloch sphere
- Plot some pure states
- Plot some mixed states

QuTiP of course supports density matrices in addition to pure states. You can convert a state to a density matrix by calling `qutip.ket2dm`. You should also try construct some density matrices directly using the expression you saw in the theory lectures earlier.

One can also plot density matrices on the Bloch. Try it out! Once you've tried a few, attempt to predict how each density matrix will be represented before calling `.show()` and see if you are right. If you're not right, read the [Wikipedia article](https://en.wikipedia.org/wiki/Bloch_sphere#Density_operators) and try again.

```{code-cell}

```

## Time-dependent operators

- Create a time-dependent operator
- What is a QobjEvo?
- Adding arguments

In QuTiP, time-dependent operators are represented by `qutip.QobjEvo` and *not* by `Qobj`. `QobjEvo` takes in a *list* of terms that are *added* together. Each term consists of a `Qobj` (the constant part) and a time-dependent coefficient function, like so:

```
sigmax_t = QobjEvo([[sx, "t"]])
```

Here `sx` is the operator and `"t"` is the time-dependent coefficient. Note that `"t"` is a string. `QobjEvo` turns this string into a very fast function for you. It's also just very convenient to be able to write the coefficient function compactly.

Try out some other operators and some other coefficients, for example, `"cos(t)"` and `"sin(t)"`.

You can evaluate a `QobjEvo` at a particular time by calling, for example, `sigmax_t(2)` or `sigmax_t(0.1)`, or any other time. Create a state object that depends on time and plot its evolution on the Bloch sphere.

Lastly you can define coefficients that depend on arguments, such as `"cos(w * t)"`. You need to supply an initial value for the arguments when creating the `QobjEvo`, like so:

```
sigmax_cos_t = QobjEvo([[sx, "cos(w * t)"]], args={"w": np.pi})
```

You can supply a different argument when calling the `QobjEvo`, for example, `sigmax_cos_t(0.1, {"w": np.pi / 4})`.

```{code-cell}

```

## Tensor products and partial traces

- Construct a state with two or three qubits
- Construct an operator on multiple qubits
- Construct an entangled state
- Take the partial trace of a product state
- Take the partial trace of an entangled state

QuTiP has two very useful methods for combining and separating sub-systems of a larger quantum state. Read about `qutip.tensor` and `Qobj.ptrace` and try them out.

```{code-cell}

```

## Hamiltonians

- Building Hamiltonians
- Plotting eigenvalues

Now use everything you've learned so far to construct a time-dependent Hamiltonian and plot the evolution of its eigenvalues over time.

```{code-cell}

```

## Solving the Schrödinger Equation using sesolve

- Plotting the states
- Plotting the states on the Bloch sphere
- Plotting expectation values


QuTiP's function for solving the Schrödinger equation is `qutip.sesolve`. `sesolve` returns a `result` object that has a `.states` attribute containing the states of the system at each requested time.

Or you can pass a list of expectation operators and then `result` will have a `.expect` attribute that contains the expectation values of each operator at each time. `.expect[0]` will contain the values of the first operator, `.expect[1]` all the values of the second operator, and so on.

Try it out and plot the results!

If you need to plot values, you can use `plt.plot(x, y, label="Label for this plot")`. You can read the help for `plt.plot` to find out more.

```{code-cell}

```

## Using functions to make your notebooks neater

- Write a function for displaying the eigenstates of an operator
- Write a function for displaying states on the Bloch sphere

```{code-cell}

```
