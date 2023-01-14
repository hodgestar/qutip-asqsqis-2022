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

# 04: Adding dissipation to the Jaynes-Cummings model

+++

In this tutorial we will add interactions with the environment to the Jaynes-Cummings model in the form of dissipation from the atom and the cavity. We will use the Wigner function to visualise the modes of the cavity certain times in the evolution to see how these evolve over time.

+++

## Tasks

- [Construct collapse operators](#Construct-collapse-operators)
- [Solve the Lindblad Master equation](#Solve-the-linblad-master-equation)
- [Visualise cavity modes as Wigner functions](#Visualise-cavity-modes-as-wigner-functions)
- [Finding the steady state](#Finding-the-steady-state)

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

## Helper functions

```{code-cell}
def jcm_h(wc, wa, g, N, atom):
    """ Construct the Jaynes-Cummings Hamiltonian (non-RWA). """
    a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
    sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())
    atom = qutip.tensor(qutip.qeye(N), atom)
    
    H = wc * a.dag() * a + wa * atom + g * (a.dag() + a) * (sm + sm.dag())
    return H
```

```{code-cell}
def jcm_rwa_h(wc, wa, g, N, atom):
    """ Construct the Jaynes-Cummings Hamiltonian (RWA). """
    a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
    sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())
    atom = qutip.tensor(qutip.qeye(N), atom)

    H = wc * a.dag() * a + wa * atom + g * (a.dag() * sm + a * sm.dag())
    return H
```

## Re-cap of solving without dissipation

Let's re-cap what we did in the previous tutorial and solve the Schrödinger equation for the Jaynes-Cummings Hamiltonian without dissipation.

```{code-cell}
# Jaynes-Cummings parameters
# system parameters
wc = 1.0 * 2 * np.pi  # cavity frequency
wa = 1.0 * 2 * np.pi  # atom frequency
g = 0.05 * 2 * np.pi  # coupling strength
N = 15  # number of cavity fock states

# Atom hamiltonian
H_atom = 0.5 * qutip.sigmaz()

H = jcm_h(wc, wa, g, N, H_atom)
```

```{code-cell}
# system operators
a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())

# relaxation operators
a  # cavity relaxation
a.dag()  # cavity excitation
sm = sm # qubit relaxation
```

```{code-cell}
# Operators to determine the expectation values of:
eop_a = a.dag() * a  # light
eop_sm = sm.dag() * sm  # matter
```

```{code-cell}
psi0 = qutip.basis([N, 2], [0, 0])  # start with an excited atom
tlist = np.linspace(0, 50, 101)

result = qutip.sesolve(H, psi0, tlist, e_ops=[eop_a, eop_sm])

plt.plot(tlist, result.expect[0], label="Light")
plt.plot(tlist, result.expect[1], label="matter")
plt.legend();
```

## Construct collapse operators

Now we will create the collapse operators for the Lindblad master equation.

These are also called jump operators or Lindblad operators.

- Define variables for the cavity dissipation rate ($\kappa$), atom dissipation rate ($\gamma$), and average number of thermal bath excitations ($N_{\mathrm{thermal}}$).

- Create a list for adding the collapse operators to.

- Create the cavity relaxation collapse operator. Consider the appropriate operator and then determine the rate based on the dissipation and thermal occupation level.

- Create the cavity excitation collapse operator, again considering the appropriate operator and rate.

- Create the atom relaxation collapse operator, again considering the appropriate operator and rate.

- Add the operators to the list if the corresponding rate is greater than zero.

Some suggested values for the initial dissipation rates and averge bath excitations:

$
  \kappa = 0.005 \\
  \gamma = 0.05 \\
  N_{\mathrm{thermal}} = 5 \\
$

```{code-cell}

```

## Solve the Linblad Master equation

- Evolve the system, starting from excited atom state, using mesolve, utilising the collapse operators from above and saving the result.

- Plot the expectation values for the atom and cavity excitations.

- What happens if you extent the evolution time?

```{code-cell}

```

## Visualise cavity modes as Wigner functions

- Note the times when the cavity excitation is at a maximum or minimum.

- Use the master equation solver to generate density matrices of the time evolution (this is the default mode when no expectation operators are given) for these times.

- For each of these density matrices, trace out (partial trace) the atom from the density matrix and create and plot the Wigner function.

```{code-cell}

```

# Finding the steady state

When you extended the time evolution, it looked like the system was heading towards a steady state. Let's find it!

You can use `qutip.steadystate` to determine it. Read the documentation by typing `qutip.steadystate?`.

Plot the expectation values and Wigner function of the steady state and compare them to the values you saw from the evolution.

```{code-cell}

```

## Links for further study

There is an excellent paper [The Jaynes-Cummings model and its descendants](https://arxiv.org/abs/2202.00330) by Larson and Mavrogordatos, that reviews the Jaynes-Cummings model and its many variations. You can try reading this paper and implementing some of the simpler variants in QuTiP.

If you do, please consider polishing your notebook, adding good explanations to it and submitting it as an example for others by opening a pull request for it at https://github.com/qutip/qutip-tutorials/.

This paper covers a *lot* of work, so don't expect to understand all of it quickly. Start with the earlier sections. We will explore the Jaynes-Cumming model more in the remaining tutorials.
