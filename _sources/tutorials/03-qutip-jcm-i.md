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

# 03: Jaynes-Cummings Model and the Rotating Wave Approximation

+++

In this tutorial we will construct the Jaynes-Cummings Hamiltonian (with and without the RWA) and see how the system evolves under the Schrodinger equation (that is, without dissipation) .
We will use this to investigate the limits of the RWA in the JCM.

+++

The Jaynes-Cumming model is the simplest possible model of quantum mechanical light-matter interaction, describing a single two-level atom interacting with a single electromagnetic cavity mode. The full Hamiltonian for the system (in dipole interaction form) is given by

$$H = H_{\rm atom} + H_{\rm cavity} + H_{\rm interact}$$

The atom Hamiltonian we use in this case is

$$\frac{1}{2} \hbar \omega_{a} \sigma_z$$

where $\omega_a$ is the system frequency.

Note that the Hamiltonian for the atom may take numerous forms. Any Hermitean operator on a two-level state is possible, but it is useful to nomalize the operator so that the difference between its eigenvalues is $1$ so that $\omega_a$ has consistent units.

The cavity Hamiltonian is given by

$$H_{\rm cavity} = \hbar \omega_c a^\dagger a$$

where $\omega_c$ is $\omega_a$ is the frequencies of the cavity and $a$ and $a^\dagger$ are the annihilation and creation operators of the cavity respectively.

The interaction Hamiltonian is given by

$$H_{\rm interact} = \hbar g(a^\dagger + a)(\sigma_- + \sigma_+)$$

or with the rotating-wave approximation

$$H_{\rm interact-RWA} = \hbar g(a^\dagger\sigma_- + a\sigma_+)$$

where $\sigma_-$ and $\sigma_+$ are the annihilation and creation operators for the atom respectively.

Note that in this notebook we will work in units where $\hbar=1$.

+++

## Tasks

- [Construct the Hamiltonian](#construct-the-hamiltonian)
- [Solve the Schrodinger equation](#solve-the-schrodinger-equation)
- [Visualise the evolution](#visualise-the-evolution)
- [Compare the RWA and non-RWA](#compare-the-rwa-and-non-rwa)

+++

## Imports

```{code-cell} ipython3
%matplotlib inline
import matplotlib.pyplot as plt
```

```{code-cell} ipython3
import qutip
import numpy as np
```

+++ {"tags": ["hide-cell"]}

## Helper functions

```{code-cell} ipython3
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

```{code-cell} ipython3
:tags: [hide-cell]

def jcm_h(wc, wa, g, N, atom):
    """ Construct the Jaynes-Cummings Hamiltonian (non-RWA). """
    a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
    sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())
    atom = qutip.tensor(qutip.qeye(N), atom)
    
    H = wc * a.dag() * a + wa * atom + g * (a.dag() + a) * (sm + sm.dag())
    return H
```

```{code-cell} ipython3
:tags: [hide-cell]

def jcm_rwa_h(wc, wa, g, N, atom):
    """ Construct the Jaynes-Cummings Hamiltonian (RWA). """
    a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
    sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())
    atom = qutip.tensor(qutip.qeye(N), atom)

    H = wc * a.dag() * a + wa * atom + g * (a.dag() * sm + a * sm.dag())
    return H
```

## Construct the Hamiltonian

- add variables for the atom and cavity parameters
- create the operators for the JCM Hamiltonian
- combine into the JCM Hamiltonian (no RWA)
- look at the energy eigenvalues and eigenstates of the Hamiltonian

Here are some example parameter values to start with:

$
  \omega_c = 2 \pi \\
  \omega_a = 2 \pi \\
  g = 0.05 \cdot 2 \pi \\
$

```{code-cell} ipython3
:tags: [hide-cell]

# system parameters
wc = 1.0 * 2 * np.pi  # cavity frequency
wa = 1.0 * 2 * np.pi  # atom frequency
g = 0.05 * 2 * np.pi # 0.05 * 2 * np.pi  # coupling strength
N = 15  # number of cavity fock states
```

```{code-cell} ipython3
:tags: [hide-cell]

# operators for the JCM
# this is the annihilation operator for the cavity (note it acts trivally on the atom) 
a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
# note the creation operator for the cavity is given by a.dag() - that is a "dagger"

# this is annihilation operator for the atom (note it acts trivally on the cavity)
sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())
```

```{code-cell} ipython3
:tags: [hide-cell]

# hamiltonian of atom
H_atom = 0.5 * qutip.sigmaz()
```

```{code-cell} ipython3
:tags: [hide-cell]

# hamiltonian (non-rwa)
H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() + a) * (sm + sm.dag())
```

```{code-cell} ipython3
:tags: [hide-cell]

display_eigenstates(H)
```

## Solve the Schrodinger equation

- create the initial state of the system (use the state with no photons and the spin system in its excited state)
- evolve the system for some time, saving the result.

If you need to remind yourself of how sesolve, remember that you can type `qutip.sesolve?` into a notebook cell to bring up the documentation.

```{code-cell} ipython3
:tags: [hide-cell]

# initial state
psi0 = qutip.basis([N, 2], [0, 0])  # start with an excited atom
psi0
```

```{code-cell} ipython3
:tags: [hide-cell]

# Solve using sesolve
tlist = np.linspace(0, 25, 101)
result = qutip.sesolve(H, psi0, tlist)
```

## Visualise the evolution

- create expectation operators for observing the state of the system.
- add these to sesolve
- plot the expectation values together on a set of axes
- change the value of g and see how it affects the period

Two good expectation operators to use are the projectors on the light and matter sub-systems.

```{code-cell} ipython3
:tags: [hide-cell]

# Operators to determine the expectation values of:
eop_a = a.dag() * a  # light
eop_sm = sm.dag() * sm  # matter
```

```{code-cell} ipython3
:tags: [hide-cell]

# Solve using sesolve
tlist = np.linspace(0, 25, 101)
result = qutip.sesolve(H, psi0, tlist, e_ops=[eop_a, eop_sm])
```

```{code-cell} ipython3
:tags: [hide-cell]

plt.plot(tlist, result.expect[0], label="Light")
plt.plot(tlist, result.expect[1], label="matter")
plt.legend();
```

```{code-cell} ipython3
:tags: [hide-cell]

# vary g
H = jcm_h(wc, wa, 0.1 * 2 * np.pi, N, H_atom)

tlist = np.linspace(0, 25, 101)

result = qutip.sesolve(H, psi0, tlist, e_ops=[eop_a, eop_sm])

plt.plot(tlist, result.expect[0], label="Light")
plt.plot(tlist, result.expect[1], label="matter")
plt.legend();
```

## Compare the RWA and non-RWA

- construct another Hamiltonian that uses the RWA
- evolve the system under the RWA Hamiltonian.
- add the results to the plot
- experiment with parameters to determine where the RWA non-RWA diverges

```{code-cell} ipython3
:tags: [hide-cell]

# Construct the RWA Hamiltonian
H = jcm_h(wc, wa, g, N, H_atom)
H_RWA = jcm_rwa_h(wc, wa, g, N, H_atom)
```

```{code-cell} ipython3
:tags: [hide-cell]

# Hmm. The results look the same.
tlist = np.linspace(0, 25, 101)

result = qutip.sesolve(H, psi0, tlist, e_ops=[eop_a, eop_sm])
plt.plot(tlist, result.expect[0], label="Light")
plt.plot(tlist, result.expect[1], label="matter")

result_rwa = qutip.sesolve(H_RWA, psi0, tlist, e_ops=[eop_a, eop_sm])
plt.plot(tlist, result_rwa.expect[0], label="Light (RWA)")
plt.plot(tlist, result_rwa.expect[1], label="matter (RWA)")

plt.legend();
```

```{code-cell} ipython3
:tags: [hide-cell]

# Try with different frequencies w_c and w_a:
tlist = np.linspace(0, 25, 101)
f = 0.9

H = jcm_h(wc, f * wa, g, N, H_atom)
H_RWA = jcm_rwa_h(wc, f * wa, g, N, H_atom)

result = qutip.sesolve(H, psi0, tlist, e_ops=[eop_a, eop_sm])
plt.plot(tlist, result.expect[0], label="Light")
plt.plot(tlist, result.expect[1], label="matter")

result_rwa = qutip.sesolve(H_RWA, psi0, tlist, e_ops=[eop_a, eop_sm])
plt.plot(tlist, result_rwa.expect[0], label="Light (RWA)")
plt.plot(tlist, result_rwa.expect[1], label="matter (RWA)")

plt.legend();
```

## Links for further study

There is an excellent paper [The Jaynes-Cummings model and its descendants](https://arxiv.org/abs/2202.00330) by Larson and Mavrogordatos, that reviews the Jaynes-Cummings model and its many variations. You can try reading this paper and implementing some of the simpler variants in QuTiP.

If you do, please consider polishing your notebook, adding good explanations to it and submitting it as an example for others by opening a pull request for it at https://github.com/qutip/qutip-tutorials/.

This paper covers a *lot* of work, so don't expect to understand all of it quickly. Start with the earlier sections. We will explore the Jaynes-Cumming model more in the remaining tutorials.
