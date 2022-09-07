---
jupytext:
  formats: ipynb,md:myst
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

# 05: Jaynes-Cummings Ultra-strong Coupling

+++

XXX Description

+++

## Tasks

- XXX

+++

# Imports

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
    sm = qutip.tensor(qutip.qeye(N), qutip.destroy(2))
    atom = qutip.tensor(qutip.qeye(N), atom)
    
    H = wc * a.dag() * a + wa * atom + g * (a.dag() + a) * (sm + sm.dag())
    return H
```

```{code-cell}
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

# XXX

```{code-cell}
# Dissipation parameters
# We use stronger dissipation to show long-term behaviour in shorter times
kappa = 0.2 # cavity dissipation rate
gamma = 0.2 # atom dissipation rate
n_th_a = 0 # avg number of thermal bath excitation
```

```{code-cell}
# system parameters
wc = 1.0 #* 2 * np.pi  # cavity frequency
wa = 1.0 #* 2 * np.pi  # atom frequency
N = 8  # 15 # number of cavity fock states
# g = 0.05 * 2 * np.pi  # coupling strength
g = 0.1

# Atom hamiltonian
H_atom = 0.5 * qutip.sigmaz()
```

```{code-cell}
def jcm_c_ops(N, n_th_a):
    """ Return basic JCM collapse operators. """
    c_ops = []

    sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())
    a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))

    # cavity relaxation
    rate = kappa * (1 + n_th_a)
    if rate > 0.0:
        c_ops.append(np.sqrt(rate) * a)

    # cavity excitation, if temperature > 0
    rate = kappa * n_th_a
    if rate > 0.0:
        c_ops.append(np.sqrt(rate) * a.dag())

    # qubit relaxation
    rate = gamma
    if rate > 0.0:
        c_ops.append(np.sqrt(rate) * sm)
        
    return c_ops
```

```{code-cell}
def matrix_element(a, x, b):
    """ Return <a|x|b>. """
    return (a * x * b).norm()
```

```{code-cell}
def jcm_c_ops_from_eigenstates(H, N, n_th_a, n_levels=None):
    """ Return full JCM collapse operators. """
    sx = qutip.tensor(qutip.qeye(N), qutip.sigmax())
    a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
    x = a + a.dag()

    energies, eigenstates = H.eigenstates()
    n_energies = len(energies)
    if n_levels is not None:
        n_energies = min(n_energies, n_levels)

    c_ops = []

    for j in range(n_energies):
        for k in range(j, n_energies):
            rate = matrix_element(eigenstates[j].dag(), x, eigenstates[k])**2 * kappa
            if rate > 0.0:
                c_ops.append(np.sqrt(rate) * eigenstates[j] * eigenstates[k].dag())

            rate = matrix_element(eigenstates[j].dag(), sx, eigenstates[k])**2 * gamma
            if rate > 0.0:
                c_ops.append(np.sqrt(rate) * eigenstates[j] * eigenstates[k].dag())

    return c_ops
```

```{code-cell}
def jcm_brmesolve(H, psi0, tlist, kappa, gamma, N, e_ops, options):
    """ Solve the given Jaynes-Cummings model use the Bloch-Redfield solver. """
    sx = qutip.tensor(qutip.qeye(N), qutip.sigmax())
    a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
    x = a + a.dag()

    cavity_spectrum = "0 if (w <= 0) else {kappa}".format(kappa=kappa)
    atom_spectrum = "0 if (w <= 0) else {gamma}".format(gamma=gamma)
    a_ops = [
        [x, cavity_spectrum],
        [sx, atom_spectrum],
    ]

    result = qutip.brmesolve(H, psi0, tlist, a_ops=a_ops, e_ops=e_ops, options=options)
    return result
```

```{code-cell}
g = 0.01  # g_values[0]
H = jcm_h(wc, wa, g, N, H_atom)

display_eigenstates(H)
```

```{code-cell}
g = 0.5

H = jcm_h(wc, wa, g, N, H_atom)

tlist = np.linspace(0, 100, 500)

gnd_energy, gnd_state = H.groundstate()
psi0 = qutip.tensor(qutip.basis(N, 0), qutip.basis(2, 0))
# psi0 = qutip.tensor(qutip.basis(N, 0), qutip.basis(2, 1))
# psi0 = gnd_state

sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())
a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
e_ops = [a.dag() * a, sm.dag() * sm]

c_ops_simple = jcm_c_ops(N, n_th_a)
c_ops_better = jcm_c_ops_from_eigenstates(H, N, n_th_a)
c_ops_better_fewer = jcm_c_ops_from_eigenstates(H, N, n_th_a, n_levels=5)

options = qutip.Options(nsteps=15000, store_states=True, rtol=1e-12, atol=1e-12)

result_br = jcm_brmesolve(H, psi0, tlist, kappa, gamma, N, e_ops=e_ops, options=options)
result_me_simple = qutip.mesolve(H, psi0, tlist, c_ops=c_ops_simple, e_ops=e_ops, options=options)
result_me_better = qutip.mesolve(H, psi0, tlist, c_ops=c_ops_better, e_ops=e_ops, options=options)
result_me_better_fewer = qutip.mesolve(H, psi0, tlist, c_ops=c_ops_better_fewer, e_ops=e_ops, options=options)

plt.plot(tlist, result_me_simple.expect[0], label="ME (simple c_ops)")
plt.plot(tlist, result_me_better.expect[0], label="ME (better c_ops)")
plt.plot(tlist, result_me_better_fewer.expect[0], label="ME (better c_ops fewer)")
plt.plot(tlist, result_br.expect[0], "-.", label="BR")
plt.xlabel("t")
plt.ylabel("Cavity occupation")
plt.legend();
```

```{code-cell}
plt.plot(tlist, result_me_simple.expect[1], label="ME (simple c_ops)")
plt.plot(tlist, result_me_better.expect[1], label="ME (better c_ops)")
plt.plot(tlist, result_me_better_fewer.expect[1], label="ME (better c_ops small)")
plt.plot(tlist, result_br.expect[1], "-.", label="BR")
plt.xlabel("t")
plt.ylabel("Atomic (bare) excited state occupation")
plt.legend();
```

## Explain the expectation values here!

```{code-cell}
a.dag() * a
```

```{code-cell}
sm.dag() * sm
```

```{code-cell}

```
