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

In this tutorial we will look at what happens when the coupling between the atom and the cavity becomes very strong. We will introduce a new model for the dissipation that more accurately describes the time evolution and compare with the standard Jaynes-Cummings dissipation model.

+++

## Introduction

If the coupling strength becomes so large that it is comparable to the system frequencies, the collapse operators from the previous tutorial do not correctly model the dissipation (for instance it does not give the ground state as the steady state at zero temperature). This is because the eigenenergies of the system are significantly altered by the coupling between the atom and the cavity.

We need to rederive the master equation, defining the system-environment coupling operator in terms of the true system eigenstates, in do so we arrive at a master equation of this form


$$\frac{\partial}{\partial t} \rho_S(t) = -i[H_S,\rho(t)] \\
+ \sum_{\Gamma=\gamma,\kappa}\left(\sum_{j,k>j}\Gamma_{j,k}n(\Delta_{j,k},T) \mathcal{D}\left[|k\rangle\langle j|\right] \rho(t) \right.\\
\left. + \sum_{j,k>j}\Gamma_{j,k}\left[1+n(\Delta_{j,k},T)\right] \mathcal{D}\left[|j\rangle\langle k|\right] \rho(t)\right)$$

Where 

$\gamma_{j,k}= \pi J(\Delta_{j,k}) |\langle j|\sigma_x | k \rangle |^2$, 

$\kappa_{j,k}= \pi J(\Delta_{j,k}) |\langle j|(a+a^{\dagger}) | k \rangle |^2$,

$\mathcal{D}\left[|j\rangle\langle k|\right]\rho(t)=2|j\rangle\langle k| \rho |k\rangle\langle j| - |k\rangle\langle k|\rho - \rho |k\rangle\langle k|$, 

$\Delta_{j,k} = E_j - E_k$.

We refer this as the ultra-strong coupling (USC) master equation

We can think of this as the collapse operators coupling to all the different possible energy jumps of the combined atom-cavity system. Note that in practice, some energy states are never occupied and so we can in fact ignore these collapse operators, as we will see later.

+++

## Tasks

Constructing the model described above is quite involved, and so that has been completed already. Our tasks will focus on analysing the output.  

- Try different values for the coupling strength ($0.01 < g < 0.9$) and observe how the eigenenergies change

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

## Helper functions

```{code-cell} ipython3
def jcm_h(wc, wa, g, N, atom):
    """ Construct the Jaynes-Cummings Hamiltonian (non-RWA). """
    a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
    sm = qutip.tensor(qutip.qeye(N), qutip.destroy(2))
    atom = qutip.tensor(qutip.qeye(N), atom)
    
    H = wc * a.dag() * a + wa * atom + g * (a.dag() + a) * (sm + sm.dag())
    return H
```

```{code-cell} ipython3
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

## Create the JCM Hamiltonian and USC collapse operators.

```{code-cell} ipython3
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

```{code-cell} ipython3
def matrix_element(a, x, b):
    """ Return <a|x|b>. """
    return (a * x * b).norm()
```

```{code-cell} ipython3
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

```{code-cell} ipython3
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

```{code-cell} ipython3
# Dissipation parameters
# We use stronger dissipation to show long-term behaviour in shorter times
kappa = 0.2 # cavity dissipation rate
gamma = 0.2 # atom dissipation rate
n_th_a = 0 # avg number of thermal bath excitation
```

```{code-cell} ipython3
# system parameters
wc = 1.0 #* 2 * np.pi  # cavity frequency
wa = 1.0 #* 2 * np.pi  # atom frequency
N = 8  # 15 # number of cavity fock states
# g = 0.05 * 2 * np.pi  # coupling strength
g = 0.1

# Atom hamiltonian
H_atom = 0.5 * qutip.sigmaz()
```

```{code-cell} ipython3
g = 0.01
H = jcm_h(wc, wa, g, N, H_atom)

display_eigenstates(H)
```

## Solving and visualising the dynamics

+++

The dynamics are solved here for the simple (from JCM 2 tutorial) and the USC dissipation model

Try changing the coupling strength and see how the two models diverge.
Try changing the initial state and see what happens.

```{code-cell} ipython3
g = 0.001

# get the JCM Hamiltonian
H = jcm_h(wc, wa, g, N, H_atom)

# set the timespace for the simulation
tlist = np.linspace(0, 100, 500)

# get the groundstate for the full system
gnd_energy, gnd_state = H.groundstate()

# set the initial state 
psi0 = qutip.tensor(qutip.basis(N, 0), qutip.basis(2, 0)) # atom excited
# psi0 = qutip.tensor(qutip.basis(N, 0), qutip.basis(2, 1)) # atom ground state
# psi0 = gnd_state

# create the expectation operators
sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())
a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
e_ops = [a.dag() * a, sm.dag() * sm] # cavity and atom excited state probabilities

# generate the collapse operators
c_ops_simple = jcm_c_ops(N, n_th_a) # basic JCM
c_ops_usc = jcm_c_ops_from_eigenstates(H, N, n_th_a) # USC 

options = qutip.Options(nsteps=15000, store_states=True, rtol=1e-12, atol=1e-12)

result_me_simple = qutip.mesolve(H, psi0, tlist, c_ops=c_ops_simple, e_ops=e_ops, options=options)
result_me_usc = qutip.mesolve(H, psi0, tlist, c_ops=c_ops_usc, e_ops=e_ops, options=options)

# set the timespace for the plot
tlist_plot = np.linspace(0, 100, 500)

plt.plot(tlist_plot, result_me_simple.expect[0], label="ME (simple)")
plt.plot(tlist_plot, result_me_usc.expect[0], label="ME (USC)")
plt.xlabel("t")
plt.ylabel("Cavity occupation")
plt.legend();
```

## Reduced level couplings (for collapse operators)

Now we investigate how many energy levels we need to consider in order for an accurate model.

You may wish to add a second plot that zooms in on the later times (t > 80) and exclude the simple model plot, to see more clearly how the USC and USC limited level plots converge.

```{code-cell} ipython3
# start by copying the previous cell
```

```{code-cell} ipython3
g = 0.5

# get the JCM Hamiltonian
H = jcm_h(wc, wa, g, N, H_atom)

# set the timespace for the simulation
tlist = np.linspace(0, 100, 500)

# get the groundstate for the full system
gnd_energy, gnd_state = H.groundstate()

# set the initial state 
psi0 = qutip.tensor(qutip.basis(N, 0), qutip.basis(2, 0)) # atom excited
# psi0 = qutip.tensor(qutip.basis(N, 0), qutip.basis(2, 1)) # atom ground state
# psi0 = gnd_state

# create the expectation operators
sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())
a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
e_ops = [a.dag() * a, sm.dag() * sm] # cavity and atom excited state probabilities

# generate the collapse operators
c_ops_simple = jcm_c_ops(N, n_th_a) # basic JCM
c_ops_usc = jcm_c_ops_from_eigenstates(H, N, n_th_a) # USC 
c_ops_usc_levlim = jcm_c_ops_from_eigenstates(H, N, n_th_a, n_levels=5) # USC with level limit

options = qutip.Options(nsteps=15000, store_states=True, rtol=1e-12, atol=1e-12)

result_me_simple = qutip.mesolve(H, psi0, tlist, c_ops=c_ops_simple, e_ops=e_ops, options=options)
result_me_usc = qutip.mesolve(H, psi0, tlist, c_ops=c_ops_usc, e_ops=e_ops, options=options)
result_me_usc_levlim = qutip.mesolve(H, psi0, tlist, c_ops=c_ops_usc_levlim, e_ops=e_ops, options=options)


# set the timespace for the plot
tlist_plot = np.linspace(0, 100, 500)

plt.plot(tlist_plot, result_me_simple.expect[0], label="ME (simple)")
plt.plot(tlist_plot, result_me_usc.expect[0], label="ME (USC)")
plt.plot(tlist_plot, result_me_usc_levlim.expect[0], label="ME (USC limited levels)")
plt.legend();

plt.figure()

tlist_plot = np.linspace(80, 100, 100)
#plt.plot(tlist_plot, result_me_simple.expect[0], label="ME (simple)")
plt.plot(tlist_plot, result_me_usc.expect[0][400:], label="ME (USC)")
plt.plot(tlist_plot, result_me_usc_levlim.expect[0][400:], label="ME (USC limited levels)")

plt.xlabel("t")
plt.ylabel("Cavity occupation")
plt.legend();
```

## Bloch-Redfield model

Our USC dissipation is actually the Bloch-Redfield model. QuTiP has a solver built in for this.

Simulate the dynamics with the Bloch-Redfield (using the helper function 'jcm_brmesolve')

```{code-cell} ipython3
# start by copying the previous cell
```

```{code-cell} ipython3
g = 0.5

# get the JCM Hamiltonian
H = jcm_h(wc, wa, g, N, H_atom)

# set the timespace for the simulation
tlist = np.linspace(0, 100, 500)

# get the groundstate for the full system
gnd_energy, gnd_state = H.groundstate()

# set the initial state 
psi0 = qutip.tensor(qutip.basis(N, 0), qutip.basis(2, 0)) # atom excited
# psi0 = qutip.tensor(qutip.basis(N, 0), qutip.basis(2, 1)) # atom ground state
# psi0 = gnd_state

# create the expectation operators
sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())
a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
e_ops = [a.dag() * a, sm.dag() * sm] # cavity and atom excited state probabilities

# generate the collapse operators
c_ops_simple = jcm_c_ops(N, n_th_a) # basic JCM
c_ops_usc = jcm_c_ops_from_eigenstates(H, N, n_th_a) # USC 
c_ops_usc_levlim = jcm_c_ops_from_eigenstates(H, N, n_th_a, n_levels=5) # USC with level limit

options = qutip.Options(nsteps=15000, store_states=True, rtol=1e-12, atol=1e-12)

result_me_simple = qutip.mesolve(H, psi0, tlist, c_ops=c_ops_simple, e_ops=e_ops, options=options)
result_me_usc = qutip.mesolve(H, psi0, tlist, c_ops=c_ops_usc, e_ops=e_ops, options=options)
result_me_usc_levlim = qutip.mesolve(H, psi0, tlist, c_ops=c_ops_usc_levlim, e_ops=e_ops, options=options)
result_br = jcm_brmesolve(H, psi0, tlist, kappa, gamma, N, e_ops=e_ops, options=options)


# set the timespace for the plot
tlist_plot = np.linspace(0, 100, 500)

plt.plot(tlist_plot, result_me_simple.expect[0], label="ME (simple)")
plt.plot(tlist_plot, result_me_usc.expect[0], label="ME (USC)")
plt.plot(tlist_plot, result_me_usc_levlim.expect[0], label="ME (USC limited levels)")
plt.plot(tlist_plot, result_br.expect[0], "--", label="BR")
plt.legend();

plt.figure()

tlist_plot = np.linspace(80, 100, 100)
#plt.plot(tlist_plot, result_me_simple.expect[0], label="ME (simple)")
plt.plot(tlist_plot, result_me_usc.expect[0][400:], label="ME (USC)")
plt.plot(tlist_plot, result_me_usc_levlim.expect[0][400:], label="ME (USC limited levels)")
plt.plot(tlist_plot, result_br.expect[0][400:], "--", label="BR")

plt.xlabel("t")
plt.ylabel("Cavity occupation")
plt.legend();
```

```{code-cell} ipython3
plt.plot(tlist, result_me_simple.expect[1], label="ME (simple c_ops)")
plt.plot(tlist, result_me_usc.expect[1], label="ME (better c_ops)")
plt.plot(tlist, result_me_usc_levlim.expect[1], label="ME (better c_ops small)")
plt.plot(tlist, result_br.expect[1], "-.", label="BR")
plt.xlabel("t")
plt.ylabel("Atomic (bare) excited state occupation")
plt.legend();
```
