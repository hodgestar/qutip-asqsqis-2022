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

# 02: Quantum Information Processing in QuTiP II

+++

In the previous tutorial, we looked at how to create and run circuits. In this tutorial we'll look at how circuits are implemented on real quantum devices and how we can use QuTiP to simulate the *quantum hardware* itself.

+++

## Tasks

- [Import qutip and qutip.qip](#imports)

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

def display_list(a, title=None):
    """ Display a list of items using Jupyter's display function. """
    if title:
        print(title)
        print("=" * len(title))
    for item in a:
        display(item)
```

```{code-cell} ipython3
:tags: [hide-cell]

def print_attrs(obj, name, attrs):
    """ Print the listed attributes of an object. """
    for attr in attrs:
        print(f"{name}.{attr}: {getattr(obj, attr)}")
```

```{code-cell} ipython3
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

## Task 1

- Sub-task 1

```{code-cell} ipython3
:tags: [hide-cell]

# sub-task 1
```

## Using functions to make your notebooks neater

- XXX

```{code-cell} ipython3
:tags: [hide-cell]

# See the helper functions at the top for examples.
```
