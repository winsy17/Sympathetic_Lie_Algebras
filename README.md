# Sympathetic Lie Algebras

## Description

Python code for the explicit construction of 
* $\mathfrak{sl}_2({\mathbb{C}})$-simple modules as submodules of a tensor product, 
* $\mathfrak{sl}_2({\mathbb{C}})$-equivariant morphisms of simple modules.

Given a Lie algebra $\mathfrak{g}$ and simple $\mathfrak{sl}_2({\mathbb{C}})$-modules $V_k, V_n, V_m$, 
* ```structure_constants.py``` verifies if candidate structure constants define a Lie algebra and obtains adjoint matrices;
* ```derivations.py``` obtains a basis for the space of derivations $\operatorname{Der}(\mathfrak{g})$;
* ```submodule_generation.py``` determines if $V_k$ is a simple $\mathfrak{sl}_2({\mathbb{C}})$-submodule of $V_n\otimes V_m$ (or $\Lambda^2 V_n$ or $S^2 (V_n))$ and, when it is so, constructs it;
* ```equivariant_morphisms.py``` constructs $\mathfrak{sl}_2(\mathbb{C})$-equivariant morphisms $V_n\otimes V_m\to V_k$, skew-symmetric $\mathfrak{sl}_2(\mathbb{C})$-equivariant morphisms $\Lambda^2(V_n)\to V_k$ and symmetric $\mathfrak{sl}_2(\mathbb{C})$-equivariant morphisms $S^2(V_n)\to V_k$.

Examples of usage are included within each module.


To reference the algorithms developed here, please cite:
> Garcia Pulido A. L; Salgado, G. *On the existence of sympathetic Lie algebras with dimension less than 25.* [arXiv:2212.01273](https://arxiv.org/abs/2212.01273)


## Requirements

Tested in [Python 3.9](https://www.python.org/) with packages:
* [Numpy](http://numpy.scipy.org/)
* [Sympy](https://www.sympy.org/)

which can be PyPI installed:
```
pip install -r requirements.txt
```

## Licensing
This package is 3-clause BSD licensed.