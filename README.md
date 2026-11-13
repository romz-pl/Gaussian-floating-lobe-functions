# Gaussian Radial Basis Functions

## Quantum Chemistry

[Application in the Hartree-Fock Method](./Gaussians-in-Hartree-Fock-Method.md)

## Methods and Theory

[Approximation by Gaussian Functions](./approximation-by-gaussians.md)

[Adaptive Widths and Multiscale Dictionaries](./approximation-rates.md)

[Applications of Adaptive Gaussian Kernels](./applications-adaptive.md)

[Radial Basis Functions for Solving PDEs](./radial-basis-functions-pde.md)

[Multiscale Representation in R^3](./multiscale-representation.md)

[Adaptive second order PDEs](./adaptive-second-order-pde.md)

[Compactly Supported Radial Basis Functions](./rbf-compact-support.md)

## References

[Books and review articles](./references-gaussian-approximation.md)

## Software and libraries

1. **treverhines/RBF**, a Python package for RBF interpolation and RBF-FD and spectral PDE solvers, available on GitHub and PyPI as treverhines-rbf. Its license still needs checking. https://github.com/treverhines/RBF

1. **SciPy's RBFInterpolator**, a Gaussian-kernel interpolation routine that is part of the SciPy library and has no PDE solver. Its license is almost certainly BSD-3-Clause, but I did not see that confirmed in the search results. https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RBFInterpolator.html

1. **SciPy's legacy Rbf class**, which is superseded by RBFInterpolator and is useful mainly for reference. https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.interpolate.Rbf.html

1. **MFREE, a MATLAB toolkit** for moving least squares, radial point interpolation, and moving kriging. It is MIT-licensed according to its paper, but it does not implement Gaussian-kernel multiscale PDE solving. https://github.com/KMountris/MFREE

1. **ferreus_rbf**, a fast global RBF interpolation package with a Rust core and Python bindings, published under the MIT license according to its PyPI page. It handles interpolation, not PDEs. https://pypi.org/project/ferreus-rbf/0.2.0/ (source: github.com/graphic-goose/ferreus_rbf_rs)

1. **RBF-Tools-PDE** by Ramstric, a Python project for Gaussian, multiquadric, and inverse multiquadric collocation for PDEs using PyTorch. Its license could not be confirmed from the search results. https://github.com/Ramstric/RBF-Tools-PDE

1. **Fasshauer's MATLAB companion code** for *Meshfree Approximation Methods with MATLAB*. The code was distributed on a CD with the book, and I found no public repository for it. https://www.mathworks.com/academia/books/meshfree-approximation-methods-with-matlab-fasshauer.html

1. **Fornberg and Flyer** book *A Primer on Radial Basis Functions* (SIAM, 2015), which is a reference rather than software, though it describes the stable evaluation methods that some codes implement. https://rbf.readthedocs.io/en/latest/pde.html

1. **CNF framework by Zhong et al. (2023)** and the Kansa collocation work by Hu et al. (2026). These are papers with Gaussian-kernel Kansa collocation for linear and nonlinear PDEs, and I did not locate a public code repository for them. https://arxiv.org/pdf/2602.07970

