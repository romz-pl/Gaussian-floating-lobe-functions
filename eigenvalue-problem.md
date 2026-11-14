# Radial Basis Functions for Eigenvalue Problems in Three Dimensions

## 1. Introduction

Eigenvalue problems posed on three-dimensional domains arise throughout computational physics, engineering and applied mathematics. Typical examples include the vibration modes of elastic solids, acoustic resonances in cavities, electromagnetic cavity modes, the stability analysis of spatially extended systems, and the stationary states of quantum mechanical particles confined to a region. In each case one seeks a scalar or vector field $u$ and a number $\lambda$ satisfying a linear operator equation

$$
\mathcal{L} u = \lambda\, \mathcal{B} u \quad \text{in } \Omega \subset \mathbb{R}^3,
$$

together with boundary conditions on $\partial\Omega$. The operator $\mathcal{L}$ is typically a second-order differential operator and $\mathcal{B}$ is either the identity (standard eigenproblem) or a weight/mass operator (generalised eigenproblem).

Radial basis function (RBF) methods offer a meshfree route to such problems. Rather than triangulating the domain into elements, one scatters nodes through $\Omega$ and approximates the unknown field as a linear combination of radially symmetric functions centred at those nodes. In three dimensions this has several attractions: complex geometries are handled by node placement alone, the approximation can be spectrally accurate for smooth solutions, and the same machinery applies to scalar Laplacian-type problems, Helmholtz-type problems, and Schrödinger-type problems.

This document surveys the mathematical formulation, the main discretisation strategies (collocation and Galerkin), the specific structure of 3D problems, the treatment of boundary conditions, the sources of spurious modes, numerical conditioning, computational cost and remedies, and worked examples with known analytic answers.

## 2. Model Problems

Several model problems illustrate the range of operators encountered.

**Laplacian (Dirichlet) problem.** Find $\lambda$ and $u \neq 0$ such that

$$
-\nabla^2 u = \lambda u \ \text{in } \Omega, \qquad u = 0 \ \text{on } \partial\Omega.
$$

**Laplacian (Neumann) problem.** Replace the Dirichlet condition with $\partial u/\partial n = 0$ on $\partial\Omega$. The constant function is then an eigenfunction with $\lambda = 0$, a fact that must be respected by any discretisation.

**Helmholtz-type problem.** Operators of the form $-\nabla^2 u + \kappa(x) u = \lambda u$ arise in waveguide and cavity problems, with $\kappa$ a spatially varying coefficient.

**Schrödinger problem.** For a particle in a potential $V(x)$,

$$
-\frac{\hbar^2}{2m}\nabla^2 \psi + V(x)\psi = E\,\psi,
$$

where in atomic units $\hbar = m = 1$ this reduces to $-\tfrac{1}{2}\nabla^2\psi + V\psi = E\psi$. The three-dimensional harmonic oscillator, with $V = \tfrac{1}{2}\omega^2 r^2$, has exact eigenvalues $E_{n} = \omega\,(n + \tfrac{3}{2})$ with degeneracies determined by the quantum numbers, which makes it a convenient test.

**Elasticity and vibration.** Coupled vector systems of the form $\mathbf{K}\mathbf{u} = \lambda\mathbf{M}\mathbf{u}$ arise from the Navier equations. RBF methods extend to these by using vector-valued interpolants, typically one scalar RBF expansion per displacement component.

## 3. Radial Basis Function Approximation

An RBF approximant takes the form

$$
u(x) \approx \sum_{j=1}^{N} c_j\, \phi\big(\|x - x_j\|\big) + p(x),
$$

where $\{x_j\}_{j=1}^{N} \subset \mathbb{R}^3$ are centres, $\phi:[0,\infty)\to\mathbb{R}$ is a radial kernel, $\|\cdot\|$ is the Euclidean norm, and $p$ is an optional low-degree polynomial whose role is discussed in Section 7. Common kernels include:

- **Gaussian:** $\phi(r) = e^{-\varepsilon^2 r^2}$, strictly positive definite in every dimension, infinitely smooth, with shape parameter $\varepsilon > 0$.
- **Multiquadric:** $\phi(r) = \sqrt{1 + \varepsilon^2 r^2}$, conditionally positive definite of order one.
- **Inverse multiquadric:** $\phi(r) = (1 + \varepsilon^2 r^2)^{-1/2}$, positive definite.
- **Polyharmonic splines:** $\phi(r) = r^{k}$ for odd $k$, or $r^{k}\log r$ for even $k$. In three dimensions $r^3$ and $r^5$ are the standard choices, usually augmented with polynomials.
- **Wendland functions:** compactly supported, positive definite functions that yield sparse interpolation matrices.

The choice of kernel controls smoothness, the existence of a unique interpolant, and the conditioning of the resulting linear systems. Smooth kernels such as the Gaussian and multiquadric give spectral (super-algebraic) convergence for sufficiently smooth solutions, while polyharmonic splines give algebraic convergence but are free of the shape-parameter sensitivity that complicates the smooth kernels.

### 3.1 Derivatives of radial kernels in three dimensions

Because the kernel depends only on $r = \|x - x_j\|$, its Laplacian in $\mathbb{R}^3$ takes a simple form. For a radial function $\phi(r)$,

$$
\nabla^2 \phi = \phi''(r) + \frac{2}{r}\phi'(r).
$$

For the Gaussian $\phi = e^{-\varepsilon^2 r^2}$ one obtains $\phi' = -2\varepsilon^2 r\,\phi$ and $\phi'' = (4\varepsilon^4 r^2 - 2\varepsilon^2)\phi$, so

$$
\nabla^2 \phi = \big(4\varepsilon^4 r^2 - 6\varepsilon^2\big)\,\phi.
$$

This closed form is what makes collocation with Gaussians particularly convenient to code. The polyharmonic kernel $r^3$ has $\nabla^2 r^3 = 12 r$, which is also elementary. Expressions of this kind should be verified symbolically or numerically before use, since sign and factor errors are common and difficult to spot in eigenvalue output.

## 4. Discretisation Strategies

### 4.1 Unsymmetric collocation (Kansa method)

Kansa's method enforces the differential equation and boundary conditions pointwise at a set of nodes. Partition the nodes into $N_I$ interior points $\{x_i\}$ and $N_B$ boundary points $\{x_b\}$, with $N = N_I + N_B$ total centres. Substituting the RBF expansion into the eigenproblem and requiring equality at the interior nodes gives

$$
\sum_{j=1}^{N} c_j\, (-\nabla^2\phi_j)(x_i) = \lambda \sum_{j=1}^{N} c_j\, \phi_j(x_i), \qquad i = 1,\dots,N_I,
$$

while the boundary condition, here Dirichlet, gives

$$
\sum_{j=1}^{N} c_j\, \phi_j(x_b) = 0, \qquad b = 1,\dots,N_B.
$$

In matrix form,

$$
\mathbf{A}\mathbf{c} = \lambda\, \mathbf{B}\mathbf{c},
$$

where $\mathbf{A} _ {ij} = -(\nabla^2\phi_j)(x_i)$ for interior rows and $\mathbf{A} _ {bj} = \phi_j(x_b)$ for boundary rows, and $\mathbf{B}_{ij} = \phi_j(x_i)$ for interior rows with zeros in the boundary rows. The matrix $\mathbf{B}$ is therefore rectangular-in-structure but square overall and singular, because its boundary rows vanish. This singularity produces infinite eigenvalues, which must be discarded.

The resulting matrices are nonsymmetric and dense. The method is simple to implement and works for arbitrary geometry and arbitrary boundary conditions, but it carries no variational guarantee: the eigenvalues may be complex, and spurious modes are common (Section 6).

### 4.2 Galerkin (Rayleigh-Ritz) formulation

The weak form of the Dirichlet Laplacian eigenproblem seeks $u \in H^1_0(\Omega)$ with

$$
\int_\Omega \nabla u \cdot \nabla v \, dx = \lambda \int_\Omega u\, v\, dx \qquad \text{for all } v \in H^1_0(\Omega).
$$

Taking $u_h = \sum_j c_j \phi_j$ with $\phi_j$ chosen to satisfy the essential boundary condition (or enforced via penalty or Lagrange multipliers) produces the symmetric generalised eigenproblem

$$
\mathbf{K}\mathbf{c} = \lambda\, \mathbf{M}\mathbf{c}, \qquad
\mathbf{K}_{ij} = \int_\Omega \nabla\phi_i \cdot \nabla\phi_j\, dx, \quad
\mathbf{M}_{ij} = \int_\Omega \phi_i \phi_j\, dx.
$$

The Galerkin approach has real, positive, and by the min-max principle correctly bounded eigenvalues. Its disadvantages are the need for accurate three-dimensional quadrature over $\Omega$, which is expensive for dense RBF matrices, and the difficulty of enforcing boundary conditions on general domains. Background cells, Monte Carlo integration, or sparse Gauss rules on a tetrahedral background mesh are used in practice. Galerkin RBF methods are closer in spirit to the finite element method but retain the meshfree character of the approximation.

### 4.3 RBF-FD (local) methods

Dense global matrices scale poorly in three dimensions. RBF-based finite differences (RBF-FD) replace the global approximation with local stencils: for each node $x_i$ one selects its $n$ nearest neighbours, builds a local RBF interpolant over that stencil (often augmented with monomials), and derives weights approximating $\nabla^2 u(x_i)$ as a linear combination of stencil values. The result is a sparse matrix with $O(n)$ nonzeros per row. The eigenproblem is then solved with sparse eigensolvers. RBF-FD sacrifices global spectral accuracy for locality and has algebraic convergence governed by stencil size and node density, but in practice it is the most scalable option for three-dimensional problems with many nodes.

## 5. Geometry, Nodes and Boundary Treatment

### 5.1 Node distributions

In three dimensions the cost of resolving a solution grows rapidly with the number of nodes, so node placement matters. Options include:

- **Regular grids** inside boxes, simple but poor for curved boundaries.
- **Quasi-random sequences** (Halton, Sobol) that fill the volume evenly while avoiding lattice artifacts.
- **Poisson-disk or advancing-front sampling** that enforces a minimum spacing and adapts to a target density.
- **Surface sampling** of the boundary using spherical Fibonacci lattices, latitude-longitude grids, or parametrisations of CAD surfaces, combined with interior sampling by repulsion or relaxation.
- **Adaptive refinement** where local node spacing follows a prescribed function $h(x)$, important when the eigenfunction has localised features such as boundary layers or potential wells.

A useful rule of thumb is that the characteristic spacing $h$ should resolve the shortest wavelength of interest, with several nodes per wavelength. For eigenfunctions of order $k$ on a domain of diameter $D$, this demands roughly $N \sim (D/h)^3$ nodes, which makes the cubic scaling of dense RBF matrices a serious constraint.

### 5.2 Boundary conditions

Dirichlet conditions are enforced by collocation at boundary nodes, as above. Neumann and Robin conditions require evaluating the normal derivative $\partial\phi_j/\partial n$ at boundary nodes, which for radial kernels is simply $\phi'(r)\,(x_b - x_j)\cdot n / r$. Corners and edges of polyhedral domains require particular care: the normal is discontinuous, so either the boundary condition is imposed on each face separately with duplicated nodes at edges, or a smoothed normal is used. Curved boundaries are handled naturally by placing nodes on the surface, but the geometric accuracy of the boundary sampling then bounds the attainable accuracy of the eigenvalues.

## 6. Spurious and Infinite Eigenvalues

Collocation eigenproblems exhibit behaviour that has no counterpart in the Galerkin setting, and any practical implementation must handle it.

**Infinite eigenvalues.** The rank deficiency of $\mathbf{B}$ introduced by zeroing the boundary rows produces infinite eigenvalues. Solvers such as LAPACK's generalised routines report these explicitly, and they must be filtered out by testing for finiteness.

**Complex eigenvalues.** Since $\mathbf{A}$ and $\mathbf{B}$ are nonsymmetric, the computed eigenvalues may be complex even when the underlying operator is self-adjoint. Their imaginary parts are usually a diagnostic of inaccuracy, and physically meaningful eigenvalues are those with negligible imaginary part.

**Spurious modes.** Collocation can yield eigenpairs that do not correspond to any true eigenfunction of the continuous problem. Typical signatures are eigenvectors concentrated near the boundary, or eigenvalues that are insensitive to refinement. A common remedy is to verify each candidate by refining the node set and retaining only eigenvalues that converge, or by checking the residual of the continuous equation at independent test points that were not used in the collocation.

**Sensitivity to the shape parameter.** For Gaussian kernels, the computed spectrum can change substantially with $\varepsilon$ even when the approximation is accurate for each $\varepsilon$. Eigenvalue stability with respect to $\varepsilon$ is therefore a useful check: a genuine eigenvalue should plateau over a range of shape parameters, while spurious ones drift.

## 7. Conditioning, Polynomial Augmentation and the Flat Limit

The Gaussian and related smooth kernels generate interpolation matrices whose condition number grows rapidly with node density and with decreasing $\varepsilon$. This is the regime of the so-called flat limit, in which the kernel becomes nearly constant over the node spacing and the basis becomes nearly linearly dependent. Schaback's uncertainty principle expresses the underlying trade-off: attainable accuracy and conditioning cannot both be improved indefinitely at fixed precision. Practical strategies include:

- Choosing $\varepsilon$ scaled to the local node spacing, for example $\varepsilon \propto 1/h$, so that the kernel width tracks resolution.
- Using contour-Padé, RBF-QR, or eigenfunction-expansion techniques to stabilise computations in the flat regime. These allow stable evaluation at small $\varepsilon$ without ill-conditioned linear solves.
- Working in extended precision for validation runs, which exposes whether the spectrum is limited by discretisation error or by round-off.

For polyharmonic kernels in three dimensions, conditional positive definiteness requires augmenting the expansion with polynomials of degree at least $m$, typically $m = 2$ for $r^3$ and $m = 3$ for $r^5$ in $\mathbb{R}^3$, under the side conditions $\sum_j c_j\, q(x_j) = 0$ for all polynomials $q$ of degree up to $m$. These constraints are appended to the linear system, and they also guarantee that polynomial solutions are reproduced exactly, a valuable property for eigenvalue problems whose eigenfunctions are polynomially structured near the boundary.

## 8. Computational Cost

Dense collocation or Galerkin matrices of size $N \times N$ require $O(N^2)$ storage and $O(N^3)$ operations for dense eigensolvers. For a three-dimensional problem resolved with $n$ nodes per direction, $N = n^3$, so the costs scale as $n^6$ and $n^9$ respectively. At $n = 30$, which is already modest for many physical problems, $N = 27{,}000$ and a dense generalised eigensolve becomes expensive on a workstation.

Several mitigations are standard:

- Compute only the eigenvalues of interest using shift-invert Arnoldi or Lanczos iterations, with a shift placed near the target spectrum. For symmetric Galerkin problems the Lanczos method is particularly efficient.
- Switch to RBF-FD to obtain sparse matrices and use sparse shift-invert factorisations.
- Use hierarchical or fast multipole techniques to accelerate matrix-vector products in dense global formulations.
- Exploit symmetry of the domain, for example by reducing a cube or sphere problem to a symmetry sector, which block-diagonalises the operator.
- Parallelise assembly, which is embarrassingly parallel over nodes, and use distributed eigensolvers such as those in SLEPc for large sparse problems.

## 9. Error Behaviour

For smooth eigenfunctions and the Gaussian kernel with appropriately chosen $\varepsilon$, global collocation and Galerkin RBF methods typically exhibit spectral convergence: the error in the eigenvalue decreases faster than any power of the node spacing until round-off or conditioning limits are reached. Eigenvalue errors are generally smaller than eigenfunction errors, consistent with standard variational arguments for the Galerkin case. For RBF-FD, convergence is algebraic and governed by the stencil order, so eigenvalues improve roughly as $O(h^{p})$ for a method of formal order $p$.

Error estimates should be checked empirically by refinement studies. A standard procedure is to compute the eigenvalues on a sequence of node sets with spacings $h_1 > h_2 > h_3$, estimate the observed order from the ratio of successive differences, and compare against known or extrapolated values. Discrepancies between the observed and expected order usually reveal boundary treatment errors, insufficient node density in regions of high curvature in the eigenfunction, or inappropriate shape-parameter scaling.

## 10. Worked Examples and Verification

Verification against exact results is essential. Several three-dimensional examples are particularly useful.

**Unit cube, Dirichlet Laplacian.** The eigenvalues are $\lambda_{pqr} = \pi^2(p^2 + q^2 + r^2)$ for positive integers $p, q, r$. The lowest is $\lambda_{111} = 3\pi^2 \approx 29.61$, followed by the threefold degenerate $\lambda_{211} = 6\pi^2 \approx 59.22$. The degeneracies make this a stringent test of symmetry handling.

**Unit ball, Dirichlet Laplacian.** The lowest eigenfunction is radial, $u = \sin(kr)/r$, which must vanish at $r = 1$, so $k = \pi$ and $\lambda = \pi^2 \approx 9.87$. Higher radial modes correspond to successive zeros of $\sin k$, and modes with angular dependence are associated with zeros of spherical Bessel functions. The ball is an excellent test of curved boundary handling, since node placement on the sphere directly affects accuracy.

**Three-dimensional harmonic oscillator.** With $V = \tfrac{1}{2}r^2$ in atomic units, $E_{n} = n + \tfrac{3}{2}$ with $n = 2n_r + \ell$ and $\ell$ the angular momentum. The ground state is $E = 1.5$, the first excited level is $E = 2.5$, threefold degenerate. Because the domain is unbounded, the test requires a sufficiently large truncation box or Gaussian-type decay of the basis, and the eigenfunctions are Gaussian in $r$, which Gaussian RBFs represent efficiently.

**Hydrogen-like problems.** The Coulomb singularity at the origin is a demanding test for RBF methods, which must resolve the cusp of the wavefunction. Node refinement near the nucleus, or a change of variables, is typically required to obtain accurate low-lying energies.

Each test should be run with several node counts, several shape parameters, and both collocation and Galerkin variants where feasible. Agreement with exact values to the expected order, together with stability of the computed spectrum across parameters, is the strongest available evidence of correctness.

## 11. Implementation Sketch

The following Python fragment illustrates the unsymmetric collocation approach for the Dirichlet Laplacian on a set of scattered nodes with Gaussian kernels. It is intended to convey structure, not to serve as production code. It omits polynomial augmentation, node generation, and filtering of spurious modes beyond finiteness.

```python
import numpy as np
from scipy.linalg import eig

def rbf_dirichlet_eigs(X_int, X_bnd, eps):
    """
    Unsymmetric collocation for -Laplacian u = lambda u with u = 0 on the boundary.
    X_int : (n_int, 3) interior nodes
    X_bnd : (n_bnd, 3) boundary nodes
    eps   : Gaussian shape parameter
    """
    X = np.vstack([X_int, X_bnd])          # all centres
    n_int = X_int.shape[0]

    # Pairwise distances between all evaluation points and all centres
    D2 = np.sum((X[:, None, :] - X[None, :, :])**2, axis=-1)
    Phi = np.exp(-eps**2 * D2)
    # Laplacian of Gaussian in 3D: (4 eps^4 r^2 - 6 eps^2) phi
    Lap = (4 * eps**4 * D2 - 6 * eps**2) * Phi

    # Interior rows: -Lap u = lambda * Phi c
    A = -Lap.copy()
    B = Phi.copy()
    B[n_int:, :] = 0.0                     # boundary rows contribute no eigen-weight

    # Boundary rows: u = 0 on the boundary
    A[n_int:, :] = Phi[n_int:, :]

    w, V = eig(A, B)
    finite = np.isfinite(w)
    w, V = w[finite], V[:, finite]
    order = np.argsort(w.real)
    return w[order], V[:, order]
```

For a unit ball, one would generate interior nodes by rejection sampling inside the sphere and boundary nodes on its surface, then compare the lowest computed eigenvalues to $\pi^2$. Refinement in $N$ and scanning in $\varepsilon$ should be performed before trusting any eigenvalue, and complex or unstable values should be discarded.

## 12. Extensions and Related Problems

**Nonlinear and parameter-dependent eigenproblems.** When $\mathcal{L}$ or $\mathcal{B}$ depends on $\lambda$, or when the operator is a polynomial or rational function of the eigenvalue, the discrete problem becomes a nonlinear eigenvalue problem. Linearisation techniques or Newton-type iterations on the characteristic determinant are used.

**Vector and coupled systems.** Elasticity and electromagnetic problems give block systems. Vector RBF expansions, with one component per field and coupling blocks assembled from the same node set, reduce the problem to a generalised eigenproblem of block size $3N$. Divergence constraints, for incompressible or divergence-free fields, require projection or augmented formulations.

**Time-dependent problems.** Eigenvalues of RBF-discretised spatial operators determine the stability and modal decomposition of time-dependent PDEs. Method-of-lines discretisations, combined with eigenanalysis of the spatial operator, give insight into stiffness and admissible time-step sizes.

**Inverse and identification problems.** Given measured eigenfrequencies, one may seek domain geometry or material coefficients. The meshfree character of RBFs makes it straightforward to vary geometry within an optimisation loop, since only node positions change between iterations.

## 13. Comparison with Other Methods

Relative to finite elements, RBF methods avoid mesh generation, handle geometric change gracefully, and can achieve higher accuracy per degree of freedom for smooth solutions. They lack the sparsity and mature error theory of FEM in their global form, and their dense matrices limit scalability unless local variants are adopted. Relative to spectral methods on tensor-product grids, RBFs are not restricted to rectangular or simple domains, but they generally require more memory and careful parameter selection. RBF-FD bridges the gap, offering sparsity and geometric flexibility at the price of lower formal order.

Relative to boundary element methods, RBF volume methods avoid the need for fundamental solutions of the operator, which may not be available for general coefficients, but BEM reduces dimensionality and can be preferable for exterior problems.

## 14. Open Challenges

Several issues remain active areas of research. Automatic and robust selection of shape parameters across scales is not fully resolved. Certified error bounds for collocation eigenvalues are less developed than those for Galerkin FEM. Efficient and reliable spurious-mode filtering without reference solutions remains difficult. Scalable solvers for large three-dimensional RBF systems, particularly in the flat regime, are an ongoing topic. Finally, extending these methods to highly anisotropic, multiscale, or singular problems, such as those with cusps or corner singularities in three dimensions, demands adaptive node strategies and enrichment that are not yet routine.

## 15. Summary

Radial basis functions provide a flexible meshfree framework for three-dimensional eigenvalue problems. Global collocation is simple and accurate for smooth solutions but produces dense nonsymmetric matrices, spurious modes and infinite eigenvalues that require careful filtering. Galerkin formulations give symmetric, variationally sound problems at the cost of demanding three-dimensional quadrature. RBF-FD offers sparsity and scalability with algebraic accuracy. Success in practice depends on node distribution, shape-parameter scaling, polynomial augmentation, and rigorous verification against analytic benchmarks such as the unit cube, the unit ball, and the three-dimensional harmonic oscillator.

## References for Further Reading

- E. J. Kansa, "Multiquadrics: a scattered data approximation scheme with applications to computational fluid dynamics," Computers and Mathematics with Applications, 1990.
- G. E. Fasshauer, *Meshfree Approximation Methods with MATLAB*, World Scientific, 2007.
- B. Fornberg and N. Flyer, *A Primer on Radial Basis Functions with Applications to the Geosciences*, SIAM, 2015.
- R. Schaback, "Error estimates and condition numbers for radial basis function interpolation," Advances in Computational Mathematics, 1995.
- G. B. Wright and B. Fornberg, "Stable computations with flat radial basis functions using vector-valued rational approximations," Journal of Computational Physics, 2017.
- V. Shankar and G. Wright on RBF-FD and stencil-based methods for elliptic eigenproblems and related PDE applications in the recent literature.
