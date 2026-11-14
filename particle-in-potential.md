# Radial Basis Function Methods for the Schrödinger Eigenvalue Problem

## 1. The problem

For a single particle of mass $m$ in a potential $V(\mathbf{x})$, the time-independent Schrödinger equation is

$$
\left[-\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf{x})\right]\psi(\mathbf{x}) = E\,\psi(\mathbf{x}),
$$

where $\psi$ is the wavefunction and $E$ the energy. Working in units where $\hbar = m = 1$ simplifies the notation to

$$
\left[-\tfrac{1}{2}\nabla^2 + V(\mathbf{x})\right]\psi = E\psi.
$$

We seek the discrete eigenvalues $E_n$ and eigenfunctions $\psi_n$ that are square-integrable over the domain $\Omega$ and satisfy the boundary conditions (usually $\psi = 0$ on $\partial\Omega$ for a confined system, or decay at infinity for a bound state). The difficulty is not conceptual but computational: for a bound state in $d$ dimensions, the eigenfunction must be resolved over a region whose size is set by the potential, and the eigenvalues must be accurate enough to be physically meaningful, often to many digits. Grid-based methods (finite differences, finite elements, spectral methods on tensor-product grids) handle one or two dimensions well but become expensive or awkward in higher dimensions or on irregular domains. Radial basis function (RBF) methods offer a meshfree alternative: the wavefunction is expanded in functions centred on a scattered set of points, and the eigenproblem is reduced to a matrix problem without any mesh connectivity.

## 2. Radial basis functions

A radial basis function is a function whose value depends only on the Euclidean distance $r = \|\mathbf{x} - \mathbf{x}_j\|$ from a centre $\mathbf{x}_j$. The most commonly used choices are:

- Gaussian: $\phi(r) = e^{-\varepsilon^2 r^2}$
- Multiquadric: $\phi(r) = \sqrt{r^2 + \varepsilon^2}$, and its inverse $1/\sqrt{r^2+\varepsilon^2}$
- Inverse quadratic: $\phi(r) = 1/(1 + \varepsilon^2 r^2)$
- Polyharmonic splines: $\phi(r) = r^{2k+1}$ or $r^{2k}\ln r$, which have no shape parameter

The parameter $\varepsilon$, the shape parameter, controls how wide each basis function is. Infinitely smooth functions such as the Gaussian and multiquadric give spectral (exponential) convergence for smooth solutions, which is the main reason they are attractive for quantum problems, where eigenfunctions are typically analytic inside the domain.

The wavefunction is approximated as

$$
\psi(\mathbf{x}) \approx \sum_{j=1}^{N} c_j\, \phi\!\left(\|\mathbf{x} - \mathbf{x}_j\|\right),
$$

with $N$ centres $\{\mathbf{x}_j\}$ scattered over $\Omega$ and its boundary. Because the basis functions are defined everywhere from the centres alone, there is no need for elements, cells, or a connectivity structure. Derivatives are obtained analytically by differentiating the basis functions. For example, for the Gaussian in $d$ dimensions,

$$
\nabla^2 \phi = \left(4\varepsilon^4 r^2 - 2d\,\varepsilon^2\right)\phi,
$$

so the kinetic-energy operator acts on the basis in closed form.

## 3. From the PDE to a matrix eigenproblem

There are two standard routes, and they lead to different matrix structures with different strengths.

### 3.1 Collocation (Kansa-type) approach

The Schrödinger equation is enforced pointwise at a set of collocation points. Split the centres into $N_i$ interior points and $N_b$ boundary points, with $N = N_i + N_b$. Substituting the RBF expansion gives, at each interior point $\mathbf{x}_i$,

$$
\sum_{j=1}^{N} c_j \left[-\tfrac{1}{2}\nabla^2\phi_j(\mathbf{x}_i) + V(\mathbf{x}_i)\phi_j(\mathbf{x}_i)\right] = E\sum_{j=1}^{N} c_j\,\phi_j(\mathbf{x}_i),
$$

and at each boundary point, the boundary condition $\sum_j c_j \phi_j(\mathbf{x}_b) = 0$. Writing these rows as a matrix pencil gives

$$
\mathbf{A}\mathbf{c} = E\,\mathbf{B}\mathbf{c},
$$

where the rows of $\mathbf{A}$ at interior points hold the Hamiltonian applied to the basis, the rows at boundary points hold the basis values, and $\mathbf{B}$ has the basis values on interior rows and zeros on boundary rows. Because $\mathbf{B}$ is rectangular or singular in the boundary rows, the pencil is treated as a generalised eigenproblem with infinite eigenvalues that must be discarded, or the problem is reduced by eliminating the boundary unknowns first.

The collocation approach is simple to implement, handles arbitrary geometry, and requires no quadrature. Its drawbacks are that the resulting matrices are nonsymmetric, which admits spurious eigenvalues that do not correspond to physical states, and that the accuracy of the eigenvalues can be sensitive to the placement of collocation points.

### 3.2 Galerkin (Rayleigh-Ritz) approach

A variational formulation is often preferable for quantum problems. Using the Rayleigh-Ritz principle, one forms the Hamiltonian and overlap matrices

$$
H_{jk} = \int_\Omega \phi_j^\ast\, \hat{H}\, \phi_k \,d\mathbf{x}, \qquad S_{jk} = \int_\Omega \phi_j^\ast\, \phi_k \,d\mathbf{x},
$$

and solves the generalised eigenproblem $\mathbf{H}\mathbf{c} = E\,\mathbf{S}\mathbf{c}$. Integrating the kinetic term by parts, which is legitimate when the boundary terms vanish, gives

$$
H_{jk} = \tfrac{1}{2}\int_\Omega \nabla\phi_j^\ast \cdot \nabla\phi_k \,d\mathbf{x} + \int_\Omega \phi_j^\ast V\,\phi_k \,d\mathbf{x}.
$$

The symmetry of $\mathbf{H}$ and $\mathbf{S}$ guarantees real eigenvalues and eliminates spurious modes, and the variational principle gives an upper bound on each eigenvalue, which is a valuable diagnostic. The cost is that the matrix elements require integrals over the domain. For Gaussian bases and polynomial or Gaussian potentials, these integrals have closed forms. For instance, the overlap of two Gaussian basis functions centred at $a$ and $b$ is

$$
S_{jk} = \sqrt{\frac{\pi}{2\varepsilon^2}}\; \exp\!\left[-\frac{\varepsilon^2}{2}(a-b)^2\right] \quad (\text{1D}),
$$

and the kinetic and potential matrix elements follow from similar Gaussian integrals. For general potentials, the integrals are evaluated by quadrature, typically on a background grid or with Gauss-Hermite rules matched to the Gaussian weight.

Galerkin is the approach most often used when accurate eigenvalues are the goal, while collocation is attractive for speed, for problems with complicated geometry, or when the potential is not easy to integrate.

## 4. The shape parameter and conditioning

The shape parameter $\varepsilon$ is the central practical difficulty of RBF methods. For a fixed set of centres, the basis becomes increasingly flat as $\varepsilon \to 0$, and the interpolation and differentiation matrices become severely ill-conditioned. Accuracy generally improves as $\varepsilon$ decreases, down to a point where round-off error dominates. This is the so-called uncertainty principle of RBF approximation: accuracy and stability pull in opposite directions.

In the eigenvalue setting this trade-off has a clear physical meaning. A large $\varepsilon$ produces narrow Gaussians that can resolve sharp features but cannot represent smooth, spread-out eigenfunctions without many terms, and the overlap matrix $\mathbf{S}$ becomes nearly the identity but the approximation of the kinetic term is poor. A small $\varepsilon$ gives broad, highly overlapping functions whose linear combinations are nearly dependent, so $\mathbf{S}$ is nearly singular.

Several strategies address this:

- Choosing $\varepsilon$ by leave-one-out cross-validation (Rippa's algorithm), or by matching the typical nearest-neighbour spacing of the centres.
- Using the RBF-QR or RBF-GA (Fornberg and coworkers) stable bases, which reformulate the expansion so that the flat limit can be computed stably.
- Contour-Padé or other rational-approximation schemes for the flat limit.
- Truncating the basis by an SVD of $\mathbf{S}$ and discarding directions with eigenvalues below a tolerance, which is simple and effective in the Galerkin setting.
- Using extended precision for the linear algebra when the target is high-accuracy spectra.

A useful practical rule is to scan $\varepsilon$ over a range and check that the low-lying eigenvalues are stable against changes in $\varepsilon$ and in the number and placement of centres. Stable eigenvalues across a plateau in $\varepsilon$ are a stronger indication of convergence than the value of any single run.

## 5. Node placement and boundary conditions

Because the method is meshfree, the centres can be placed freely, which is a significant advantage for quantum problems. Physically relevant regions, such as the classically allowed interval near the turning points of a bound state or the neighbourhood of a potential minimum, can be densified, while the far field where the wavefunction decays can be sparsely sampled. Quasi-random sequences such as Halton or Sobol points, or Chebyshev-like clustering near the boundary, are common choices.

Boundary conditions for confined systems are imposed either collocationally, as described above, or in the Galerkin setting by including the boundary points in the basis and enforcing the Dirichlet condition weakly or through a penalty. For unbounded problems one usually truncates to a large domain where the bound-state wavefunction is negligible. The truncation radius must be chosen so that the eigenvalues of interest are insensitive to it, which is an important convergence check that is easy to overlook.

## 6. A worked example: the harmonic oscillator

The one-dimensional harmonic oscillator, $V(x) = \tfrac{1}{2}\omega^2 x^2$ with exact eigenvalues $E_n = \omega\,(n + \tfrac{1}{2})$, is a standard benchmark. It illustrates the strengths of RBFs particularly well. The ground state is $\psi_0 \propto e^{-\omega x^2/2}$, which is itself a single Gaussian. With a Gaussian basis of matching width, the exact ground state lies in the span of one basis function, so the method reproduces it to machine precision in this trivial case, and the eigenvalue error is governed by how well the higher states, which require Hermite-polynomial multiplications of the Gaussian, are approximated by the rest of the basis.

Using $N$ equally or quasi-uniformly spaced Gaussian centres on a domain $[-L, L]$ with $L$ a few times the classical turning point, one typically finds that the lowest several eigenvalues converge exponentially in $N$ until the conditioning limit is reached. Higher states require more centres and a more careful choice of $\varepsilon$, because their wavefunctions have more oscillations and extend further into the classically forbidden region. This example is a good test of any implementation: the eigenvalues should approach $\omega(n+\tfrac{1}{2})$ with an error that drops rapidly as the basis is enriched, and the eigenvectors should show the expected node counts.

## 7. Potentials with singularities

The method's spectral accuracy relies on smoothness of the solution and, usually, of the potential. The hydrogen atom, with the Coulomb potential $V(r) = -1/r$, is a demanding case because the potential is singular at the origin and the wavefunction has a cusp there (the Kato cusp condition). Infinitely smooth RBFs struggle to represent the cusp, and convergence degrades to algebraic. Remedies include:

- Using polyharmonic splines or multiquadrics, which are less smooth but better suited to nonsmooth features.
- Adding centres densely near the singularity.
- Working in transformed coordinates, for example in a regularised radial variable, so that the cusp is removed before the expansion.
- Subtracting the singular part analytically and expanding only the smooth remainder.

These are the same concerns that arise for any high-order method applied to quantum singularities, and they are not unique to RBFs, but they must be addressed explicitly.

## 8. Spurious eigenvalues and verification

Collocation-based eigenproblems are notorious for producing spurious eigenvalues, particularly when the boundary conditions are imposed in a way that does not match the underlying operator's spectrum. Several checks are essential:

- Compare the collocation spectrum against a Galerkin computation with the same centres.
- Verify that low-lying eigenvalues are stable under changes in $\varepsilon$, in $N$, and in the centre distribution.
- Check the eigenfunctions for the correct number of nodes, since the $n$-th bound state must have $n$ sign changes in one dimension.
- Test the variational bound in the Galerkin case, since each computed eigenvalue should lie above the corresponding exact value.
- Test the method on a problem with a known exact solution before applying it to an unknown potential.

## 9. Higher dimensions and many-body systems

One advantage of RBFs over tensor-product grids is that the dimension enters only through the distance function and the Laplacian formula, not through an exponential growth in the number of grid nodes for a fixed resolution. For a fixed number of centres, the cost is dominated by dense linear algebra of size $N \times N$, which scales as $O(N^3)$. The curse of dimensionality is therefore not removed; it reappears through the number of centres needed to resolve the solution in $d$ dimensions. Still, meshfree placement allows the centres to concentrate where the wavefunction has structure, which can make two- and three-dimensional problems tractable on irregular domains, such as potential wells with complicated shapes or molecular-scale geometries. For many-particle problems, the wavefunction lives in $3n$ dimensions and the approach is generally impractical without additional structure, such as a separable or low-rank representation of the wavefunction.

## 10. Time-dependent extension

The same spatial discretisation can be used for the time-dependent Schrödinger equation $i\,\partial_t \psi = \hat{H}\psi$. Using the method of lines, one writes $\psi(\mathbf{x},t) \approx \sum_j c_j(t)\phi_j(\mathbf{x})$ and obtains a linear system of ODEs for the coefficients, $\mathbf{S}\dot{\mathbf{c}} = -i\,\mathbf{H}\mathbf{c}$. Time stepping is then performed by a unitary scheme such as Crank-Nicolson or a Magnus-type exponential integrator. Unitarity is important for quantum dynamics, because it preserves the norm of the wavefunction, and a norm drift is a useful diagnostic of the spatial discretisation's accuracy.

## 11. Implementation notes

For a C++ implementation, the dense generalised eigenproblem $\mathbf{H}\mathbf{c} = E\,\mathbf{S}\mathbf{c}$ can be solved with LAPACK's `sygvd` or `ggev`, or with Eigen's `GeneralizedSelfAdjointEigenSolver`, which is appropriate for the Galerkin symmetric case. Only the lowest few eigenvalues are usually needed, so iterative eigensolvers for the smallest eigenvalues are attractive for large $N$. The distance matrix and basis evaluations are embarrassingly parallel and benefit from vectorisation, and the Gaussian kernels can be evaluated with a fused exponential. Care is needed with ill-conditioned $\mathbf{S}$: a pivoted Cholesky factorisation or an SVD-based truncation should be used rather than naive inversion. Tolerances should be set relative to the largest eigenvalue of $\mathbf{S}$, and results should be reported together with the condition number, the value of $\varepsilon$, and the number of centres, so that the computation can be reproduced and its accuracy judged.

## 12. Summary

RBF methods reduce the Schrödinger eigenvalue problem to a matrix problem on a scattered set of centres, with no mesh and with derivatives computed analytically. Galerkin formulations give symmetric matrices, real spectra, and variational bounds, and are the preferred route for accurate bound-state energies. Collocation formulations are simpler and more flexible but require care with spurious modes. The shape parameter governs a fundamental trade-off between accuracy and conditioning, and it must be managed with stable algorithms or by checking the plateau of the computed eigenvalues. For smooth potentials and analytic eigenfunctions, convergence is spectral. Singular potentials such as the Coulomb potential require special treatment. The method is most effective for low-dimensional problems on irregular domains, and it extends naturally to time-dependent dynamics. Rigorous verification against benchmark problems, such as the harmonic oscillator, is essential before any method of this kind is trusted on a new problem.
