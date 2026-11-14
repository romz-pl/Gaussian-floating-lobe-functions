# Solving the 3D Schrödinger Eigenvalue Problem with Gaussian RBFs (Galerkin / Rayleigh-Ritz)

## 1. The problem

For a particle of mass $m$ in a potential $V(\mathbf{r})$, the stationary states satisfy

$$
\hat{H}\psi(\mathbf{r}) = \left[ -\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf{r}) \right]\psi(\mathbf{r}) = E\,\psi(\mathbf{r}),
$$

with $\psi$ square-integrable over $\mathbb{R}^3$. Solving this exactly is possible only for special potentials, so we seek an approximate solution in a finite-dimensional subspace.

## 2. The variational (Rayleigh-Ritz) principle

For any normalisable trial function $\psi$, the Rayleigh quotient

$$
E[\psi] = \frac{\langle \psi | \hat{H} | \psi \rangle}{\langle \psi | \psi \rangle}
$$

satisfies $E[\psi] \geq E_0$, where $E_0$ is the ground-state energy, with equality only for the exact ground state. Stationary points of $E[\psi]$ are the eigenstates of $\hat{H}$.

The Galerkin (Rayleigh-Ritz) approach restricts the search to a finite basis $\{\varphi_1, \dots, \varphi_N\}$:

$$
\psi(\mathbf{r}) \approx \sum_{j=1}^{N} c_j\, \varphi_j(\mathbf{r}).
$$

Requiring the residual $(\hat{H} - E)\psi$ to be orthogonal to every basis function, $\langle \varphi_i | \hat{H} - E | \psi \rangle = 0$ for all $i$, gives the matrix equation

$$
\sum_{j} H_{ij} c_j = E \sum_{j} S_{ij} c_j, \qquad \text{i.e.} \qquad \mathbf{H}\mathbf{c} = E\,\mathbf{S}\mathbf{c},
$$

with

$$
H_{ij} = \langle \varphi_i | \hat{H} | \varphi_j \rangle, \qquad S_{ij} = \langle \varphi_i | \varphi_j \rangle.
$$

This is a generalised Hermitian eigenvalue problem. Because the basis is nonorthogonal, $\mathbf{S}$ is not the identity. Since the subspace is finite, the eigenvalues are upper bounds to the true levels, and adding basis functions can only lower them (nested subspaces). This monotone convergence is a useful diagnostic.

## 3. Gaussian radial basis functions

Take isotropic Gaussians centred at points $\mathbf{R}_i$ with exponents $\alpha_i > 0$:

$$
\varphi_i(\mathbf{r}) = \exp\!\left( -\alpha_i \, |\mathbf{r} - \mathbf{R}_i|^2 \right).
$$

Each function is smooth, decays rapidly, and needs no explicit boundary conditions for bound states. The centres can lie on a Cartesian grid, a Halton or Sobol sequence, or be optimised. The exponents are often chosen as an even-tempered series $\alpha_k = \alpha_0 \beta^{k}$, which spreads resolution over several length scales.

The key advantage is that all the integrals needed for $\mathbf{S}$, $\mathbf{T}$ and many potentials reduce to closed forms through the Gaussian product rule.

### 3.1 The Gaussian product rule

For two Gaussians with exponents $a$ and $b$ and centres $\mathbf{A}$ and $\mathbf{B}$, define

$$
p = a + b, \qquad \mathbf{P} = \frac{a\mathbf{A} + b\mathbf{B}}{p}, \qquad K_{AB} = \exp\!\left( -\frac{ab}{p}|\mathbf{A} - \mathbf{B}|^2 \right).
$$

Then the product is a single Gaussian times a constant:

$$
e^{-a|\mathbf{r}-\mathbf{A}|^2}\, e^{-b|\mathbf{r}-\mathbf{B}|^2} = K_{AB}\, e^{-p|\mathbf{r}-\mathbf{P}|^2}.
$$

Integrals of $e^{-p|\mathbf{r}-\mathbf{P}|^2}$ against polynomials in $\mathbf{r}$ are straightforward because $\int e^{-p|\mathbf{r}-\mathbf{P}|^2} d^3r = (\pi/p)^{3/2}$ and the Gaussian moments are $\langle |\mathbf{r}-\mathbf{P}|^2 \rangle = 3/(2p)$ under the normalised weight.

## 4. Matrix elements

### 4.1 Overlap

$$
S_{ij} = \left( \frac{\pi}{p_{ij}} \right)^{3/2} K_{ij}, \qquad p_{ij} = \alpha_i + \alpha_j, \quad K_{ij} = \exp\!\left( -\frac{\alpha_i\alpha_j}{p_{ij}} d_{ij}^2 \right), \quad d_{ij} = |\mathbf{R}_i - \mathbf{R}_j|.
$$

### 4.2 Kinetic energy

The kinetic matrix element is $T_{ij} = \frac{\hbar^2}{2m}\int \nabla\varphi_i \cdot \nabla\varphi_j \, d^3r$ after one integration by parts. Using $\nabla\varphi_i = -2\alpha_i(\mathbf{r}-\mathbf{R}_i)\varphi_i$ and the product rule above gives

$$
T_{ij} = \frac{\hbar^2}{2m} \left( \frac{\pi}{p_{ij}} \right)^{3/2} K_{ij} \left( \frac{6\,\alpha_i\alpha_j}{p_{ij}} - \frac{4\,\alpha_i^2\alpha_j^2\, d_{ij}^2}{p_{ij}^2} \right).
$$

Check: for $i = j$ with $d = 0$ this gives $T = \frac{\hbar^2}{2m}\left(\frac{\pi}{2\alpha}\right)^{3/2} 3\alpha$, which agrees with a direct calculation using $-\nabla^2 e^{-\alpha r^2} = (6\alpha - 4\alpha^2 r^2)e^{-\alpha r^2}$.

### 4.3 Potential energy

For a general $V(\mathbf{r})$ there is no universal closed form, so the matrix elements are computed by quadrature (see Section 6). Several common cases are analytic.

**Polynomial potentials.** For $V = \tfrac{1}{2}m\omega^2 r^2$, using $\langle r^2 \rangle = 3/(2p) + |\mathbf{P}|^2$ for the product Gaussian,

$$
V_{ij} = \frac{1}{2}m\omega^2 \left( \frac{\pi}{p_{ij}} \right)^{3/2} K_{ij} \left( \frac{3}{2p_{ij}} + |\mathbf{P}_{ij}|^2 \right).
$$

**Gaussian potentials.** For $V = -V_0\, e^{-\beta|\mathbf{r}-\mathbf{C}|^2}$, the product of three Gaussians combines with exponent $q = p + \beta$, giving

$$
V_{ij} = -V_0 \left( \frac{\pi}{q} \right)^{3/2} K_{ij}\, \exp\!\left( -\frac{p\beta}{q}|\mathbf{P}_{ij} - \mathbf{C}|^2 \right).
$$

**Coulomb potentials.** For $V = -Z/r$ (with the nucleus at the origin), the Boys-type integral gives



$$
\left\langle \varphi_i \left| \frac{1}{r} \right| \varphi_j \right\rangle = 
\frac{2\pi}{p_{ij}} K_{ij} F_0 \left(p_{ij}|\mathbf{P}_{ij}|^2\right),
$$

$$
\qquad F_0(t) = \int_0^1 e^{-t u^2} du = \frac{\sqrt{\pi}}{2}\frac{\{erf}(\sqrt{t})}{\sqrt{t}}.
$$

Setting $\mathbf{R}_i = \mathbf{R}_j = 0$ recovers $2\pi/p$, which matches $\int e^{-pr^2}/r\, d^3r = 2\pi/p$.

## 5. Assembling and solving

The full procedure is:

1. Choose $N$ centres $\mathbf{R}_i$ and exponents $\alpha_i$.
2. Build $\mathbf{S}$, $\mathbf{T}$ and $\mathbf{V}$ (or $\mathbf{H} = \mathbf{T} + \mathbf{V}$).
3. Remove near-linear dependence in $\mathbf{S}$ (see Section 7).
4. Solve $\mathbf{H}\mathbf{c} = E\,\mathbf{S}\mathbf{c}$ for the lowest eigenpairs.
5. Reconstruct $\psi(\mathbf{r}) = \sum_i c_i \varphi_i(\mathbf{r})$ and normalise with $\mathbf{c}^{\ast}\mathbf{S}\mathbf{c} = 1$.

Because $\mathbf{H}$ and $\mathbf{S}$ are Hermitian, and $\mathbf{S}$ is positive definite, the problem can be solved stably with a Cholesky factorisation $\mathbf{S} = \mathbf{L}\mathbf{L}^{\ast}$, which reduces it to a standard Hermitian eigenproblem for $\mathbf{L}^{-1}\mathbf{H}\mathbf{L}^{-\ast}$.

Expectation values follow directly: $\langle \hat{O} \rangle = \mathbf{c}^{\ast}\mathbf{O}\mathbf{c}$ for any operator whose matrix $\mathbf{O}$ is available.

## 6. Quadrature for general potentials

When $V$ has no closed form, use a product quadrature. Two common options:

- **Gauss-Hermite in each Cartesian direction** centred on each product point $\mathbf{P}_{ij}$. Because the integrand is a Gaussian times a smooth function, a modest number of nodes gives high accuracy.
- **A fixed 3D grid** with the basis evaluated on the grid, so that $\mathbf{V}$ is a weighted sum $\sum_g w_g \varphi_i(\mathbf{r} _ g)\, V(\mathbf{r} _ g)\, \varphi _ j(\mathbf{r} _ g)$. This is simple and scales well with $N$ when the basis is stored as a matrix $\Phi_{gi}$, since $\mathbf{V} = \Phi^{\mathsf{T}} \mathrm{diag}(w \cdot V) \Phi$.

## 7. Numerical issues

**Ill-conditioning.** Gaussians with nearby centres and similar exponents are almost linearly dependent, so $\mathbf{S}$ becomes nearly singular and the generalised eigenproblem produces spurious, very negative energies. Remedies:

- Diagonalise $\mathbf{S} = \mathbf{U}\boldsymbol{\Lambda}\mathbf{U}^{\ast}$ and discard eigenvectors with $\lambda < \epsilon$ (canonical orthogonalisation), or use Löwdin symmetric orthogonalisation $\mathbf{S}^{-1/2}$ on the retained subspace.
- Monitor $\kappa(\mathbf{S})$ and choose the exponent series so that adjacent functions overlap moderately.

**Convergence and the cusp problem.** Exact Coulomb eigenfunctions have a cusp at the nucleus, $\psi \sim 1 - Zr$ near $r = 0$, and Gaussians are smooth there. Energies therefore converge algebraically rather than exponentially with $N$. Even-tempered exponents and a few functions tuned to the nuclear region help considerably.

**Optimising the basis.** The exponents and centres can be treated as nonlinear variational parameters. Minimising the lowest eigenvalue over them is a valid, if costly, refinement, since the result remains an upper bound for every choice.

## 8. Worked checks

**Harmonic oscillator.** The ground state is exactly Gaussian, $\psi_0 \propto e^{-m\omega r^2/2\hbar}$. A single s-type Gaussian at the origin, $\varphi = e^{-\alpha r^2}$, gives the energy function

$$
E(\alpha) = \frac{3\hbar^2\alpha}{2m} + \frac{3m\omega^2}{8\alpha}.
$$

Minimising over $\alpha$ yields $\alpha = m\omega/(2\hbar)$ and $E = \tfrac{3}{2}\hbar\omega$, which is exact. This is a useful unit test for any implementation.

**Hydrogen atom.** In Hartree atomic units ($\hbar = m = e = 1$) with $V = -1/r$, a single Gaussian gives $E(\alpha) = \tfrac{3}{2}\alpha - 2\sqrt{2\alpha/\pi}$. The minimum is at $\alpha = 8/(9\pi)$ with $E = -4/(3\pi) \approx -0.4244$ Ha, compared with the exact $-0.5$ Ha. Adding a few more Gaussians at different exponents (and an even-tempered series) brings the ground-state energy steadily toward $-0.5$ Ha, and the excited states are captured as well.

## 9. Sketch of an implementation

```python
import numpy as np
from scipy.linalg import eigh

def gaussian_matrices(alpha, R, hbar=1.0, m=1.0, omega=None):
    """Overlap, kinetic and (harmonic) potential matrices for isotropic s-Gaussians."""
    a = alpha[:, None]
    b = alpha[None, :]
    p = a + b
    diff = R[:, None, :] - R[None, :, :]
    d2 = np.sum(diff**2, axis=-1)
    K = np.exp(-a * b * d2 / p)
    pref = (np.pi / p) ** 1.5 * K

    S = pref
    T = (hbar**2 / (2 * m)) * pref * (6 * a * b / p - 4 * (a * b) ** 2 * d2 / p**2)

    V = None
    if omega is not None:
        P = (a[..., None] * R[:, None, :] + b[..., None] * R[None, :, :]) / p[..., None]
        P2 = np.sum(P**2, axis=-1)
        V = 0.5 * m * omega**2 * pref * (3 / (2 * p) + P2)
    return S, T, V

# Harmonic oscillator test (hbar = m = omega = 1)
alpha = np.array([0.5])               # alpha = m*omega/(2*hbar) = 0.5
R = np.zeros((1, 3))
S, T, V = gaussian_matrices(alpha, R, hbar=1.0, m=1.0, omega=1.0)
E, c = eigh(T + V, S)
print(E)                              # expect [1.5]

# Larger basis: Gaussians on a grid with an even-tempered set of exponents
alphas = 0.1 * 2.0 ** np.arange(6)
grid = np.linspace(-1.5, 1.5, 3)
Rg = np.array([[x, y, z] for x in grid for y in grid for z in grid])
A = np.repeat(alphas, len(Rg))
Rr = np.tile(Rg, (len(alphas), 1))
S, T, V = gaussian_matrices(A, Rr, omega=1.0)
# Canonical orthogonalisation to guard against near-linear dependence
lam, U = np.linalg.eigh(S)
keep = lam > 1e-8 * lam.max()
X = U[:, keep] / np.sqrt(lam[keep])
E_all = np.linalg.eigvalsh(X.T @ (T + V) @ X)
print(E_all[:4])                      # lowest levels, should approach 1.5, 2.5, ...
```

This is a sketch rather than production code. For large problems, one would exploit symmetry, use sparse quadrature for generic potentials, and replace the dense eigensolver with an iterative one for the lowest few states.

## 10. Summary

The Galerkin method turns the Schrödinger equation into the generalised eigenproblem $\mathbf{H}\mathbf{c} = E\mathbf{S}\mathbf{c}$. Gaussian radial basis functions make this practical in three dimensions because the overlap and kinetic matrices are closed-form, the product rule handles many potentials analytically, and the variational principle gives monotone upper bounds that can be used to check convergence. The main difficulties are ill-conditioning of the nonorthogonal basis and slow convergence for Coulomb-type cusps, both of which are manageable with orthogonalisation, even-tempered exponents, and careful basis design.


# Hierarchies of Gaussian Functions for High-Accuracy Eigenvalues and Eigenfunctions

This document extends *Solving the 3D Schrödinger Eigenvalue Problem with Gaussian RBFs (Galerkin / Rayleigh-Ritz)*. It assumes the notation and matrix elements established there: the basis $\varphi_i(\mathbf{r}) = \exp(-\alpha_i |\mathbf{r} - \mathbf{R}_i|^2)$, the matrices $\mathbf{S}$, $\mathbf{T}$ and $\mathbf{V}$, and the generalised eigenproblem $\mathbf{H}\mathbf{c} = E\,\mathbf{S}\mathbf{c}$.

The central idea is simple. A single Gaussian can only approximate a wavefunction well over a limited range of length scales. A hierarchy of Gaussians with exponents spanning many orders of magnitude, combined with optimised exponents, can represent the nucleus-region cusp, the intermediate shell structure, and the far-field decay at once. The result is a sequence of nested variational spaces whose eigenvalues converge monotonically from above.

## 1. The hierarchy

### 1.1 Nested subspaces

Let $\mathcal{V}_1 \subset \mathcal{V}_2 \subset \cdots$ be a sequence of basis sets, where $\mathcal{V}_k$ is spanned by $N_k$ Gaussians and each level contains all functions from the previous level. Denote by $E_n^{(k)}$ the $n$-th eigenvalue at level $k$. The Hylleraas-Undheim-MacDonald theorem gives

$$
E_n^{(k+1)} \leq E_n^{(k)}, \qquad E_n^{(k)} \geq E_n \quad \text{for all } k,
$$

where $E_n$ is the exact level. Two consequences matter in practice. First, every level is a rigorous upper bound, so the sequence can be monitored without knowing the exact answer. Second, the differences $E_n^{(k)} - E_n^{(k+1)}$ provide an empirical error estimate that is trustworthy as long as the sequence is still decreasing.

### 1.2 Tempering families

The usual way to generate the exponents of a level is a geometric progression. Three families are common:

**Even-tempered.** $\alpha_j = \alpha_0\, \beta^{\,j}$ for $j = 0, \dots, N-1$. There are only two free parameters, $\alpha_0$ and $\beta$, which keeps optimisation cheap. The exponents cover a range of about $\beta^{N}$ in scale, so increasing $N$ at fixed $\beta$ extends both the compact and the diffuse regions.

**Well-tempered.** A two-parameter family, usually with a different growth law in the tail, designed so that the exponent spacing is finer where the wavefunction varies most. It is a practical modification of the even-tempered form when a single ratio $\beta$ is too coarse near the nucleus.

**Generalised or hybrid hierarchy.** Exponents are grouped in blocks, each block an even-tempered series with its own ratio. For example, a compact block for the core, a shell block for the valence region, and a diffuse block for the tail. Each block can be optimised separately, which is the natural choice when the eigenfunction has several distinct length scales.

For high accuracy the even-tempered form is a good starting point, and a hybrid hierarchy is the usual refinement once the ground state has been converged to a few digits.

### 1.3 Refinement within a hierarchy

Given a converged level $k$, the next level can be built in three ways:

1. **Add exponents between existing ones.** Insert $\alpha_{j+1/2} = \sqrt{\alpha_j\alpha_{j+1}}$ to halve the logarithmic spacing. This is the most conservative refinement and preserves the existing optimum.
2. **Extend the range.** Add exponents above the largest or below the smallest, to capture the tail or the cusp more fully.
3. **Add new centres.** For molecular or non-central problems, add functions at additional centres. For a central potential this is usually unnecessary.

Each refinement strictly enlarges the subspace, so the new eigenvalues are no higher than the old ones. Re-optimising the exponents after each refinement often gives a larger improvement than adding functions alone.

## 2. Angular structure and excited states

### 2.1 Why s-functions are not enough

An s-type Gaussian $e^{-\alpha r^2}$ is spherically symmetric. It can describe $l=0$ states but cannot represent a $p$, $d$ or higher state, whose wavefunction has a nodal angular dependence. The minimum requirement for a $p$ state is a set of functions proportional to $z\, e^{-\alpha r^2}$ and the analogous $x$- and $y$-type functions. Similarly, $d$ states need $x^2 - y^2$-type or $z^2 - r^2/3$-type functions.

### 2.2 Generating $p$, $d$ and higher functions by differentiation

A convenient route to higher angular momentum uses the centre-derivative identity for a displaced s-Gaussian. Since

$$
\frac{\partial}{\partial A_z}\, e^{-a|\mathbf{r}-\mathbf{A}|^2} = 2a\,(z - A_z)\, e^{-a|\mathbf{r}-\mathbf{A}|^2},
$$

a $p_z$-type function at the origin is obtained as $z\,e^{-a r^2} = (2a)^{-1}\,\partial_{A_z} e^{-a|\mathbf{r}-\mathbf{A}|^2}\big|_{\mathbf{A}=0}$. Matrix elements between $p$-type functions then follow by differentiating the s-type closed forms with respect to the two centres. Higher-order derivatives produce $d$, $f$ and higher Cartesian Gaussians. Libraries typically use recursion relations (Obara-Saika or McMurchie-Davidson) instead of explicit differentiation, because the recursions scale better and avoid cancellation.

### 2.3 Restricting to the symmetry sector

For a central potential, $\mathbf{L}^2$ and $L_z$ commute with $\hat{H}$. The basis can therefore be split into angular-momentum blocks, and each block is diagonalised separately. An $l$-block contains only functions of angular type $l$, which reduces the matrix size and avoids spurious mixing. This is the standard way to obtain the lowest state in each channel to high accuracy.

## 3. Optimising the exponents

### 3.1 The variational objective

For a fixed set of centres, the energy of the lowest eigenvalue is a function of the exponents $\boldsymbol{\alpha} = (\alpha_1, \dots, \alpha_N)$. Minimising $E(\boldsymbol{\alpha})$ is a nonlinear variational problem. The result is always an upper bound, so there is no risk of overshooting the exact energy.

### 3.2 Analytic gradient

Let $\mathbf{c}$ be the normalised eigenvector, $\mathbf{c}^{\ast}\mathbf{S}\mathbf{c} = 1$. The Hellmann-Feynman theorem for the generalised problem gives

$$
\frac{\partial E}{\partial \alpha_k} = \mathbf{c}^{\ast}\!\left( \frac{\partial \mathbf{H}}{\partial \alpha_k} - E\,\frac{\partial \mathbf{S}}{\partial \alpha_k} \right)\mathbf{c}.
$$

Because $\partial\mathbf{H}/\partial\alpha_k$ and $\partial\mathbf{S}/\partial\alpha_k$ are obtained by differentiating the closed-form matrix elements, the gradient costs roughly the same as one matrix assembly. This makes quasi-Newton optimisers (BFGS, L-BFGS) practical for bases with tens of exponents.

### 3.3 Parametrisation and conditioning

Optimise $\log\alpha_k$ rather than $\alpha_k$. This enforces positivity and balances the wide dynamic range of the exponents. For the even-tempered family, only $\log\alpha_0$ and $\log\beta$ are needed. When a basis becomes nearly linearly dependent, the gradient is dominated by noise, so the optimisation should stop, or the offending function should be removed.

### 3.4 Stochastic and greedy strategies

For large or irregular bases, a stochastic variational method (SVM) is often more effective than gradient descent alone. New exponents are drawn at random over a suitable range, and a candidate is accepted only if it lowers the energy. This "add the best random Gaussian" procedure is greedy, but it explores the exponent landscape better than a fixed geometric grid when the eigenfunction has nonuniform structure. Explicitly correlated Gaussians, which include inter-particle distances in the exponent, extend the same idea to few-body problems.

## 4. Convergence behaviour and extrapolation

### 4.1 Rates

Gaussian bases converge more slowly than polynomial or Chebyshev bases for smooth problems, but they are far more practical in three dimensions. The convergence rate depends on the singularity structure:

- For smooth potentials and eigenfunctions with no cusp, the energy error falls quickly with basis size.
- For Coulomb problems, the cusp at $r = 0$ limits convergence. The energy error decreases algebraically with the number of functions in a naive basis, but tempered expansions with well-optimised exponents convert this into a much faster, approximately exponential-type decay in the number of functions. The exact rate is problem-dependent, so it should be measured for each system rather than assumed.

### 4.2 Practical error estimates

For a sequence of levels $E^{(k)}$ with basis sizes $N_k$, use:

- **Successive differences.** $\Delta_k = E^{(k)} - E^{(k+1)}$. A trustworthy estimate of the remaining error requires $\Delta_k$ to decrease steadily; a sudden plateau often indicates a basis deficiency rather than convergence.
- **Extrapolation.** If the errors fall geometrically, $\Delta_k \approx C\,\rho^{k}$, then the limit can be estimated as $E_\infty \approx E^{(k+1)} - \Delta_{k+1}\,\rho/(1-\rho)$ with $\rho$ fitted from the last few differences. This is a heuristic; it should be checked against an independently converged result when possible.
- **Comparison across families.** Converge the same quantity with two different tempering families. Agreement to the expected number of digits is strong evidence of convergence.

### 4.3 Expected accuracy

Well-optimised even-tempered or hybrid Gaussian expansions with a few tens of functions have been reported to reproduce light-atom ground-state energies to high relative accuracy. The precise accuracy depends on the system, the number of functions, and the optimisation effort, and should be established for each calculation through the convergence tests above rather than taken from a table.

## 5. Accuracy of eigenfunctions

### 5.1 Energy versus wavefunction error

The energy error is quadratic in the wavefunction error. If the trial function differs from the exact eigenfunction by $\delta\psi$ with $\|\delta\psi\| = \epsilon$, then $E - E_{\text{exact}} = O(\epsilon^2)$. The consequence is that an energy converged to $10^{-8}$ can be accompanied by a wavefunction that is accurate to only $10^{-4}$ in the $L^2$ norm. Eigenfunction accuracy must therefore be checked separately.

### 5.2 Residual-based diagnostics

Compute the residual in the Galerkin sense:

$$
\mathbf{r} = (\mathbf{H} - E\,\mathbf{S})\mathbf{c}.
$$

For an exact solution in the basis, $\mathbf{r} = \mathbf{0}$. A small residual indicates a good eigenvector within the subspace, but it does not prove that the subspace contains the exact eigenfunction. Use it together with the convergence tests of Section 4.

### 5.3 The cusp condition

For a Coulomb potential $-Z/r$ with a nucleus at the origin, the exact s-wave eigenfunction satisfies the Kato cusp condition

$$
\left.\frac{1}{\psi}\frac{\partial \psi}{\partial r}\right|_{r=0} = -Z.
$$

Evaluate the left-hand side from the expansion, using $\partial_r \varphi_i = -2\alpha_i r\,\varphi_i$ for the s-functions, which gives $\partial_r \psi\,|_{r=0} = 0$ term by term. This means the expansion cannot satisfy the cusp condition for any finite set of Gaussians with the standard form, and the discrepancy is a direct measure of basis incompleteness near the nucleus. Tracking how the cusp mismatch decreases with the number of tight (large-$\alpha$) functions is a useful eigenfunction convergence test. Remedies include adding functions that vanish at the origin, such as $r\,e^{-\alpha r^2}$-type functions for odd-parity pieces, or explicitly correlated functions that include $r$ linearly, $e^{-\alpha r + \dots}$-type behaviour in the trial space.

### 5.4 Virial test

For a Coulomb problem with $\hat{H} = \hat{T} - Z/r$, the exact eigenstates satisfy the virial theorem $2\langle \hat{T} \rangle + \langle \hat{V} \rangle = 0$, which gives $\langle \hat{T}\rangle = -E$ and $\langle \hat{V}\rangle = 2E$. Evaluate both expectation values from $\mathbf{c}$ and compute

$$
\eta = \frac{2\langle \hat{T}\rangle + \langle \hat{V}\rangle}{|\langle \hat{V}\rangle|}.
$$

The quantity $\eta$ vanishes for the exact eigenfunction. Its size in the approximate solution is an independent check on both energy and wavefunction quality, and it is sensitive to the same basis deficiencies as the cusp test.

### 5.5 Expectation values

Any observable with a matrix representation gives an expectation value $\langle \hat{O}\rangle = \mathbf{c}^{\ast}\mathbf{O}\mathbf{c}$. Expectation values that weight the nucleus region, such as $\langle \delta(\mathbf{r}) \rangle$ or $\langle 1/r \rangle$, converge more slowly than the energy and are the most demanding test of eigenfunction quality. Report them with the convergence data alongside the energy.

## 6. Numerical stability at high accuracy

### 6.1 Canonical orthogonalisation

For large hierarchies the overlap matrix becomes nearly singular, and the naive generalised eigenproblem loses digits rapidly. Use the eigendecomposition $\mathbf{S} = \mathbf{U}\boldsymbol{\Lambda}\mathbf{U}^{\ast}$, discard the eigenvalues below a threshold $\lambda_{\min}$, and form $\mathbf{X} = \mathbf{U}_{\text{kept}}\boldsymbol{\Lambda}_{\text{kept}}^{-1/2}$. The reduced problem $\mathbf{X}^{\ast}\mathbf{H}\mathbf{X}\,\mathbf{y} = E\,\mathbf{y}$ is a standard eigenproblem, and the eigenvectors in the original basis are $\mathbf{c} = \mathbf{X}\mathbf{y}$. The threshold should be chosen by checking that the lowest eigenvalues are insensitive to it.

### 6.2 Precision

Double precision is often insufficient once the condition number of $\mathbf{S}$ exceeds about $10^{12}$. In that regime, cancellation in $\mathbf{H}$ and $\mathbf{S}$ limits the attainable accuracy regardless of the eigensolver. Options include quadruple or arbitrary-precision arithmetic for the matrix assembly and eigensolve, or reducing the condition number by a better choice of exponents. Precision requirements should be established by repeating the calculation at two precisions and comparing.

### 6.3 Incremental construction

When levels are nested, the matrices of the previous level can be reused. Only the new rows and columns need to be computed, and the Cholesky factor or the canonical transformation can be updated rather than recomputed. This matters when the hierarchy has many levels.

## 7. A protocol for high-accuracy calculations

1. **Choose the problem class.** Determine the symmetry sector (for example $l=0$ or $l=1$) and the potential form.
2. **Build level 1.** A small even-tempered s-set, optimised for the ground state, gives a baseline energy and a rough length scale.
3. **Refine the hierarchy.** Double the number of exponents by inserting geometric midpoints, and re-optimise. Repeat until successive energy differences fall below the target.
4. **Add angular functions.** For excited or non-s states, add the appropriate $p$, $d$ or higher functions in the relevant symmetry block.
5. **Monitor convergence.** Track $\Delta_k$, extrapolate, and compare two tempering families.
6. **Check the eigenfunction.** Compute the residual, the cusp mismatch, the virial ratio $\eta$, and the expectation values of interest.
7. **Verify numerical stability.** Repeat the final level with a different canonical threshold and, if needed, a higher precision. The eigenvalues should agree to the claimed number of digits.
8. **Report.** Give the basis size, tempering parameters, final energies with their convergence estimates, and the eigenfunction diagnostics.

## 8. Code: even-tempered hierarchy for an s-state

The following code applies to the hydrogen-like Coulomb problem $\hat{H} = -\tfrac{1}{2}\nabla^2 - Z/r$ in atomic units, using s-type Gaussians centred at the origin. The matrix elements are

$$
S_{ij} = \left(\frac{\pi}{p}\right)^{3/2}, \qquad T_{ij} = \frac{3\alpha_i\alpha_j}{p}\left(\frac{\pi}{p}\right)^{3/2}, \qquad V_{ij} = -Z\,\frac{2\pi}{p}, \qquad p = \alpha_i + \alpha_j.
$$

```python
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize

def s_matrices(alpha, Z=1.0):
    """Overlap, kinetic and Coulomb matrices for s-Gaussians at the origin (a.u.)."""
    a = alpha[:, None]
    b = alpha[None, :]
    p = a + b
    S = (np.pi / p) ** 1.5
    T = 3.0 * a * b / p * (np.pi / p) ** 1.5
    V = -Z * 2.0 * np.pi / p
    return S, T, V

def solve(alpha, Z=1.0, n_states=1, thresh=1e-12):
    """Generalised eigenproblem with canonical orthogonalisation."""
    S, T, V = s_matrices(alpha, Z)
    lam, U = np.linalg.eigh(S)
    keep = lam > thresh * lam.max()
    X = U[:, keep] / np.sqrt(lam[keep])
    E, Y = np.linalg.eigh(X.T @ (T + V) @ X)
    C = X @ Y                      # eigenvectors in the original basis
    return E[:n_states], C[:, :n_states], alpha

def even_tempered(params, N):
    a0, beta = np.exp(params)
    return a0 * beta ** np.arange(N)

def energy_et(params, N, Z):
    alpha = even_tempered(params, N)
    E, _, _ = solve(alpha, Z)
    return E[0]

def optimise_et(N, Z=1.0, x0=(np.log(0.1), np.log(2.5))):
    res = minimize(energy_et, x0=np.array(x0), args=(N, Z), method="Nelder-Mead",
                   options={"xatol": 1e-10, "fatol": 1e-14, "maxiter": 20000})
    return res.fun, even_tempered(res.x, N)

if __name__ == "__main__":
    Z = 1.0
    exact = -0.5 * Z**2
    prev = None
    for N in range(2, 21, 2):
        E, alpha = optimise_et(N, Z)
        diff = None if prev is None else prev - E
        print(f"N={N:2d}  E={E:.10f}  err={E - exact:.3e}  "
              f"dE={'-' if diff is None else f'{diff:.2e}'}")
        prev = E

    # Eigenfunction diagnostics for the largest basis
    E, C, alpha = solve(alpha, Z)
    c = C[:, 0]
    S, T, V = s_matrices(alpha, Z)
    Ekin = c @ T @ c
    Epot = c @ V @ c
    print("virial ratio eta =", (2 * Ekin + Epot) / abs(Epot))
```

Notes on the code:

- The energy and eigenvector are obtained through canonical orthogonalisation, so the near-linear dependence of tightly spaced exponents does not produce spurious energies.
- The Nelder-Mead optimisation is simple and robust for two parameters. For larger, non-tempered exponent sets, replace it with BFGS using the gradient in Section 3.2.
- Double precision limits the attainable accuracy for large $N$. Repeat the final calculation with higher precision before reporting digits beyond the third or fourth.
- For $p$-states, extend the basis with $z\,e^{-\alpha r^2}$-type functions using the derivative construction of Section 2.2, and solve in the $l=1$ block.

## 9. Limitations and alternatives

The hierarchy of Gaussians is a practical and rigorous method, but it has limits:

- Gaussians do not reproduce the cusp, so high accuracy for nucleus-sensitive quantities requires either very large bases or cusp-adapted functions.
- For few-body systems with strong inter-particle correlation, explicitly correlated Gaussians (ECGs), which include inter-particle distances in the exponent, converge far faster than single-centre bases. Stochastic variational methods are the standard approach for these systems.
- For many-body systems, Gaussian bases become expensive because the number of integrals grows rapidly with the number of centres, and grid or spectral methods can be more efficient.
- Hylleraas-type bases, which include explicit powers of the inter-particle distance, remain the benchmark for two-electron atoms, and Gaussian results should be validated against them where available.

The method is best suited to systems where a few digits of accuracy are needed for many states, where the potential is smooth or Coulombic, and where the eigenfunction is needed as an analytic expansion for further matrix elements.

## 10. Summary

A hierarchy of Gaussian functions gives a rigorous upper-bound sequence of approximations to the eigenvalues. High accuracy comes from tempered exponent sets spanning the relevant length scales, nonlinear optimisation of the exponents using analytic gradients, angular functions for non-s states, and careful numerical treatment of the nearly singular overlap matrix. Eigenfunction accuracy must be checked separately from energy accuracy, using residuals, the Kato cusp condition, the virial ratio, and expectation values of nucleus-sensitive observables. With these checks in place, the hierarchy delivers eigenvalues and eigenfunctions whose accuracy is established by the convergence data themselves rather than assumed.

# Bibliography

## 1. Variational and Galerkin foundations

- Ritz, W. (1909). Über eine neue Methode zur Lösung gewisser Variationsprobleme der mathematischen Physik. *Journal für die reine und angewandte Mathematik*, 135, 1–61.
  The original statement of the Ritz method, which underlies the Rayleigh-Ritz approach used in Section 2.

- Galerkin, B. G. (1915). Rods and plates: Series occurring in various questions concerning the elastic equilibrium of rods and plates. *Vestnik Inzhenerov i Tekhnikov*, 19, 897–908.
  The source of the projection method that bears Galerkin's name. The paper is in Russian and is often cited through secondary sources.

- Strang, G., & Fix, G. (1973). *An Analysis of the Finite Element Method*. Prentice-Hall.
  A clear treatment of Galerkin and Ritz error bounds, useful for the convergence discussion.

- Reddy, J. N. (2013). *An Introduction to the Finite Element Method* (4th ed.). McGraw-Hill.
  Accessible coverage of weighted-residual methods and their equivalence to variational formulations.

## 2. Quantum mechanics textbooks

- Griffiths, D. J., & Schroeter, D. F. (2018). *Introduction to Quantum Mechanics* (3rd ed.). Cambridge University Press.
  Covers the variational principle, the harmonic oscillator, and the hydrogen atom, which serve as test cases here.

- Bransden, B. H., & Joachain, C. J. (2003). *Physics of Atoms and Molecules* (2nd ed.). Pearson.
  Standard reference for the Ritz variational method applied to atoms, including helium and hydrogen trial functions.

- Sakurai, J. J., & Napolitano, J. (2017). *Modern Quantum Mechanics* (2nd ed.). Cambridge University Press.
  Good background on the Rayleigh-Ritz method in the context of time-independent perturbation theory.

## 3. Gaussian integrals and Gaussian-type orbitals

- Boys, S. F. (1950). Electronic wave functions I. A general method of calculation for the stationary states of any molecular system. *Proceedings of the Royal Society A*, 200, 542–554.
  Introduces Gaussian-type functions for molecular integrals and the Boys function $F_0$ used for Coulomb matrix elements in Section 4.3.

- Shavitt, I. (1963). The Gaussian function in calculations of statistical mechanics and quantum mechanics. In B. Alder, S. Fernbach, & M. Rotenberg (Eds.), *Methods in Computational Physics*, Vol. 2. Academic Press.
  An early comprehensive review of Gaussian integral techniques, including the product rule.

- Hehre, W. J., Stewart, R. F., & Pople, J. A. (1969). Self-consistent molecular-orbital methods. I. Use of Gaussian expansions of Slater-type atomic orbitals. *The Journal of Chemical Physics*, 51, 2657–2664.
  The STO-nG basis-set construction, showing how contracted Gaussians approximate Slater orbitals.

- Hehre, W. J., Ditchfield, R., & Pople, J. A. (1972). Self-consistent molecular orbital methods. XII. Further extensions of Gaussian-type basis sets for use in molecular orbital studies of organic molecules. *The Journal of Chemical Physics*, 56, 2257–2261.
  Extends the Gaussian basis-set approach used for the even-tempered and split-valence ideas discussed in Section 7.

- Dunning, T. H., Jr. (1989). Gaussian basis sets for use in correlated molecular calculations. I. The atoms boron through neon and hydrogen. *The Journal of Chemical Physics*, 90, 1007–1023.
  Introduces correlation-consistent basis sets, relevant to the convergence and basis-design questions.

- Frost, A. A. (1967). Gaussian orbitals for atoms and molecules. *The Journal of Chemical Physics*, 47, 3707–3713.
  Early work on optimising Gaussian exponents for atomic and molecular calculations.

## 4. Quantum chemistry methods and basis-set practice

- Szabo, A., & Ostlund, N. S. (1996). *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory* (revised ed.). Dover.
  Practical treatment of Gaussian integrals, the Galerkin-type matrix eigenproblem $\mathbf{F}\mathbf{C} = \mathbf{S}\mathbf{C}\boldsymbol{\varepsilon}$, and Löwdin orthogonalisation.

- Helgaker, T., Jørgensen, P., & Olsen, J. (2000). *Molecular Electronic-Structure Theory*. Wiley.
  The most thorough modern reference on Gaussian integral evaluation, the product rule, and basis-set theory.

- Löwdin, P.-O. (1950). On the non-orthogonality problem connected with the use of atomic wave functions in the theory of molecules and crystals. *The Journal of Chemical Physics*, 18, 365–375.
  The symmetric orthogonalisation $\mathbf{S}^{-1/2}$ used to handle ill-conditioned overlap matrices in Section 7.

- Pekeris, C. L. (1958). 1¹S, 2¹S, 2³S, 2¹P, and 3¹S states of helium. *Physical Review*, 112, 1649–1658.
  Demonstrates very high accuracy from correlated trial functions, providing a benchmark for basis-set convergence.

- Hylleraas, E. A. (1929). Neue Berechnung der Energie des Heliums im Grundzustande, sowie des tiefsten Terms von Wasserstoff. *Zeitschrift für Physik*, 54, 347–366.
  The original correlated-function calculation for helium, which motivated the inclusion of interparticle coordinates.

## 5. Cusp conditions and convergence

- Kato, T. (1957). On the eigenfunctions of many-particle systems in quantum mechanics. *Communications on Pure and Applied Mathematics*, 10, 151–177.
  Derives the electron-nucleus cusp condition that explains the algebraic convergence of Gaussian bases discussed in Section 7.

- Schwartz, C. (1962). Importance of angular correlations between atomic electrons. *Physical Review*, 126, 1015–1019.
  Shows how slow basis-set convergence arises from the cusp and electron correlation, and how to accelerate it.

- Bingel, W. A. (1963). Die Berechnung von Mehrelektronen-Eigenfunktionen mit Hilfe der Methode von Ritz/Hylleraas. *Zeitschrift für Naturforschung A*, 18, 1249–1250. (Note: the Bingel cusp-based approach. Verify exact details before citing.)

## 6. Radial basis function theory

- Buhmann, M. D. (2003). *Radial Basis Functions: Theory and Implementations*. Cambridge University Press.
  Standard monograph on RBF approximation theory, including Gaussian kernels and their native spaces.

- Wendland, H. (2005). *Scattered Data Approximation*. Cambridge University Press.
  Covers positive-definite kernels, approximation orders, and the stability of RBF interpolation.

- Fasshauer, G. E. (2007). *Meshfree Approximation Methods with MATLAB*. World Scientific.
  Practical guide to Gaussian RBFs, including the shape-parameter trade-off and ill-conditioning.

- Fornberg, B., & Flyer, N. (2015). *A Primer on Radial Basis Functions with Applications to the Geosciences*. SIAM.
  Accessible introduction covering Gaussian kernels, shape parameters, and the stable-computation methods relevant to Section 7.

- Schaback, R. (1995). Error estimates and condition numbers for radial basis function interpolation. *Advances in Computational Mathematics*, 3, 251–264.
  Establishes the link between kernel smoothness, condition numbers, and interpolation error, which underlies the conditioning discussion.

- Hardy, R. L. (1971). Multiquadric equations of topography and other irregular surfaces. *Journal of Geophysical Research*, 76, 1905–1915.
  The original RBF paper, introducing multiquadrics and motivating the broader RBF literature.

- Kansa, E. J. (1990). Multiquadrics: A scattered data approximation scheme with applications to computational fluid dynamics. *Computers & Mathematics with Applications*, 19, 127–145.
  The first RBF collocation method for PDEs, the precursor to using RBFs as trial spaces in Galerkin formulations.

## 7. RBFs and Galerkin methods for PDEs

- Kansa, E. J. (1990). See Section 6 above. The collocation variant is the closest relative of the Galerkin approach used here.

- Fasshauer, G. E. (2002). Newton iteration with multiquadrics for the solution of nonlinear PDEs. *Computers & Mathematics with Applications*, 43, 423–438.
  Demonstrates RBF trial spaces for nonlinear PDEs, an extension relevant to self-consistent potentials.

- Chen, W., Fu, Z.-J., & Chen, C. S. (2014). *Recent Advances in Radial Basis Function Collocation Methods*. Springer.
  Survey of RBF collocation and Galerkin-type schemes, including boundary-value and eigenvalue problems.

## 8. Numerical linear algebra and eigenproblems

- Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins University Press.
  The definitive reference for the generalised symmetric eigenproblem $\mathbf{H}\mathbf{c} = E\mathbf{S}\mathbf{c}$, Cholesky reduction, and conditioning.

- Trefethen, L. N., & Bau, D. (1997). *Numerical Linear Algebra*. SIAM.
  Clear discussion of conditioning, eigenvalue computation, and stability.

- Parlett, B. N. (1998). *The Symmetric Eigenvalue Problem*. SIAM.
  Detailed treatment of symmetric eigensolvers, useful for choosing a method for the lowest eigenpairs.

- Saad, Y. (2011). *Numerical Methods for Large Eigenvalue Problems* (revised ed.). SIAM.
  Covers iterative eigensolvers for large sparse or dense problems, relevant when $N$ is large.

## 9. Special functions and integrals

- Abramowitz, M., & Stegun, I. A. (Eds.). (1964). *Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables*. National Bureau of Standards.
  Source for the error function and related special functions used in the Coulomb matrix element.

- Gradshteyn, I. S., & Ryzhik, I. M. (2015). *Table of Integrals, Series, and Products* (8th ed.). Academic Press.
  Reference for Gaussian moment and product integrals.

- Olver, F. W. J., Lozier, D. W., Boisvert, R. F., & Clark, C. W. (Eds.). (2010). *NIST Handbook of Mathematical Functions*. Cambridge University Press.
  Modern successor to Abramowitz and Stegun, with current definitions of $F_0$-type functions.

## 10. Numerical methods and software

- Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.
  General numerical methods, including quadrature and eigensolvers.

- Virtanen, P., et al. (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272.
  Citation for `scipy.linalg.eigh` and related routines used in the implementation sketch.

- Harris, C. R., et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362.
  Citation for the NumPy operations used in the implementation sketch.

## 11. Further reading on spectral and high-accuracy methods

- Boyd, J. P. (2001). *Chebyshev and Fourier Spectral Methods* (2nd ed.). Dover.
  Comparison point for spectral basis expansions versus Gaussian bases, especially regarding exponential convergence.

- Thijssen, J. M. (2007). *Computational Physics* (2nd ed.). Cambridge University Press.
  Covers variational eigenvalue computations and basis-set methods for quantum problems.
