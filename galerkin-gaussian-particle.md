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
\left\langle \varphi_i \left| \frac{1}{r} \right| \varphi_j \right\rangle = \frac{2\pi}{p_{ij}}\, K_{ij}\, F_0\!\left(p_{ij}|\mathbf{P}_{ij}|^2\right), \qquad F_0(t) = \int_0^1 e^{-t u^2} du = \frac{\sqrt{\pi}}{2}\frac{\operatorname{erf}(\sqrt{t})}{\sqrt{t}}.
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
- **A fixed 3D grid** with the basis evaluated on the grid, so that $\mathbf{V}$ is a weighted sum $\sum_g w_g \varphi_i(\mathbf{r}_g)\, V(\mathbf{r}_g)\, \varphi_j(\mathbf{r}_g)$. This is simple and scales well with $N$ when the basis is stored as a matrix $\Phi_{gi}$, since $\mathbf{V} = \Phi^{\mathsf{T}} \mathrm{diag}(w \cdot V) \Phi$.

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
