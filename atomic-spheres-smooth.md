# Hartree-Fock with Confined Numerical Atomic Orbitals with Zero Condition and Spherical Gaussians

## 1. Scope and motivation

This document describes a mixed-basis Hartree-Fock (HF) method. The one-electron space is spanned by two kinds of functions:

1. **Confined numerical atomic orbitals (NAOs)**, which live inside non-overlapping atomic spheres and vanish at each sphere boundary.
2. **s-type spherical Gaussians**, which are nonzero everywhere, both in the interstitial region and inside the spheres.

The motivation is a division of labour. Inside the spheres, the wavefunction has sharp, nodal structure inherited from the atomic core and valence regions, and numerical radial functions represent it compactly and accurately. In the interstitial region, the wavefunction is smooth and slowly varying, and a modest number of Gaussians suffices. The Gaussians also allow the orbitals to be represented in a region where no atomic-centred basis is natural, and they give the basis enough flexibility to be variationally complete in the interstitial region.

The method is variational. The total energy is an upper bound to the exact Hartree-Fock energy for the given basis, and the orbitals are determined by the Roothaan-Hall equations solved self-consistently. The principal technical difficulties are (i) the enforcement of the sphere-boundary conditions, (ii) the evaluation of mixed NAO-Gaussian integrals, and (iii) the cost of the two-electron integrals and exchange build.

## 2. Geometry and partition of space

### 2.1 Atomic spheres

Let the system contain $N_{\text{at}}$ atoms at positions $\mathbf{R} _ A$, $A = 1, \dots, N _ {\text{at}}$. Each atom is assigned a sphere $\Omega_A$ of radius $R_A$ centred on $\mathbf{R}_A$. The spheres are chosen to be non-overlapping:

$$
|\mathbf{R}_A - \mathbf{R}_B| > R_A + R_B \qquad (A \neq B).
$$

The interstitial region $\Omega_{\text{int}}$ is the complement of the union of all spheres within the computational domain (a molecule embedded in vacuum, or a unit cell of a periodic solid):

$$
\Omega_{\text{int}} = \Omega \setminus \bigcup_A \Omega_A.
$$

The sphere radii are chosen so that the spheres fill a substantial fraction of the cell where the electron density is concentrated, while leaving a nonvanishing interstitial region. Typical choices set $R_A$ near the covalent or ionic radius of the element, scaled so that adjacent spheres do not touch. For the method to be well defined, the spheres must not overlap, because the NAO radial functions are defined separately in each sphere.

### 2.2 Why non-overlapping spheres

Overlapping spheres would require either a prescription for the value of each NAO in the overlap region, or a partition-of-unity scheme, and this would introduce additional ambiguities. Non-overlapping spheres make the NAO on atom $A$ a well-defined function that is zero outside $\Omega_A$. The entire problem is then specified by the sphere radii, the radial functions, and the Gaussian set.

## 3. Basis set construction

### 3.1 Numerical atomic orbitals

For atom $A$, the NAO basis is

$$
\varphi_{n l m}^{A}(\mathbf{r}) =
\begin{cases}
R_{n l}^{A}(r)\, Y_{l m}(\hat{\mathbf{r}}), & r = |\mathbf{r} - \mathbf{R}_A| \le R_A, \\[4pt]
0, & r > R_A,
\end{cases}
$$

where $Y_{lm}$ are real spherical harmonics (or complex ones, with appropriate bookkeeping), and $n$ labels radial functions within a given angular momentum channel $l$. The radial functions $R_{nl}^{A}$ are given numerically on a radial grid $r_1 < r_2 < \dots < r_{N_r} = R_A$.

Three properties of the radial functions are important.

**Boundary zero condition.** Every radial function satisfies

$$
R_{n l}^{A}(R_A) = 0.
$$

This is the condition that makes the NAO vanish smoothly at the sphere edge in value.

**Origin of the radial functions.** The radial functions are typically obtained from one of the following:

- Solutions of the radial Schrödinger or Kohn-Sham equations for the free atom, cut off at $R_A$ (a "confined" atom, where the atomic potential is supplemented by a confining potential so that the eigenfunctions vanish at $R_A$ naturally).
- Numerically optimised radial functions, obtained by minimising the energy of a small molecular or solid-state test system subject to $R_{nl}^{A}(R_A) = 0$.
- Pseudo-atomic orbitals or Slater-type radial functions multiplied by a smooth window function that brings the value to zero at $R_A$.

The confined-atom route is the most physically grounded. A radial equation of the form

$$
\left[-\frac{1}{2}\frac{d^2}{dr^2} + \frac{l(l+1)}{2r^2} + V_{\text{at}}(r) + V_{\text{conf}}(r)\right] u_{nl}(r) = \varepsilon_{nl}\, u_{nl}(r), \qquad u_{nl}(r) = r\, R_{nl}(r),
$$

with a confining potential $V_{\text{conf}}$ that diverges or becomes large at $R_A$, produces radial functions that are naturally small at $R_A$. The zero condition is then enforced by requiring $u_{nl}(R_A) = 0$, which is a Dirichlet boundary condition. The confining potential has the effect of compressing the orbitals inside the sphere, so the zero condition is satisfied with a smooth, physically motivated shape rather than an ad hoc cut.

**Normalisation.** Each radial function is normalised within its sphere:

$$
\int_0^{R_A} \left[R_{nl}^{A}(r)\right]^2 r^2\, dr = 1.
$$

The angular normalisation of $Y_{lm}$ then gives $\langle \varphi_{nlm}^{A} | \varphi_{nlm}^{A}\rangle = 1$, up to the nonorthogonality of different $n$ and $A$.

### 3.2 Spherical Gaussian set

The interstitial basis is a set of s-type Gaussians

$$
g_\nu(\mathbf{r}) = \exp\!\left(-\alpha_\nu\, |\mathbf{r} - \mathbf{R}_\nu|^2\right),
$$

with exponents $\alpha_\nu > 0$ and centres $\mathbf{R}_\nu$. The centres may coincide with atomic positions, or may be placed at bond midpoints, in the interstitial region, or at other points chosen to describe the diffuse tail of the density. The exponents form a set spanning several decades, from diffuse functions (small $\alpha$) to moderately compact functions (large $\alpha$), so that both the long-range interstitial behaviour and the short-range behaviour near sphere boundaries are representable.

Important features of the Gaussian set:

- Each $g_\nu$ is smooth and nonzero everywhere, including inside the spheres.
- Because Gaussians are not confined, they can represent the portion of the orbital that extends into the interstitial region.
- The Gaussians are not required to vanish at sphere boundaries. Their role is to supply the smooth background that the confined NAOs cannot provide.

The Gaussians are not required to be orthogonal to each other or to the NAOs. This is handled in the overlap matrix, as described in Section 6.

### 3.3 The combined basis

The total basis is the union of the two sets:

$$
\mathcal{B} = \{\varphi_{n l m}^{A}\}_{A, n, l, m} \;\cup\; \{g_\nu\}_\nu.
$$

Each orbital is expanded as

$$
\psi_i(\mathbf{r}) = \sum_{\mu \in \mathcal{B}} c_{i\mu}\, \chi_\mu(\mathbf{r}),
$$

where $\chi_\mu$ ranges over both NAOs and Gaussians. The total dimension of the basis is $N_{\text{NAO}} + N_{G}$, where $N_{\text{NAO}}$ counts all radial-angular NAO functions on all atoms and $N_G$ counts the Gaussians.

## 4. Boundary conditions and continuity

This section requires particular care, because the conditions at the sphere boundary are the point where the method most easily goes wrong.

### 4.1 What the zero condition guarantees

Consider a point $\mathbf{r}_s$ on the boundary $\partial\Omega_A$ of sphere $A$. Approached from inside the sphere, the orbital is

$$
\psi_i(\mathbf{r}_s^{-}) = \sum_{\mu \in \text{NAO}(A)} c_{i\mu}\, \varphi_\mu(\mathbf{r}_s) + \sum_\nu c_{i\nu}\, g_\nu(\mathbf{r}_s).
$$

The NAO part vanishes at $\mathbf{r} _ s$, because every radial function satisfies $R _ {nl}^{A}(R_A) = 0$, and the real spherical harmonics are bounded. So

$$
\psi_i(\mathbf{r} _ s^{-}) = \sum_\nu c_{i\nu}\, g_\nu(\mathbf{r}_s) = \psi_i(\mathbf{r}_s^{+}),
$$

where the last equality holds because, outside the sphere, the NAO part is zero and the Gaussian part is the same function. Therefore the **value** of the orbital is continuous across the sphere boundary, exactly and by construction, for every choice of coefficients.

### 4.2 What the zero condition does not guarantee

The derivative is a different matter. The normal derivative of an NAO at the boundary is

$$
\frac{\partial}{\partial r}\varphi_{nlm}^{A}\Big|_{r = R_A} = \frac{d R_{nl}^{A}}{dr}\bigg|_{R_A} Y_{lm}(\hat{\mathbf{r}}),
$$

which is generally nonzero, because $R_{nl}^{A}$ is zero at $R_A$ but its slope is not. Outside the sphere, the NAO contribution to the derivative is zero. So a generic linear combination of NAOs has a normal derivative that jumps across the boundary, unless the combined coefficients are chosen to cancel the derivative.

The Gaussians have a normal derivative that is continuous across the boundary, since they are smooth everywhere. Therefore the derivative mismatch is entirely due to the NAO part, and the continuity of the derivative of the full orbital is not automatic.

This means that the statement "the zero condition enforces continuity of the orbital and its derivative" must be read with care. The zero condition enforces value continuity exactly. For derivative continuity, the method relies on the variational principle: the expansion coefficients are determined by minimising the energy, and the minimisation drives the combined orbital toward one whose derivative matches across the boundary, to the extent that the basis can represent such a function. This is an approximate, basis-dependent enforcement rather than an exact constraint.

This distinction has practical consequences. The derivative mismatch can be reduced by:

- Enlarging the NAO radial basis, so that more radial functions are available to describe the boundary behaviour.
- Enlarging the Gaussian set, so that the interstitial part of the orbital can compensate for mismatches at the boundary.
- Choosing radial functions whose derivatives at $R_A$ are small, for example by using a confining potential that makes the radial functions flatten smoothly toward the boundary.

If an exact derivative constraint is required, one option is to add explicit Lagrange multipliers or a penalty term to the variational functional, enforcing $\partial \psi/\partial n$ continuity on each sphere boundary. This is more involved, and most implementations of the basic mixed-basis method accept the variational approximation.

### 4.3 Summary of boundary treatment

- Value continuity: exact, guaranteed by the zero condition on the NAOs together with smoothness of the Gaussians.
- Derivative continuity: approximate, variational, and improving with basis size.
- No explicit matching equations appear in the working equations. The Rayleigh-Ritz procedure applies directly to the full expansion.

## 5. Hartree-Fock theory in this basis

### 5.1 Energy functional

For a closed-shell system with $N_e$ electrons and $N_{\text{occ}} = N_e/2$ doubly occupied orbitals, the Hartree-Fock energy is

$$
E[P] = \sum_{\mu\nu} P_{\mu\nu}\, H^{\text{core}}_{\mu\nu} + \frac{1}{2}\sum_{\mu\nu\lambda\sigma} P_{\mu\nu} P_{\lambda\sigma}\left[(\mu\nu|\lambda\sigma) - \frac{1}{2}(\mu\lambda|\nu\sigma)\right] + E_{\text{nuc}}.
$$

Here $P_{\mu\nu} = 2\sum_{i \in \text{occ}} c_{i\mu} c_{i\nu}^{\ast}$ is the total density matrix, $H^{\text{core}} _ {\mu\nu} = T_{\mu\nu} + V_{\mu\nu}$ collects kinetic and nuclear attraction integrals, and $E_{\text{nuc}}$ is the nuclear repulsion energy. The two-electron integrals are

$$
(\mu\nu|\lambda\sigma) = \iint \frac{\chi_\mu^{\ast}(\mathbf{r}_1)\chi_\nu(\mathbf{r}_1)\chi_\lambda^{\ast}(\mathbf{r}_2)\chi_\sigma(\mathbf{r}_2)}{|\mathbf{r}_1 - \mathbf{r}_2|}\, d\mathbf{r}_1\, d\mathbf{r}_2.
$$

The basis functions $\chi_\mu$ are real in the present treatment, so the complex conjugations can be dropped.

### 5.2 Fock matrix

Define the Coulomb and exchange matrices

$$
J_{\mu\nu}[P] = \sum_{\lambda\sigma} (\mu\nu|\lambda\sigma)\, P_{\lambda\sigma}, \qquad K_{\mu\nu}[P] = \sum_{\lambda\sigma} (\mu\lambda|\nu\sigma)\, P_{\lambda\sigma}.
$$

With $G[P] = J[P] - \tfrac{1}{2}K[P]$, the Fock matrix is

$$
F_{\mu\nu} = H^{\text{core}}_{\mu\nu} + G_{\mu\nu}[P].
$$

### 5.3 Roothaan-Hall equations

The orbitals satisfy the generalised eigenvalue problem

$$
F\, C = S\, C\, \varepsilon,
$$

where $S_{\mu\nu} = \langle \chi_\mu | \chi_\nu \rangle$ is the overlap matrix, $C$ is the matrix of orbital coefficients, and $\varepsilon$ is diagonal. The basis is nonorthogonal, so $S \neq 1$ in general. The energy is stationary with respect to variations of $C$ subject to the orthonormality condition $C^{\dagger} S C = 1$.

The Roothaan-Hall equations are the working equations. The mixed basis does not change their form; it changes only the matrices $S$, $H^{\text{core}}$, and the two-electron integrals that enter $G$.

### 5.4 Spin-unrestricted and open-shell cases

For open-shell systems, one uses separate alpha and beta density matrices $P^{\alpha}$ and $P^{\beta}$, with

$$
F^{\sigma} = H^{\text{core}} + J[P^{\alpha} + P^{\beta}] - K[P^{\sigma}], \qquad \sigma \in \{\alpha, \beta\}.
$$

The unrestricted equations are solved as two coupled eigenproblems. The same integrals serve both spin channels, so the additional cost is mainly in the diagonalisation and the DIIS extrapolation.

## 6. Overlap, orthogonalisation, and conditioning

### 6.1 One-electron matrices

The one-electron matrices $S$, $T$, and $V$ are built in the nonorthogonal basis. Their evaluation is described in Section 7.

### 6.2 Conditioning of the basis

Because the Gaussians can be very diffuse and because NAOs on neighbouring atoms can overlap in the interstitial region (their spheres do not overlap, but the Gaussians do), the overlap matrix can become nearly singular. A near-singular $S$ leads to large, unphysical coefficients and numerical noise in the eigenproblem.

The standard remedy is to transform the problem to an orthonormal basis. Two options are common.

**Symmetric orthogonalisation.** Form $X = S^{-1/2}$ by diagonalising $S = U s U^{\dagger}$ and setting $X = U s^{-1/2} U^{\dagger}$. This is symmetric in the basis functions but is sensitive to near-linear dependence.

**Canonical orthogonalisation.** Diagonalise $S = U s U^{\dagger}$ and form $X = U_{\text{kept}}\, s_{\text{kept}}^{-1/2}$, keeping only eigenvalues $s_k$ above a threshold $\tau$ (typically $10^{-6}$ to $10^{-8}$). This removes linearly dependent combinations and is the robust default.

After the transformation, the generalised problem $FC = SC\varepsilon$ becomes the standard problem $F' C' = C' \varepsilon$, with $F' = X^{\dagger} F X$ and $C = X C'$.

### 6.3 Removing near-linear dependence

Removing eigenvectors of $S$ reduces the number of independent orbitals. The effective basis dimension is then $N_{\text{eff}} \le N_{\text{NAO}} + N_G$. The threshold must be chosen carefully: too aggressive a cut discards physically important functions, while too conservative a cut leaves ill-conditioning. A practical procedure monitors the total energy as a function of $\tau$ and selects a value at which the energy is stable.

## 7. Integral evaluation

This is the most technically demanding part of the method. The integrals fall into several classes, depending on whether the basis functions are NAOs or Gaussians, and where they are centred.

### 7.1 Classes of integrals

For one-electron integrals (overlap, kinetic, nuclear attraction) and for two-electron integrals, the basis function pairs fall into three categories:

1. **Gaussian-Gaussian.** These have closed analytic forms, given by the Gaussian product theorem. The overlap of two s-type Gaussians is

$$
\langle g_\mu | g_\nu \rangle = \left(\frac{\pi}{\alpha_\mu + \alpha_\nu}\right)^{3/2} \exp\!\left(-\frac{\alpha_\mu\alpha_\nu}{\alpha_\mu + \alpha_\nu}|\mathbf{R}_\mu - \mathbf{R}_\nu|^2\right),
$$

and the nuclear attraction and electron repulsion integrals reduce to the Boys function $F_0$ in the standard way.

2. **NAO-Gaussian.** These are the mixed integrals. Because NAOs are numerical, the analytic Gaussian machinery does not apply directly. The integrand is a product of a numerical radial function, a spherical harmonic, and a Gaussian, integrated over space.

3. **NAO-NAO.** These are numerical integrals over one or two spheres. For two NAOs on the same atom, the integrals reduce to radial quadrature over $[0, R_A]$ combined with angular integration over the sphere. For NAOs on different atoms, the spheres are separated, and the integral is a double sphere integral, which is evaluated using a one-centre expansion or a numerical grid.

### 7.2 Numerical quadrature inside spheres

For integrals that involve an NAO on atom $A$, the natural quadrature is a radial-angular grid centred on $\mathbf{R}_A$:

$$
\int_{\Omega_A} f(\mathbf{r})\, d\mathbf{r} \approx \sum_{p} w_p^{\text{rad}} \sum_{q} w_q^{\text{ang}}\, r_p^2\, f(r_p, \hat{\mathbf{r}}_q).
$$

The radial grid uses Gauss-Legendre or Gauss-Chebyshev nodes mapped onto $[0, R_A]$, with sufficient density near the nucleus and near the sphere boundary. The angular grid is typically a Lebedev quadrature, which integrates spherical harmonics exactly up to a chosen degree. The required degree depends on the highest angular momentum appearing in the integrand.

For integrals involving Gaussians, the same grid can be used when the Gaussian is centred on the same atom. When the Gaussian is centred elsewhere, its spherical-harmonic expansion about $\mathbf{R}_A$ is needed. This is a standard procedure: the Gaussian is expanded in regular solid harmonics about $\mathbf{R}_A$, and the coefficients are computed from the Gaussian product and the appropriate addition theorem.

### 7.3 Interstitial integrals

Integrals over the interstitial region are more awkward, because the region is not a simple geometric shape. Two approaches are common.

**Full-cell grid.** Integrate over the entire computational cell with a uniform or adaptive grid, and mask the contributions inside the spheres (or use the spherical grids inside the spheres as above). This is straightforward but expensive for large cells.

**Atom-centred partition.** Use a partition of unity, such as Becke weights, to split the cell into atom-centred regions. Each region is integrated on its own radial-angular grid. This is accurate and scales well, but requires care with the partition near sphere boundaries.

In both approaches, the Gaussian contributions in the interstitial region are typically the easiest to integrate, because the Gaussians are smooth and decay rapidly outside their centres.

### 7.4 Two-electron integrals

The two-electron integrals $(\mu\nu|\lambda\sigma)$ are the dominant cost. They are classified by the types of the four basis functions:

- Four Gaussians: analytic, using the Obara-Saika or McMurchie-Davidson recursion for s-type functions, which reduces to Boys functions.
- Mixed NAO-Gaussian: a numerical integral over one sphere, with the Gaussian contribution handled by its expansion about the sphere centre.
- Four NAOs on two or more spheres: the most expensive class, requiring a double sphere integration, but also the class most amenable to sparsity, because NAOs on distant atoms produce negligible contributions.

### 7.5 Integral screening

Because NAOs vanish outside their spheres, many integrals involving NAOs on distant atoms are exactly zero. Integral screening exploits this. For each pair of basis functions $(\mu\nu)$, the Schwarz inequality gives a bound

$$
|(\mu\nu|\lambda\sigma)| \le \sqrt{(\mu\nu|\mu\nu)\,(\lambda\sigma|\lambda\sigma)},
$$

and integrals whose bound falls below a threshold are skipped. Shell-pair screening and spatial screening based on sphere overlap are both effective. For NAOs, the sphere geometry gives an exact cutoff: if two NAO supports are disjoint, the integral is zero.

### 7.6 Storage and direct evaluation

The two-electron integrals can be stored in memory or on disk, or computed on the fly ("direct SCF"). Storage is fast per iteration but scales as $N^4$ in the worst case. Direct evaluation avoids storage but repeats the integral cost every SCF cycle. Most practical implementations use a hybrid: store integrals that are expensive to recompute and screen the rest.

## 8. Initial guess

A good initial guess reduces the number of SCF cycles. Options include:

- **Core Hamiltonian guess.** Diagonalise $F' = X^{\dagger} H^{\text{core}} X$ and fill the lowest $N_{\text{occ}}$ orbitals. This is simple but often converges slowly or to a poor state for systems with significant electron correlation in the density.
- **Superposition of atomic densities.** Build the initial density as a sum of atomic densities, each computed from the confined-atom radial functions in its own sphere. This is physically motivated because the NAOs already come from atomic calculations.
- **Previous geometry.** For geometry optimisations or molecular dynamics, reuse the converged density from the previous step, transformed to the new basis if the basis moved with the atoms.

The superposition-of-atomic-densities guess is particularly natural here and is recommended as the default.

## 9. Self-consistent field procedure

### 9.1 Algorithm

Given the basis, integrals, and an initial guess $P^{(0)}$, the SCF loop proceeds as follows at iteration $k$:

1. Build $J[P^{(k)}]$ and $K[P^{(k)}]$ from the two-electron integrals, using screening.
2. Form the Fock matrix $F^{(k)} = H^{\text{core}} + J[P^{(k)}] - \tfrac{1}{2}K[P^{(k)}]$.
3. Compute the DIIS error vector $e^{(k)} = F^{(k)}P^{(k)}S - SP^{(k)}F^{(k)}$, which vanishes at convergence.
4. Apply DIIS extrapolation to the Fock matrix, using the stored history of Fock matrices and error vectors. This produces an extrapolated $\tilde{F}^{(k)}$.
5. Transform to the orthonormal basis, $F' = X^{\dagger}\tilde{F}^{(k)}X$, and diagonalise to get $C'$ and $\varepsilon$.
6. Back-transform, $C = XC'$, and form the new density $P^{(k+1)}_{\mu\nu} = 2\sum_{i \in \text{occ}} C_{\mu i} C_{\nu i}^{\ast}$.
7. Compute the energy $E^{(k+1)}$ from the density, and check convergence.

### 9.2 DIIS

Direct inversion in the iterative subspace (DIIS) extrapolates the Fock matrix by minimising the norm of a linear combination of past error vectors. Given a history of Fock matrices $\{F_j\}$ and error vectors $\{e_j\}$, one solves

$$
\min_{\{c_j\}} \left\| \sum_j c_j e_j \right\|^2 \quad \text{subject to} \quad \sum_j c_j = 1,
$$

and sets $\tilde{F} = \sum_j c_j F_j$. The subspace typically holds 6 to 10 vectors.

DIIS accelerates convergence dramatically when the iterations are in the basin of a stable solution. It can fail or oscillate in difficult cases, which is handled by damping, level shifting, or restarting the subspace.

### 9.3 Damping and level shifting

Two common stabilisers are used when DIIS alone is insufficient.

**Density damping.** Mix the new density with the old one, $P^{(k+1)} \leftarrow (1-\lambda) P^{(k+1)} + \lambda P^{(k)}$, with $\lambda$ between 0.2 and 0.5 in early iterations.

**Level shifting.** Add a positive shift to the virtual orbital energies, $F \leftarrow F + \mu\, S P_{\text{virt}} S$, which increases the HOMO-LUMO gap and suppresses oscillations between occupied and virtual spaces. The shift is removed once the density is stable.

### 9.4 Convergence criteria

Convergence is declared when all of the following hold:

- The change in total energy between iterations is below $\Delta E_{\text{tol}}$, typically $10^{-6}$ to $10^{-8}$ Hartree.
- The DIIS error norm $\| e^{(k)} \|$ is below $\varepsilon_{\text{tol}}$, typically $10^{-5}$ to $10^{-7}$.
- The density matrix change $\| P^{(k+1)} - P^{(k)} \|$ is below a threshold.

Tighter thresholds are needed for properties such as forces and dipole moments.

## 10. Total energy and observables

At convergence the total energy is

$$
E = \tfrac{1}{2}\,\mathrm{Tr}\!\left[P\,(H^{\text{core}} + F)\right] + E_{\text{nuc}}.
$$

The orbital energies are the eigenvalues $\varepsilon_i$ of $F$. The electron density is

$$
\rho(\mathbf{r}) = \sum_{\mu\nu} P_{\mu\nu}\, \chi_\mu(\mathbf{r})\chi_\nu(\mathbf{r}),
$$

which is continuous across sphere boundaries because the basis functions are (the NAO contributions vanish at the boundary, so $\rho$ is continuous in value even though its derivative is only approximately continuous). Dipole moments, population analyses, and other one-electron properties follow directly from $P$ and the corresponding one-electron operator matrices.

## 11. Forces and geometry optimisation

### 11.1 Hellmann-Feynman plus Pulay forces

The force on atom $A$ is the negative derivative of the energy with respect to its position. Because the NAOs are attached to atoms and move with them, and because the basis is not complete, the basis-set (Pulay) contribution is not zero. The total force has the form

$$
\mathbf{F}_A = -\frac{\partial E}{\partial \mathbf{R}_A} = -\mathrm{Tr}\!\left[P \frac{\partial H^{\text{core}}}{\partial \mathbf{R}_A}\right] - \mathrm{Tr}\!\left[\,\cdots\right] + \mathbf{F}_A^{\text{Pulay}},
$$

where the Pulay term accounts for the change of the basis functions with atomic position and includes overlap-derivative contributions weighted by the energy-weighted density matrix $W = \sum_i \varepsilon_i C_i C_i^{\dagger}$.

For NAOs, the derivative of a basis function with respect to its atom's position is the negative of its derivative with respect to the electron coordinate, moved to the sphere-centred frame. For Gaussians, the derivatives are analytic. The NAO derivative terms require numerical differentiation of the radial functions, or a closed-form expression for the radial derivative of the basis functions inside the sphere.

### 11.2 Geometry optimisation

Forces drive a geometry optimiser, typically a quasi-Newton scheme such as BFGS. When atoms move, the spheres move with them, and the NAO and Gaussian basis functions move accordingly. The sphere radii are usually held fixed. Care is needed to ensure that the spheres do not overlap as atoms move, which can be enforced by a constraint or by restricting the step size.

## 12. Periodic systems

If the system is a periodic solid, the basis must be expressed in a Bloch-compatible form. The NAOs in each cell are replicated by lattice translations, and the Fock and overlap matrices are built at each $k$-point. The Coulomb and exchange sums become sums over lattice vectors, with the long-range Coulomb interaction handled by an Ewald or similar method.

Two points are specific to the mixed basis:

- The sphere-boundary condition applies to every replicated NAO, so the replicas do not interact through overlapping spheres.
- The Gaussian set must be chosen so that the Bloch sums converge over the $k$-mesh. Diffuse Gaussians can cause convergence problems and may need to be removed or modified in the periodic setting.

The exchange term in a periodic system is the most expensive part and typically requires careful truncation or density fitting.

## 13. Density fitting and resolution of the identity

The two-electron integral cost can be reduced by density fitting (also called resolution of the identity). The product of two basis functions is expanded in an auxiliary basis $\{\gamma_P\}$:

$$
\chi_\mu(\mathbf{r})\chi_\nu(\mathbf{r}) \approx \sum_P c^P_{\mu\nu}\, \gamma_P(\mathbf{r}),
$$

and the four-index integrals are approximated as

$$
(\mu\nu|\lambda\sigma) \approx \sum_{PQ} (\mu\nu|P)\,[V^{-1}]_{PQ}\,(Q|\lambda\sigma),
$$

where $V_{PQ} = (P|Q)$ is the auxiliary metric. This reduces the cost from $N^4$ to $N^3$ in the number of basis functions.

Adapting density fitting to the mixed NAO-Gaussian basis requires an auxiliary basis that can represent the NAO-NAO and NAO-Gaussian products. Gaussian auxiliary functions are the natural choice for the interstitial part of the product. For the sphere-confined part, a numerical auxiliary basis on each sphere is needed. This adaptation is nontrivial, and the accuracy depends on the quality of the auxiliary set on each sphere.

## 14. Computational cost

The leading costs per SCF iteration are:

- Coulomb and exchange builds: formally $O(N^4)$ in the basis size without screening or fitting, reducible to $O(N^2)$ to $O(N^3)$ with effective screening and density fitting.
- Diagonalisation: $O(N^3)$ in the effective basis size.
- DIIS: $O(N_{\text{DIIS}}^2 N^2)$ for the subspace, which is small in practice.
- Integral evaluation (one-time, or per iteration for direct SCF): dominated by the NAO-NAO and NAO-Gaussian classes.

The NAO part of the basis has a locality advantage: each NAO is zero outside its sphere, so the number of nonzero overlap, kinetic, and nuclear attraction matrix elements grows linearly with the number of atoms in a large system. The Gaussian part has no such locality in principle, but its diffuse components can be cut off in practice, which restores locality at the cost of some accuracy.

## 15. Parallelisation

The method parallelises naturally over several dimensions:

- **Atoms.** Each sphere's NAO integrals are independent of other spheres, so they can be distributed over processes.
- **Integral batches.** Two-electron integrals can be partitioned over shell pairs, with each process building its share of $J$ and $K$, followed by a reduction.
- **Grid points.** The interstitial and sphere quadrature grids can be distributed, with the density and potential assembled through reductions.
- **$k$-points** (periodic case). Each $k$-point is an independent eigenproblem, once the integrals are built.

The main communication cost is the reduction of $J$ and $K$ and of the density across processes. Memory scaling is dominated by the storage of integrals or of the density and Fock matrices. Replicating the full matrices on each process is simple but limits system size; block distribution of matrices is necessary for large calculations.

## 16. Validation and testing

A mixed-basis method of this kind requires systematic validation:

- **Atomic limits.** For a single atom, the energy should converge to the exact Hartree-Fock value for that atom as the NAO and Gaussian sets are enlarged.
- **Diatomic molecules.** Standard test molecules (H$_2$, LiH, N$_2$, H$_2$O) can be compared against converged Gaussian-basis or plane-wave reference values.
- **Boundary behaviour.** The derivative mismatch at sphere boundaries should be measured directly, by evaluating the normal derivative of the full orbital on both sides of the boundary, and its decrease with basis size should be documented.
- **Basis-set convergence.** Energies and properties should be checked for stability as the Gaussian exponent set and the NAO radial basis are enlarged.
- **Conditioning.** The eigenvalue spectrum of $S$ and the energy sensitivity to the orthogonalisation threshold $\tau$ should be reported.
- **Gradient check.** Analytic forces should be compared against finite differences of the energy, to verify the Pulay terms.

## 17. Limitations and open questions

Several issues remain for any implementation of this approach.

**Derivative continuity is approximate.** As discussed in Section 4, the zero condition enforces value continuity exactly, but the derivative continuity is only variational. This is the main conceptual limitation of the boundary treatment, and it should be stated clearly whenever the method is described or compared to others.

**Sphere geometry.** Non-overlapping spheres leave gaps where the NAO description is absent. The Gaussians fill these gaps, but the quality of the description in the gaps depends entirely on the Gaussian set. Systems with unusual bonding geometries, such as those with very short or very long interatomic distances, may need careful choice of radii.

**Basis redundancy.** The union of NAOs and Gaussians is generally linearly dependent, especially where a Gaussian and an NAO represent similar functions near the sphere edge. Canonical orthogonalisation handles this, but the threshold choice affects accuracy and must be validated.

**Integral cost.** The numerical NAO-NAO and NAO-Gaussian integrals are the main computational bottleneck. Efficient implementation requires careful quadrature design and a good spherical-harmonic machinery.

**Density fitting adaptation.** Extending density fitting to the mixed basis is possible but nontrivial, and the auxiliary basis design is a significant implementation choice.

## 18. Summary

The method describes a Hartree-Fock calculation in a mixed basis of confined numerical atomic orbitals inside non-overlapping atomic spheres and s-type Gaussians covering the whole cell. The NAOs are numerical radial functions times spherical harmonics, constrained to vanish at their sphere boundaries. The Gaussians are smooth and nonzero everywhere. The zero condition guarantees exact continuity of the orbital value across each sphere boundary, while continuity of the normal derivative is achieved only approximately through the variational principle. The Hartree-Fock equations take the standard Roothaan-Hall form in this nonorthogonal basis, solved by a self-consistent field procedure with DIIS acceleration. The two-electron integrals are the dominant cost. They are evaluated by a combination of analytic Gaussian techniques and numerical sphere quadrature, and are made tractable by the locality of the NAOs, integral screening, and, where adapted, density fitting. Forces require Pulay corrections because the basis moves with the atoms.
