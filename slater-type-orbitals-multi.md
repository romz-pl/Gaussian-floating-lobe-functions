# Multi-centre integrals for S-type Slater-type orbitals

This note extends the one-centre algorithm to molecular settings, where the basis functions sit on several atomic centres and the potential is a sum of nuclear attractions. It assumes the one-centre formulas for S-type STOs, which are set out in the companion note on spherically symmetric STO integrals.

---

## 1. Setting and notation

Basis functions are centred at positions $\mathbf R_a$ and have the form

$$
\chi_a(\mathbf r) = N_a\, r_a^{\,k_a}\, e^{-\zeta_a r_a}\, Y_{00}, \qquad r_a = |\mathbf r - \mathbf R_a|, \qquad k_a = n_a - 1 .
$$

The potential is a sum of nuclear attractions,

$$
V(\mathbf r) = -\sum_{C} \frac{Z_C}{|\mathbf r - \mathbf R_C|},
$$

with nuclei $C$ at positions $\mathbf R_C$ and charges $Z_C$. A nucleus may coincide with a basis centre. The matrices to be built are

$$
S_{ab} = \langle a|b\rangle, \qquad T_{ab} = \langle a|-\tfrac12\nabla^2|b\rangle, \qquad V_{ab} = \langle a|V|b\rangle, \qquad H = T + V .
$$

Each matrix element is a sum over pairs of basis functions and, for $V$, over nuclei. The integrals fall into three classes.

1. **Two-centre** integrals, where both basis functions lie on different centres and the operator is either the identity, the Laplacian, or a nuclear attraction at one of the two centres.
2. **Three-centre** integrals, where the nuclear attraction sits on a third centre.
3. **Multi-centre** assembly, where all of the above are combined into the molecular matrices.

Two-centre integrals have a closed form in prolate spheroidal coordinates. Three-centre integrals do not, and need one of the methods in Section 5.

---

## 2. Two-centre integrals in prolate spheroidal coordinates

### 2.1 Coordinates

Let $R = |\mathbf R_a - \mathbf R_b|$ be the internuclear distance. Define prolate spheroidal coordinates $(\xi,\eta,\varphi)$ with foci at the two centres:

$$
\xi = \frac{r_a + r_b}{R} \in [1,\infty), \qquad \eta = \frac{r_a - r_b}{R} \in [-1,1], \qquad \varphi \in [0,2\pi).
$$

Inverting gives $r_a = \tfrac{R}{2}(\xi+\eta)$ and $r_b = \tfrac{R}{2}(\xi-\eta)$. The volume element is

$$
d^3\mathbf r = \frac{R^3}{8}\,(\xi^2-\eta^2)\,d\xi\,d\eta\,d\varphi .
$$

The exponential factor becomes separable:

$$
\zeta_a r_a + \zeta_b r_b = \rho\,\xi + \sigma\,\eta, \qquad \rho = \frac{R(\zeta_a+\zeta_b)}{2}, \qquad \sigma = \frac{R(\zeta_a-\zeta_b)}{2},
$$

with $\rho > |\sigma|$ for positive exponents.

### 2.2 Master integral

Any two-centre integral of the type

$$
\int r_a^{\alpha}\, r_b^{\beta}\, e^{-\zeta_a r_a - \zeta_b r_b}\, d^3\mathbf r
$$

can be written, using the Jacobian above and $(\xi^2-\eta^2) = (\xi+\eta)(\xi-\eta)$, as

$$
M(\alpha,\beta) = \frac{\pi R^3}{4}\left(\frac{R}{2}\right)^{\alpha+\beta} \int_1^\infty\!\!\int_{-1}^{1} (\xi+\eta)^{\alpha+1}\,(\xi-\eta)^{\beta+1}\, e^{-\rho\xi - \sigma\eta}\, d\eta\, d\xi .
$$

The integrand is a polynomial whenever $\alpha+1 \ge 0$ and $\beta+1 \ge 0$ are integers. Expanding the polynomial into monomials $\xi^i\eta^j$ reduces the double integral to a finite sum of products of two one-dimensional auxiliary functions, which are the Löwdin $\alpha$- and $\beta$-type functions:

$$
A_i(\rho) = \int_1^\infty \xi^i e^{-\rho\xi}\,d\xi, \qquad B_j(\sigma) = \int_{-1}^{1} \eta^j e^{-\sigma\eta}\,d\eta .
$$

Thus

$$
M(\alpha,\beta) = \frac{\pi R^3}{4}\left(\frac{R}{2}\right)^{\alpha+\beta}\sum_{i,j} c_{ij}\,A_i(\rho)\,B_j(\sigma),
$$

where $c_{ij}$ are the coefficients of $(\xi+\eta)^{\alpha+1}(\xi-\eta)^{\beta+1}$.

### 2.3 Evaluating the auxiliary functions

The auxiliary functions satisfy simple recurrences that follow from integration by parts.

For $A_i$, with $A_0(\rho) = e^{-\rho}/\rho$,

$$
A_i(\rho) = \frac{e^{-\rho}}{\rho} + \frac{i}{\rho}\,A_{i-1}(\rho).
$$

For $B_j$, with $B_0(\sigma) = (e^{\sigma}-e^{-\sigma})/\sigma$,

$$
B_j(\sigma) = \frac{(-1)^j e^{\sigma} - e^{-\sigma}}{\sigma} + \frac{j}{\sigma}\,B_{j-1}(\sigma).
$$

Implementation notes.

- The $A$ recurrence is a forward recurrence with positive terms and is stable.
- The $B$ recurrence suffers cancellation when $|\sigma|$ is small, which happens for nearly equal exponents. For $|\sigma|$ below a threshold of order $0.5$, evaluate $B_j$ from its Taylor series in $\sigma$, $B_j(\sigma) = \sum_{m\ge 0} \frac{(-\sigma)^m}{m!}\cdot\frac{1-(-1)^{m+j+1}}{m+j+1}$, truncated when terms fall below machine precision.
- When $\rho$ is large, $A_i$ involves $e^{-\rho}$ and the recurrence loses nothing, but the overall integrals become small. Scale the contribution with the prefactor $N_aN_b$ and the exponential, computed in log form, to avoid underflow.

### 2.4 Overlap

With $\alpha = k_a$ and $\beta = k_b$, the overlap integrand is $r_a^{k_a}r_b^{k_b}$ times the Jacobian, so

$$
S_{ab} = N_aN_b\, M(k_a,k_b)\quad\text{(with the Jacobian factor } \xi^2-\eta^2 \text{ already inside } M).
$$

For two 1s orbitals with equal exponents, $\sigma = 0$, $\rho = \zeta R$, and the result reduces to Mulliken's formula

$$
S = e^{-\rho}\left(1+\rho+\tfrac{\rho^2}{3}\right),
$$

which is the standard check for this part of the code.

### 2.5 Nuclear attraction to a centre of the pair

If the nucleus $C$ coincides with centre $b$, the integrand has $r_b^{-1}$, so $\beta = k_b - 1$. Similarly, a nucleus at $a$ corresponds to $\alpha = k_a - 1$. Therefore

$$
V_{ab}^{(b)} = -Z_b\, N_aN_b\, M(k_a,\,k_b-1), \qquad V_{ab}^{(a)} = -Z_a\, N_aN_b\, M(k_a-1,\,k_b).
$$

The polynomial condition is met because for $k_b = 0$ the exponent $\beta+1 = 0$ is still a valid non-negative power.

**Check.** As $R \to 0$ the two-centre attraction must tend to the one-centre value $-Z\,N_aN_b (s+1)!/p^{s+2}$ from the companion note.

### 2.6 Kinetic energy through the Laplacian

For a basis function centred at $b$, the Laplacian satisfies

$$
\nabla^2\!\left(r_b^{k} e^{-\zeta r_b}\right) = \left[\,\zeta^2 r_b^{k} - 2\zeta(k+1)\,r_b^{k-1} + k(k+1)\,r_b^{k-2}\,\right]e^{-\zeta r_b}.
$$

Using $T_{ab} = -\tfrac12\langle a|\nabla^2|b\rangle$ and expanding the bra and ket, every term is again a master integral $M(\alpha,\beta)$ with $\alpha = k_a$ and $\beta \in \{k_b, k_b-1, k_b-2\}$. The factor $k_b(k_b+1)$ removes the one case where $\beta+1$ would be negative, so the polynomial condition is always satisfied.

In practice the kinetic element is assembled as

$$
T_{ab} = -\tfrac12 N_aN_b\Big[\,\zeta_b^2\,M(k_a,k_b) - 2\zeta_b(k_b+1)\,M(k_a,k_b-1) + k_b(k_b+1)\,M(k_a,k_b-2)\,\Big]
$$

with the convention that a term whose coefficient vanishes (for $k_b = 0$ the last term) is simply omitted. The one-centre limit provides a check, as does Hermiticity: $T_{ab}$ must equal $T_{ba}$, which is a useful test because the formula is written with the Laplacian acting on $b$.

For 1s functions the same result is compactly $T_{ab} = -\tfrac12\zeta_b^2 S_{ab} + \zeta_b\,\langle a|r_b^{-1}|b\rangle$ with unit coefficient on the attraction-type integral.

---

## 3. Auxiliary expansion for general pairs

The method above requires both basis functions to be S-type and centred on the two foci. For basis functions with $l>0$, or for a molecular frame that is not aligned with the internuclear axis, the usual procedure is to rotate to a frame where the internuclear axis is the $z$-axis, then apply the same prolate construction with the angular factors handled by Gaunt coefficients. This extends the method to the full basis, but the algebra grows quickly with angular momentum and is normally done with generated code.

---

## 4. Three-centre attraction integrals

For an attraction to a nucleus $C$ that is not one of the two basis centres, the prolate construction does not apply directly, because $r_c$ is not separable in $(\xi,\eta,\varphi)$. Four families of method are in use.

### 4.1 Expansion methods (Löwdin, Ruedenberg, Guseinov)

The function $1/r_c$, or equivalently the STO on centre $b$ expressed about centre $a$, is expanded in spherical harmonics about a convenient origin using the addition theorem. The STO on centre $b$ is re-expressed as a series of functions on centre $a$ (translation formulas), and the three-centre integral then becomes a sum of two-centre integrals of the form already discussed.

Algorithm:

1. Choose an origin, usually the centre of the bra function or the nucleus.
2. Expand each off-centre factor $r_x^{\,k}e^{-\zeta r_x}$ in solid harmonics about the origin, using the Guseinov or Löwdin translation coefficients.
3. Use Gaunt coefficients to restrict the angular sums.
4. Evaluate the remaining radial integrals in closed form via the Laplace master integral or via the $\alpha$-type auxiliary functions.

This is the standard analytic route for STO codes. Its accuracy is limited by the truncation of the angular-momentum sums, which is controlled by a chosen maximum $L$ and must be tested for convergence.

### 4.2 Fourier-space method

The Fourier transform of $1/r$ is $4\pi/k^2$, and that of an S-type STO is $8\pi\zeta/(k^2+\zeta^2)^2$. Three-centre attraction integrals can be written as a momentum-space convolution with phase factors $e^{i\mathbf k\cdot(\mathbf R_a - \mathbf R_C)}$ and similar terms. The phases are absorbed by the Fourier shift theorem, and the angular integral reduces to a Bessel-function average.

Algorithm:

1. Write the three-centre integral as $\frac{1}{(2\pi)^3}\int d^3\mathbf k\, \tilde f_a(\mathbf k)^{\ast}\tilde V(\mathbf k)\tilde f_b(\mathbf k)\,e^{i\mathbf k\cdot(\mathbf R_b-\mathbf R_a)}$, with $\ast$ denoting complex conjugation.
2. Use the spherical symmetry of the S-type transforms to reduce the angular integration to $\sin(kR)/(kR)$, where $R$ is the separation of the two relevant centres.
3. Integrate over $k$ by a one-dimensional quadrature, typically Gauss–Laguerre or a mapped rule on $[0,\infty)$ with the decay scale set by $\zeta$.

This method is robust and does not need angular truncation, but it becomes costly when many centres and high accuracy are required.

### 4.3 Gaussian-fitting method

Each STO is approximated by a linear combination of Gaussian primitives, $\chi \approx \sum_k d_k\, g_k$, with $g_k \propto e^{-\alpha_k r^2}$. This is the STO-$nG$ approach. Every integral then reduces to Gaussian integrals, which have closed forms for overlap, kinetic energy and nuclear attraction, so the three-centre integral becomes a sum over Gaussian pairs. The formulas for a pair of primitives $g_a$ at $\mathbf R_A$ with exponent $\alpha$ and $g_b$ at $\mathbf R_B$ with exponent $\beta$ are as follows.

Let $p = \alpha+\beta$, $\mu = \alpha\beta/p$, $\mathbf P = (\alpha\mathbf R_A+\beta\mathbf R_B)/p$, and $K_{AB} = \exp(-\mu R_{AB}^2)$. Then

$$
S = \left(\frac{\pi}{p}\right)^{3/2} K_{AB}, \qquad T = \mu\,(3-2\mu R_{AB}^2)\left(\frac{\pi}{p}\right)^{3/2}K_{AB},
$$

and the nuclear attraction to a point $\mathbf R_C$ with charge $Z_C$ is

$$
V = -Z_C\,\frac{2\pi}{p}\,K_{AB}\,F_0\!\left(p\,|\mathbf P-\mathbf R_C|^2\right),
$$

where $F_0$ is the Boys function

$$
F_m(T) = \int_0^1 t^{2m}\,e^{-Tt^2}\,dt .
$$

$F_0$ is evaluated with the error function, $F_0(T) = \tfrac12\sqrt{\pi/T}\,\mathrm{erf}(\sqrt T)$, or from a precomputed table with Taylor interpolation for small $T$. For a contracted STO, the matrix element is the sum over primitive pairs weighted by the contraction coefficients. The method is exact for the fitted functions and its accuracy is controlled by the fit quality, which should be checked against the STO itself.

Points on this method:

- The Gaussian product rule makes the cost scale with the number of primitive pairs, which is usually modest.
- The Boys function evaluation dominates the run time for attraction integrals.
- Large $R$ produces tiny $K_{AB}$, so screening can skip pairs entirely.

### 4.4 Numerical multi-centre quadrature

When the analytic routes are not wanted, integrals can be evaluated directly on a grid built from atom-centred pieces. Becke partitioning is the standard choice.

Each centre $A$ gets a weight function

$$
w_A(\mathbf r) = \frac{P_A(\mathbf r)}{\sum_{B} P_B(\mathbf r)}, \qquad P_A(\mathbf r) = \prod_{B\neq A} s(\mu_{AB}), \qquad \mu_{AB} = \frac{r_A - r_B}{R_{AB}},
$$

with $s(\mu) = \tfrac12\left[1 - f(f(f(\mu)))\right]$ and $f(\mu) = \tfrac32\mu - \tfrac12\mu^3$. The weights sum to one at every point, so the molecular integral becomes a sum over atoms of atomic integrals:

$$
\int F(\mathbf r)\,d^3\mathbf r = \sum_{A}\sum_{i} w_A(\mathbf r_{Ai})\,F(\mathbf r_{Ai})\,\omega_{Ai}.
$$

Each atomic grid is a product of a radial rule and an angular rule. The radial rule is a mapped Gauss–Chebyshev or Gauss–Laguerre scheme with $r^2$ included in the weight, and the angular rule is a Lebedev grid. The standard refinement (Stratmann, Scuseria and Frisch) replaces the Becke smoothing with a cell-function that reduces grid artefacts.

Numerical quadrature is simple and general, handles any potential and any basis, and is the natural choice for exchange-correlation-type work. Its accuracy is limited by grid density and it is slower than analytic methods for a fixed accuracy target.

---

## 5. Assembling molecular matrices

The overall procedure loops over pairs of basis functions and, for the potential, over nuclei.

```
input : basis {(centre R_i, n_i, zeta_i)}, nuclei {(R_C, Z_C)}, method choices, tolerances
output: S, T, V (symmetric)

build neighbour list: for each pair (i,j) compute R_ij
screen: if exp(-min(zeta_i, zeta_j) * R_ij) < tol_screen, set S, T, V entries to zero
          (skip the pair entirely, subject to a check that the cutoff is safe)

for each retained pair (i,j), i >= j:
  if R_ij == 0:
    use one-centre closed forms (companion note)
  else:
    S[i][j] = N_i N_j M(k_i, k_j)                       # Section 2.4
    T[i][j] = kinetic via Laplacian identity             # Section 2.6
    for each nucleus C:
      if C == centre j:  V[i][j] += -Z_C N_i N_j M(k_i, k_j - 1)
      if C == centre i:  V[i][j] += -Z_C N_i N_j M(k_i - 1, k_j)
      else:              V[i][j] += three_centre_attraction(i, j, C)   # Section 4
  S[j][i] = S[i][j]; T[j][i] = T[i][j]; V[j][i] = V[i][j]

H = T + V
```

Points on assembly:

- Screening is essential for large molecules. STO overlap decays as $e^{-\zeta R}$, so distant pairs contribute negligibly. A cutoff on $\zeta R$ is safe for overlap; for the attraction term, screening must account for the fact that the nucleus is independent of the pair distance.
- Three-centre integrals are the main cost. Each is independent, which makes them natural for parallel execution over pairs and nuclei.
- Symmetry (point-group or translation) reduces the number of distinct integrals. For molecules with equivalent centres, many pairs are related by symmetry and need computing only once.
- Caching is effective when several basis functions share the same centre and exponent, since the auxiliary functions $A_i$ and $B_j$ depend only on $\rho$ and $\sigma$, which are shared by all pairs with the same exponent pair and separation.
- Numerical stability of the auxiliary functions is handled at the level of $\sigma$ (Section 2.3), so the assembly loop itself does not need special cases.

After assembly, the generalized eigenproblem $H c = E S c$ is solved exactly as in the one-centre case, using Cholesky factorization or canonical orthogonalization. Because molecular bases with many centres can be strongly linearly dependent, canonical orthogonalization with a threshold on the eigenvalues of $S$ is the safer default.

---

## 6. Validation

1. **One-centre limit.** Set $R \to 0$ for two centres with equal exponents and check that the two-centre formulas reproduce the one-centre closed forms.
2. **Mulliken overlap.** For two 1s STOs with equal exponents, check $S = e^{-\rho}(1+\rho+\rho^2/3)$.
3. **Hermiticity.** $T_{ab} = T_{ba}$ and $V_{ab} = V_{ba}$ for every pair, which tests the Laplacian-based kinetic formula.
4. **Method agreement.** The analytic two-centre results, the Gaussian-fitted results, and the Becke quadrature results must agree within their stated tolerances for a set of test geometries.
5. **Angular-momentum convergence.** For the expansion methods, the three-centre result should converge as the truncation $L_{\max}$ increases.
6. **Invariance.** Rigidly translating or rotating the whole molecule, with the basis and nuclei moved together, must leave all eigenvalues unchanged to machine precision.
7. **H$_2^+$ benchmark.** A one-electron, two-proton system with a 1s STO basis on each proton is a standard test. The LCAO energy curve should show a minimum at an internuclear distance near 2.5 bohr, in agreement with published LCAO-STO values for this molecule, and should approach the one-centre hydrogen-atom energy as $R$ grows.
8. **Variational bound.** Adding basis functions on any centre must never raise the lowest eigenvalue.

---

## Summary

1. Classify each integral as two-centre or three-centre, and screen distant pairs.
2. Compute two-centre overlap, kinetic and nuclear-attraction elements in closed form using prolate spheroidal coordinates, the master integral $M(\alpha,\beta)$, and the Löwdin auxiliary functions $A_i$ and $B_j$.
3. For nuclei not on either centre of a pair, use an expansion method, a Fourier-space integral, Gaussian fitting with the Boys function, or Becke-partitioned numerical quadrature, depending on the accuracy and cost target.
4. Assemble $S$, $T$ and $V$ over all pairs and nuclei, form $H = T + V$, and solve the generalized eigenproblem with canonical orthogonalization.
5. Validate with the one-centre limit, Mulliken's overlap, Hermiticity, method cross-checks, rotational and translational invariance, and the H$_2^+$ benchmark.
