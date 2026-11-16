# Matrix elements of spherically symmetric Slater-type orbitals

This note describes an algorithm for evaluating the overlap, kinetic-energy and potential-energy matrix elements of spherically symmetric (S-type, $l=0$) Slater-type orbitals (STOs) for the one-electron Hamiltonian

$$
\hat H = -\tfrac{1}{2}\nabla^2 + V(\mathbf r),
$$

and for assembling and solving the generalized eigenvalue problem $\mathbf H\mathbf c = E\,\mathbf S\,\mathbf c$. All orbitals are taken as real-valued and centred on a single origin, which is the standard one-centre setting.

---

## 1. Basis functions

An STO with principal quantum number $n$ and orbital exponent $\zeta > 0$ is

$$
\chi_{n\zeta}(\mathbf r) = N_{n\zeta}\, r^{\,n-1}\, e^{-\zeta r}\, Y_{00}(\theta,\varphi), \qquad Y_{00} = \frac{1}{\sqrt{4\pi}}.
$$

For S-type functions only the radial part carries information, and it is convenient to introduce the shorthand

$$
k \equiv n-1 \ge 0 .
$$

The normalization constant follows from the radial integral $\int_0^\infty r^{2n}e^{-2\zeta r}\,dr = (2n)!/(2\zeta)^{2n+1}$:

$$
N_{n\zeta} = \frac{(2\zeta)^{\,n+\frac12}}{\sqrt{(2n)!}} .
$$

For the 1s STO this reduces to $N = 2\zeta^{3/2}$, which gives the familiar $\psi_{1s} = (\zeta^3/\pi)^{1/2}e^{-\zeta r}$.

Because $(2n)!$ and $(2\zeta)^{n}$ overflow quickly for large $n$ or $\zeta$, the normalization should be evaluated in logarithmic form:

$$
\ln N_{n\zeta} = \left(n+\tfrac12\right)\ln(2\zeta) - \tfrac12\,\ln\Gamma(2n+1),
$$

where $\ln\Gamma(2n+1)$ is the log-gamma function (`lgamma` in C/C++ `<cmath>`).

---

## 2. Reduction of all integrals to one master form

Consider two S-type STOs $a = (n_a,\zeta_a)$ and $b = (n_b,\zeta_b)$. Define

$$
p = \zeta_a + \zeta_b, \qquad k_a = n_a - 1, \qquad k_b = n_b - 1, \qquad s = k_a + k_b = n_a + n_b - 2,
$$

and the combined normalization $N_{ab} = N_a N_b$. Every radial integrand that appears has the form $r^{q}e^{-pr}$ multiplied by a polynomial prefactor, and the only integral needed is the Laplace-type master integral

$$
J_q(p) = \int_0^\infty r^{q}\, e^{-pr}\, dr = \frac{\Gamma(q+1)}{p^{\,q+1}} = \frac{q!}{p^{\,q+1}}, \qquad \text{Re}\,p > 0,\ q > -1 .
$$

The volume element for an S-type integrand is $d^3\mathbf r = r^2\,dr\,d\Omega$, and the angular integral of $Y_{00}^2$ equals one. Therefore all radial integrals carry an extra factor $r^2$.

---

## 3. Overlap integral

The overlap is

$$
S_{ab} = \int \chi_a\,\chi_b\, d^3\mathbf r = N_{ab}\int_0^\infty r^{\,k_a+k_b}\, e^{-pr}\, r^2\,dr = N_{ab}\, J_{s+2}(p).
$$

Hence

$$
\boxed{\,S_{ab} = N_{ab}\,\frac{(s+2)!}{p^{\,s+3}} = N_{ab}\,\frac{(n_a+n_b)!}{(\zeta_a+\zeta_b)^{\,n_a+n_b+1}}\,}
$$

**Check.** For $a=b$ with $n=1$, $\zeta$ equal: $N^2 = 4\zeta^3$, $p=2\zeta$, and $S = 4\zeta^3\cdot 2/(8\zeta^3) = 1$.

---

## 4. Kinetic-energy integral

Because the orbitals are real and the operator is Hermitian, integration by parts moves the Laplacian onto the gradient:

$$
\langle a|-\tfrac12\nabla^2|b\rangle = \tfrac12\int \nabla\chi_a\cdot\nabla\chi_b\, d^3\mathbf r = \tfrac12 \int_0^\infty R_a'(r)\,R_b'(r)\, r^2\,dr,
$$

where the boundary terms vanish because each STO is finite at the origin and decays exponentially. Here $R(r)=N r^{k} e^{-\zeta r}$ is the radial part, and

$$
R_a' = N_a\left(k_a r^{k_a-1} - \zeta_a r^{k_a}\right)e^{-\zeta_a r}.
$$

Multiplying the two derivatives gives three terms, with powers $r^{s-2}$, $r^{s-1}$ and $r^{s}$. After multiplication by $r^2$ each becomes a master integral. The result is

$$
\boxed{\;
T_{ab} = \tfrac12 N_{ab}\left[\,k_a k_b\,\frac{s!}{p^{\,s+1}} \;-\; (k_a\zeta_b + k_b\zeta_a)\,\frac{(s+1)!}{p^{\,s+2}} \;+\; \zeta_a\zeta_b\,\frac{(s+2)!}{p^{\,s+3}}\,\right]\;}
$$

**Checks.**

- For $n=1$ (so $k=0$, $s=0$) the first two terms vanish and $T = N^2\zeta_a\zeta_b/p^3$. For identical 1s orbitals this gives $T = \zeta^2/2$, the known result.
- The first term has a prefactor $k_ak_b$, so it vanishes automatically for 1s functions. Its $s!$ is never evaluated with a negative argument because $s \ge 0$.

---

## 5. Potential-energy integrals

### 5.1 Coulomb potential $V = -Z/r$

The integrand reduces to $r^{s+1}e^{-pr}$, so

$$
\boxed{\,V_{ab} = -Z\,N_{ab}\,\frac{(s+1)!}{p^{\,s+2}}\,}
$$

**Check.** For a single 1s orbital, $V = -Z\zeta$ in Hartree units. The variational energy of the normalized 1s STO is therefore $E(\zeta) = \zeta^2/2 - Z\zeta$, which is minimized at $\zeta = Z$ with $E = -Z^2/2$, the exact hydrogen-like ground-state energy.

### 5.2 Power-law potentials $V = c\, r^{-m}$ (any real $m$)

The same reduction works with a shifted factorial:

$$
V_{ab} = c\,N_{ab}\,\frac{\Gamma(s+3-m)}{p^{\,s+3-m}},
$$

which is finite when $s+3-m > 0$. For $m$ non-integer or $m > s+3$ the Gamma function handles the analytic continuation directly.

### 5.3 Sums of exponential-type terms

If the potential can be written as a finite sum $V(r) = \sum_j c_j\, r^{q_j} e^{-\lambda_j r}$ (Yukawa, screened Coulomb, Slater-type potentials, and so on), then each term is again a master integral with the exponent shifted:

$$
V_{ab} = N_{ab}\sum_j c_j\, \frac{(s+2+q_j)!}{(p+\lambda_j)^{\,s+3+q_j}}.
$$

This gives exact results with no quadrature error.

### 5.4 General radial potential $V(r)$

For an arbitrary smooth function the potential integral is

$$
V_{ab} = N_{ab}\int_0^\infty r^{\,s+2}\,V(r)\,e^{-pr}\,dr.
$$

The natural quadrature is Gauss–Laguerre. Substituting $x = pr$ gives

$$
V_{ab} = N_{ab}\,p^{-(s+3)}\int_0^\infty e^{-x}\,\big[\,x^{\,s+2}\,V(x/p)\,\big]\,dx \;\approx\; N_{ab}\,p^{-(s+3)}\sum_{i=1}^{N_q} w_i\,x_i^{\,s+2}\,V(x_i/p),
$$

where $\{x_i, w_i\}$ are the nodes and weights of the $N_q$-point Gauss–Laguerre rule for the weight $e^{-x}$.

Points on this method:

- The rule is exact for polynomial integrands of degree up to $2N_q - 1$. The Coulomb case is therefore reproduced exactly when $s+1 \le 2N_q - 1$, which serves as a built-in test.
- The exponential decay of the STO is absorbed by the weight, so no truncation of the radial domain is needed.
- If $V$ has an integrable singularity $r^{-\mu}$ at the origin, use generalized Gauss–Laguerre with weight $x^{\alpha}e^{-x}$ and $\alpha = s+2-\mu$ (this requires $\alpha > -1$), or subtract the singular part and treat it analytically with Section 5.2.
- The number of nodes needed grows with $s$ and with how rapidly $V$ varies on the length scale $1/\zeta$. In practice, increase $N_q$ until the matrix element changes by less than a chosen tolerance, typically $10^{-12}$ relative.
- For potentials with structure at large $r$ (long-range tails) or highly oscillatory character, a composite rule (split $[0,\infty)$ into segments) or tanh–sinh/adaptive quadrature may be preferable.

### 5.5 Non-spherical potentials

If $V$ depends on angles, expand it in spherical harmonics, $V(\mathbf r) = \sum_{LM} V_{LM}(r)\,Y_{LM}(\theta,\varphi)$. Because $\int Y_{00}\,Y_{LM}\,Y_{00}\,d\Omega = \delta_{L0}\,\delta_{M0}/\sqrt{4\pi}$, only the spherical average $V_{00}(r)$ contributes to S-type matrix elements. The procedure above then applies with $V(r)$ replaced by the angular average $\bar V(r) = \frac{1}{4\pi}\int V\,d\Omega$. This is a useful simplification for the $ss$ block of the Hamiltonian.

---

## 6. Assembling the matrices

For a basis of $M$ functions $\{(n_i,\zeta_i)\}_{i=1}^{M}$ the matrices $\mathbf S$, $\mathbf T$ and $\mathbf V$ are symmetric, so only the lower triangle needs computing. The full Hamiltonian is $\mathbf H = \mathbf T + \mathbf V$.

Pseudocode:

```
input : basis {(n_i, zeta_i)}, potential V(r) or its parameters, quadrature order N_q
output: S, T, V (symmetric M x M)

for i = 0 .. M-1
  for j = 0 .. i
    k_i = n_i - 1;  k_j = n_j - 1
    s   = k_i + k_j
    p   = zeta_i + zeta_j
    logN = logNorm(n_i, zeta_i) + logNorm(n_j, zeta_j)

    # all factorial/power combinations evaluated in log space
    S[i][j] = exp( logN + lgamma(s+3) - (s+3) log p )

    T[i][j] = 0.5 * ( k_i*k_j * exp(logN + lgamma(s+1)   - (s+1) log p)
                    - (k_i*zeta_j + k_j*zeta_i) * exp(logN + lgamma(s+2) - (s+2) log p)
                    +  zeta_i*zeta_j * exp(logN + lgamma(s+3) - (s+3) log p) )

    V[i][j] = potential_integral(s, p, logN, V)   # Sections 5.1 to 5.5

    S[j][i] = S[i][j];  T[j][i] = T[i][j];  V[j][i] = V[i][j]

return S, T, V
```

Two implementation details matter for accuracy.

- The kinetic expression contains a difference of positive terms. When orbitals have similar exponents and high $n$, the terms can nearly cancel. Use extended precision or a compensated sum when the magnitude of the result is much smaller than the magnitude of the individual terms.
- Working in log space keeps every intermediate value in range. The subtraction is then performed on exponentials of moderate size.

---

## 7. Solving the generalized eigenproblem

The energies and coefficients satisfy $\mathbf H\mathbf c = E\,\mathbf S\,\mathbf c$. Two standard routes exist.

**Cholesky route.** Factor $\mathbf S = \mathbf L\mathbf L^{\mathsf T}$, form $\tilde{\mathbf H} = \mathbf L^{-1}\mathbf H\,\mathbf L^{-\mathsf T}$, diagonalize it with a symmetric eigensolver, and back-transform the eigenvectors as $\mathbf c = \mathbf L^{-\mathsf T}\tilde{\mathbf c}$. This is fast and stable when $\mathbf S$ is well conditioned.

**Canonical orthogonalization.** STO bases with nearly equal exponents or many high-$n$ functions produce nearly linearly dependent sets, and $\mathbf S$ becomes almost singular. Diagonalize $\mathbf S = \mathbf U\,\mathrm{diag}(\sigma)\,\mathbf U^{\mathsf T}$, discard eigenvalues below a threshold $\tau$ (typically $10^{-8}$ to $10^{-10}$), form $\mathbf X = \mathbf U_{\text{kept}}\,\mathrm{diag}(\sigma_{\text{kept}}^{-1/2})$, and solve the reduced problem $\mathbf X^{\mathsf T}\mathbf H\mathbf X\,\mathbf d = E\,\mathbf d$.

Canonical orthogonalization is the safer default for basis sets that are chosen automatically.

---

## 8. Validation procedure

1. **Normalization.** Each diagonal overlap must equal $1$ to machine precision. This tests $N_{n\zeta}$ and the master integral.
2. **Symmetry and positivity.** $\mathbf S$ must be symmetric and positive definite, and $\mathbf T$ and $\mathbf V$ must be symmetric.
3. **Analytic benchmark.** For $V = -Z/r$ with a single 1s STO, the energy must satisfy $E(\zeta) = \zeta^2/2 - Z\zeta$, with minimum $-Z^2/2$ at $\zeta = Z$.
4. **Quadrature convergence.** For a general $V$, the matrix elements should converge as $N_q$ increases. Compare against the closed-form result of Section 5.3 for a Yukawa or exponential test potential.
5. **Variational bound.** Enlarging the basis (adding functions or extra exponents) must never raise the lowest eigenvalue.
6. **Cross-check with the radial ODE.** For a central potential, the lowest eigenvalue of the STO basis should approach the value obtained from a high-accuracy numerical solution of the radial Schrödinger equation as the basis grows.

---

## 9. Alternative and complementary approaches

- **Gaussian expansion of STOs (STO-$nG$).** Each STO is fitted by a linear combination of Gaussians, $e^{-\zeta r} \approx \sum_k d_k e^{-\alpha_k r^2}$. All integrals then become analytic through the Gaussian product rule. This is the traditional route in quantum chemistry and avoids numerical quadrature for Gaussian-friendly potentials, at the cost of a fitting error that is controlled by the number of Gaussians.
- **Direct radial discretization.** Finite-difference, finite-element or spectral methods on $[0,R_{\max}]$ avoid basis-set issues altogether, but lose the compact analytic structure that the STO approach provides.
- **Non-origin potentials.** If the potential is centred elsewhere, or if several nuclei are present, the origin-centred formulas no longer apply directly. One route uses the Barnett–Coulson or Guseinov translation formulas, which express an STO on one centre in terms of functions on another. A second route expands the potential in a multipole series about the basis centre and combines it with the Gaunt coefficients.

---

## 10. Extensions to $l > 0$

The framework extends to higher angular momentum with three changes.

- The angular part of each basis function is $Y_{lm}$, and the radial part is $N\,r^{n-1}e^{-\zeta r}$ with the same master-integral structure.
- The kinetic operator gains the centrifugal term $l(l+1)/r^2$, which adds one more master integral $J_{s}$ to the kinetic expression.
- Cross-matrix elements between different $l$ values are scaled by Gaunt coefficients $\langle Y_{l_1 m_1}|Y_{LM}|Y_{l_2 m_2}\rangle$, and the potential contributes only through the multipole components $V_{LM}(r)$ that satisfy the triangle condition.

---

## Summary of the algorithm

1. Choose the basis $\{(n_i,\zeta_i)\}$ and compute the log-normalizations.
2. For each pair, compute $p = \zeta_a+\zeta_b$ and $s = n_a+n_b-2$.
3. Evaluate $S_{ab}$ and $T_{ab}$ from the closed forms in Sections 3 and 4.
4. Evaluate $V_{ab}$ in closed form when the potential is Coulombic, power-law, or a finite exponential sum; otherwise use Gauss–Laguerre quadrature (Section 5.4) with convergence testing.
5. Form $\mathbf H = \mathbf T + \mathbf V$, then solve $\mathbf H\mathbf c = E\,\mathbf S\,\mathbf c$ by Cholesky or canonical orthogonalization.
6. Validate with the normalization, symmetry, hydrogenic and variational checks of Section 8.
