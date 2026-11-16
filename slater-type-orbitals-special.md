# Multi-centre integrals for pure exponential S-type orbitals

This note specialises the multi-centre algorithms to the simplest STO, the pure exponential orbital with no radial power prefactor,

$$
\chi_a(\mathbf r) = N_a\, e^{-\zeta_a r_a}, \qquad r_a = |\mathbf r - \mathbf R_a| .
$$

This is the $n=1$ case, so $k_a = n_a - 1 = 0$ throughout. The normalization for a spatially normalized orbital is

$$
N_a = \left(\frac{\zeta_a^3}{\pi}\right)^{1/2}, \qquad \int |\chi_a|^2\, d^3\mathbf r = N_a^2\,\frac{\pi}{\zeta_a^3} = 1 .
$$

If the spherical harmonic $Y_{00}$ is kept explicitly, $N_a$ changes by the factor $\sqrt{4\pi}$. The formulas below use the spatially normalized convention.

The $k=0$ simplification removes every power-of-$r$ factor. The kinetic operator acting on the ket becomes especially simple:

$$
\nabla^2 e^{-\zeta r} = \left(\zeta^2 - \frac{2\zeta}{r}\right)e^{-\zeta r}.
$$

All three matrix types then reduce to two kinds of building block: overlap-type integrals and single-centre-weighted integrals $\langle a|1/r_x|b\rangle$. The first is available in closed form. The second is also closed form for the two-centre case, because the nucleus sits at one of the two basis centres. A third centre needs a separate treatment, which Section 4 gives.

---

## 1. Notation and the two-centre setting

Let $R = |\mathbf R_a - \mathbf R_b|$, and define

$$
\rho = \frac{R(\zeta_a+\zeta_b)}{2}, \qquad \sigma = \frac{R(\zeta_a-\zeta_b)}{2}, \qquad \rho > |\sigma| .
$$

Use prolate spheroidal coordinates $\xi = (r_a+r_b)/R$ and $\eta = (r_a-r_b)/R$. The volume element is $d^3\mathbf r = \tfrac{R^3}{8}(\xi^2-\eta^2)\,d\xi\,d\eta\,d\varphi$, and the exponent separates as $\zeta_a r_a + \zeta_b r_b = \rho\xi + \sigma\eta$.

The auxiliary functions are

$$
A_0(\rho) = \int_1^\infty e^{-\rho\xi}d\xi = \frac{e^{-\rho}}{\rho}, \qquad A_1(\rho) = \int_1^\infty \xi\, e^{-\rho\xi}d\xi = e^{-\rho}\left(\frac1\rho+\frac1{\rho^2}\right),
$$

$$
B_0(\sigma) = \int_{-1}^1 e^{-\sigma\eta}d\eta = \frac{2\sinh\sigma}{\sigma}, \qquad B_1(\sigma) = \int_{-1}^1 \eta\, e^{-\sigma\eta}d\eta = \frac{2\sinh\sigma}{\sigma^2} - \frac{2\cosh\sigma}{\sigma},
$$

and

$$
A_2(\rho) = e^{-\rho}\left(\frac1\rho+\frac2{\rho^2}+\frac2{\rho^3}\right), \qquad B_2(\sigma) = \frac{2\sinh\sigma}{\sigma} - \frac{4\cosh\sigma}{\sigma^2} + \frac{4\sinh\sigma}{\sigma^3}.
$$

For $|\sigma|$ small, evaluate $B_1$ and $B_2$ from their Taylor series rather than from these expressions, because the terms cancel. The series is

$$
B_j(\sigma) = \sum_{m\ge 0}\frac{(-\sigma)^m}{m!}\,\frac{1-(-1)^{m+j+1}}{m+j+1}.
$$

---

## 2. Closed-form two-centre integrals

### 2.1 Overlap

The overlap integrand in prolate coordinates is $(\xi^2-\eta^2)$ times the exponential, which gives

$$
S_{ab} = N_aN_b\,\frac{\pi R^3}{4}\left[A_2(\rho)B_0(\sigma) - A_0(\rho)B_2(\sigma)\right].
$$

Substituting the auxiliary functions, the result simplifies to

$$
\boxed{\;S_{ab} = N_aN_b\,\pi R^3\,e^{-\rho}\left[\left(\frac1{\rho^2}+\frac1{\rho^3}\right)\frac{\sinh\sigma}{\sigma} + \frac1\rho\left(\frac{\cosh\sigma}{\sigma^2}-\frac{\sinh\sigma}{\sigma^3}\right)\right]\;}
$$

The second term in the bracket is finite as $\sigma\to 0$, because $\cosh\sigma/\sigma^2 - \sinh\sigma/\sigma^3 \to 1/3$.

**Check, equal exponents.** At $\sigma = 0$ the bracket tends to $\tfrac{1}{\rho^2}+\tfrac{1}{\rho^3}+\tfrac{1}{3\rho}$, and with $N^2 = \zeta^3/\pi$, $\rho = \zeta R$, this gives $S = e^{-\rho}(1+\rho+\rho^2/3)$, Mulliken's result.

### 2.2 Nuclear attraction to a basis centre

The attraction integral $J^{(B)}_{ab} = \langle a|1/r_b|b\rangle$ has the integrand $r_b^{-1}$, which corresponds to $\beta = -1$ in the master integral. In prolate coordinates it is

$$
J^{(B)}_{ab} = N_aN_b\,\frac{\pi R^2}{2}\left[A_1(\rho)B_0(\sigma) + A_0(\rho)B_1(\sigma)\right].
$$

Substituting the explicit functions,

$$
\boxed{\;J^{(B)}_{ab} = N_aN_b\,\pi R^2\, e^{-\rho}\left[\left(\frac1\rho+\frac1{\rho^2}\right)\frac{\sinh\sigma}{\sigma} + \frac1\rho\left(\frac{\sinh\sigma}{\sigma^2}-\frac{\cosh\sigma}{\sigma}\right)\right]\;}
$$

The attraction to the nucleus on centre $b$ is $V^{(B)}_{ab} = -Z_B\,J^{(B)}_{ab}$.

The attraction to the nucleus on centre $a$ is obtained by exchanging the roles of the two centres. This changes $\sigma \to -\sigma$, which flips the sign of the odd-in-$\sigma$ bracket term:

$$
J^{(A)}_{ab} = N_aN_b\,\pi R^2\, e^{-\rho}\left[\left(\frac1\rho+\frac1{\rho^2}\right)\frac{\sinh\sigma}{\sigma} - \frac1\rho\left(\frac{\sinh\sigma}{\sigma^2}-\frac{\cosh\sigma}{\sigma}\right)\right].
$$

**Check, equal exponents.** At $\sigma = 0$ the second term vanishes, and with $N^2 = \zeta^3/\pi$ and $\rho = \zeta R$,

$$
J^{(B)}_{ab} = \zeta\,(1+\zeta R)\,e^{-\zeta R}.
$$

**Check, one-centre limit.** As $R\to 0$ the value tends to $N^2\pi/\zeta^2 = \zeta$, the one-centre result $\langle 1/r\rangle = \zeta$ for a normalized pure exponential.

### 2.3 Kinetic energy

Because $\nabla^2\chi_b = (\zeta_b^2 - 2\zeta_b/r_b)\chi_b$, the kinetic matrix element is

$$
T_{ab} = -\tfrac12\langle a|\nabla^2|b\rangle = -\tfrac12\zeta_b^2\,S_{ab} + \zeta_b\,J^{(B)}_{ab}.
$$

Hermiticity allows the Laplacian to act on the bra instead, which gives the same value:

$$
T_{ab} = -\tfrac12\zeta_a^2\,S_{ab} + \zeta_a\,J^{(A)}_{ab}.
$$

The two forms are equal by the Hermiticity of the kinetic operator, so the difference between them is a useful numerical test. In a code, the average of the two is a convenient symmetric choice.

**Check, equal exponents.** Substituting the closed forms gives

$$
T_{ab} = \frac{\zeta^2}{2}\,e^{-\zeta R}\left(1+\zeta R-\frac{\zeta^2R^2}{3}\right),
$$

which reduces to $\zeta^2/2$ at $R = 0$, the one-centre value.

### 2.4 Summary of the closed-form two-centre set

The complete two-centre matrix for a pair on centres $a$ and $b$ is therefore

$$
S_{ab},\qquad T_{ab} = -\tfrac12\zeta_b^2S_{ab}+\zeta_bJ^{(B)}_{ab},\qquad V_{ab} = -Z_A J^{(A)}_{ab} - Z_B J^{(B)}_{ab},
$$

where the last term includes every nucleus that happens to sit on one of the two basis centres. Each of $S$, $J^{(A)}$ and $J^{(B)}$ depends only on $\rho$, $\sigma$, $R$ and the normalizations, so these can be cached for pairs that share an exponent pair and separation.

---

## 3. Numerical stability of the closed forms

- For large $\rho$, $e^{-\rho}$ underflows while the prefactors grow. Work with $\ln N_a + \ln N_b - \rho$ and combine the logarithms before exponentiating, or scale the whole pair by $e^{-\rho}$ and keep that factor separately.
- For $|\sigma|$ small, use the Taylor series for $\sinh\sigma/\sigma^k$ and $\cosh\sigma/\sigma^k$ combinations. The differences $\sinh\sigma/\sigma^3 - \cosh\sigma/\sigma^2$ and $\sinh\sigma/\sigma^2 - \cosh\sigma/\sigma$ cancel severely near $\sigma=0$. A threshold near $|\sigma| \approx 10^{-2}$ with a few series terms keeps full precision.
- For very small $\rho$ (nearly coincident centres), the terms $1/\rho^2$ and $1/\rho^3$ are large and cancel against the $\sinh$ and $\cosh$ combinations. Switch to the one-centre formulas when $R$ falls below a small tolerance, and check continuity at the switch.

---

## 4. Three-centre attraction: Gaussian superposition

For a nucleus $C$ that is not on either centre of the pair, the prolate construction does not apply. For the pure exponential there is an exact and very convenient route through the Gaussian representation

$$
e^{-\zeta r} = \frac{1}{\sqrt\pi}\int_0^\infty x^{-1/2}\,e^{-x}\,\exp\!\left(-\frac{\zeta^2 r^2}{4x}\right)dx .
$$

This identity holds for all $r \ge 0$ and is checked at $r = 0$, where the integral gives $\Gamma(\tfrac12)/\sqrt\pi = 1$.

### 4.1 Discrete expansion

The integral has the form of a generalized Laguerre weight $x^{\alpha}e^{-x}$ with $\alpha = -\tfrac12$. Applying the generalized Gauss–Laguerre rule with $N$ nodes $x_i$ and weights $w_i$ for this weight gives

$$
\boxed{\;e^{-\zeta r} \approx \sum_{i=1}^{N} d_i\, e^{-t_i r^2}, \qquad d_i = \frac{w_i}{\sqrt\pi}, \qquad t_i = \frac{\zeta^2}{4x_i}\;}
$$

Each term is a Gaussian primitive with exponent $t_i$ and contraction coefficient $d_i$. The expansion is exact at $r = 0$ for every $N$, because the rule integrates the constant $1$ exactly.

The accuracy depends on the distance range. The Gaussian $e^{-\zeta^2r^2/(4x)}$ is only appreciable for $x$ of order $\zeta^2r^2/4$, so a radial distance $r$ requires nodes reaching that $x$. As a rough guide, the rule resolves distances up to about $r \approx 4\sqrt{N}/\zeta$ to reasonable accuracy, and matrix elements over a larger range need more nodes or a tail correction. This is a property of the discretization, not of the method, and convergence should be checked by increasing $N$.

### 4.2 Multi-centre matrix elements from primitives

Each pure exponential is now a contracted Gaussian set $\{(d_i^{a}, t_i^{a})\}$ on its centre. For a pair of primitives $t$ on $\mathbf R_a$ and $u$ on $\mathbf R_b$, define

$$
p = t+u, \qquad \mu = \frac{tu}{t+u}, \qquad \mathbf P = \frac{t\mathbf R_a+u\mathbf R_b}{p}, \qquad K = \exp(-\mu R^2).
$$

The primitive integrals are standard:

$$
s(t,u) = \left(\frac{\pi}{p}\right)^{3/2}K, \qquad \tau(t,u) = \mu\,(3-2\mu R^2)\,s(t,u),
$$

$$
v_C(t,u) = -Z_C\,\frac{2\pi}{p}\,K\,F_0\!\left(p\,|\mathbf P-\mathbf R_C|^2\right),
$$

where $F_0$ is the Boys function, $F_0(T) = \tfrac12\sqrt{\pi/T}\,\mathrm{erf}(\sqrt T)$, with $F_0(0) = 1$. The contracted matrix elements are then

$$
\boxed{\;S_{ab} = N_aN_b\sum_{i,j} d^{a}_i d^{b}_j\, s(t_i^a,t_j^b),\qquad T_{ab} = N_aN_b\sum_{i,j} d^{a}_i d^{b}_j\, \tau(t_i^a,t_j^b),\qquad V_{ab} = \sum_C N_aN_b\sum_{i,j} d^{a}_i d^{b}_j\, v_C(t_i^a,t_j^b)\;}
$$

For the general case of any nucleus, this is a sum of $N^2$ closed-form terms per matrix element and per nucleus. Because the Gaussian overlap decays as $\exp(-\mu R^2)$, the sum is very amenable to screening.

### 4.3 Why this covers the three-centre case

The formulas in Section 4.2 do not distinguish between the nucleus being on a basis centre and being elsewhere. They therefore handle the two-centre attraction and the three-centre attraction uniformly, and they can be used as an independent check on the closed forms of Section 2.

### 4.4 Alternative three-centre route in Fourier space

The Fourier transform of the pure exponential is $\int e^{-\zeta r}e^{-i\mathbf k\cdot\mathbf r}\,d^3\mathbf r = 8\pi\zeta/(k^2+\zeta^2)^2$, and the transform of $1/r$ is $4\pi/k^2$. The three-centre attraction can then be written as a momentum-space integral with phase factors $e^{i\mathbf k\cdot(\mathbf R_b-\mathbf R_a)}$ and $e^{i\mathbf k\cdot(\mathbf R_C - \mathbf R_a)}$. This route avoids the Boys function but requires a multidimensional oscillatory integral, which is less efficient than the Gaussian route for most purposes. It is mentioned here because it gives an independent numerical check.

### 4.5 Numerical multi-centre quadrature

The Becke-partitioned atomic-grid quadrature described in the multi-centre note applies unchanged to the pure exponential basis. It is the most general fallback and the most convenient for checking the other routes, but it is slower for a given accuracy.

---

## 5. Assembly algorithm

The overall procedure is the same as in the general multi-centre note, with the two-centre and three-centre blocks replaced by the forms above.

```
input : basis {(R_i, zeta_i)}, nuclei {(R_C, Z_C)}, N_gauss (nodes), tolerances
output: S, T, V (symmetric)

precompute for each basis function i:
  normalization N_i = sqrt(zeta_i^3 / pi)
  Gauss-Laguerre (alpha = -1/2, N_gauss nodes): d_i^k, t_i^k  for k = 1..N_gauss
  (node/weight table computed once and rescaled by zeta_i)

for each pair (i,j), i >= j:
  R_ij = |R_i - R_j|

  if R_ij < tol_coincident:
    p = zeta_i + zeta_j      (spatially normalized: N = sqrt(zeta^3/pi))
    S = 8 pi N_i N_j / p^3
    T = 4 pi zeta_i zeta_j N_i N_j / p^3
    V = -sum_C Z_C * 4 pi N_i N_j / p^2      (nuclei at the common centre)
  else:
    S[i][j] = closed form (Section 2.1)
    J_B = closed form (Section 2.2) for nuclei on centre j
    J_A = closed form (Section 2.2) for nuclei on centre i
    T[i][j] = average of -1/2 zeta_j^2 S + zeta_j J_B   and   -1/2 zeta_i^2 S + zeta_i J_A
    V[i][j] = -Z_j J_B - Z_i J_A
    for each nucleus C not at R_i or R_j:
      V[i][j] += sum_{k,l} N_i N_j d_k^i d_l^j v_C(t_k^i, t_l^j)     # Section 4.2
  mirror into (j,i)

H = T + V
```

Practical points:

- Screening: the closed-form overlap decays as $e^{-\rho}$ and the Gaussian-route terms decay as $e^{-\mu R^2}$. A pair whose estimated magnitude falls below the tolerance can be skipped. A pair with $R$ beyond a few divided by the smaller exponent is negligible.
- Three-centre cost: each nucleus not on the pair requires $N^2$ Boys-function evaluations per pair. The Boys function can be tabulated, and its argument $p|\mathbf P-\mathbf R_C|^2$ is cheap to compute, so the dominant cost is the sum over primitive pairs. Reducing $N$ while still meeting the accuracy target is the main optimization.
- Caching: the closed forms depend only on $\rho$, $\sigma$ and $R$, so pairs sharing an exponent pair and distance reuse the same auxiliary values.
- Solve $\mathbf H\mathbf c = E\mathbf S\mathbf c$ by Cholesky or canonical orthogonalization, as in the general note.

---

## 6. Validation

1. **Equal-exponent closed forms.** For two identical exponents, check $S = e^{-\rho}(1+\rho+\rho^2/3)$, $J^{(B)} = \zeta(1+\rho)e^{-\rho}$, and $T = \tfrac{\zeta^2}{2}e^{-\rho}(1+\rho-\rho^2/3)$ with $\rho = \zeta R$.
2. **One-centre limit.** As $R\to 0$, check $S\to 1$, $T\to\zeta^2/2$, and $J\to\zeta$.
3. **Gaussian-route normalization.** For a single centre, $\sum_{i,j} d_id_j\,(\pi/(t_i+t_j))^{3/2}$ times $N^2$ must approach $1$ as $N$ grows. This tests both the Gauss–Laguerre rule and the weight $x^{-1/2}$.
4. **Route agreement.** For a set of geometries, the closed two-centre forms and the Gaussian-superposition forms must agree, with the discrepancy shrinking as $N$ increases.
5. **Hermiticity.** The two forms of $T_{ab}$ in Section 2.3 must agree, and $V_{ab} = V_{ba}$.
6. **Invariance.** Rigid translations and rotations of the whole molecule must leave the eigenvalues unchanged.
7. **Variational behaviour.** Adding basis functions must not raise the lowest eigenvalue.
8. **H$_2^+$ benchmark.** One-electron, two-proton system with one pure exponential per proton. The LCAO energy should show a minimum at an internuclear distance near 2.5 bohr, and the curve should approach the one-centre hydrogen-atom value at large separation, with the energy-curve shape compared against published LCAO-STO results for this molecule.

---

## Summary

1. For $\chi = N e^{-\zeta r}$, set $\rho = R(\zeta_a+\zeta_b)/2$ and $\sigma = R(\zeta_a-\zeta_b)/2$.
2. Compute $S_{ab}$ and the nuclear-attraction integrals $J^{(A)}_{ab}$, $J^{(B)}_{ab}$ in closed form using the sinh/cosh expressions of Sections 2.1 and 2.2, with Taylor series for small $\sigma$.
3. Obtain $T_{ab}$ from the Laplacian identity $T = -\tfrac12\zeta^2 S + \zeta J$, using either centre and averaging.
4. For any nucleus not on the pair, expand each orbital in a Gaussian superposition using generalized Gauss–Laguerre with $\alpha=-\tfrac12$, then evaluate the primitive overlap, kinetic and Boys-function attraction terms.
5. Assemble, solve the generalized eigenproblem, and validate against the equal-exponent closed forms, the Gaussian-route normalization, Hermiticity, and the H$_2^+$ benchmark.
