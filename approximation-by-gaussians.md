# Approximation of Functions on $\mathbb{R}^3$ by Gaussian Functions

## 1. Problem statement

Let $f:\mathbb{R}^3 \to \mathbb{R}$ (or $\mathbb{C}$) be a function, typically defined on a bounded domain $\Omega \subset \mathbb{R}^3$ and decaying at infinity, or known only through samples $f(x_i)$ at points $x_i$. The goal is to find a finite *Gaussian expansion*

$$
g_N(x) = \sum_{k=1}^{N} c_k \, \exp\!\left( -\alpha_k \, \lVert x - \mu_k \rVert^2 \right),
\qquad x \in \mathbb{R}^3,
$$

where the parameters are the weights $c_k$, the centres $\mu_k \in \mathbb{R}^3$, and the inverse widths $\alpha_k > 0$, such that $g_N \approx f$ in a chosen norm.

Several generalisations are common:

**Anisotropic Gaussians.** Replacing the isotropic kernel by a full covariance matrix $\Sigma_k$ gives

$$
g_N(x) = \sum_{k=1}^{N} c_k \exp\!\left( -\tfrac{1}{2}(x-\mu_k)^{\mathsf{T}} \Sigma_k^{-1} (x-\mu_k) \right),
$$

which is the standard building block in Gaussian mixture models and in Gaussian splatting for graphics. Anisotropy allows fewer terms for elongated features but makes fitting harder.

**Gaussians with polynomial prefactors.** Gaussian-type orbitals in quantum chemistry use $x^a y^b z^c \exp(-\alpha r^2)$, which are needed to represent angular structure. These are not treated further here, but the same fitting ideas apply.

**Error norms.** The appropriate norm depends on the application:

- $L^2(\mathbb{R}^3)$ error for energy-like or projection-based problems,
- $L^\infty$ (sup) error for pointwise accuracy and graphics,
- weighted $L^2$ norms, for example $L^2(\mathbb{R}^3, e^{-|x|^2}\,dx)$, when a Gaussian weight is natural.

## 2. Why Gaussians?

Gaussians are unusually convenient, which explains their widespread use.

**Closed-form algebra.** Products of Gaussians are Gaussians:

$$
\exp(-\alpha|x-\mu|^2)\,\exp(-\beta|x-\nu|^2) = K \exp\!\left(-\gamma |x-\rho|^2\right),
$$

with $\gamma = \alpha+\beta$, $\rho = (\alpha\mu+\beta\nu)/(\alpha+\beta)$, and $K = \exp\!\left(-\tfrac{\alpha\beta}{\alpha+\beta}|\mu-\nu|^2\right)$. Integrals of Gaussian products, overlap and Coulomb-type matrix elements, therefore have analytic forms.

**Separability.** In three dimensions,

$$
\exp(-\alpha|x|^2) = \exp(-\alpha x_1^2)\,\exp(-\alpha x_2^2)\,\exp(-\alpha x_3^2),
$$

so a single isotropic Gaussian is a tensor product of one-dimensional Gaussians. This is the foundation for tensor-product and low-rank methods discussed in Section 5.

**Fourier duality.** The Fourier transform of a Gaussian is a Gaussian. With the convention $\hat f(k) = \int f(x) e^{-i k\cdot x}\,dx$,

$$
\widehat{\exp(-\alpha|x|^2)}(k) = \left(\frac{\pi}{\alpha}\right)^{3/2} \exp\!\left(-\frac{|k|^2}{4\alpha}\right).
$$

Thus a Gaussian expansion is simultaneously a localised spatial object and a localised spectral object.

**Convolution semigroup.** Gaussians are closed under convolution: normalised Gaussians of variances $\sigma^2$ and $\tau^2$ convolve into a Gaussian of variance $\sigma^2+\tau^2$. This makes Gaussian smoothing (the heat semigroup) a natural approximation tool.

**Smoothness and positivity.** Gaussians are $C^\infty$, strictly positive, and rapidly decaying, which matches many physical densities and potentials.

**Machine learning connection.** Gaussian kernels are reproducing kernels of an infinitely smooth RKHS and underpin radial basis function (RBF) interpolation, kernel density estimation and Gaussian processes.

## 3. Theoretical foundations

### 3.1 Density (existence of approximations)

A basic question is whether Gaussian sums can approximate arbitrary functions at all. The answer is yes, under mild conditions. Two classical routes are:

1. **Stone–Weierstrass type arguments.** Sums of Gaussians with varying centres and widths form an algebra-like family that separates points and does not vanish, so it is dense in the continuous functions on compact sets. Density in sup norm on compacts follows from this type of argument, with extensions to $L^p$ on $\mathbb{R}^3$ when $f$ and the family are suitably integrable.

2. **Fourier argument.** If $\hat f$ is sufficiently integrable, the inverse Fourier representation of $f$ can be approximated by Gaussian-weighted superpositions, because Gaussians in $k$-space are the Fourier images of Gaussians in $x$-space.

These results guarantee existence but say nothing about how many terms are needed. Quantitative rates are the central practical concern.

### 3.2 Gaussian mollification

The simplest constructive approximation is convolution with a Gaussian kernel. Let

$$
G_\sigma(x) = (2\pi\sigma^2)^{-3/2} \exp\!\left(-\frac{|x|^2}{2\sigma^2}\right),
$$

and define $f_\sigma = f \ast G_\sigma$. Because $G_\sigma$ integrates to one and is concentrated near the origin, a Taylor expansion of $f$ gives

$$
f_\sigma(x) - f(x) = \frac{\sigma^2}{2}\,\Delta f(x) + O(\sigma^4) \quad \text{as } \sigma \to 0,
$$

so the error is $O(\sigma^2)$ for $f \in C^2$ with bounded Laplacian. The convolution $f_\sigma$ is not itself a finite Gaussian sum, but it can be discretised: replacing the integral by a quadrature over a grid of centres $\mu_k$ with weights $w_k$ gives

$$
f_\sigma(x) \approx \sum_k w_k\, G_\sigma(x-\mu_k).
$$

This is the standard *Gaussian blur plus quadrature* approximation. Its accuracy is controlled by both the smoothing error $O(\sigma^2)$ and the quadrature error, which becomes small when the grid spacing is at most a fraction of $\sigma$.

The price is that the number of centres scales as $(L/h)^3$ on a domain of size $L$ with spacing $h$, which is the curse of dimensionality in its simplest form. Adaptive or sparse quadrature mitigates this.

### 3.3 Approximation rates

Rates depend strongly on the regularity of $f$ and on how the width is chosen.

- **Fixed width, many centres.** With a fixed width $\alpha$, the span of Gaussians is dense, but for generic $f \in H^s(\mathbb{R}^3)$ the error decays only algebraically in $N$, with exponent depending on $s$ and the dimension $d=3$.

- **Adaptive width.** Allowing $\alpha_k$ to vary, especially across scales, improves convergence. Multiscale dictionaries combining widths $\alpha_j = 2^{j}$ behave like wavelet frames and yield better rates for functions with localised features.

- **Analytic functions.** For functions analytic in a neighbourhood of a compact set, exponential convergence in $N$ is achievable with well-chosen widths and centres. The precise rate depends on the strip of analyticity and the chosen grid, and is a current research topic in the context of Gaussian and exponential sums.

- **Singular kernels.** A key success story is the approximation of singular radial kernels such as $1/r$ on a bounded interval $[\delta, R]$ by exponential sums (Section 4.4). Here the number of Gaussians needed grows only logarithmically in $1/\varepsilon$ and in $R/\delta$.

### 3.4 Ill-conditioning

Gaussians with widely separated widths are nearly linearly dependent. Consequently, the design matrix (Section 4.1) often has a condition number that grows exponentially with the number of terms. Practical consequences:

- Coefficients may become very large with alternating signs, leading to catastrophic cancellation.
- Regularisation (Tikhonov, truncated SVD) is usually essential.
- Numerical evaluation of the expansion may need higher working precision than the target accuracy.

## 4. Fitting methods

### 4.1 Linear least squares with fixed centres and widths

If the centres $\mu_k$ and widths $\alpha_k$ are fixed, the problem is linear in the weights. Given sample points $x_i$ ($i=1,\dots,M$) and values $f_i = f(x_i)$, form the design matrix

$$
A_{ik} = \exp\!\left(-\alpha_k\,\lVert x_i - \mu_k\rVert^2\right),
$$

and solve $\min_c \lVert A c - f \rVert_2^2$. Ridge regularisation gives

$$
c_\lambda = (A^{\mathsf{T}}A + \lambda I)^{-1} A^{\mathsf{T}} f.
$$

Choosing $\lambda$ by generalised cross-validation, the L-curve, or a known noise level is standard practice.

Choosing centres and widths matters greatly:

- **Centres** are often placed on a regular grid, on a random subset of the sample points, or by k-means clustering of the data.
- **Widths** are commonly chosen to match the local feature size, for example via the distance to the nearest neighbouring centre, or from a geometric sequence over a range of scales.

### 4.2 Nonlinear least squares

When centres, widths and weights are all unknown, the problem becomes nonlinear. Write the parameter vector as $\theta = (c_k, \mu_k, \alpha_k)_{k=1}^N$ and minimise

$$
E(\theta) = \sum_{i=1}^{M} \bigl( g_N(x_i;\theta) - f_i \bigr)^2.
$$

Gradients are analytic and inexpensive:

$$
\frac{\partial g_N}{\partial c_k} = \varphi_k(x), \qquad
\frac{\partial g_N}{\partial \alpha_k} = -c_k\,\lVert x-\mu_k\rVert^2\,\varphi_k(x), \qquad
\frac{\partial g_N}{\partial \mu_k} = 2\,\alpha_k\,c_k\,(x-\mu_k)\,\varphi_k(x),
$$

where $\varphi_k(x) = \exp(-\alpha_k\lVert x-\mu_k\rVert^2)$. The problem is non-convex, so good initialisation matters. Typical strategies are:

- initialise centres from a linear fit on a coarse dictionary, then refine;
- use continuation: start with few, wide Gaussians and add terms where the residual is large;
- use Levenberg–Marquardt, which handles the least-squares structure well.

A minimal Python sketch of the nonlinear fit using SciPy is given in Section 7.

### 4.3 Greedy and sparse methods

To avoid nonlinear optimisation, one can select Gaussians from a large dictionary $\mathcal{D}$ of candidate triples $(\mu, \alpha)$, for example centres on a grid and widths on a geometric sequence.

- **Orthogonal matching pursuit (OMP).** At each step, add the dictionary atom most correlated with the current residual, then re-solve the least-squares problem on the selected support. This is simple, fast and interpretable.
- **$\ell_1$ regularisation (LASSO / basis pursuit denoising).** Solve $\min_c \lVert Ac - f\rVert_2^2 + \lambda \lVert c\rVert_1$. This is convex and promotes sparsity. The selected support yields a compact Gaussian expansion.
- **Thresholding.** Solve on the full dictionary, then prune coefficients below a tolerance and refit.

Sparse methods are particularly effective when $f$ has a few localised features, as in density fitting and splatting.

### 4.4 Exponential sum approximation of radial kernels

A powerful result (Beylkin and Monzón, and subsequent work by Hackbusch and coworkers) concerns approximating radial functions $\phi(r)$ with a sum of exponentials in $r^2$:

$$
\phi(r) \approx \sum_{k=1}^{M} w_k \exp(-\alpha_k r^2),
$$

valid on an interval $r \in [\delta, R]$. For the Coulomb kernel $\phi(r) = 1/r$ the approximation error $\varepsilon$ requires a number of terms $M$ that grows only logarithmically in $1/\varepsilon$ and in the ratio $R/\delta$. Because the approximation is an exponential sum in $r^2$, each term is a three-dimensional isotropic Gaussian.

This has major consequences: long-range $1/r$ interactions can be represented by a small number of Gaussians, and the resulting sums can be evaluated by separable one-dimensional operations. This is exploited in fast Gaussian-based Coulomb solvers and in Gaussian-based density fitting.

Practical notes:

- The nodes $\alpha_k$ are typically spaced geometrically, and the weights $w_k$ are obtained by a least-squares or Remez-type procedure on a logarithmic grid.
- The approximation is accurate only on $[\delta,R]$. Behaviour outside this range is not controlled, which is why short-range singularities are often handled separately.

### 4.5 Gaussian mixture fitting via EM

When $f$ is a probability density, or can be normalised as one, the natural tool is the expectation–maximisation (EM) algorithm for Gaussian mixture models:

$$
p(x) \approx \sum_{k=1}^{N} \pi_k\, \mathcal{N}(x; \mu_k, \Sigma_k), \qquad \sum_k \pi_k = 1.
$$

EM alternates between computing responsibilities and updating the weights, means and covariances. It converges to a local maximum of the likelihood. Model order $N$ is selected by information criteria (BIC, AIC) or by cross-validation. Variants include variational Bayes mixtures, which prune unneeded components automatically.

### 4.6 Sampling-based and particle methods

If $f$ is given as a set of particles or a point cloud, a Gaussian expansion can be built directly: place one Gaussian at each point with a width tied to the local spacing, then optimise weights. This is essentially how 3D Gaussian splatting represents radiance fields. Optimisation of positions, covariances, opacities and colours is done with differentiable rendering and gradient descent, together with adaptive densification (splitting and cloning) and pruning.

## 5. Tensor and separable approaches

When $f$ has low-rank structure, the three-dimensional problem reduces to one-dimensional problems. Suppose $f$ is sampled on a tensor grid and admits a CP or Tucker decomposition,

$$
f(x_1,x_2,x_3) \approx \sum_{r=1}^{R} u^{(1)}_r(x_1)\,u^{(2)}_r(x_2)\,u^{(3)}_r(x_3).
$$

Each one-dimensional factor $u^{(j)}_r$ can be approximated by a one-dimensional exponential sum $\sum_m a_m e^{-\beta_m t^2}$. The product of one-dimensional Gaussians is a three-dimensional Gaussian, so the full approximation becomes

$$
f(x) \approx \sum_{r=1}^{R}\sum_{m_1,m_2,m_3} a^{(1)}_{r,m_1}a^{(2)}_{r,m_2}a^{(3)}_{r,m_3}\,\exp\!\left(-\bigl(\beta^{(1)}_{m_1}x_1^2 + \beta^{(2)}_{m_2}x_2^2 + \beta^{(3)}_{m_3}x_3^2\bigr)\right).
$$

Only the axis-aligned, anisotropic Gaussians appear here, and the number of terms is $R \cdot M^3$ in the naive form. Compression by subsequent sum-of-exponentials truncation or by a second tensor decomposition reduces this considerably. This route is attractive for high-resolution fields and for separable kernels such as those arising in Ewald-type splittings.

Tensor-train formats extend the idea to higher dimensions, but for $d=3$ the CP and Tucker forms are usually sufficient.

## 6. Error analysis and diagnostics

For a fitted expansion, practical diagnostics include:

- **Residual on held-out points.** Partition the samples into training and validation sets, and track both errors as $N$ grows. Overfitting shows as a falling training error and a rising validation error.
- **Spectral check.** Compute the Fourier transform of $g_N - f$ on a grid, or equivalently check the residual in a Gaussian-weighted norm. Large high-frequency components indicate insufficient resolution.
- **Coefficient magnitudes.** Very large $|c_k|$ with alternating signs indicate ill-conditioning. Report $\lVert c \rVert$ alongside the error.
- **Tail behaviour.** Because Gaussians decay super-exponentially, a Gaussian expansion can badly misrepresent slowly decaying tails of $f$. Check the far-field separately.
- **Positivity.** If $f \ge 0$ is required, enforce it through nonnegative least squares, through a logarithmic parametrisation, or by squaring a Gaussian sum, though the last option changes the function class.

A posteriori error estimators for Gaussian expansions are less developed than for finite elements, so empirical validation on test points remains the norm.

## 7. Practical algorithm and code sketch

A typical workflow is:

1. Sample $f$ on a set of points $\{x_i\}$ that covers the region of interest, with extra density near features.
2. Choose a dictionary of widths, geometric in $\alpha$ across the expected scale range.
3. Choose centres, either on a grid or at the sample points.
4. Solve a regularised linear least-squares problem for the weights.
5. Optionally, refine centres and widths by nonlinear least squares.
6. Validate on held-out points and check tail behaviour and coefficient size.
7. Prune small coefficients and refit if a compact representation is needed.

The following sketch illustrates steps 2 to 5 for a radial function:

```python
import numpy as np
from scipy.optimize import least_squares

def design_matrix(X, centers, alphas):
    # X: (M, 3), centers: (N, 3), alphas: (N,)
    d2 = np.sum((X[:, None, :] - centers[None, :, :]) ** 2, axis=-1)
    return np.exp(-alphas[None, :] * d2)

def fit_linear(X, f, centers, alphas, lam=1e-8):
    A = design_matrix(X, centers, alphas)
    # Tikhonov-regularised normal equations
    c = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ f)
    return c

def fit_nonlinear(X, f, centers0, alphas0, c0):
    N = len(alphas0)

    def unpack(p):
        c = p[:N]
        alphas = np.exp(p[N:2 * N])          # log-parametrise to keep alpha > 0
        centers = p[2 * N:].reshape(N, 3)
        return c, alphas, centers

    def residual(p):
        c, alphas, centers = unpack(p)
        return design_matrix(X, centers, alphas) @ c - f

    p0 = np.concatenate([c0, np.log(alphas0), centers0.ravel()])
    sol = least_squares(residual, p0, method="lm")
    return unpack(sol.x)

# Example: approximate f(x) = exp(-|x|) (a cusp, not Gaussian) on a sample cloud
rng = np.random.default_rng(0)
X = rng.uniform(-3, 3, size=(2000, 3))
f = np.exp(-np.linalg.norm(X, axis=1))

centers0 = np.zeros((1, 3))
alphas0 = np.array([1.0])
c0 = np.array([1.0])
c, alphas, centers = fit_nonlinear(X, f, centers0, alphas0, c0)
```

For a cusp such as $e^{-|x|}$, a single Gaussian cannot match the non-differentiable point at the origin, so many terms with small widths concentrate near $x=0$. This illustrates an important general point: the Gaussian basis is smooth, so non-smooth features require a cluster of terms.

## 8. Pitfalls

- **Curse of dimensionality on grids.** Full tensor grids scale as $n^3$. Use adaptive sampling or sparse grids for large domains.
- **Width selection.** Too small a width gives spiky, poorly generalising fits. Too large a width blurs features. A multiscale dictionary avoids committing to one scale.
- **Ill-conditioning.** Near-linear dependence of Gaussians with similar centres and widths can destroy numerical accuracy. Regularise and monitor the condition number.
- **Local minima.** Nonlinear fits depend on initialisation. Use continuation or multiple random starts.
- **Tails and far-field.** Gaussians vanish too fast for slowly decaying functions, such as $1/r$ far from the source. Handle long-range behaviour with the exponential-sum construction or separately.
- **Non-smooth features.** Cusps, discontinuities and sharp edges need many terms. Consider other bases or adaptive refinement.
- **Scaling of coefficients.** Working with raw coefficients over widely varying $\alpha$ can produce overflow or underflow. Normalise the dictionary columns.

## 9. Applications

- **Quantum chemistry.** Gaussian-type orbitals are the standard basis. Slater-type orbitals are approximated by fixed linear combinations of Gaussians (STO-$nG$ expansions, as in Hehre, Stewart and Pople, 1969). Density fitting or resolution-of-the-identity methods expand electron densities in auxiliary Gaussian sets, and Gaussian-based Coulomb kernels make the integrals tractable.
- **Electrostatics and Poisson solvers.** Representing $1/r$ by Gaussians allows fast evaluation of Coulomb potentials, including fast Gauss transforms and multilevel summation methods.
- **Molecular dynamics and Ewald-type methods.** Gaussian charge smearing splits long-range and short-range interactions cleanly.
- **Kernel methods and interpolation.** Gaussian RBF interpolation is a standard approach for scattered data in three dimensions. Shape parameter selection is equivalent to choosing the width.
- **Density estimation and statistics.** Gaussian mixture models approximate probability densities in three dimensions for clustering and generative modelling.
- **Computer graphics.** 3D Gaussian splatting represents scenes and volumetric data as collections of anisotropic Gaussians rendered by projection and alpha-blending.
- **Scientific computing with convolution.** Gaussian smoothing and Gaussian-based image or volume filters regularise noisy 3D data.

## 10. Summary

Approximating a function on $\mathbb{R}^3$ by Gaussians is attractive because Gaussians are closed under products, convolution and Fourier transform, and separate across coordinates. Existence of approximations follows from density arguments, and constructive methods range from Gaussian mollification to linear least squares with fixed dictionaries, nonlinear refinement, sparse greedy selection, EM for densities, and exponential-sum constructions for radial kernels. The main difficulties are ill-conditioning, the choice of widths and centres, non-smooth features and slowly decaying tails. Exponential-sum approximations of singular kernels such as $1/r$ are the most striking success, because they convert long-range three-dimensional interactions into compact sums of Gaussians. For practical work, the choice of dictionary and regularisation usually matters more than the choice of solver.
