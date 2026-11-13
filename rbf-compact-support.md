# Compactly Supported Radial Basis Functions in $\mathbb{R}^3$

## 1. Radial basis functions in brief

A radial basis function (RBF) is a function whose value depends only on the Euclidean distance from a fixed center. For a univariate profile $\varphi : [0,\infty) \to \mathbb{R}$, the associated RBF centered at $\mathbf{c} \in \mathbb{R}^3$ is

$$
\Phi(\mathbf{x}) = \varphi\big(\lVert \mathbf{x} - \mathbf{c} \rVert_2\big), \qquad \mathbf{x} \in \mathbb{R}^3 .
$$

Given data sites $X = \{\mathbf{x}_1, \dots, \mathbf{x}_N\} \subset \mathbb{R}^3$ and values $f_1, \dots, f_N$, the standard interpolant is

$$
s(\mathbf{x}) = \sum_{j=1}^{N} \lambda_j \, \varphi\big(\lVert \mathbf{x} - \mathbf{x}_j \rVert_2\big),
$$

where the coefficients $\boldsymbol{\lambda}$ solve the linear system

$$
A \boldsymbol{\lambda} = \mathbf{f}, \qquad A_{ij} = \varphi\big(\lVert \mathbf{x}_i - \mathbf{x}_j \rVert_2\big).
$$

The matrix $A$ is symmetric. Whether the interpolation problem is uniquely solvable depends on the positive definiteness of $\varphi$ (or, for conditionally positive definite profiles, on augmentation with low-order polynomials).

The appeal of RBFs is that they are meshless, dimension-independent in formulation, and can interpolate scattered data of arbitrary geometry. Their main drawback is that the global profiles, such as Gaussians, multiquadrics, and thin-plate splines, produce dense matrices $A$, which makes large problems expensive. Compact support addresses exactly this drawback.

## 2. What compact support means

An RBF has compact support of radius $\delta > 0$ if

$$
\varphi(r) = 0 \quad \text{for all } r \ge \delta .
$$

By rescaling, one usually works with a fixed unit-radius profile $\varphi$ with support $[0,1)$ and defines the shape-adjusted function

$$
\varphi_\delta(r) = \varphi\!\left(\frac{r}{\delta}\right).
$$

Consequences of compact support:

1. **Sparsity.** $A_{ij} = 0$ whenever $\lVert \mathbf{x}_i - \mathbf{x}_j \rVert \ge \delta$. For quasi-uniform data with fill distance $h$, each row contains only $O\big((\delta/h)^3\big)$ nonzero entries in $\mathbb{R}^3$. The system matrix is therefore sparse and can be factored with sparse Cholesky or solved with sparse iterative methods.

2. **Locality.** Changing one data value affects the interpolant only within distance $\delta$ of that data site. This is useful for local updates, streaming data, and domain decomposition.

3. **Finite-range covariance.** In geostatistics and Gaussian processes, a compactly supported covariance function means that spatial correlations vanish beyond a finite distance, which can be exploited in sparse linear algebra for kriging and likelihood evaluation.

4. **Smoothness trade-off.** A compactly supported function is not analytic at the boundary of its support unless it is identically zero. Its smoothness is therefore determined by how gracefully it decays to zero at $r = 1$. Achieving high smoothness requires a high-order vanishing of $\varphi$ and its derivatives at $r = 1$.

## 3. Positive definiteness in $\mathbb{R}^3$

A radial function $\Phi(\mathbf{x}) = \varphi(\lVert \mathbf{x} \rVert)$ is positive definite on $\mathbb{R}^d$ if and only if, by Schoenberg's theorem, $\varphi$ is a completely monotone–type object in the sense of the Schoenberg representation, in which case $\varphi$ admits a Fourier representation with a nonnegative spectral measure. For compactly supported functions the practical test is Fourier-based:

$$
\widehat{\Phi}(\boldsymbol{\omega}) = \int_{\mathbb{R}^3} \Phi(\mathbf{x}) e^{-i \boldsymbol{\omega} \cdot \mathbf{x}} \, d\mathbf{x} \;\ge 0 \quad \text{with } \widehat{\Phi} \not\equiv 0 .
$$

Compactly supported functions that satisfy this are called strictly positive definite. Askey's family of truncated powers gives a simple sufficient condition:

$$
\varphi_\nu(r) = (1 - r)_+^{\nu}
$$

is positive definite on $\mathbb{R}^d$ whenever $\nu \ge \lfloor d/2 \rfloor + 1$. For $d = 3$ this means $\nu \ge 2$. The exponent $\nu = 2$ is the lowest admissible power, so $(1-r)_+^2$ is the smallest truncated power that is positive definite in $\mathbb{R}^3$, while $\nu = 1$ fails in this dimension.

## 4. The Wendland functions in $\mathbb{R}^3$

Holger Wendland's construction produces a family of compactly supported, strictly positive definite functions $\varphi_{d,k}$ on $\mathbb{R}^d$ with prescribed smoothness $C^{2k}$ and a corresponding native space. For $d = 3$ the functions are:

| Smoothness | Profile $\varphi_{3,k}(r)$, $0 \le r < 1$ | Minimal polynomial degree |
|---|---|---|
| $C^0$ | $(1-r)^2$ | 2 |
| $C^2$ | $(1-r)^4 (4r + 1)$ | 5 |
| $C^4$ | $\dfrac{(1-r)^6}{3} \left(35 r^2 + 18 r + 3\right)$ | 8 |
| $C^6$ | $(1-r)^8 \left(32 r^3 + 25 r^2 + 8 r + 1\right)$ | 11 |

Here "$(1-r)^n$" is shorthand for $(1-r)^n_+$, meaning the function is zero for $r \ge 1$. The $C^0$ member coincides with Askey's truncated power, which is the minimal admissible choice.

The $C^2$ function is the one most often used in practice. Writing it out explicitly,

$$
\varphi(r) =
\begin{cases}
(1-r)^4 (4r + 1), & 0 \le r < 1, \\[4pt]
0, & r \ge 1 .
\end{cases}
$$

Its first and second derivatives vanish at $r = 1$, so $\Phi$ is $C^2$ as a function on $\mathbb{R}^3$. The third derivative has a jump at $r=1$, which is the source of the restriction to $C^2$.

### 4.1 Derivatives

Because the profile is radial, derivatives of $\Phi$ at $\mathbf{x} \ne \mathbf{c}$ follow from the chain rule:

$$
\nabla \Phi(\mathbf{x}) = \varphi'(r) \, \frac{\mathbf{x} - \mathbf{c}}{r}, \qquad r = \lVert \mathbf{x} - \mathbf{c} \rVert .
$$

For the $C^2$ profile,

$$
\varphi'(r) = -20 r (1 - r)^3 .
$$

This makes gradients, Laplacians, and other differential operators easy to evaluate, which matters for RBF collocation and RBF-finite-difference methods. The support of every derivative up to order $2k$ is the same ball of radius $\delta$.

### 4.2 Native spaces

Each Wendland function has an associated reproducing kernel Hilbert space, the native space, which is a Sobolev space. For $\varphi_{3,k}$ the native space is

$$
\mathcal{N}_{\varphi_{3,k}}(\mathbb{R}^3) = H^{\tau}(\mathbb{R}^3), \qquad \tau = k + 2 .
$$

Thus the $C^0$ function corresponds to $H^2$, the $C^2$ function to $H^3$, and so on. Higher smoothness of the kernel means functions in the reconstructed space are smoother, and interpolation error bounds improve accordingly.

### 4.3 Fourier transform

For Wendland functions the Fourier transform is strictly positive on $\mathbb{R}^3$, and it decays algebraically as $\lVert \boldsymbol{\omega} \rVert^{-(3+2\tau)}$ in the sense of the Sobolev class. This algebraic decay is the spectral counterpart of the compact support in physical space: a function that is exactly zero outside a ball cannot be band-limited, so its spectrum spreads across all frequencies, but the spreading is controlled by the smoothness index.

## 5. Other compactly supported families

Wendland's functions are the most widely used, but several other families are useful in $\mathbb{R}^3$.

**Wu's functions.** Zongmin Wu proposed a family of compactly supported radial functions built from integrals of truncated powers. Several members are positive definite on $\mathbb{R}^3$ and attain specified smoothness. Some are only conditionally positive definite and therefore require polynomial augmentation.

**Buhmann's functions.** Martin Buhmann constructed compactly supported functions with an explicit closed form based on the logarithm of the radius in certain parameter ranges, and in $\mathbb{R}^3$ these give interpolation schemes with good error behavior. They are less common in production code.

**Gneiting's functions.** Tilmann Gneiting developed compactly supported covariance functions for spatiotemporal statistics, including families derived from generalized Wendland constructions with a free shape parameter. They are useful when a correlation range must be tuned separately from smoothness.

**Piecewise polynomial kernels from SPH.** Smoothed particle hydrodynamics uses compactly supported kernels such as the cubic spline and the Wendland kernels, often without requiring positive definiteness, because the goal is local averaging rather than exact interpolation. These kernels are usually normalized so that their integral over $\mathbb{R}^3$ equals one.

## 6. Choosing the support radius

The support radius $\delta$ is the central design parameter. It controls both accuracy and cost.

**Accuracy side.** A larger $\delta$ gives more overlap between basis functions, which generally improves the conditioning of the interpolation problem and the stability of the reconstructed field. Error estimates for Wendland functions take the form

$$
\lVert f - s \rVert_{L_\infty(\Omega)} \le C \, h^{\,k + 2} \lVert f \rVert_{H^{\tau}(\Omega)},
$$

where $h$ is the fill distance and $C$ depends on the domain and the kernel but not on $h$, provided $\delta$ is a fixed fixed multiple of $h$.

**Cost side.** The number of nonzero entries per row scales like $(\delta/h)^3$. Doubling $\delta$ multiplies the work per row by roughly eight.

A common practical rule is to choose $\delta$ such that each center has between 20 and 60 neighbors within its support, which typically means $\delta \approx 2.5$ to $3$ times the average spacing between data sites in $\mathbb{R}^3$. The exact value should be tuned by checking condition numbers and reconstruction error on a benchmark.

### 6.1 Conditioning

Even though $A$ is sparse, it can be badly conditioned. For Wendland functions in $\mathbb{R}^3$ with $\delta$ fixed relative to $h$, the condition number grows only mildly with $N$, because the smallest eigenvalue is bounded below by a quantity proportional to the Sobolev-type lower bound of the native space. When $\delta$ is taken much larger than $h$, the matrix becomes denser and more ill-conditioned in the flat-limit sense, so the trade-off above is also a conditioning trade-off.

## 7. Computational implementation

Implementing compactly supported RBF methods efficiently in $\mathbb{R}^3$ involves several components.

**Neighbor search.** Building $A$ requires identifying all pairs of data sites within distance $\delta$. A k-d tree or an octree reduces this from $O(N^2)$ to roughly $O(N \log N)$ or $O(N)$ for quasi-uniform data. Uniform grid bucketing with cell size $\delta$ is a simple and fast alternative.

**Assembly.** The matrix is assembled in compressed sparse row (CSR) or compressed sparse column (CSC) format. Symmetry allows storing only the upper triangle.

**Factorization.** Sparse Cholesky factorization with a nested-dissection or approximate-minimum-degree ordering is effective. The fill-in depends on the geometry of the point cloud. For data arranged on surfaces embedded in $\mathbb{R}^3$, the factor behaves more like a 2D problem and fills in less than a fully volumetric configuration.

**Evaluation.** Evaluating $s(\mathbf{x})$ at a point requires only the centers within distance $\delta$ of $\mathbf{x}$. Again a spatial index gives $O(1)$ neighbors on average per query.

**Preconditioning.** When direct factorization is too expensive, conjugate gradients with an incomplete Cholesky or multigrid-inspired preconditioner converge quickly, since the sparsity pattern is fixed and the matrix is well structured.

**Numerical care.** Because the profile is piecewise polynomial, evaluate it in a numerically stable form. Writing $(1-r)^n$ directly is safe for $r \in [0,1]$, but the test $r < 1$ should be done on the squared distance when possible to avoid an unnecessary square root.

## 8. Applications in $\mathbb{R}^3$

**Scattered-data interpolation and approximation.** Reconstructing scalar or vector fields from point clouds, such as medical imaging, geophysical sampling, or LiDAR data, is a standard use. Compact support makes large datasets tractable.

**Surface reconstruction.** Implicit surfaces are built by fitting an RBF interpolant to signed-distance or inside/outside values at surface samples and extracting the zero level set. Compact support keeps the implicit function local and reduces artifacts far from the data.

**Meshless PDE solvers.** RBF collocation and RBF-generated finite differences (RBF-FD) solve partial differential equations on irregular three-dimensional domains. Compactly supported kernels yield sparse differentiation matrices, which dramatically reduce memory compared with global RBFs. RBF-FD in particular uses local stencils formed from the nearest neighbors of each node, which is closely aligned with compact support.

**Geostatistics and kriging.** Covariance models with finite range permit sparse kriging systems. This is valuable for mineral exploration, hydrology, and climate data on large 3D grids.

**Machine learning.** Compactly supported kernels define Gaussian process priors with sparse Gram matrices. They also appear in local regression and kernel density models in three-dimensional settings such as point-cloud classification.

**Physics simulation.** SPH and related particle methods use compactly supported kernels to compute densities and pressure gradients from neighboring particles in three dimensions.

**Shape deformation and registration.** Compactly supported RBFs define smooth displacement fields that move only nearby geometry, which is useful for mesh morphing and deformation of volumetric models.

## 9. Advantages and limitations

**Advantages.**

- Sparse, scalable linear systems
- Local control and local updates
- Explicit smoothness and native-space theory
- Closed-form derivatives
- Easy parallelization, since computations are localized

**Limitations.**

- The support radius must be chosen carefully, and the best value depends on data density and smoothness.
- Accuracy for smooth globally defined functions can be lower than that of global RBFs with spectral convergence, since compactly supported kernels have algebraic rather than exponential convergence in the flat limit.
- Very small supports produce ill-conditioned or nearly diagonal systems that fail to propagate information across gaps in the data.
- Non-uniform point distributions require adaptive support radii, which complicates the theory and implementation.

## 10. Summary

A compactly supported RBF in $\mathbb{R}^3$ is a radial function that vanishes beyond a finite radius, typically chosen as one of the Wendland functions with profile $(1-r)^{n}_+$ times a polynomial. These functions are strictly positive definite, have explicit smoothness and native-space characterizations, and produce sparse interpolation matrices whose sparsity pattern is controlled by the support radius relative to the data spacing. Their combination of local support, algebraic convergence, and closed-form derivatives makes them the standard choice when meshless methods must scale to large three-dimensional problems, in applications ranging from surface reconstruction and geostatistics to meshless PDE solvers and particle methods.
