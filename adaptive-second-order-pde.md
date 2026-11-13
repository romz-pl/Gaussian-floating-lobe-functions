# Gaussian-Kernel Multiscale Representations for Second-Order PDEs in $\mathbb{R}^3$

## 1. Introduction

Meshless kernel methods approximate a solution by a finite combination of translates of a positive definite kernel centred at scattered points, so no mesh is needed. The Gaussian kernel

$$
K_\sigma(x,y) = \exp\!\left(-\frac{\|x-y\|_2^2}{2\sigma^2}\right), \qquad x,y \in \mathbb{R}^3,
$$

is the most widely studied choice. It is strictly positive definite, infinitely smooth, and its native space consists of functions with rapidly decaying Fourier transforms. The parameter $\sigma$ controls the width of the kernel: small $\sigma$ gives localised basis functions and well-conditioned but low-accuracy local interpolation, while large $\sigma$ gives flat, highly accurate but severely ill-conditioned systems. A single-scale method must compromise between these regimes. Multiscale representations resolve the compromise by using a hierarchy of widths, each responsible for a different band of frequencies in the solution.

This document describes how Gaussian multiscale representations are applied to second-order elliptic and parabolic PDEs in three dimensions, how adaptive refinement is organised, and how a posteriori and a priori error estimates are obtained.

## 2. Model problem

Consider a second-order linear operator on a bounded domain $\Omega \subset \mathbb{R}^3$:

$$
\mathcal{L}u = -\nabla\cdot\big(A(x)\nabla u\big) + \mathbf{b}(x)\cdot\nabla u + c(x)\,u = f(x) \quad \text{in } \Omega,
$$

with boundary conditions

$$
u = g_D \ \text{on } \Gamma_D, \qquad \partial_n u = g_N \ \text{on } \Gamma_N,
$$

where $A(x)$ is a symmetric positive definite $3\times 3$ matrix field. The Laplacian is the special case $A = I$, $\mathbf{b}=0$, $c=0$. For the Poisson equation $-\Delta u = f$, the kernel has to be differentiated twice, which is the essential ingredient of collocation methods.

For parabolic problems $\partial_t u = \Delta u + f$ the same spatial machinery applies after a time discretisation, and the heat kernel, itself a Gaussian in $x$, provides a natural multiscale structure in time-stepping schemes.

## 3. The Gaussian kernel as a basis for PDE approximation

### 3.1 Trial space and collocation

Choose a set of centres $X = \{x_1,\dots,x_N\} \subset \overline{\Omega}$ and form the trial space

$$
s(x) = \sum_{j=1}^N c_j\, K_\sigma(x,x_j).
$$

Because $\mathcal{L}$ is a differential operator, the collocation approach (Kansa's method) requires applying $\mathcal{L}$ to each basis function. Writing $\mathcal{L}_x$ for the operator acting in the first argument, the discrete problem is

$$
\sum_{j=1}^N c_j\, \mathcal{L}_x K_\sigma(x_i,x_j) = f(x_i), \qquad x_i \in X_\Omega,
$$

$$
\sum_{j=1}^N c_j\, \mathcal{B}_x K_\sigma(x_i,x_j) = g(x_i), \qquad x_i \in X_{\partial\Omega},
$$

where $\mathcal{B}$ denotes the boundary operator. The resulting matrix $M \in \mathbb{R}^{N\times N}$ has entries $M_{ij} = \mathcal{L}_x K_\sigma(x_i,x_j)$ or $\mathcal{B}_x K_\sigma(x_i,x_j)$ depending on the row. For the Gaussian these derivatives are explicit. For instance, for the Laplacian, with $r = \|x-y\|_2$,

$$
\Delta_x K_\sigma(x,y) = \frac{r^2 - 3\sigma^2}{\sigma^4}\,K_\sigma(x,y).
$$

The collocation matrix is asymmetric in general, which is the price of the simplicity of Kansa's method. A symmetric variant, the Hermite–Birkhoff collocation of Fasshauer, uses derivatives of the kernel at both arguments and yields a symmetric positive definite system when the operator is self-adjoint and the boundary conditions are of the same type.

### 3.2 Native space and stability

The Gaussian kernel generates a reproducing kernel Hilbert space $\mathcal{N}_{K_\sigma}(\mathbb{R}^3)$. For a bounded domain, restrictions of this space carry a norm that equivalently controls Sobolev norms of order $\le$ a large integer, so the approximation is effectively in $C^\infty$ for smooth solutions. Stability is governed by the smallest eigenvalue of the interpolation matrix $K = [K_\sigma(x_i,x_j)]$, which decays super-exponentially in $N$ for fixed $\sigma$ and fixed spacing. Schaback's uncertainty principle states that one cannot simultaneously achieve high approximation accuracy and a well-conditioned matrix with a fixed kernel; the Gaussian is the extreme case where accuracy improves dramatically as $\sigma$ grows while stability deteriorates.

### 3.3 Choosing the shape parameter

A standard heuristic is Rippa's leave-one-out cross-validation. For the interpolant with coefficients $c = K^{-1}f$, the LOO error at point $i$ is

$$
e_i = \frac{c_i}{\big[K^{-1}\big]_{ii}},
$$

and $\sigma$ is chosen to minimise the root mean square of these values. For PDE problems the same idea is applied to the collocation residual, but the lack of a cheap closed form for asymmetric systems makes it expensive, so it is typically used offline or in a coarse pre-computation.

## 4. Multiscale representation

### 4.1 Hierarchy of widths

Introduce a sequence of widths $\sigma_0 > \sigma_1 > \cdots > \sigma_L > 0$, usually geometric, $\sigma_\ell = \sigma_0\, 2^{-\ell}$, paired with nested or refined centre sets $X_0 \subset X_1 \subset \cdots \subset X_L$ such that the fill distance satisfies $h_\ell \approx \sigma_\ell$. The multiscale approximant is

$$
s_L(x) = \sum_{\ell=0}^{L} \sum_{x_j \in X_\ell} c_{\ell,j}\, K_{\sigma_\ell}(x,x_j).
$$

Each level contributes a band of frequencies: coarse levels capture the low-frequency, smooth part of $u$ with broad kernels, while fine levels resolve local features with narrow kernels. This is the Gaussian analogue of a multiresolution decomposition. Compared with a single-level expansion at the finest scale, the coarse levels keep the global system well conditioned because each level is solved for a smaller correction.

### 4.2 Successive residual correction

The practical construction is a level-by-level residual correction. Set $s_{-1} = 0$ and for $\ell = 0,1,\dots,L$:

1. Compute the residual of the current approximation on the level's collocation points, $r_\ell(x) = f(x) - \mathcal{L}s_{\ell-1}(x)$ in $\Omega$, with boundary residuals analogously.
2. Solve, at width $\sigma_\ell$ and centres $X_\ell$, the collocation problem $\mathcal{L}\,t_\ell = r_\ell$ with the boundary conditions adapted to the residual.
3. Update $s_\ell = s_{\ell-1} + t_\ell$.

Each subproblem involves a matrix with a moderate condition number because the width is matched to the spacing of its centres. The procedure is related to multilevel interpolation and to greedy and wavelet-type decompositions. Convergence of the sequence follows from the contraction properties of the residual correction when the coarse level has already reduced the residual below a fixed fraction of its previous size; in practice one observes a geometric decay of the residual norm per level, whose rate depends on the ratio $\sigma_{\ell+1}/\sigma_\ell$ and on the smoothness of $u$.

### 4.3 Scale-selection and conditioning

The condition number of each level's matrix can be controlled by choosing $\sigma_\ell$ relative to $h_\ell$ so that $\sigma_\ell / h_\ell$ stays in a fixed window, typically between 1 and 3 for three-dimensional problems. Using stable evaluation algorithms, such as the contour–Padé method or RBF-QR factorisations, one can push the usable range of $\sigma/h$ to smaller values without loss of accuracy. These algorithms reparametrise the kernel expansion in a well-conditioned orthonormal basis, so the ill-conditioning is decoupled from the representation of the approximant.

## 5. Adaptive procedures

### 5.1 Adaptive loop

The adaptive method follows the classical SOLVE → ESTIMATE → MARK → REFINE cycle familiar from finite element adaptivity, adapted to scattered centres:

- SOLVE: compute the multiscale approximant $s$ on the current level structure.
- ESTIMATE: compute local error indicators $\eta_\xi$ over a partition of $\Omega$ into cells or over a set of test points $\{\xi\}$.
- MARK: select a set $\mathcal{M}$ of subregions satisfying a Dörfler criterion,

$$
\sum_{\xi \in \mathcal{M}} \eta_\xi^2 \ \ge\ \theta \sum_{\xi} \eta_\xi^2, \qquad 0 < \theta < 1.
$$

- REFINE: in each marked region, add centres and possibly a new level with a smaller width, so that $h$ and $\sigma$ decrease locally.

Because the approximant is global in its support, refinement is local in the coefficient space: new centres only add rows and columns, and previously computed levels remain usable as preconditioners or as initial guesses for iterative solvers.

### 5.2 Choosing where and how to refine

Two refinement strategies are common. Centre refinement inserts new points at the same width, improving the approximation in the marked region through a smaller fill distance. Width refinement adds a finer level with a smaller $\sigma$ restricted to the marked region, capturing features such as boundary layers or steep gradients of $A$ or $c$ without enlarging the whole system. The choice is guided by the estimator: if the residual is dominated by high-frequency content, which can be diagnosed from the decay of the level coefficients $\|c_{\ell}\|$, width refinement is preferred; otherwise centre refinement suffices.

A greedy variant chooses new centres at maximal values of the power function or the residual. The power function, defined below, is maximised where the current centre set provides the least information, so greedy selection drives the centre set toward a quasi-uniform distribution that is adapted to the solution's behaviour in regions where the residual is large.

### 5.3 Adaptive treatment of anisotropy and boundaries

When $A(x)$ is strongly anisotropic, the isotropic Gaussian can be generalised to a Gaussian with a symmetric positive definite shape matrix $\Sigma(x)$, which is locally adapted to the diffusion tensor. This anisotropic multiscale kernel reduces the number of centres needed along directions of weak diffusion. Near the boundary, centres can be concentrated by using a graded point distribution, and the width at boundary-adjacent levels can be reduced to resolve Neumann or Robin layers. Adaptive refinement is driven by the boundary residual in the same way as the interior residual.

## 6. Error estimation

### 6.1 Power function and pointwise bounds

The power function of a kernel with centres $X$ is

$$
P_{K,X}(x) = \sqrt{K(x,x) - \mathbf{k}(x)^T K^{-1}\mathbf{k}(x)}, \qquad \mathbf{k}(x) = \big(K(x,x_j)\big)_{j=1}^N.
$$

For the function $u$ in the native space, it bounds the pointwise interpolation error:

$$
|u(x) - s(x)| \ \le\ P_{K,X}(x)\,\|u - s\|_{\mathcal{N}_K}.
$$

For PDE approximation, the analogous bound replaces the point evaluation by the operator. If $\mathcal{L}$ is bounded from the native space into a Sobolev space $H^{k}(\Omega)$, then

$$
\|u - s\|_{L^2(\Omega)} \ \le\ C_{\text{stab}}\,\|\mathcal{L}(u-s)\|_{L^2(\Omega)} \ \le\ C_{\text{stab}}\,\|f - \mathcal{L}s\|_{L^2(\Omega)} + \text{(boundary terms)},
$$

where $C_{\text{stab}}$ is a stability constant related to the Green's function of $\mathcal{L}$ and to the smallest eigenvalue of the collocation operator. The residual $f - \mathcal{L}s$ is computable, which makes this a posteriori estimate attractive.

### 6.2 A posteriori estimator

Define the local residual indicator on a cell $T$ (or the Voronoi cell of a centre $x_j$) as

$$
\eta_T^2 = h_T^2\,\|f - \mathcal{L}s\|_{L^2(T)}^2 + h_T\,\|g - \mathcal{B}s\|_{L^2(\partial T\cap\Gamma)}^2 + \sum_{F \subset \partial T} h_F\,\|[\![A\nabla s\cdot n]\!]\|_{L^2(F)}^2,
$$

where $h_T$ is the local size, $h_F$ the size of an interface, and $[\![\cdot]\!]$ the jump. For smooth Gaussian approximants the interface jump terms vanish identically in the interior, because $s$ is $C^\infty$, so the estimator reduces to the interior and boundary residual terms. The global estimator $\eta = (\sum_T \eta_T^2)^{1/2}$ is reliable and efficient up to a constant that depends on the shape regularity of the cells and on the conditioning of the local collocation problems. In a meshless setting the Voronoi cells of the centres serve as the partition; their diameters are of the order of the fill distance.

The residual is evaluated by applying $\mathcal{L}$ analytically to the Gaussian expansion, which is cheap because the derivatives of $K_\sigma$ are explicit polynomials times Gaussians. Quadrature over each cell uses a moderate number of points, since the integrand is smooth.

### 6.3 Level-wise and cumulative estimates

In a multilevel scheme, the residual decreases level by level. Let $\rho_\ell = \|f - \mathcal{L}s_\ell\|$. The per-level increments $\|t_\ell\|_{\mathcal{N}}$ estimate the contribution of each scale, and a geometric decay $\rho_{\ell+1} \le \gamma\,\rho_\ell$ with $\gamma < 1$ signals convergence. Stopping when $\rho_L \le \varepsilon\,\rho_0$, or when the increment $\|t_L\|$ falls below a tolerance, gives a computable termination criterion that does not require the exact solution. The residual at the finest level also indicates whether additional levels are needed in particular regions, which feeds directly into the marking step.

### 6.4 A priori estimates

For $u$ in a Sobolev space of sufficient regularity, classical results for Gaussian kernels give exponential convergence in the fill distance when $\sigma$ is held fixed, of the form

$$
\|u - s\|_{L^\infty(\Omega)} \le C\, e^{-c/h}\,\|u\|_{\mathcal{N}},
$$

valid for analytic $u$ and for a suitable choice of $\sigma$ scaling with $h$. For finite smoothness, algebraic rates of the form $h^{k}$ are obtained by using the Wendland-type kernels instead; the Gaussian yields spectral accuracy at the cost of the stability issues discussed above. The multiscale construction removes part of this tension: each level uses a width matched to its spacing, so the a priori rates per level combine into a bound on the total error that is governed by the finest spacing and by the smoothness of the residual rather than by $\sigma_0$ alone.

### 6.5 Practical error control

In practice, one monitors three quantities during an adaptive run: the residual estimator $\eta$, the level increments $\|t_\ell\|$, and an independent check against a reference solution on a sparse validation set of points not used in collocation. Agreement between the estimator and the validation error increases confidence in the result. Discrepancies usually indicate either insufficient resolution of boundary layers, which calls for boundary refinement, or excessive $\sigma$ relative to $h$, which calls for width refinement or a stabler evaluation algorithm.

## 7. Implementation notes

The principal computational cost is the assembly and solution of the collocation system at each level. For $N$ centres the dense system costs $O(N^3)$ to factorise, which limits direct methods to a few thousand centres. Two strategies extend the reach. Domain decomposition with overlapping Gaussian patches, each with its own multiscale hierarchy, reduces the global system to a sequence of local ones coupled by interface conditions. Iterative solvers, such as GMRES preconditioned by the coarse-level solution, exploit the fact that coarse levels are smooth and well conditioned. Fast multipole methods accelerate matrix-vector products for large centre sets, since the Gaussian has rapidly decaying far-field interactions that can be truncated or compressed.

Floating-point stability requires care when $\sigma/h$ is small or when $N$ is large. Using extended precision for the Cholesky or QR factorisation of the interpolation matrix, or switching to the stable basis reparametrisation mentioned in Section 4.3, avoids catastrophic cancellation. The derivatives of the kernel should be evaluated using the closed-form expressions rather than automatic differentiation of the exponential, to preserve accuracy for large $r/\sigma$.

## 8. Summary

Gaussian-kernel multiscale representations offer a meshless, spectrally accurate framework for second-order PDEs in $\mathbb{R}^3$. A hierarchy of widths matched to a hierarchy of centre sets produces a well-conditioned, level-by-level residual correction. Adaptive refinement uses residual-based indicators and Dörfler marking to concentrate centres and finer levels where the solution is poorly resolved. Error control rests on the power function for pointwise bounds, on residual-based a posteriori estimators for the operator, on geometrically decaying level increments as a termination criterion, and on validation against independent points. The main limitations are the cubic cost of dense factorisations, the sensitivity of stability to $\sigma/h$, and the need for stable evaluation algorithms when accuracy requirements are high. Addressing these through domain decomposition, stable bases, and fast summation is what makes the approach practical for realistic three-dimensional problems.
