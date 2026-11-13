# Radial Basis Functions for Solving PDEs

## 1. Overview

Radial basis function (RBF) methods are meshless numerical techniques for approximating functions and solving partial differential equations (PDEs). Instead of partitioning the domain into elements (finite elements) or a structured grid (finite differences), an RBF method works with a scattered set of nodes $\Big( \mathbf{x} _ j \Big)_{j=1}^{N}$ in the domain $\Omega \subset \mathbb{R}^d$. The only geometric information required is the node locations and a notion of distance between them.

The central idea is to approximate the unknown solution as a linear combination of translates of a single radially symmetric function, centred at the nodes:

$$
u(\mathbf{x}) \approx \sum_{j=1}^{N} \lambda_j \, \varphi\!\left(\|\mathbf{x} - \mathbf{x}_j\|\right),
$$

where $\varphi : [0,\infty) \to \mathbb{R}$ is the RBF, $\|\cdot\|$ is typically the Euclidean norm, and $\lambda_j$ are coefficients to be determined. Because the basis depends only on distance from a centre, the same function form works in any dimension and on irregular geometries, which is the primary source of the method's flexibility.

The problem is then reduced to finding the coefficients $\lambda_j$ so that the PDE and its boundary conditions are satisfied, typically at the nodes themselves. This transforms the PDE into a linear (or nonlinear) algebraic system.

---

## 2. RBF Interpolation: The Foundation

Before discussing PDEs, it helps to understand RBF interpolation. Given data values $f_i = f(\mathbf{x}_i)$ at nodes $\mathbf{x}_i$, one seeks

$$
s(\mathbf{x}) = \sum_{j=1}^{N} \lambda_j \, \varphi(\|\mathbf{x} - \mathbf{x}_j\|)
$$

such that $s(\mathbf{x}_i) = f_i$ for all $i$. This leads to the linear system

$$
A \boldsymbol{\lambda} = \mathbf{f}, \qquad A_{ij} = \varphi(\|\mathbf{x}_i - \mathbf{x}_j\|).
$$

For many common kernels the matrix $A$ is symmetric positive definite for distinct nodes, guaranteeing a unique solution. Some kernels (for example polyharmonic splines with odd powers) are only conditionally positive definite and require an augmentation with low-degree polynomials to ensure solvability.

### 2.1 Common Kernels

| Name | $\varphi(r)$ | Notes |
|---|---|---|
| Gaussian (GA) | $e^{-(\varepsilon r)^2}$ | Infinitely smooth, strictly positive definite, very ill-conditioned for small $\varepsilon$ |
| Multiquadric (MQ) | $\sqrt{r^2 + c^2}$ | Infinitely smooth, globally supported, often very accurate |
| Inverse multiquadric (IMQ) | $(r^2 + c^2)^{-1/2}$ | Decays with distance, positive definite |
| Polyharmonic spline (PHS) | $r^{k}$ for odd $k$, $r^{k}\log r$ for even $k$ | Finitely smooth, shape-parameter free, requires polynomial augmentation |
| Wendland $C^2$ | $(1 - r/\delta)_+^4 (4r/\delta + 1)$ | Compactly supported, sparse matrices, finite smoothness |

The **shape parameter** $\varepsilon$ (or $c$) controls how "flat" the basis is. Flat kernels (small $\varepsilon$) tend to give higher accuracy but drive the condition number of the collocation matrix up dramatically. This tension between accuracy and stability is one of the central practical issues in RBF methods.

---

## 3. Solving PDEs with RBFs: Collocation

Consider a general boundary value problem

$$
\mathcal{L} u(\mathbf{x}) = f(\mathbf{x}), \quad \mathbf{x} \in \Omega, \qquad \mathcal{B} u(\mathbf{x}) = g(\mathbf{x}), \quad \mathbf{x} \in \partial\Omega,
$$

where $\mathcal{L}$ is a differential operator (for example $\nabla^2$ for Poisson's equation) and $\mathcal{B}$ encodes boundary conditions (Dirichlet, Neumann, Robin, or mixed).

### 3.1 Node Sets

The nodes are split into interior nodes $\{\mathbf{x}_i\}_{i=1}^{N_I}$ and boundary nodes $\{\mathbf{x}_i\}_{i=N_I+1}^{N}$. Both sets are used as collocation points. A common practice is to place more nodes near boundaries or regions of steep gradients.

### 3.2 Unsymmetric Collocation (Kansa's Method)

Introduced by E. J. Kansa in 1990, this is the most widely known RBF PDE solver. Substituting the expansion into the PDE and boundary conditions and enforcing them at the nodes yields

$$
\sum_{j=1}^{N} \lambda_j \, \mathcal{L}\varphi(\|\mathbf{x}_i - \mathbf{x}_j\|) = f(\mathbf{x}_i), \qquad i = 1,\dots,N_I,
$$

$$
\sum_{j=1}^{N} \lambda_j \, \mathcal{B}\varphi(\|\mathbf{x}_i - \mathbf{x}_j\|) = g(\mathbf{x}_i), \qquad i = N_I+1,\dots,N.
$$

In matrix form,

$$
A\boldsymbol{\lambda} = \mathbf{b}, \qquad A_{ij} = \begin{cases} \mathcal{L}\varphi(\|\mathbf{x}_i - \mathbf{x}_j\|), & i \le N_I, \\ \mathcal{B}\varphi(\|\mathbf{x}_i - \mathbf{x}_j\|), & i > N_I. \end{cases}
$$

The operator is applied to the kernel analytically. For the Laplacian in $\mathbb{R}^2$ with $r = \|\mathbf{x} - \mathbf{x}_j\|$ and a radial kernel $\varphi(r)$,

$$
\nabla^2 \varphi(r) = \varphi''(r) + \frac{1}{r}\varphi'(r).
$$

Kansa's method is very easy to implement because it needs only derivatives of the kernel, and it works for arbitrary linear and nonlinear operators. Its drawback is that the matrix is nonsymmetric and, for some operator and kernel combinations, the system can be singular. Kansa's original formulation is known to be unsolvable in certain cases, which motivated the symmetric variant below.

### 3.3 Symmetric (Hermitian) Collocation

Due to Fasshauer, this approach uses the kernel $\varphi(\|\mathbf{x}-\mathbf{y}\|)$ as a reproducing kernel and applies the operator to both arguments. For a symmetric problem where interior collocation uses $\mathcal{L}$ and boundary collocation uses $\mathcal{B}$, the blocks become

$$
A_{11} = \mathcal{L}_{\mathbf{x}}\mathcal{L}_{\mathbf{y}} \varphi(\|\mathbf{x}-\mathbf{y}\|)\big|_{\mathbf{x}=\mathbf{x}_i,\,\mathbf{y}=\mathbf{x}_j}, \quad \text{and so on for other blocks.}
$$

The resulting matrix is symmetric positive definite whenever the kernel is positive definite, which guarantees solvability for well-posed problems. The price is roughly twice the number of kernel derivative evaluations.

---

## 4. A Worked Example: Poisson's Equation

Consider

$$
-\nabla^2 u = f \ \text{in } \Omega, \qquad u = g \ \text{on } \partial\Omega.
$$

Using the Gaussian kernel $\varphi(r) = e^{-\varepsilon^2 r^2}$ in two dimensions, with $r^2 = (x-x_j)^2 + (y-y_j)^2$, the Laplacian of the kernel is

$$
\nabla^2 \varphi = \left(4\varepsilon^2 r^2 - 4\varepsilon^2 \right)\varepsilon^0 e^{-\varepsilon^2 r^2}\cdot\frac{1}{1} = 4\varepsilon^2\left(\varepsilon^2 r^2 - 1\right)e^{-\varepsilon^2 r^2}.
$$

Let $N_I$ interior and $N_B$ boundary nodes be given. The collocation matrix has entries

$$
A_{ij} = -4\varepsilon^2\left(\varepsilon^2 r_{ij}^2 - 1\right)e^{-\varepsilon^2 r_{ij}^2}\quad (i \le N_I), \qquad A_{ij} = e^{-\varepsilon^2 r_{ij}^2} \quad (i > N_I),
$$

with $r_{ij} = \|\mathbf{x}_i - \mathbf{x}_j\|$. The right-hand side is $b_i = f(\mathbf{x}_i)$ for interior rows and $b_i = g(\mathbf{x}_i)$ for boundary rows. Solving $A\boldsymbol{\lambda} = \mathbf{b}$ gives coefficients, and the approximate solution is evaluated anywhere by summing the kernel expansion.

The method's accuracy typically converges spectrally (exponentially fast) with node refinement for smooth solutions when the shape parameter is chosen appropriately, though the conditioning deteriorates as the nodes become denser or the kernel becomes flatter.

---

## 5. Time-Dependent Problems

For evolution equations such as the heat equation

$$
\frac{\partial u}{\partial t} = \nabla^2 u + f(\mathbf{x},t),
$$

there are two main strategies.

### 5.1 Method of Lines

Here the spatial discretisation is performed with RBF collocation, and the time derivative is kept continuous. Writing $u(\mathbf{x},t) \approx \sum_j \lambda_j(t)\,\varphi(\|\mathbf{x}-\mathbf{x}_j\|)$ and collocating at the nodes produces the semi-discrete system

$$
A\,\dot{\boldsymbol{\lambda}}(t) = B\,\boldsymbol{\lambda}(t) + \mathbf{f}(t),
$$

where $A$ is the interpolation matrix and $B$ is the matrix of the spatial operator applied to the kernel. This is an ODE system that can be advanced with standard integrators, such as backward Euler, Crank–Nicolson, or Runge–Kutta methods, with an appropriate stability analysis.

### 5.2 Rothe's Method (Time Stepping First)

Here the time derivative is discretised first, for example with backward Euler,

$$
\frac{u^{n+1} - u^n}{\Delta t} = \nabla^2 u^{n+1} + f^{n+1},
$$

giving a sequence of stationary problems of the form $(I - \Delta t\,\nabla^2)u^{n+1} = u^n + \Delta t\, f^{n+1}$. Each step is then solved as an RBF collocation problem. This is conceptually simple and allows the same solver to be reused at every step, though it requires a new linear solve at each time level unless the matrix is factored once.

---

## 6. Handling Boundary Conditions

Boundary conditions are enforced through the operator $\mathcal{B}$ applied to the kernel. Examples:

- **Dirichlet**: $\mathcal{B} = I$, so $\mathcal{B}\varphi = \varphi$.
- **Neumann**: $\mathcal{B} = \partial/\partial \mathbf{n}$, so the boundary row uses $\mathbf{n}\cdot\nabla\varphi$.
- **Robin**: $\mathcal{B} = \alpha I + \beta\,\partial/\partial\mathbf{n}$.

Because boundary nodes are included directly in the global system, no special boundary treatment is needed beyond selecting the right operator. This is a practical advantage over methods that require explicit boundary integration or ghost-node schemes.

---

## 7. Choosing the Shape Parameter and the Conditioning Problem

The shape parameter $\varepsilon$ is the most important user-controlled setting. Its effect is a trade-off:

- **Large $\varepsilon$** (narrow kernels): the matrix is well-conditioned and close to diagonal, but accuracy is limited because the basis resolves only local features.
- **Small $\varepsilon$** (flat kernels): accuracy can improve dramatically, sometimes to machine precision on the scale of the problem, but the condition number grows as roughly $\mathcal{O}(\varepsilon^{-2d})$ or worse.

This is known as the **uncertainty principle** of RBF approximation. Several approaches address it:

- **Trial-and-error and optimisation**: choose $\varepsilon$ by minimising a leave-one-out cross-validation error or a residual-based criterion.
- **Contour–Padé and RBF-QR methods**: reformulate the basis so that the flat limit can be computed stably, avoiding the explicit ill-conditioned matrix.
- **Variable shape parameters**: let $\varepsilon$ vary with node position to adapt to local solution features.
- **Preconditioning**: use local or sparse preconditioners to speed up iterative solvers.

Stable computation at small shape parameters is an active research area, and the choice of method often depends on whether the problem is steady or time-dependent.

---

## 8. Global Versus Local RBF Methods

### 8.1 Global Methods

The methods described above are global: each node contributes to every entry in the matrix. This gives high accuracy and spectral convergence for smooth problems but produces dense matrices. For a problem with $N$ nodes, storage and solve costs scale as $\mathcal{O}(N^2)$ and $\mathcal{O}(N^3)$ respectively for direct solvers.

### 8.2 Local Methods: RBF-FD

The RBF-generated finite difference (RBF-FD) method computes stencil weights for derivatives by fitting an RBF interpolant to a small set of neighbouring nodes. Each derivative at a node is then a weighted sum over its stencil, similar to classical finite differences but applicable to scattered nodes in arbitrary dimension. The resulting matrices are sparse, making the method attractive for large-scale problems, and it retains much of the geometric flexibility of global RBF methods. Typical stencils contain between about 20 and 50 nodes in two dimensions.

### 8.3 Compactly Supported RBFs

Compact kernels such as Wendland functions vanish beyond a support radius, which directly produces sparse collocation matrices. They sacrifice some smoothness and sometimes accuracy but can be very efficient for large problems.

---

## 9. Nonlinear Problems

Nonlinear PDEs are handled by substituting the RBF expansion into the nonlinear residual and applying Newton's method (or a quasi-Newton variant). If the residual is $\mathbf{R}(\boldsymbol{\lambda}) = 0$, then each iteration solves

$$
J(\boldsymbol{\lambda}^{(k)})\,\Delta\boldsymbol{\lambda} = -\mathbf{R}(\boldsymbol{\lambda}^{(k)}), \qquad \boldsymbol{\lambda}^{(k+1)} = \boldsymbol{\lambda}^{(k)} + \Delta\boldsymbol{\lambda},
$$

where $J$ is the Jacobian of the residual with respect to the coefficients. Because the kernel derivatives are known analytically, the Jacobian is available in closed form. Examples include reaction–diffusion systems, Navier–Stokes flows (via streamfunction or velocity–pressure formulations), and nonlinear elasticity.

---

## 10. Theoretical Properties

RBF collocation has a developing but substantial error theory:

- **Interpolation error** for sufficiently smooth functions decays exponentially with the fill distance $h$ (maximum distance from any point in the domain to the nearest node) when the kernel is infinitely smooth, provided the shape parameter is fixed.
- **Native space** framework: functions in the reproducing kernel Hilbert space of the RBF define the regularity needed for convergence. Smoothness requirements on the solution affect the attainable order.
- **Stability** of collocation is related to the **separation distance** and the kernel's smoothness, and to the behaviour of the local Lebesgue constant.

For practical purposes, the key takeaways are that accuracy is excellent for smooth solutions on well-distributed nodes, and that convergence in the regime of small shape parameter requires special care.

---

## 11. Advantages and Limitations

### Advantages

- **Geometric flexibility**: works on scattered nodes and irregular domains without meshing.
- **Dimension independence**: the same formulation applies in 1D, 2D, 3D, and higher dimensions.
- **High accuracy**: spectral convergence for smooth solutions with infinitely smooth kernels.
- **Simple implementation**: the basic collocation method requires only kernel derivatives.
- **Adaptivity**: nodes can be refined locally, and the shape parameter can be varied.
- **Meshless moving-boundary capability**: useful in fluid–structure interaction and free-surface problems.

### Limitations

- **Dense matrices** for global methods, limiting problem size without sparse or local variants.
- **Ill-conditioning** for flat kernels, requiring stabilisation techniques.
- **Choice of shape parameter** is not always intuitive and affects both accuracy and stability.
- **Nonsymmetric collocation** (Kansa) may be singular for some problems.
- **Less mature software ecosystem** than finite element or finite difference libraries, although this is improving.
- **Boundary and node placement** still requires care for accuracy near corners and sharp features.

---

## 12. Implementation Sketch

The following Python code illustrates Kansa's method for the Poisson problem on the unit square with a Gaussian kernel:

```python
import numpy as np

def gaussian_kernel_and_laplacian(X, Y, eps):
    """Return the Gaussian kernel matrix and its Laplacian w.r.t. the first argument."""
    dx = X[:, None, 0] - Y[None, :, 0]
    dy = X[:, None, 1] - Y[None, :, 1]
    r2 = dx**2 + dy**2
    phi = np.exp(-eps**2 * r2)
    lap = 4 * eps**2 * (eps**2 * r2 - 1) * phi
    return phi, lap

def solve_poisson(nodes, interior_mask, f, g, eps):
    """Solve -Laplacian(u) = f in the domain, u = g on the boundary, by Kansa collocation."""
    phi, lap = gaussian_kernel_and_laplacian(nodes, nodes, eps)
    A = np.where(interior_mask[:, None], -lap, phi)
    b = np.where(interior_mask, f(nodes), g(nodes))
    lam = np.linalg.solve(A, b)
    return lam

def evaluate(x, nodes, lam, eps):
    """Evaluate the RBF approximation at points x."""
    r2 = ((x[:, None, :] - nodes[None, :, :])**2).sum(axis=2)
    return np.exp(-eps**2 * r2) @ lam

# Example usage
rng = np.random.default_rng(0)
nodes = rng.random((400, 2))
interior = (nodes > 0).all(axis=1) & (nodes < 1).all(axis=1)
interior &= (nodes > 0.01).all(axis=1) & (nodes < 0.99).all(axis=1)
f = lambda X: 2 * np.pi**2 * np.sin(np.pi * X[:, 0]) * np.sin(np.pi * X[:, 1])
g = lambda X: np.zeros(len(X))
lam = solve_poisson(nodes, interior, f, g, eps=3.0)
u_val = evaluate(np.array([[0.5, 0.5]]), nodes, lam, eps=3.0)
print("u(0.5, 0.5) ≈", u_val[0], "exact =", 1.0)
```

This sketch is intentionally minimal. A production implementation would use a structured node distribution, a boundary mask derived from the domain geometry, and a stable method for flat kernels.

---

## 13. Summary

Radial basis function methods solve PDEs by representing the solution as a sum of kernels centred at scattered nodes, then enforcing the governing equations and boundary conditions at those nodes. The approach is meshless, dimension-independent, and capable of spectral accuracy for smooth problems, making it attractive for complex geometries and high-accuracy applications. Its main challenges are the dense or ill-conditioned matrices that arise from global kernels, the sensitivity to the shape parameter, and the need for careful node and boundary treatment. Local variants such as RBF-FD and compactly supported kernels address scalability, while stabilised flat-kernel techniques address accuracy. Together these ideas make RBF methods a powerful and active area of numerical analysis.
