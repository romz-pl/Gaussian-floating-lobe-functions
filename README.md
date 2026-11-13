# Off-Center Gaussian Basis Expansions in the Hartree-Fock Method

## Scope and Caveats

This document combines two overlapping descriptions: the **Gaussian lobe function expansion** and the **Gaussian floating function expansion**. Both are families of off-center Gaussian basis functions used in the LCAO-Hartree-Fock framework. The term "Gaussian lobe function" is not tied to a single canonical reference in my knowledge, so the lobe material describes the concept as it is generally understood. Historical attributions for the floating-Gaussian methods are given at a general level and should be checked against primary sources before citation.

## Motivation

Standard Hartree-Fock calculations expand molecular orbitals in atom-centered Gaussian-type orbitals (GTOs),

$$
\phi_\mu(\mathbf r) = (x-A_x)^{l}(y-A_y)^{m}(z-A_z)^{n}\, e^{-\alpha |\mathbf r - \mathbf A|^2}.
$$

Atom-centered sets are systematic and well understood, but they have two practical limitations:

1. **Shape representation.** A single s-type Gaussian is spherically symmetric. Describing p, d, or bonding-directed density requires either higher angular momentum functions or many contracted primitives.
2. **Density between nuclei.** Features such as bond midpoints, lone pairs, and diffuse inter-atomic density do not coincide with nuclei. Capturing them with atom-centered functions can require many additional functions.

Both families of methods address these problems by moving Gaussian centers away from the nuclei. The lobe approach uses displacement to build directional shape. The floating approach treats the centers as variational parameters that the energy determines.

## The Basic Functions

### Floating spherical Gaussians

A primitive s-type Gaussian centered at an arbitrary point $\mathbf R$ is

$$
g(\mathbf r; \alpha, \mathbf R) = e^{-\alpha |\mathbf r - \mathbf R|^2}.
$$

In the floating-Gaussian approach, both the exponent $\alpha$ and the center $\mathbf R$ are free parameters. The center is not tied to any nucleus. Because the function is spherically symmetric, it is best suited to density that is roughly isotropic about its center, such as electron pairs in bonds or core-like regions. Some implementations use floating functions with p-type or higher angular momentum, but the spherical form is the classic starting point because its integrals are the simplest.

### Gaussian lobe functions

A lobe function is an s-type Gaussian whose center is displaced from a nucleus along a chosen direction by a distance $a$, which is either fixed or optimized. For a p-like lobe pair along the $z$ axis, a common construction is the antisymmetric combination

$$
\chi_{p_z}(\mathbf r) = g(\mathbf r; \alpha, \mathbf A + a\hat{z}) - g(\mathbf r; \alpha, \mathbf A - a\hat{z}).
$$

This function has a nodal plane through the nucleus and two lobes of opposite sign, approximating a $p_z$ orbital. Symmetric pairs or varied displacement magnitudes approximate d-like or more diffuse angular shapes. Each lobe is still a single s-type Gaussian, so the integral machinery is unchanged.

The two constructions share a single principle. The lobe method fixes the displacement direction from the desired orbital shape, while the floating method lets the displacement be optimized freely. A lobe basis is therefore a directed, constrained form of a floating basis.

## Expansion and the Roothaan Equations

In both methods, spatial orbitals are expanded in the basis:

$$
\psi_i(\mathbf r) = \sum_{\mu} c_{\mu i}\, \chi_\mu(\mathbf r).
$$

Substituting into the Hartree-Fock equations and projecting onto the basis yields the Roothaan-Hall equations

$$
\mathbf F \mathbf C = \mathbf S \mathbf C \boldsymbol\varepsilon,
$$

with overlap matrix

$$
S_{\mu\nu} = \int \chi_\mu^\ast(\mathbf r)\, \chi_\nu(\mathbf r)\, d\mathbf r,
$$

and Fock matrix

$$
F_{\mu\nu} = h_{\mu\nu} + \sum_{\lambda\sigma} P_{\lambda\sigma}\left[(\mu\nu|\lambda\sigma) - \tfrac{1}{2}(\mu\lambda|\nu\sigma)\right].
$$

For a closed-shell system, $P_{\lambda\sigma} = 2\sum_{i\,\text{occ}} c_{\lambda i} c_{\sigma i}$. The core Hamiltonian $h_{\mu\nu}$ contains kinetic and nuclear attraction terms, and $(\mu\nu|\lambda\sigma)$ denotes two-electron repulsion integrals.

The SCF cycle is the usual one: build $\mathbf F$ from the current density, solve the generalized eigenvalue problem, form a new density, and iterate to convergence. The differences lie in the basis, the integrals, and, for the floating approach, the outer optimization of nonlinear parameters.

For the lobe method, the nonlinear parameters are typically the centers and exponents chosen by a geometric rule or optimized directly. For the floating method, the energy is minimized with respect to both the linear coefficients $c_{\mu i}$ and the nonlinear parameters $\{\alpha_\mu, \mathbf R_\mu\}$. At fixed nonlinear parameters, the coefficients satisfy the Roothaan-Hall equations above.

## Integral Evaluation

Because every basis function is an s-type Gaussian, the central tool is the Gaussian product theorem. For two primitives centered at $\mathbf A$ and $\mathbf B$,

$$
\exp(-\alpha |\mathbf r - \mathbf A|^2) \exp(-\beta |\mathbf r - \mathbf B|^2) = 
K_{AB} \exp(-\gamma|\mathbf r - \mathbf P|^2),
$$

where

$$
\gamma = \alpha + \beta
$$

$$
\mathbf P = (\alpha \mathbf A + \beta \mathbf B) / \gamma
$$

$$
K_{AB} = \exp \left( -\frac{\alpha\beta}{\gamma}|\mathbf A - \mathbf B|^2 \right)
$$

Overlap and kinetic integrals have closed forms. Nuclear attraction and two-electron repulsion integrals reduce to the Boys function $F_m(T)$ after the product theorem is applied. Because all functions are s-type, the recursion relations are short, and the integral code is simpler and faster than for general Cartesian or spherical Gaussians. The cost grows with the number of functions and with the number of distinct centers in multi-center integrals.

## Gradients with Respect to Centers

For the floating approach, the energy gradient with respect to a center $\mathbf R_\mu$ is needed for efficient optimization. For an s-type Gaussian,

$$
\frac{\partial g}{\partial \mathbf R_\mu} = 2\alpha_\mu (\mathbf r - \mathbf R_\mu) g
$$

and the resulting derivative integrals can be evaluated with the same product machinery using shifted indices. The gradient also includes contributions from the non-orthogonality of the basis, which enter through the overlap matrix and the energy-weighted density matrix. These terms make analytic gradients harder to implement than for fixed-center methods, and they are a main reason floating-Gaussian methods were historically more difficult to program.

The lobe method, with centers tied to a geometric rule, needs gradients only if the displacements are optimized, and then the same expressions apply.

## Parameter Choice and Optimization

Both families involve a few choices that define a specific implementation:

- **Placement.** Centers may be fixed by geometry, such as bond midpoints, bond directions, or lone-pair axes, or optimized as nonlinear parameters. Good starting guesses at chemically meaningful points are common practice.
- **Displacement and exponents.** These control how compact or diffuse each function is. Their optimization is nonlinear.
- **Linear coefficients.** Once nonlinear parameters are fixed, the coefficients follow from the SCF procedure.

Floating-Gaussian calculations are typically organized in two nested loops. The inner loop is the SCF procedure at fixed nonlinear parameters. The outer loop adjusts exponents and centers to lower the total energy, using gradient-based or direct-search methods. Some workflows add or remove functions during optimization, guided by energy changes or by the size of the coefficients.

## Advantages

- The basis can be compact and directional (lobes) or compact and adaptive (floating) for a given accuracy.
- Only s-type integrals are needed in the basic form, which simplifies implementation.
- Density between atoms and anisotropic density can be captured without high angular momentum functions.
- Floating functions can be placed where the energy prefers them, rather than where a rule dictates.

## Difficulties

- **Linear dependence.** Off-center functions are not orthogonal to each other or to atom-centered functions. The overlap matrix can become ill-conditioned when functions overlap strongly, requiring canonical orthogonalization or pruning.
- **Nonlinear optimization.** Optimization is expensive and can converge to local minima. Results depend on starting configurations.
- **Gradient complexity.** Analytic gradients for floating centers are more involved than for fixed-center methods.
- **Systematic convergence.** Neither family has a hierarchy comparable to correlation-consistent basis sets, so convergence toward the basis-set limit is harder to control.
- **Integral cost.** Multi-center integrals involving many displaced functions can become costly.

## Historical Context

The floating-Gaussian idea grew out of early work on Gaussian functions in molecular integrals, with contributions associated with Boys and Bingel in the 1950s and 1960s. The floating spherical Gaussian orbital (FSGO) approach, developed by Frost and coworkers in the 1960s, is the best-known application, using spherical floating functions for valence description in semi-empirical and minimal-basis contexts. The lobe approach is less standardized, and I cannot give a reliable single attribution for it.

## Relation to Other Approaches

Off-center Gaussian methods sit between two extremes. Atom-centered GTO sets are systematic and standard for routine quantum chemistry. Fully variational schemes that treat all positions as free parameters are more flexible but harder to control. The lobe method can be viewed as a directed floating basis, in which displacement follows a chosen orbital shape. The classic FSGO method uses spherical floating functions with no preferred direction.

In practice, off-center Gaussian methods have been applied mainly in specialized studies and in minimal or semi-empirical contexts, rather than as general replacements for standard Gaussian basis sets. Their practical value depends on whether the reduction in function count outweighs the cost of nonlinear optimization and the numerical difficulties of non-orthogonal, off-center functions.

## Bibliography

### 1950

+ S. F. Boys, **Electronic Wave Functions. I. A General Method of Calculation for the Stationary States of Any Molecular System**, [Proc. R. Soc. Lond. A. 200 542–554 (1950)](https://doi.org/10.1098/rspa.1950.0036)

### 1956

+ H. Preuss, **Über die Verwendung von Gauß-Funktionen in der Wellenmechanik (On the Use of Gaussian Functions in Wave Mechanics)**, [Zeitschrift für Naturforschung A, Volume 11a, Issue 10, Pages 823–831 (1956)](https://doi.org/0.1515/zna-1956-1007)


### 1963

+ J.L. Whitten, **Gaussian Expansion of Hydrogen Atom Wavefunctions**, [J. Chem. Phys. 39, 349–352 (1963)](https://doi.org/10.1063/1.1734251)


### 1965

+ H. Sambe **Use of 1s Gaussian Wavefunctions for Molecular Calculations. I. The Hydrogen Atom and the Hydrogen Molecule Ion**, [J. Chem. Phys. 42, 1732–1738 (1965)](https://doi.org/10.1063/1.1696185)


### 1966

+ J. L. Whitten, **Gaussian Lobe Function Expansions of Hartree—Fock Solutions for the FirstRow Atoms and Ethylene**, [J. Chem. Phys. 44, 359 (1966)](https://doi.org/10.1063/1.1726470)

+ S. D. Peyerimhoff, R. J. Buenker, L. C. Allen, **Geometry of Molecules. I. Wavefunctions for Some Six‐ and Eight‐Electron Polyhydrides**, [J. Chem. Phys. 45, 734–749 (1966)](https://doi.org/10.1063/1.1727634)



### 1967

+ M. E. Schwartz, L. J. Schaad, **Ab Initio Studies of Small Molecules Using 1s Gaussian Basis Functions. I. Exploratory Calculations**, [J. Chem. Phys. 46, 4112–4120 (1967)](https://doi.org/10.1063/1.1840493)


+ J. F. Harrison **On the Gaussian‐Lobe Representation of Atomic Orbitals**, [J. Chem. Phys. 46, 1115–1118 (1967)](https://doi.org/10.1063/1.1840776)

+ A. A. Frost, **Floating Spherical Gaussian Orbital Model of Molecular Structure. I. Computational Procedure. LiH as an Example**, [J. Chem. Phys. 47, 3707 (1967)](https://doi.org/10.1063/1.1701524)

+ A. A. Frost, **Floating Spherical Gaussian Orbital Model of Molecular Structure. II. One and Two Electron Pair Systems**, [J. Chem. Phys. 47, 3714 (1967)](https://doi.org/10.1063/1.1701525)


### 1969

+ J. D. Petke, J. L. Whitten, A. W. Douglas, **Gaussian lobe function expansions of Hartree–Fock solutions for the second-row atoms**, [J. Chem. Phys. 51, 256–262 (1969)]( https://doi.org/10.1063/1.1671716)


### 1997

+ A.H. Pakiari, M. Oftadeh, **First derivative of molecular integrals for electronic energy using spherical gaussian type orbitals**, [Journal of Molecular Structure (Theochem) 389, 169-177 (1997)](https://doi.org/10.1016/S0166-1280(96)04727-6)


### 2001

+ S. Garashchuk, J. C. Light, **Quasirandom distributed Gaussian bases for bound problems**, [J. Chem. Phys. 114, 3929–3939 (2001)](https://doi.org/10.1063/1.1348022)




