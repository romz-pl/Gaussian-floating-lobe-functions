# Approximation Rates of Functions on $\mathbb{R}^3$ by Gaussian Functions: Adaptive Widths and Multiscale Dictionaries

## 1. Setting and notation

Let $f:\mathbb{R}^3\to\mathbb{R}$ be a function in some smoothness or integrability class (for example $L_2(\mathbb{R}^3)$, or a Sobolev or Besov space on a bounded domain $\Omega\subset\mathbb{R}^3$). For a width parameter $\alpha>0$ and a center $y\in\mathbb{R}^3$, define the Gaussian

$$
\phi_{\alpha,y}(x) = e^{-\alpha\|x-y\|^2}, \qquad x\in\mathbb{R}^3,
$$

whose spatial scale is $\sigma=(2\alpha)^{-1/2}$. The fixed-width dictionary is

$$
\mathcal{G}_\alpha = \{\phi_{\alpha,y} : y\in\mathbb{R}^3\},
$$

and a multiscale dictionary is the union over a set of widths $A\subset(0,\infty)$:

$$
\mathcal{D}_A = \bigcup_{\alpha\in A}\mathcal{G}_\alpha.
$$

An $n$-term approximant is a linear combination of $n$ dictionary elements,

$$
f_n(x) = \sum_{k=1}^{n} c_k\,\phi_{\alpha_k,y_k}(x),
$$

where in the adaptive setting both the centers $y_k$ and the widths $\alpha_k$ are chosen depending on $f$. The best $n$-term error is

$$
\sigma_n (f)_ {X} = \inf \Big( \|f-f_n\|_{X} : f_n \text{ is an } n\text{-term combination from } \mathcal{D} \Big),
$$

with $X$ the chosen error norm, typically $L_2$ or $L_\infty$ on $\Omega$. The rate question is how fast $\sigma_n(f)_X$ decays as $n\to\infty$, and how this depends on the smoothness of $f$ and on the dictionary.

Throughout, the Fourier transform convention is $\hat g(\omega)=\int_{\mathbb{R}^3} g(x)e^{-i\omega\cdot x}\,dx$. The Fourier transform of a Gaussian is again a Gaussian:

$$
\widehat{\phi_{\alpha,0}}(\omega) = \Big(\frac{\pi}{\alpha}\Big)^{3/2} e^{-\|\omega\|^2/(4\alpha)}.
$$

This identity is the source of most of the multiscale behaviour discussed below.

---

## 2. Why a single fixed width is a poor dictionary for multiscale functions

Consider first a fixed width $\alpha_0$. The dictionary $\mathcal{G}_{\alpha_0}$ contains functions that are all smooth at one spatial scale $\sigma_0=(2\alpha_0)^{-1/2}$. Its Fourier support is effectively the ball $\|\omega\|\lesssim \sqrt{\alpha_0}$, with Gaussian decay beyond that, so the dictionary is essentially a low-pass family with cutoff frequency tied to $\alpha_0$.

Two consequences follow.

**(a) Features finer than $\sigma_0$ are expensive.** Suppose $f$ contains a localized bump of width $h \ll \sigma_0$. Representing it with Gaussians of width $\sigma_0$ requires either many terms, or coefficients of very large magnitude that cancel one another. The second route is ill-conditioned: deconvolving a Gaussian of width $\sigma_0$ amplifies high frequencies by a factor $\exp(c\,\sigma_0^{-2}h^{-2}\cdot\ldots)$ in the Fourier picture, which is numerically unusable. Hence the $n$-term error for such features stagnates at a level determined by $\sigma_0$ rather than decaying with $n$.

**(b) Features coarser than $\sigma_0$ are cheap but waste resolution.** A broad feature is represented well, but the dictionary has many nearly redundant elements at that scale if $\alpha_0$ is small relative to the feature.

The conclusion is that a dictionary with a single width has an approximation error with a saturation phenomenon: for smoothness beyond what the fixed width can represent, the rate does not improve no matter how many terms are used. Allowing the width to vary removes this obstruction.

---

## 3. Known background: variation-norm and Sobolev-type results for Gaussian networks

Two classical lines of results frame the discussion.

**Variation-norm (Barron-type) bounds.** For a fixed Gaussian dictionary, define the variation norm

$$
\|f\|_{\mathcal{G}_\alpha} = \inf\Big( \sum_k |c_k| : f = \sum_k c_k\,\phi_{\alpha,y_k} \Big)
$$

(with the infimum over finite or convergent representations). A Maurey-type sparsification argument then gives, for $f$ with finite variation norm, an $n$-term approximant with

$$
\|f - f_n\|_{L_2(\Omega)} \le \frac{C\,\|f\|_{\mathcal{G}_\alpha}}{\sqrt{n}}.
$$

This rate $n^{-1/2}$ is independent of the dimension $d=3$, which is the central attraction of this framework. The price is that the function class $\{f:\|f\|_{\mathcal{G}_\alpha}<\infty\}$ is defined through a fixed scale and is not a standard smoothness class. It contains functions that are smooth at scale $\sigma$ but says nothing about finer structure.

**Classical RBF and Gaussian-network results.** Girosi and Anzellotti studied Gaussian radial basis approximation through Fourier-based smoothness conditions, and Maiorov and Meir analysed rates for Gaussian networks in Sobolev-type classes. These works establish that fixed-width Gaussian approximation can achieve algebraic rates, but the exponents are constrained by the fixed width and by how the smoothness of $f$ compares with the dictionary's implicit smoothness. Interpreted through Section 2, the saturation effect is what these analyses reveal: smoothness of $f$ beyond the dictionary's native class does not translate into faster rates.

These results motivate the question of what changes when the width is treated as a free parameter.

---

## 4. Adaptive widths

### 4.1 The adaptive dictionary

Define the adaptive $n$-term error

$$
\sigma_n^{\mathrm{ad}}(f) = \inf\Big( \|f-f_n\|_{L_2(\Omega)} : f_n = \sum_{k=1}^n c_k\,\phi_{\alpha_k,y_k},\ \alpha_k>0,\ y_k\in\mathbb{R}^3\Big).
$$

Clearly $\sigma_n^{\mathrm{ad}}(f)\le\sigma_n(f)$ for any fixed-width dictionary, since fixed-width approximants are a special case. The adaptive problem is a nonlinear approximation problem in the sense of DeVore: the $n$-term approximants do not form a linear space, because both the centers and the widths are selected.

### 4.2 Why widths should adapt

Consider a function with features at several spatial scales, for instance

$$
f(x) = \sum_{j=0}^{J} a_j\, g\big(2^{j/2}(x-y_j)\big)
$$

where $g$ is a fixed localized bump and each $y_j$ is a center. Each term has its own width. A fixed-width dictionary needs roughly one family of terms per scale and cannot represent a term at scale $j$ with a single element of a different width. An adaptive dictionary represents each term with one element, so the $n$-term error is governed by the number $J+1$ of scales rather than by the number of spatial points needed to resolve the finest scale.

This is the basic mechanism by which varying width improves convergence: the dictionary is matched to the local scale of the function.

### 4.3 Rate heuristics for localized features

A useful benchmark is the nonlinear wavelet approximation theory on $\mathbb{R}^3$. For a Besov space $B^s_{\tau,\tau}$ with $\tau$ determined by

$$
\frac1\tau = \frac{s}{3} + \frac12,
$$

best $n$-term wavelet approximation in $L_2$ achieves

$$
\sigma_n(f)_{L_2} \lesssim n^{-s/3}.
$$

The exponent $s/3$ reflects the spatial dimension: a feature of scale $h$ in $\mathbb{R}^3$ contributes $h^{-3}$ coefficients at that resolution, and sparsity comes from the decay of coefficients across scales.

The adaptive Gaussian heuristic is that a Gaussian of width $\sigma\sim h$ centered on a localized feature can capture it with one term, so features are "sparse in the dictionary" to the same degree that they are sparse in a wavelet basis. For functions whose significant features are localized, this suggests that adaptive Gaussians can reach rates comparable to nonlinear wavelet approximation, with the potential for better constants because a Gaussian is a single smooth bump rather than a wavelet with compact support and vanishing moments.

This should be read as a heuristic. A complete rigorous statement requires specifying the admissible class of $f$, the error norm, and the selection strategy, and the constants and possible logarithmic factors need to be tracked carefully.

### 4.4 Why adaptivity helps: summary of mechanisms

Three mechanisms explain the improvement:

1. **Scale matching.** Each term can take the width of the feature it represents.
2. **Reduced conditioning problems.** Narrow features are not forced through wide Gaussians, so the coefficients do not need to cancel in an ill-conditioned way.
3. **Escape from saturation.** The effective smoothness class is no longer tied to one fixed scale, so higher smoothness can translate into faster decay of the error.

---

## 5. Multiscale dictionaries with dyadic widths

### 5.1 The dictionary

Take widths on a geometric progression, $\alpha_j = 2^j$ for $j\in\mathbb{Z}$, and let

$$
\mathcal{D}_{\mathrm{dy}} = \bigcup_{j\in\mathbb{Z}} \mathcal{G}_{2^j} = \{\,e^{-2^j\|x-y\|^2} : j\in\mathbb{Z},\ y\in\mathbb{R}^3\,\}.
$$

The spatial scale of the $j$-th family is $\sigma_j = 2^{-(j+1)/2}$, so consecutive families differ by a factor $\sqrt2$ in scale. This is a half-octave spacing in spatial scale. A dyadic-in-scale dictionary in the wavelet sense would instead use $\alpha_j = 4^j$, which gives scale ratio $2$. The choice $\alpha_j=2^j$ is therefore a finer ladder, and it gives more redundancy across scales, which can help approximation at the cost of coherence (Section 6).

### 5.2 Fourier picture and bandpass behaviour

The Fourier transform of the $j$-th Gaussian is, up to a constant,

$$
\widehat{\phi_{2^j,0}}(\omega) \propto 2^{-3j/2}\, e^{-\|\omega\|^2/(4\cdot 2^j)}.
$$

Its frequency profile is concentrated at $\|\omega\|\lesssim 2^{j/2}$. Differences of consecutive-scale Gaussians are therefore bandpass in frequency. For example, the difference

$$
\psi_j(x) = e^{-2^j\|x\|^2} - e^{-2^{j+1}\|x\|^2}
$$

has Fourier transform supported essentially in an annulus $\|\omega\|\sim 2^{j/2}$. A sum over $j$ of such differences resembles a Littlewood–Paley decomposition of $f$: each band captures the content at a dyadic frequency range, and the bands together cover the frequency axis.

This is the main analytic reason the dyadic Gaussian family behaves like a multiscale system. Band-by-band decomposition is the defining feature of wavelet frames, and Gaussians inherit it through the explicit Fourier formula above.

### 5.3 Frame-like behaviour: what holds and what does not

The statement that dyadic Gaussians "behave like wavelet frames" is best understood as an analogy with precise limits.

**What supports the analogy.**

- The scale ladder is geometric, and the bands of differences are localized in frequency, as in Littlewood–Paley theory.
- If translates are sampled on a lattice whose spacing is proportional to the scale, $\{e^{-2^j\|x-2^{-j/2}k\|^2}: k\in\mathbb{Z}^3, j\in\mathbb{Z}\}$, one expects upper and lower frame-type bounds for appropriate coefficient classes, because the lattice density matches the scale at every level. Establishing the lower bound requires controlling the Gram matrix across scales and positions.
- Sparse representations of functions with localized features exist: a function consisting of a few localized bumps at different scales has a sparse expansion.

**What does not hold in the classical sense.**

- Gaussians are not compactly supported, so they do not form a local basis in the usual sense.
- Gaussians have no vanishing moments. A wavelet with vanishing moments annihilates polynomials, and this property is essential for the sparsity of wavelet expansions of smooth functions. A positive Gaussian does not annihilate anything. Smooth functions therefore do not benefit from the same cancellation that makes wavelets efficient for them. Differences of Gaussians recover some of this (their integrals against low-order polynomials are small), but they are not an exact replacement.
- The dictionary is continuous in translation and therefore highly redundant. A frame in the strict sense is a discrete object, and the continuous family must be discretized before frame theory applies.

The honest summary is that dyadic Gaussians provide a multiscale, approximately bandpass, redundant family with frame-like stability under lattice sampling, and this is sufficient to explain good approximation for localized features. They are not a wavelet frame in the strict construction sense.

### 5.4 Rates for localized features with dyadic dictionaries

Consider a function that is a finite or summable sum of localized features,

$$
f = \sum_{m} b_m\, \theta\big(2^{j_m/2}(x-z_m)\big),
$$

where $\theta$ is a fixed smooth localized bump, $j_m$ are scale indices and $z_m$ are locations. If the dictionary contains $\theta$-like elements at the matching scales and locations, each feature is captured by one term. The $n$-term error for such a function is then the tail of its coefficient sequence, and the decay rate of that tail is governed by the smoothness and localization of $f$ in the same way as in wavelet theory.

For $f$ in a Besov-type class with the same scaling as the wavelet benchmark of Section 4.3, the heuristic rate is $\sigma_n(f)_{L_2}\lesssim n^{-s/3}$, possibly with a logarithmic factor from the redundancy of the dictionary in scale. The careful statement requires an analysis of the best $n$-term selection in the presence of coherence, which is discussed below.

---

## 6. Coherence and conditioning across scales

A central practical and theoretical issue is the Gram structure of the dictionary. For two Gaussians with the same center and widths $\alpha$ and $\beta$, the normalized inner product is

$$
\frac{\langle \phi_{\alpha,y},\phi_{\beta,y}\rangle}{\|\phi_{\alpha,y}\|\,\|\phi_{\beta,y}\|} = \Big(\frac{2\sqrt{\alpha\beta}}{\alpha+\beta}\Big)^{3/2}.
$$

For $\beta = 2\alpha$ (widths in a dyadic ladder with ratio 2), this equals $(2\sqrt2/3)^{3/2}\approx 0.915$. For the half-octave ladder $\alpha_j = 2^j$, the ratio between neighbouring widths is $2$ in $\alpha$ but $\sqrt2$ in scale, and the same formula with $\beta=2\alpha$ applies. Neighbouring-scale Gaussians at the same center are therefore highly correlated. This is the coherence problem.

Consequences:

- **Greedy selection** (matching pursuit or orthogonal matching pursuit) can select similar atoms across adjacent scales, which slows convergence and can cause the algorithm to choose a suboptimal scale.
- **Coefficient instability**: near-collinear atoms produce large, cancelling coefficients unless the least-squares step is well regularized.
- **Frame constants** deteriorate: the lower frame bound decreases as the density of scales increases, unless the translation lattice is refined correspondingly.

These effects are partly offset by the fact that the correlation decays rapidly with separation in position or in $\log\alpha$ beyond one or two steps. Practical methods therefore use:

1. Sparse scale grids (a few widths per octave) combined with fine positional sampling near features.
2. Orthogonalization of coarse-to-fine residuals, so each scale contributes its innovation.
3. Ridge or Tikhonov regularization in the least-squares step.

---

## 7. Algorithms for adaptive Gaussian approximation

The theory above suggests the following algorithmic structure.

**Greedy residual scheme.** Start with $f_0=0$ and residual $r_0=f$. At step $k$, solve

$$
(y_k,\alpha_k) = \arg\max_{y,\alpha}\ \frac{|\langle r_{k-1},\phi_{\alpha,y}\rangle|}{\|\phi_{\alpha,y}\|},
$$

then refit all coefficients by least squares on $\{\phi_{\alpha_i,y_i}\}_{i\le k}$, and update $r_k$. This is orthogonal matching pursuit over a continuous dictionary. The inner product maximization is a nonconvex problem in $(y,\alpha)$, usually handled by a coarse grid followed by local refinement.

**Coarse-to-fine multiscale scheme.** Compute the residual at scale $j$ by convolution with the Gaussian $\phi_{2^j,0}$, which can be done with fast Fourier transforms. Use the band $\psi_j$ of Section 5.2 to isolate the content at each dyadic frequency range, and select centers where the band coefficient is large. This mirrors the thresholding step in wavelet compression.

**Scale-space interpretation.** Since $\phi_{2^j,0}$ is the heat kernel at time $\propto 2^{-j}$ (up to normalization), $f\ast\phi_{2^j,0}$ is a scale-space representation of $f$. The bands $f\ast\psi_j$ are differences of scale-space levels, which is the classical multiscale edge and blob detection construction. This gives a concrete numerical interpretation of the dictionary.

---

## 8. Summary of rates and regimes

| Regime | Dictionary | Expected $n$-term behaviour in $L_2$ | Comment |
|---|---|---|---|
| Fixed width, smooth $f$ at scale $\sigma_0$ | $\mathcal{G}_{\alpha_0}$ | Barron-type $n^{-1/2}$ for variation-norm-bounded $f$ | Dimension-free, but saturates for features finer than $\sigma_0$ |
| Fixed width, finer features | $\mathcal{G}_{\alpha_0}$ | Stagnation at a level set by $\sigma_0$ | Deconvolution is ill-conditioned |
| Adaptive width, localized features | $\bigcup_\alpha \mathcal{G}_\alpha$ | Heuristically $n^{-s/3}$ for Besov-type $f$ | Matches nonlinear wavelet benchmark up to logs; rigorous constants need care |
| Dyadic multiscale, $\alpha_j=2^j$ | $\mathcal{D}_{\mathrm{dy}}$ | Same exponent, possibly with a log factor from redundancy | Frame-like behaviour under lattice sampling; coherence across adjacent scales |

---

## 9. Caveats and open points

- The $n^{-s/3}$ benchmark is a heuristic transferred from wavelet theory. A full proof requires specifying the function class, the selection rule, and handling the nonconvexity of joint $(y,\alpha)$ selection.
- The frame analogy is qualitative. Precise frame bounds depend on the lattice density chosen for each scale and on the decay of the Gram matrix.
- Gaussians lack vanishing moments, so smooth functions without localized features do not automatically enjoy wavelet-like sparsity. Differences of Gaussians partially recover this but do not fully replace it.
- Coherence between adjacent scales is the main obstacle to stable adaptive selection, and it motivates sparser scale grids in practice.
- The choice $\alpha_j=2^j$ versus $\alpha_j=4^j$ changes the scale ratio from $\sqrt2$ to $2$. Which ladder gives the best trade-off between redundancy and coherence is an empirical and theoretical question that depends on the target class.

---

## 10. Conclusion

Allowing the Gaussian width to vary is the key step that moves approximation by Gaussians out of the saturation regime of a fixed-scale dictionary. Each term can be matched to the local scale of a feature, which turns multiscale structure into sparsity. Dyadic families $\alpha_j=2^j$ make this multiscale structure explicit: their Fourier transforms are Gaussians, their differences are approximately bandpass, and their lattice-sampled versions behave like frames in the Littlewood–Paley sense. These properties heuristically give rates comparable to nonlinear wavelet approximation for functions with localized features on $\mathbb{R}^3$. The rigorous details, including constants, logarithmic factors and coherence control, remain the places where careful analysis is required.
