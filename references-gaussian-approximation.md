# References: Gaussian Approximation Rates, Adaptive Widths and Multiscale Dictionaries

Entries marked **[verified]** were checked against the source text or abstract during preparation. Entries marked **[from memory]** are standard citations I could not re-check in this session; confirm details (year, volume, pages) before citing them in print.

---

## 1. Universal approximation and Barron-type rates

**Barron, A. R. (1993).** Universal approximation bounds for superpositions of a sigmoidal function. *IEEE Transactions on Information Theory*, 39(3), 930–945. DOI: 10.1109/18.256500. **[verified]**
The foundational dimension-independent result: networks with one hidden layer of sigmoidal units achieve integrated squared error of order $O(1/n)$ (so $O(n^{-1/2})$ in norm) when the target has a bounded first moment of its Fourier magnitude distribution. This is the origin of the variation-norm argument used in Section 3 of the write-up.

**Barron, A. R. (1994).** Approximation and estimation bounds for artificial neural networks. *Machine Learning*, 14, 115–133. **[verified]**
Extends the approximation bound to a statistical estimation setting and computes the index of resolvability for minimum-complexity estimation over these networks.

**Jones, L. K. (1992).** A simple lemma on greedy approximation in Hilbert space and convergence rates for projection pursuit regression and neural network training. *Annals of Statistics*, 20(1), 608–613. **[from memory]**
The greedy Hilbert-space lemma underlying Barron's rate; the basis for the greedy algorithms in Section 7.

**Maurey, B.** (unpublished, cited by Pisier, 1981). **[verified as cited in Šámalová 2008, see §5]**
The sparsification idea behind the $1/\sqrt{n}$ bound for convex hulls in Banach spaces.

**Klusowski, J. M., & Barron, A. R. (2016).** Approximation by combinations of ReLU and squared ReLU ridge functions with $\ell^1$ and $\ell^0$ controls. *IEEE Transactions on Information Theory*, 64(12), 7649–7656. **[from memory]**
Modern sharpening of Barron-type bounds; useful when comparing sigmoidal and ReLU-type ridge dictionaries.

**Siegel, J. W., & Xu, J. (2020/2022).** Approximation rates for neural networks with general activation functions. *Neural Networks* (2020) and related arXiv preprints. **[from memory]**
Rates for shallow networks under general smoothness assumptions.

**Siegel, J. W. (2020 preprint, revised later).** High-order approximation rates for shallow neural networks with cosine and ReLU$^k$ activation functions. arXiv:2012.07205. **[verified]**
Extends rates in the spectral Barron space to higher smoothness $s$, with improved $n^{-(k+1)}$-type rates for ReLU$^k$ and lower bounds showing optimality under the stated assumptions. Good for context on how smoothness beyond Barron's class translates into rates.

---

## 2. Gaussian and radial basis function approximation with fixed width

**Girosi, F., & Anzellotti, G. (1993).** Rates of convergence for radial basis functions and neural networks. In *Artificial Neural Networks for Speech and Vision* (R. J. Mammone, ed.), Chapman & Hall, 97–113. **[from memory]**
Also circulated as MIT AI Memo 1288 ("Rates of convergence of approximation by translates", 1992). Gives $O(n^{-1/2})$ bounds for approximation by translates of a fixed kernel on band-limited or Sobolev-type classes. This is the fixed-width baseline against which adaptive widths are compared. Verified indirectly via citations in Kainen, Kůrková & Sanguineti (2008) and the MIT course notes listed below.

**Girosi, F., & Anzellotti, G. (1995).** Convergence rates of approximation by translates. *Proceedings of ICANN 1995*, Vol. 1, 295–302, Paris. **[verified as cited]**
The conference version of the AI Memo above. Cited by Recht (MIT 9.520 notes) and Kainen et al.

**Maiorov, V., & Meir, R. (1998).** On the near optimality of the stochastic approximation of smooth and analytic functions by neural networks. *IEEE Transactions on Neural Networks*, 9(5)? **[from memory — verify volume and pages]**
Establishes near-optimal rates for Gaussian and sigmoidal networks in Sobolev-type classes. The title and journal are from memory and should be confirmed.

**Mhaskar, H. N. (1996).** Neural networks for optimal approximation of smooth and analytic functions. *Neural Computation*, 8(1), 164–177. **[verified as cited]**
Rates for smooth and analytic functions by networks with fixed-width units; contrast point for the adaptive discussion.

**Kainen, P. C., Kůrková, V., & Sanguineti, M. (2006).** Rates of approximation of smooth functions by Gaussian radial-basis networks. Technical Report 976, Institute of Computer Science, Academy of Sciences of the Czech Republic. **[verified]**
Earlier technical report version of the 2009 paper below.

**Kainen, P. C., Kůrková, V., & Sanguineti, M. (2009).** Complexity of Gaussian-radial-basis networks approximating smooth functions. *Journal of Complexity*, 25(1), 63–74. **[verified]**
Central reference for varying-width Gaussian networks. Gives upper bounds on $L_2$ approximation rates for smooth functions by networks whose units have varying widths, with bounds expressed in Sobolev and Bessel-potential norms, and discusses tractability. Directly supports Section 4 of the write-up.

**Kůrková, V., & Sanguineti, M. (2002).** Comparison of worst case errors in linear and neural network approximation. *IEEE Transactions on Information Theory*, 48(1), 264–275. **[from memory]**
Shows that for some function classes nonlinear (neural) approximation outperforms linear approximation in worst-case error.

**Krzyżak, A., & Linder, T. (1998).** Radial basis function networks and complexity regularization in function learning. *IEEE Transactions on Neural Networks*, 9(2), 247–256. **[verified]**
Covers RBF networks in estimation settings, including approximation error rates for RBF networks with complexity regularization. Useful for the statistical side of the story.

**Lin, S., Liu, X., Rong, Y., & Xu, Z. (2014).** Almost optimal estimates for approximation and learning by radial basis function networks. *Machine Learning*, 2014. doi:10.1007/S10994-013-5406-Z. **[verified]**
Near-optimal approximation and learning estimates for RBF networks, including Gaussian kernels. Good for checking the current state of fixed-kernel rates.

---

## 3. Nonlinear approximation by Gaussians with adaptive widths (core references)

**Hangelbroek, T., & Ron, A. (2010).** Nonlinear approximation using Gaussian kernels. *Journal of Functional Analysis*, 259(1), 203–219. DOI: 10.1016/j.jfa.2010.02.001. arXiv:0911.2803. **[verified]**
The most directly relevant reference for this topic. Introduces an $N$-term scheme using translates of Gaussians with spatially varying tension (width) parameters, chosen from local information about the target. Shows approximation rates comparable to spline and wavelet nonlinear approximation, with rates $N^{-s/d}$ in Besov and Triebel–Lizorkin settings, "as high as needed," saturated only by the smoothness of the approximand. Their scheme resolves local singularities non-locally. This paper supports the adaptive-width claims in Sections 4.3 and 8 of the write-up and corrects my earlier heuristic: the $n^{-s/3}$ rate in $\mathbb{R}^3$ is the right benchmark and has been established in this setting, not only heuristically. Read this paper before anything else in the list.

**DeVore, R. A., & Ron, A. (2010).** Approximation using scattered shifts of a multivariate function. *Transactions of the American Mathematical Society*, 362(12), 6205–6229. **[from memory — verify]**
Nonlinear approximation with homogeneous radial kernels (surface splines), whose $N$-term spaces are invariant under rescaling. Hangelbroek and Ron note that the DeVore–Ron kernels saturate, which motivates the Gaussian result.

**DeVore, R. A. (1998).** Nonlinear approximation. *Acta Numerica*, 7, 51–150. **[from memory]**
The standard survey of nonlinear approximation theory, including $N$-term approximation, Besov-space characterizations and the wavelet benchmark $N^{-s/d}$ used in Section 4.3.

**Cohen, A., Dahmen, W., & DeVore, R. (2001).** Adaptive wavelet methods for elliptic operator equations: convergence rates. *Mathematics of Computation*, 70(233), 27–75. **[from memory]**
Adaptive wavelet frame/basis schemes with rate guarantees; the wavelet-side comparison class.

---

## 4. Multiscale dictionaries, Gaussians, and frames

**Hangelbroek, T., & Ron, A.** (see §3). The construction uses multiple scales and selects centers and tensions locally; it is the closest published analogue to the dyadic dictionary $\alpha_j = 2^j$ discussed in Section 5 of the write-up.

**Author list not confirmed (2012 preprint, Friedrich-Schiller-Universität Jena, Novak group).** Nonlinear frame widths for elliptic equations. Available at users.fmi.uni-jena.de/~novak/DNS12.pdf. **[content verified from the excerpt; read the PDF's title page for authors and publication details before citing]**
Shows that for nonlinear best $n$-term approximation with frames, the approximation order matches the nonlinear widths associated with Riesz bases, and that frames give rates better than linear widths when $p<2$. Useful to make precise the sense in which "frame-like" dictionaries inherit wavelet-type rates.

**Daubechies, I. (1992).** *Ten Lectures on Wavelets*. SIAM, CBMS-NSF Series, 61. **[from memory]**
Standard reference for wavelet frames, Littlewood–Paley decomposition and the frame bounds needed to make the analogy in Section 5.3 precise.

**Christensen, O. (2016).** *An Introduction to Frames and Riesz Bases* (2nd ed.). Birkhäuser. **[from memory]**
Frame theory in Hilbert spaces, including the discrete-versus-continuous frame issue raised in the write-up.

**Mallat, S. (2009).** *A Wavelet Tour of Signal Processing: The Sparse Way* (3rd ed.). Academic Press. **[from memory]**
Standard text on sparse approximation and multiscale representations; a good reference for the Littlewood–Paley and difference-of-Gaussians reasoning.

**Lindeberg, T. (1994).** *Scale-Space Theory in Computer Vision*. Kluwer. **[from memory]**
Gaussian scale-space and differences of Gaussians as bandpass filters; supports the scale-space interpretation in Section 7.

---

## 5. Background on coherence and sparse approximation

**Tropp, J. A. (2004).** Greed is good: algorithmic results for sparse approximation. *IEEE Transactions on Information Theory*, 50(10), 2231–2242. **[from memory]**
Greedy selection and coherence conditions for sparse recovery; relevant to Section 6 on coherence across scales.

**Pati, Y. C., Rezaiifar, R., & Krishnaprasad, P. S. (1993).** Orthogonal matching pursuit: recursive function approximation with applications to wavelet decomposition. *Proceedings of the 27th Asilomar Conference on Signals, Systems and Computers*, 40–44. **[from memory]**
The orthogonal matching pursuit algorithm referenced in Section 7.

**Mallat, S. G., & Zhang, Z. (1993).** Matching pursuits with time-frequency dictionaries. *IEEE Transactions on Signal Processing*, 41(12), 3397–3415. **[from memory]**
The original matching pursuit paper, relevant to greedy selection over redundant dictionaries.

---

## 6. Surveys and books

**Pinkus, A. (1999).** Approximation theory of the MLP model in neural networks. *Acta Numerica*, 8, 143–195. **[from memory]**
Survey covering density, rates and the characterization of activation functions; referenced as the source for the Leshno et al. theorem in the approximation-theory survey below.

**Leshno, M., Lin, V. Y., Pinkus, A., & Schocken, S. (1993).** Multilayer feedforward networks with a nonpolynomial activation function can approximate any function. *Neural Networks*, 6(6), 861–867. **[from memory]**
Universal approximation characterization; background for the density results in Section 1.

**Cybenko, G. (1989).** Approximation by superpositions of a sigmoidal function. *Mathematics of Control, Signals and Systems*, 2(4), 303–314. **[from memory]**
The original universal approximation theorem for sigmoidal networks.

**Hornik, K., Stinchcombe, M., & White, H. (1989).** Multilayer feedforward networks are universal approximators. *Neural Networks*, 2(5), 359–366. **[from memory]**
Extension of universal approximation to $L^p(\mu)$ and general activations.

**Wendland, H. (2005).** *Scattered Data Approximation*. Cambridge University Press. **[from memory]**
Standard reference on radial basis functions, native spaces and approximation orders; useful background for Gaussian kernels in the interpolation setting.

**Buhmann, M. D. (2003).** *Radial Basis Functions: Theory and Implementations*. Cambridge University Press. **[from memory]**
Comprehensive RBF text covering Gaussian kernels, shape parameters, and the trade-off between width and conditioning.

**Fasshauer, G. E. (2007).** *Meshfree Approximation Methods with MATLAB*. World Scientific. **[from memory]**
Practical treatment of Gaussian kernels with variable shape parameters, relevant to the stability issues in Section 6.

**Rasmussen, C. E., & Williams, C. K. I. (2006).** *Gaussian Processes for Machine Learning*. MIT Press. **[from memory]**
Gaussian kernels with length-scale hyperparameters; background for the connection between width selection and kernel learning.

---

## 7. Course notes and supplementary sources

**Recht, B. (2008).** Lecture notes, MIT 9.520 "Statistical Learning Theory and Applications", lecture on approximation and rates of convergence (dated 4 July 2008). Available at mit.edu/~9.520/spring08/Classes/recht_040708.pdf. **[verified]**
Covers Maurey–Barron–Jones approximation schemes and "dimension-free" convergence, with Girosi–Anzellotti cited for the Gaussian case.

**Niyogi, P., & Girosi, F. (1996).** On the relationship between generalization error, hypothesis complexity, and sample complexity for radial basis functions. *Neural Computation*, 8(4), 819–842. **[verified as cited]**
Joint approximation and estimation trade-offs for RBF networks.

**Niyogi, P., & Girosi, F. (1999).** Generalization bounds for function approximation from scattered noisy data. *Advances in Computational Mathematics*, 10, 51–80. **[verified as cited]**

**Šámalová, T. (2008).** Rates of approximation by neural networks with sigmoidal activation functions and $G$-variation. Technical Report 1021, Institute of Computer Science, Czech Academy of Sciences. **[verified]**
Discusses $G$-variation as a sufficient condition for good approximation rates and the limits of variation-based estimates.

**Lin, S., et al. (2014).** (see §2.)

---

## Notes for the write-up

1. **Correction to the main document.** Section 4.3 of the write-up describes the $n^{-s/3}$ rate for adaptive Gaussians as a heuristic. Hangelbroek and Ron (2010) establish rates of this type for Gaussians with varying tension, in Besov and Triebel–Lizorkin settings, so the statement can now be cited rather than conjectured. The write-up should be revised to cite that paper and to state the result with its hypotheses.

2. **Possible error to check.** Section 3 of the write-up cites Maiorov and Meir for rates of Gaussian networks in Sobolev-type classes. The title and venue in §2 above are from memory and should be verified before the write-up is circulated.

3. **Coverage gaps.** I did not locate a published result that establishes the "frame-like" behaviour of dyadic Gaussian dictionaries $\alpha_j = 2^j$ in the strict sense. The Hangelbroek–Ron construction and the Novak frame-widths paper are the closest supporting results. If a precise frame statement is needed, search the arXiv listings for "Gaussian frames" and "Gaussian multiresolution" after 2010.
