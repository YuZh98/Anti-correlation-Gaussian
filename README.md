# Anti-correlation-Gaussian
Codes for the paper "Gibbs Sampling using Anti-correlation Gaussian Data Augmentation, with Applications to L1-ball-type Models" https://arxiv.org/abs/2309.09371

```SourceCode.R``` and ``SliceSampler.R``` are two source files. The former defines a function using the anti-correlation Gaussian data augmentation trick for linear regression. We update the threshold kappa using a slice sampler, so the latter source file is necessary.

```LinearReg.R``` shows a toy example: linear regression. It will also plot a comparison between the posterior estimation and ground truth, trace plots of the parameters involved, and ACF plots.

```TruncatedMVN.R``` shows another toy example which is an extension of the anti-correlation Gaussian trick. We aim to sample from a truncated multivariate normal distribution with a separable truncation. The file contains a 2-dimensional toy example and it will plot the marginal densities and the joint density.
