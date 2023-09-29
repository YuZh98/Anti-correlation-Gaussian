# Anti-correlation-Gaussian
Codes for the paper [Gibbs Sampling using Anti-correlation Gaussian Data Augmentation, with Applications to L1-ball-type Models](https://arxiv.org/abs/2309.09371).

```SourceCode.R``` and ```SliceSampler.R``` are two **source files**. The former defines a function using the anti-correlation Gaussian data augmentation trick for linear regression. The latter uses the stepping-out and the shrinkage techniques in [Slice Sampling](https://arxiv.org/pdf/physics/0009028.pdf) by Radford M. Neal. You can source it if a slice sampler is needed in your project. It works both for scalar and multi-dimensional cases with a given direction.

```LinearReg.R``` contains a toy example: variable selection in linear regression. It will plot a comparison between the posterior estimation and ground truth, trace plots of the parameters involved, and ACF plots.

```TruncatedMVN.R``` contains an extension of the anti-correlation Gaussian trick. We aim to sample from a truncated multivariate normal distribution with a separable truncation (meaning the constraint on the random vector is a high-dimensional box). The file contains a 2-dimensional toy example and it will plot the marginal densities and the joint density. You can increase the dimension ```p``` to suit your needs.



