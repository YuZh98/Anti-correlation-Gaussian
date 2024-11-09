# Anti-correlation-Gaussian

This repository contains the code and the funnctional magnetic resonance
imaging (fMRI) data for the paper [*Gibbs Sampling using Anti-correlation Gaussian Data Augmentation, with Applications to L1-ball-type Models*](https://arxiv.org/abs/2309.09371).

## Overview of Files

- **`SourceCode.R`**: Defines a function that implements the anti-correlation Gaussian data augmentation technique for linear regression models. This file provides the core functionality for the methods described in the paper.

- **`SliceSampler.R`**: Provides an implementation of slice sampling based on the stepping-out and shrinkage techniques detailed in [*Slice Sampling*](https://arxiv.org/pdf/physics/0009028.pdf) by Radford M. Neal. This code can be sourced if a slice sampler is required for your project and supports both scalar and multi-dimensional cases with a specified direction.

- **`LinearReg.R`**: Demonstrates a toy example of variable selection in linear regression. It generates plots comparing posterior estimates with ground truth values, along with trace plots and autocorrelation function (ACF) plots for the parameters.

- **`TruncatedMVN.R`**: Extends the anti-correlation Gaussian technique to sampling from a truncated multivariate normal distribution with separable truncation (i.e., the constraint on the random vector is a high-dimensional box). The file includes a 2-dimensional example with plots of the marginal densities and joint density. The dimension `p` can be modified to accommodate higher dimensions as needed.
