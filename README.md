# Exploration-of-Mixture-g-Priors-BMA-Variable-Selection-in-Generalized-Linear-Model-with-Microarray-D

Bayesian model average is a well-known method that provides a coherent mechanism for accounting for the model uncertainty. In BMA, we place a prior probability $\pi$ on the model in the model space such that $\pi(\mathcal{M}m) > 0$ and $\sum^M{m = 1}\pi(\mathcal{M}_m) = 1$. The BMA framework can be written as:

$$\begin{array}{lr} \begin{aligned} X_{i} \mid \boldsymbol{\theta}, \mathcal{M}{m} & \stackrel{\text { ind }}{\sim} f{m}\left(x_{i} \mid \boldsymbol{\theta}{m}\right), & i=1, \ldots, n \ \boldsymbol{\theta}{m} & \stackrel{\text { ind }}{\sim} g_{m}\left(\boldsymbol{\theta}{m}\right), & m=1, \ldots, M \ \mathcal{M}{m} & \sim \pi\left(\mathcal{M}_{m}\right) & \end{aligned}&(1) \end{array}$$

