# Details of the parameter estimation

## The probability model

Maximum likelihood estimates are based on the probability model for the observed responses.
In the probability model the distribution of the responses is expressed as a function of one or more *parameters*.

For a continuous distribution the probability density is a function of the responses, given the parameters.
The *likelihood* function is the same expression as the probability density but regarding the observed values as fixed and the parameters as varying.

In general a mixed-effects model incorporates two random variables: $\mathcal{B}$, the $q$-dimensional vector of random effects, and $\mathcal{Y}$, the $n$-dimensional response vector.
The value, $\bf y$, of $\mathcal{Y}$ is observed; the value, $\bf b$, of $\mathcal{B}$ is not.

## Linear Mixed-Effects Models

In a linear mixed model the unconditional distribution of $\mathcal{B}$ and the conditional distribution, $(\mathcal{Y} | \mathcal{B}=\bf{b})$, are both multivariate Gaussian distributions,
```math
\begin{equation}
\begin{aligned}
  (\mathcal{Y} | \mathcal{B}=\bf{b}) &\sim\mathcal{N}(\bf{ X\beta + Z b},\sigma^2\bf{I})\\\\
  \mathcal{B}&\sim\mathcal{N}(\bf{0},\Sigma_\theta) .
\end{aligned}
\end{equation}
```
The *conditional mean* of $\mathcal Y$, given $\mathcal B=\bf b$, is the *linear predictor*, $\bf X\bf\beta+\bf Z\bf b$, which depends on the $p$-dimensional *fixed-effects parameter*, $\bf \beta$, and on $\bf b$.
The *model matrices*, $\bf X$ and $\bf Z$, of dimension $n\times p$ and $n\times q$, respectively, are determined from the formula for the model and the values of covariates. 
Although the matrix $\bf Z$ can be large (i.e. both $n$ and $q$ can be large), it is sparse (i.e. most of the elements in the matrix are zero).

The *relative covariance factor*, $\Lambda_\theta$, is a $q\times q$ lower-triangular matrix, depending on the *variance-component parameter*, $\bf\theta$, and generating the symmetric $q\times q$ variance-covariance matrix, $\Sigma_\theta$, as
```math
\begin{equation}
\Sigma_\theta=\sigma^2\Lambda_\theta\Lambda_\theta'
\end{equation}
```

The *spherical random effects*, $\mathcal{U}\sim\mathcal{N}(\bf{0},\sigma^2\bf{I}_q)$, determine $\mathcal B$ according to
```math
\begin{equation}
\mathcal{B}=\Lambda_\theta\mathcal{U}.
\end{equation}
```
The *penalized residual sum of squares* (PRSS),
```math
\begin{equation}
r^2(\theta,\beta,\bf{u})=\|\bf{y} - \bf{X}\beta -\bf{Z}\Lambda_\theta\bf{u}\|^2+\|\bf{u}\|^2,
\end{equation}
```
is the sum of the residual sum of squares, measuring fidelity of the model to the data, and a penalty on the size of $\bf u$, measuring the complexity of the model.
Minimizing $r^2$ with respect to $\bf u$,
```math
\begin{equation}
r^2_{\beta,\theta} =\min_{\bf{u}}\left(\|\bf{y} -\bf{X}{\beta} -\bf{Z}\Lambda_\theta\bf{u}\|^2+\|\bf{u}\|^2\right)
\end{equation}
```
is a direct (i.e. non-iterative) computation.
The particular method used to solve this generates a *blocked Choleksy factor*, $\bf{L}_\theta$, which is an lower triangular $q\times q$ matrix satisfying
```math
\begin{equation}
\bf{L}_\theta\bf{L}_\theta'=\Lambda_\theta'\bf{Z}'\bf{Z}\Lambda_\theta+\bf{I}_q .
\end{equation}
```
where ${\bf I}_q$ is the $q\times q$ *identity matrix*.

Negative twice the log-likelihood of the parameters, given the data, $\bf y$, is
```math
\begin{equation}
d({\bf\theta},{\bf\beta},\sigma|{\bf y})
=n\log(2\pi\sigma^2)+\log(|{\bf L}_\theta|^2)+\frac{r^2_{\beta,\theta}}{\sigma^2}.
\end{equation}
```
where $|{\bf L}_\theta|$ denotes the *determinant* of ${\bf L}_\theta$.
Because ${\bf L}_\theta$ is triangular, its determinant is the product of its diagonal elements.

Because the conditional mean, $\bf\mu_{\mathcal Y|\mathcal B=\bf b}=\bf
X\bf\beta+\bf Z\Lambda_\theta\bf u$, is a linear function of both $\bf\beta$ and $\bf u$, minimization of the PRSS with respect to both $\bf\beta$ and $\bf u$ to produce
```math
\begin{equation}
r^2_\theta =\min_{{\bf\beta},{\bf u}}\left(\|{\bf y} -{\bf X}{\bf\beta} -{\bf Z}\Lambda_\theta{\bf u}\|^2+\|{\bf u}\|^2\right)
\end{equation}
```
is also a direct calculation.
The values of $\bf u$ and $\bf\beta$ that provide this minimum are called, respectively, the *conditional mode*, $\tilde{\bf u}_\theta$, of the spherical random effects and the conditional estimate, $\widehat{\bf\beta}_\theta$, of the fixed effects.
At the conditional estimate of the fixed effects the objective is
```math
\begin{equation}
d({\bf\theta},\widehat{\beta}_\theta,\sigma|{\bf y})
=n\log(2\pi\sigma^2)+\log(|{\bf L}_\theta|^2)+\frac{r^2_\theta}{\sigma^2}.
\end{equation}
```
Minimizing this expression with respect to $\sigma^2$ produces the conditional estimate
```math
\begin{equation}
\widehat{\sigma^2}_\theta=\frac{r^2_\theta}{n}
\end{equation}
```
which provides the *profiled log-likelihood* on the deviance scale as
```math
\begin{equation}
\tilde{d}(\theta|{\bf y})=d(\theta,\widehat{\beta}_\theta,\widehat{\sigma}_\theta|{\bf y})
=\log(|{\bf L}_\theta|^2)+n\left[1+\log\left(\frac{2\pi r^2_\theta}{n}\right)\right],
\end{equation}
```
a function of $\bf\theta$ alone.

The MLE of $\bf\theta$, written $\widehat{\bf\theta}$, is the value that minimizes this profiled objective.
We determine this value by numerical optimization.
In the process of evaluating $\tilde{d}(\widehat{\theta}|{\bf y})$ we determine $\widehat{\beta}=\widehat{\beta}_{\widehat\theta}$, $\tilde{\bf u}_{\widehat{\theta}}$ and $r^2_{\widehat{\theta}}$, from which we can evaluate $\widehat{\sigma}=\sqrt{r^2_{\widehat{\theta}}/n}$.

The elements of the conditional mode of $\mathcal B$, evaluated at the parameter estimates,
```math
\begin{equation}
\tilde{\bf b}_{\widehat{\theta}}=\Lambda_{\widehat{\theta}}\tilde{\bf u}_{\widehat{\theta}}
\end{equation}
```
are sometimes called the *best linear unbiased predictors* or BLUPs of the random effects.
Although BLUPs an appealing acronym, I don’t find the term particularly instructive (what is a “linear unbiased predictor” and in what sense are these the “best”?) and prefer the term “conditional modes”, because these are the values of $\bf b$ that maximize the density of the conditional distribution $\mathcal{B} | \mathcal{Y} = {\bf y}$.
For a linear mixed model, where all the conditional and unconditional distributions are Gaussian, these values are also the *conditional means*.

## Internal structure of $\Lambda_\theta$ and $\bf Z$

In the types of `LinearMixedModel` available through the `MixedModels` package, groups of random effects and the corresponding columns of the model matrix, $\bf Z$, are associated with *random-effects terms* in the model formula.

For the simple example

````julia
julia> fm1 = fit!(LinearMixedModel(@formula(Y ~ 1 + (1|G)), dat[:Dyestuff]))
Linear mixed model fit by maximum likelihood
 Y ~ 1 + (1 | G)
   logLik   -2 logLik     AIC        BIC    
 -163.66353  327.32706  333.32706  337.53065

Variance components:
              Column    Variance  Std.Dev. 
 G        (Intercept)  1388.3334 37.260347
 Residual              2451.2500 49.510100
 Number of obs: 30; levels of grouping factors: 6

  Fixed-effects parameters:
──────────────────────────────────────────────────
             Estimate  Std.Error  z value  P(>|z|)
──────────────────────────────────────────────────
(Intercept)    1527.5    17.6946   86.326   <1e-99
──────────────────────────────────────────────────

````




the only random effects term in the formula is `(1|G)`, a simple, scalar random-effects term.
````julia
julia> t1 = first(fm1.reterms)
30×6 MixedModels.ReMat{Float64,1}:
 1.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 ⋮                        ⋮  
 0.0  0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0
 0.0  0.0  0.0  0.0  0.0  1.0
 0.0  0.0  0.0  0.0  0.0  1.0
 0.0  0.0  0.0  0.0  0.0  1.0
 0.0  0.0  0.0  0.0  0.0  1.0

````




```@docs
ReMat
```

This `RandomEffectsTerm` contributes a block of columns to the model matrix $\bf Z$ and a diagonal block to $\Lambda_\theta$.
In this case the diagonal block of $\Lambda_\theta$ (which is also the only block) is a multiple of the $6\times6$
identity matrix where the multiple is 
````julia
julia> t1.λ
1×1 LinearAlgebra.LowerTriangular{Float64,Array{Float64,2}}:
 0.75258072175308

````





Because there is only one random-effects term in the model, the matrix $\bf Z$ is the indicators matrix shown as the result of `Matrix(t1)`, but stored in a special sparse format.
Furthermore, there is only one block in $\Lambda_\theta$.


For a vector-valued random-effects term, as in
````julia
julia> fm2 = fit!(LinearMixedModel(@formula(Y ~ 1 + U + (1+U|G)), dat[:sleepstudy]))
Linear mixed model fit by maximum likelihood
 Y ~ 1 + U + (1 + U | G)
   logLik   -2 logLik     AIC        BIC    
 -875.96967 1751.93934 1763.93934 1783.09709

Variance components:
              Column    Variance  Std.Dev.   Corr.
 G        (Intercept)  565.51068 23.780468
          U             32.68212  5.716828  0.08
 Residual              654.94145 25.591824
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
───────────────────────────────────────────────────
             Estimate  Std.Error   z value  P(>|z|)
───────────────────────────────────────────────────
(Intercept)  251.405     6.63226  37.9064    <1e-99
U             10.4673    1.50224   6.96781   <1e-11
───────────────────────────────────────────────────

````




the model matrix $\bf Z$ for is of the form
````julia
julia> t21 = first(fm2.reterms);

julia> convert(Array{Int}, Matrix(t21))  # convert to integers for more compact printing
180×36 Array{Int64,2}:
 1  0  0  0  0  0  0  0  0  0  0  0  0  …  0  0  0  0  0  0  0  0  0  0  0  0
 1  1  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0
 1  2  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0
 1  3  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0
 1  4  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0
 1  5  0  0  0  0  0  0  0  0  0  0  0  …  0  0  0  0  0  0  0  0  0  0  0  0
 1  6  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0
 1  7  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0
 1  8  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0
 1  9  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  0  0
 ⋮              ⋮              ⋮        ⋱     ⋮              ⋮              ⋮
 0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  1  1
 0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  1  2
 0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  1  3
 0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  1  4
 0  0  0  0  0  0  0  0  0  0  0  0  0  …  0  0  0  0  0  0  0  0  0  0  1  5
 0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  1  6
 0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  1  7
 0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  1  8
 0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  0  0  0  0  0  0  0  1  9

````




and $\Lambda_\theta$ is a $36\times36$ block diagonal matrix with $18$ diagonal blocks, all of the form
````julia
julia> t21.λ
2×2 LinearAlgebra.LowerTriangular{Float64,Array{Float64,2}}:
 0.929221    ⋅      
 0.0181684  0.222645

````




The $\theta$ vector is
````julia
julia> MixedModels.getθ(t21)
3-element Array{Float64,1}:
 0.9292213194756659 
 0.01816837678765467
 0.2226448725352735 

````





__This model cannot be fit in v"2.0.0".  Watch for updates.__
Random-effects terms in the model formula that have the same grouping factor are amagamated into a single `VectorFactorReTerm` object.
````julia
julia> #fm3 = fit(LinearMixedModel(@formula(Y ~ 1 + U + (1|G) + (0+U|G)), dat[:sleepstudy]))
#t31 = first(fm3.retrms)

````




For this model the matrix $\bf Z$ is the same as that of model `fm2` but the diagonal blocks of $\Lambda_\theta$ are themselves diagonal.
````julia
julia> #getΛ(t31)
#getθ(t31)

````




Random-effects terms with distinct grouping factors generate distinct elements of the `trms` member of the `LinearMixedModel` object.
Multiple `AbstractFactorReTerm` (i.e. either a `ScalarFactorReTerm` or a `VectorFactorReTerm`) objects are sorted by decreasing numbers of random effects.
````julia
julia> fm4 = fit!(LinearMixedModel(@formula(Y ~ 1 + (1|H) + (1|G)), dat[:Penicillin]))
Linear mixed model fit by maximum likelihood
 Y ~ 1 + (1 | H) + (1 | G)
   logLik   -2 logLik     AIC        BIC    
 -166.09417  332.18835  340.18835  352.06760

Variance components:
              Column    Variance  Std.Dev. 
 G        (Intercept)  0.7149795 0.8455646
 H        (Intercept)  3.1351931 1.7706477
 Residual              0.3024264 0.5499331
 Number of obs: 144; levels of grouping factors: 24, 6

  Fixed-effects parameters:
──────────────────────────────────────────────────
             Estimate  Std.Error  z value  P(>|z|)
──────────────────────────────────────────────────
(Intercept)   22.9722   0.744596  30.8519   <1e-99
──────────────────────────────────────────────────

julia> t41 = first(fm4.reterms)
144×24 MixedModels.ReMat{Float64,1}:
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
 ⋮                        ⋮              ⋱                 ⋮                 
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  1.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  1.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  1.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  1.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  1.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  1.0

julia> t42 = last(fm4.reterms)
144×6 MixedModels.ReMat{Float64,1}:
 1.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0
 1.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  0.0  0.0
 ⋮                        ⋮  
 0.0  0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0
 1.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0  1.0

````




Note that the first `ReMat` in `fm4.terms` corresponds to grouping factor `G` even though the term `(1|G)` occurs in the formula after `(1|H)`.

### Progress of the optimization

An optional named argument, `verbose=true`, in the call to `fit!` of a `LinearMixedModel` causes printing of the objective and the $\theta$ parameter at each evaluation during the optimization.
````julia
julia> fit!(LinearMixedModel(@formula(Y ~ 1 + (1|G)), dat[:Dyestuff]), verbose=true);
f_1: 327.76702 [1.0]
f_2: 331.03619 [1.75]
f_3: 330.64583 [0.25]
f_4: 327.69511 [0.9761896354666064]
f_5: 327.56631 [0.9285689063998191]
f_6: 327.3826 [0.8333274482662446]
f_7: 327.35315 [0.8071883308459074]
f_8: 327.34663 [0.7996883308459073]
f_9: 327.341 [0.7921883308459072]
f_10: 327.33253 [0.7771883308459073]
f_11: 327.32733 [0.7471883308459073]
f_12: 327.32862 [0.7396883308459073]
f_13: 327.32706 [0.7527765100341628]
f_14: 327.32707 [0.7535265100341627]
f_15: 327.32706 [0.7525837540989767]
f_16: 327.32706 [0.7525087540989767]
f_17: 327.32706 [0.7525912540989768]
f_18: 327.32706 [0.75258072175308]

julia> fit!(LinearMixedModel(@formula(Y ~ 1 + U + (1+U|G)), dat[:sleepstudy]), verbose=true);
f_1: 1784.6423 [1.0, 0.0, 1.0]
f_2: 1790.12564 [1.75, 0.0, 1.0]
f_3: 1798.99962 [1.0, 1.0, 1.0]
f_4: 1803.8532 [1.0, 0.0, 1.75]
f_5: 1800.61398 [0.25, 0.0, 1.0]
f_6: 1798.60463 [1.0, -1.0, 1.0]
f_7: 1752.26074 [1.0, 0.0, 0.25]
f_8: 1797.58769 [1.183261296536725, -0.008661887958066734, 0.0]
f_9: 1754.95411 [1.075, 0.0, 0.32499999999999996]
f_10: 1753.69568 [0.8166315695342622, 0.011167254457032028, 0.28823768689684415]
f_11: 1754.817 [1.0, -0.07071067811865475, 0.19696699141100893]
f_12: 1753.10673 [0.9436827046395598, 0.06383542916411844, 0.26269630296448143]
f_13: 1752.93938 [0.9801419885634035, -0.026656844944251402, 0.2747427560953546]
f_14: 1752.25688 [0.9843428851818266, -0.013234749183508296, 0.24719098754480387]
f_15: 1752.05745 [0.9731403970870448, 0.0025378492296018644, 0.23791031400557064]
f_16: 1752.02239 [0.9545259030378296, 0.0038642106167098463, 0.23589201227662854]
f_17: 1752.02273 [0.935928530095792, 0.001331797324013128, 0.2344453466986732]
f_18: 1751.97169 [0.9549646039741277, 0.00790664248220189, 0.22904616789238916]
f_19: 1751.9526 [0.9533132639128648, 0.016627370637358883, 0.22576831302956576]
f_20: 1751.94852 [0.9469287318320833, 0.013076079966332387, 0.2228711267557771]
f_21: 1751.98718 [0.9334175303177042, 0.006137673805892002, 0.2189509416694232]
f_22: 1751.98321 [0.9515444328070808, 0.005788999134550992, 0.2206181988176504]
f_23: 1751.95197 [0.9528093408771897, 0.019033192065953525, 0.224177609032511]
f_24: 1751.94628 [0.946321530415129, 0.015373858744692083, 0.2250881772580525]
f_25: 1751.9467 [0.9471235457951435, 0.014889405832371933, 0.22489234774227196]
f_26: 1751.94757 [0.9464970169223144, 0.015464270391508158, 0.22581419823656373]
f_27: 1751.94531 [0.9460858413030404, 0.01579336991653446, 0.2244494625491751]
f_28: 1751.94418 [0.9453036925782654, 0.016690245709506217, 0.22336052906104453]
f_29: 1751.94353 [0.9440720728911274, 0.017210606451832967, 0.22271587721759273]
f_30: 1751.94244 [0.94127109818191, 0.016309946412355933, 0.2225226318161082]
f_31: 1751.94217 [0.9390004000469168, 0.01589901699551205, 0.22213197694735162]
f_32: 1751.94237 [0.9389790831398233, 0.016547964693256263, 0.22156175049004037]
f_33: 1751.94228 [0.9388628187235076, 0.015246587891193314, 0.22268346178654286]
f_34: 1751.9422 [0.9382687959892757, 0.015732966953942926, 0.22202359841547592]
f_35: 1751.94131 [0.9388391783099166, 0.016637330527596, 0.22261144012238268]
f_36: 1751.94093 [0.9383965533006099, 0.017396535205100986, 0.22281726224477705]
f_37: 1751.94057 [0.9370059167075939, 0.018044488667373536, 0.22253447647612865]
f_38: 1751.94018 [0.934109475252778, 0.01873542053799496, 0.22194958606264223]
f_39: 1751.94008 [0.9326416035982451, 0.01892417257690312, 0.2217257570735155]
f_40: 1751.94027 [0.9313571393417835, 0.019008175269225207, 0.22130945573424496]
f_41: 1751.9415 [0.9328207247273719, 0.020645436215409395, 0.221367296381061]
f_42: 1751.93949 [0.9318674733476721, 0.017957359912079057, 0.2225636491046244]
f_43: 1751.93939 [0.929167399868579, 0.017782430245528478, 0.22253384253403247]
f_44: 1751.9394 [0.9296587373859203, 0.01777209182293817, 0.22250843991222882]
f_45: 1751.93943 [0.9291934142211301, 0.018780637888612995, 0.2225704191376973]
f_46: 1751.93935 [0.9289856340300151, 0.018236600857577914, 0.22248440302631292]
f_47: 1751.93949 [0.9286970520601006, 0.018293686503017284, 0.2231753352863793]
f_48: 1751.93936 [0.9282426480418671, 0.018269516670966327, 0.22258371254427078]
f_49: 1751.93934 [0.9291127805145617, 0.018179125962924505, 0.22262389260385707]
f_50: 1751.93934 [0.9291906403862651, 0.01816575833837419, 0.22264320646294164]
f_51: 1751.93935 [0.9292543460617592, 0.018209275110801674, 0.22262081505106296]
f_52: 1751.93935 [0.9291892648085425, 0.01812980794781097, 0.22257323421427294]
f_53: 1751.93934 [0.9292535261081944, 0.018167626594423336, 0.22264990005039495]
f_54: 1751.93934 [0.9292145071088831, 0.01817173496298417, 0.22264674287059041]
f_55: 1751.93934 [0.9292083688509699, 0.018171468877843655, 0.22264619423897955]
f_56: 1751.93934 [0.9292093002228291, 0.018172962099089134, 0.2226520618262842]
f_57: 1751.93934 [0.9292213194756659, 0.01816837678765467, 0.2226448725352735]

````





A shorter summary of the optimization process is always available as an
```@docs
OptSummary
```
object, which is the `optsum` member of the `LinearMixedModel`.
````julia
julia> fm2.optsum
Initial parameter vector: [1.0, 0.0, 1.0]
Initial objective value:  1784.6422961924686

Optimizer (from NLopt):   LN_BOBYQA
Lower bounds:             [0.0, -Inf, 0.0]
ftol_rel:                 1.0e-12
ftol_abs:                 1.0e-8
xtol_rel:                 0.0
xtol_abs:                 [1.0e-10, 1.0e-10, 1.0e-10]
initial_step:             [0.75, 1.0, 0.75]
maxfeval:                 -1

Function evaluations:     57
Final parameter vector:   [0.9292213194756659, 0.01816837678765467, 0.2226448725352735]
Final objective value:    1751.939344464691
Return code:              FTOL_REACHED


````





### Modifying the optimization process

The `OptSummary` object contains both input and output fields for the optimizer.
To modify the optimization process the input fields can be changed after constructing the model but before fitting it.

Suppose, for example, that the user wishes to try a [Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) optimization method instead of the default [`BOBYQA`](https://en.wikipedia.org/wiki/BOBYQA) (Bounded Optimization BY Quadratic Approximation) method.
````julia
julia> fm2 = LinearMixedModel(@formula(Y ~ 1 + U + (1+U|G)), dat[:sleepstudy]);

julia> fm2.optsum.optimizer = :LN_NELDERMEAD;

julia> fit!(fm2)
Linear mixed model fit by maximum likelihood
 Y ~ 1 + U + (1 + U | G)
   logLik   -2 logLik     AIC        BIC    
 -875.96967 1751.93934 1763.93934 1783.09709

Variance components:
              Column    Variance   Std.Dev.   Corr.
 G        (Intercept)  565.528831 23.780850
          U             32.681047  5.716734  0.08
 Residual              654.941678 25.591828
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
──────────────────────────────────────────────────
             Estimate  Std.Error  z value  P(>|z|)
──────────────────────────────────────────────────
(Intercept)  251.405     6.63233  37.906    <1e-99
U             10.4673    1.50222   6.9679   <1e-11
──────────────────────────────────────────────────

julia> fm2.optsum
Initial parameter vector: [1.0, 0.0, 1.0]
Initial objective value:  1784.6422961924686

Optimizer (from NLopt):   LN_NELDERMEAD
Lower bounds:             [0.0, -Inf, 0.0]
ftol_rel:                 1.0e-12
ftol_abs:                 1.0e-8
xtol_rel:                 0.0
xtol_abs:                 [1.0e-10, 1.0e-10, 1.0e-10]
initial_step:             [0.75, 1.0, 0.75]
maxfeval:                 -1

Function evaluations:     140
Final parameter vector:   [0.9292360739538559, 0.018168794976407835, 0.22264111430139058]
Final objective value:    1751.9393444750278
Return code:              FTOL_REACHED


````





The parameter estimates are quite similar to those using `:LN_BOBYQA` but at the expense of 140 functions evaluations for `:LN_NELDERMEAD` versus 57 for `:LN_BOBYQA`.

See the documentation for the [`NLopt`](https://github.com/JuliaOpt/NLopt.jl) package for details about the various settings.

### Convergence to singular covariance matrices

To ensure identifiability of $\Sigma_\theta=\sigma^2\Lambda_\theta \Lambda_\theta$, the elements of $\theta$ corresponding to diagonal elements of $\Lambda_\theta$ are constrained to be non-negative.
For example, in a trivial case of a single, simple, scalar, random-effects term as in `fm1`, the one-dimensional $\theta$ vector is the ratio of the standard deviation of the random effects to the standard deviation of the response.
It happens that $-\theta$ produces the same log-likelihood but, by convention, we define the standard deviation to be the positive square root of the variance.
Requiring the diagonal elements of $\Lambda_\theta$ to be non-negative is a generalization of using this positive square root.

If the optimization converges on the boundary of the feasible region, that is if one or more of the diagonal elements of $\Lambda_\theta$ is zero at convergence, the covariance matrix $\Sigma_\theta$ will be *singular*.
This means that there will be linear combinations of random effects that are constant.
Usually convergence to a singular covariance matrix is a sign of an over-specified model.

## Generalized Linear Mixed-Effects Models

In a [*generalized linear model*](https://en.wikipedia.org/wiki/Generalized_linear_model) the responses are modelled as coming from a particular distribution, such as `Bernoulli` for binary responses or `Poisson` for responses that represent counts.
The scalar distributions of individual responses differ only in their means, which are determined by a *linear predictor* expression $\eta=\bf X\beta$, where, as before, $\bf X$ is a model matrix derived from the values of covariates and $\beta$ is a vector of coefficients.

The unconstrained components of $\eta$ are mapped to the, possiby constrained, components of the mean response, $\mu$, via a scalar function, $g^{-1}$, applied to each component of $\eta$.
For historical reasons, the inverse of this function, taking components of $\mu$ to the corresponding component of $\eta$ is called the *link function* and more frequently used map from $\eta$ to $\mu$ is the *inverse link*.

A *generalized linear mixed-effects model* (GLMM) is defined, for the purposes of this package, by
```math
\begin{equation}
\begin{aligned}
  (\mathcal{Y} | \mathcal{B}=\bf{b}) &\sim\mathcal{D}(\bf{g^{-1}(X\beta + Z b)},\phi)\\\\
  \mathcal{B}&\sim\mathcal{N}(\bf{0},\Sigma_\theta) .
\end{aligned}
\end{equation}
```
where $\mathcal{D}$ indicates the distribution family parameterized by the mean and, when needed, a common scale parameter, $\phi$.
(There is no scale parameter for `Bernoulli` or for `Poisson`.
Specifying the mean completely determines the distribution.)
```@docs
Bernoulli
Poisson
```

A `GeneralizedLinearMixedModel` object is generated from a formula, data frame and distribution family.

````julia
julia> mdl = GeneralizedLinearMixedModel(@formula(r2 ~ 1 + a + g + b + s + (1|id) + (1|item)),
           dat[:VerbAgg], Bernoulli());

julia> typeof(mdl)
MixedModels.GeneralizedLinearMixedModel{Float64}

````





A separate call to `fit!` is required to fit the model.
This involves optimizing an objective function, the Laplace approximation to the deviance, with respect to the parameters, which are $\beta$, the fixed-effects coefficients, and $\theta$, the covariance parameters.
The starting estimate for $\beta$ is determined by fitting a GLM to the fixed-effects part of the formula

````julia
julia> mdl.β
6-element Array{Float64,1}:
  0.20605302210323584 
  0.039940376051149744
  0.2313166767498431  
 -0.7941857249205405  
 -1.5391882085456963  
 -0.7766556048305944  

````





and the starting estimate for $\theta$, which is a vector of the two standard deviations of the random effects, is chosen to be

````julia
julia> mdl.θ
2-element Array{Float64,1}:
 1.0
 1.0

````





The Laplace approximation to the deviance requires determining the conditional modes of the random effects.
These are the values that maximize the conditional density of the random effects, given the model parameters and the data.
This is done using Penalized Iteratively Reweighted Least Squares (PIRLS).
In most cases PIRLS is fast and stable.
It is simply a penalized version of the IRLS algorithm used in fitting GLMs.

The distinction between the "fast" and "slow" algorithms in the `MixedModels` package (`nAGQ=0` or `nAGQ=1` in `lme4`) is whether the fixed-effects parameters, $\beta$, are optimized in PIRLS or in the nonlinear optimizer.
In a call to the `pirls!` function the first argument is a `GeneralizedLinearMixedModel`, which is modified during the function call.
(By convention, the names of such *mutating functions* end in `!` as a warning to the user that they can modify an argument, usually the first argument.)
The second and third arguments are optional logical values indicating if $\beta$ is to be varied and if verbose output is to be printed.

````julia
julia> pirls!(mdl, true, true)
varyβ = true, obj₀ = 10210.853438905404, β = [0.20605302210323584, 0.039940376051149744, 0.2313166767498431, -0.7941857249205405, -1.5391882085456963, -0.7766556048305944]
   1: 8301.483049027265
   2: 8205.604285133919
   3: 8201.89659746689
   4: 8201.848598910705
   5: 8201.848559060705
   6: 8201.848559060621
Generalized Linear Mixed Model fit by maximum likelihood (nAGQ = 1)
  r2 ~ 1 + a + g + b + s + (1 | id) + (1 | item)
  Distribution: Distributions.Bernoulli{Float64}
  Link: GLM.LogitLink()

  Deviance: 8201.8486

Variance components:
          Column   Variance Std.Dev. 
 id   (Intercept)         1        1
 item (Intercept)         1        1
 Number of obs: 7584; levels of grouping factors: 316, 24

Fixed-effects parameters:
─────────────────────────────────────────────────────
               Estimate  Std.Error   z value  P(>|z|)
─────────────────────────────────────────────────────
(Intercept)   0.218535    0.464651   0.47032   0.6381
a             0.0514385   0.012319   4.17556   <1e-4 
g: M          0.290225    0.140555   2.06485   0.0389
b: scold     -0.979124    0.476395  -2.05528   0.0399
b: shout     -1.95402     0.477182  -4.09491   <1e-4 
s: self      -0.979493    0.389283  -2.51615   0.0119
─────────────────────────────────────────────────────

````



````julia
julia> deviance(mdl)
8201.848559060621

````



````julia
julia> mdl.β
6-element Array{Float64,1}:
  0.21853493716521263
  0.05143854258081151
  0.29022454166301326
 -0.9791237061900561 
 -1.9540167628140472 
 -0.97949257180371   

````



````julia
julia> mdl.θ # current values of the standard deviations of the random effects
2-element Array{Float64,1}:
 1.0
 1.0

````





If the optimization with respect to $\beta$ is performed within PIRLS then the nonlinear optimization of the Laplace approximation to the deviance requires optimization with respect to $\theta$ only.
This is the "fast" algorithm.
Given a value of $\theta$, PIRLS is used to determine the conditional estimate of $\beta$ and the conditional mode of the random effects, **b**.

````julia
julia> mdl.b # conditional modes of b
2-element Array{Array{Float64,2},1}:
 [-0.6007716038488825 -1.9322680866219635 … -0.1445537397533577 -0.5752238433557008]
 [-0.18636418747905698 0.18055184071670155 … 0.2820923275093635 -0.2219744597240774]

````



````julia
julia> fit!(mdl, fast=true, verbose=true);
varyβ = true, obj₀ = 10251.003116042964, β = [0.21853493716521263, 0.05143854258081151, 0.29022454166301326, -0.9791237061900561, -1.9540167628140472, -0.97949257180371]
   1: 8292.390783437771
   2: 8204.692089323944
   3: 8201.87681054392
   4: 8201.848569551961
   5: 8201.848559060627
   6: 8201.848559060623
f_1: 8201.848559060623 [1.0, 1.0]
varyβ = true, obj₀ = 10565.66000371903, β = [0.21853493716513464, 0.05143854258081691, 0.2902245416630119, -0.9791237061900924, -1.9540167628141278, -0.9794925718036764]
   1: 8356.488185490087
   2: 8200.932708793818
   3: 8190.472561852921
   4: 8190.11991835326
   5: 8190.117817016412
   6: 8190.117816802472
f_2: 8190.117816802472 [1.75, 1.0]
varyβ = true, obj₀ = 10317.246304689326, β = [0.1934205873511403, 0.057387763683876505, 0.3179960395260651, -1.0485464263732174, -2.096031981211341, -1.051128801683702]
   1: 8322.599130522533
   2: 8227.326059488892
   3: 8224.479329416501
   4: 8224.450989053388
   5: 8224.450978980776
   6: 8224.45097898077
f_3: 8224.45097898077 [1.0, 1.75]
varyβ = true, obj₀ = 9776.19623800526, β = [0.21876438895177894, 0.05148884100230788, 0.29046070478860203, -0.9797643406721465, -1.957296279410692, -0.9812079030632912]
   1: 9035.725240344302
   2: 9026.060558389585
   3: 9026.003998095694
   4: 9026.00390640421
   5: 9026.003906403308
f_4: 9026.003906403308 [0.25, 1.0]
varyβ = true, obj₀ = 10149.608772216949, β = [0.19893695106103745, 0.04180670725991035, 0.2418986481305357, -0.8153959340434698, -1.6334718386669982, -0.8201772064235269]
   1: 8286.346044727727
   2: 8208.252412394584
   3: 8205.81282801732
   4: 8205.793778239711
   5: 8205.793775487931
f_5: 8205.793775487931 [1.0, 0.25]
varyβ = true, obj₀ = 10406.590456413764, β = [0.21728094120285796, 0.050731092268406296, 0.28679620590417065, -0.9699881548970194, -1.9127335969466073, -0.9591823894790004]
   1: 8290.678720764576
   2: 8163.613865769221
   3: 8157.167433289948
   4: 8157.041222813374
   5: 8157.041032394054
   6: 8157.041032392803
f_6: 8157.041032392803 [1.385829382367978, 0.7364567143287161]
varyβ = true, obj₀ = 10334.285609061462, β = [0.20721523412589943, 0.05490998576682537, 0.30652969414239317, -1.0227806700669033, -2.0401155431546196, -1.022735054186784]
   1: 8461.553656168606
   2: 8371.201859251125
   3: 8367.771392346047
   4: 8367.72426348817
   5: 8367.72422310067
   6: 8367.724223100573
f_7: 8367.724223100573 [1.337152446554522, 0.0]
varyβ = true, obj₀ = 10441.99288890887, β = [0.21310775427581224, 0.05262369623724979, 0.29594851619262597, -0.9987562912090301, -1.933937463915685, -0.9735353380195877]
   1: 8308.813446103555
   2: 8177.318105403585
   3: 8170.433438731543
   4: 8170.289084465008
   5: 8170.288828193681
   6: 8170.288828191517
f_8: 8170.288828191517 [1.4136494266701973, 1.1104233507216201]
varyβ = true, obj₀ = 10394.123754834098, β = [0.2065845371410874, 0.05521337318328264, 0.3079133180419655, -1.0262965361212943, -2.0504566216499214, -1.0280362485643955]
   1: 8283.5403452246
   2: 8163.918696581621
   3: 8158.90579582785
   4: 8158.829375849112
   5: 8158.829317592085
   6: 8158.8293175919935
f_9: 8158.8293175919935 [1.2722464441311803, 0.7628110428958444]
varyβ = true, obj₀ = 10449.122512319096, β = [0.21122282489395758, 0.0540054840572109, 0.30230959849208483, -1.0122995687448604, -2.0190655839378895, -1.0121221456081149]
   1: 8299.982309744812
   2: 8168.539378419905
   3: 8162.06111254569
   4: 8161.933600781366
   5: 8161.933408558582
   6: 8161.933408557296
f_10: 8161.933408557296 [1.4093623667537525, 0.8680844411698196]
varyβ = true, obj₀ = 10414.107933396324, β = [0.20653729580020003, 0.05513368643868309, 0.3075576744964962, -1.0253459434061463, -2.0468884488247925, -1.0261895534474195]
   1: 8286.138316663313
   2: 8161.902045658136
   3: 8156.392925946188
   4: 8156.301068837265
   5: 8156.300980018685
   6: 8156.300980018441
f_11: 8156.300980018441 [1.3269393938015637, 0.7210153433276605]
varyβ = true, obj₀ = 10407.14590574727, β = [0.20929363010827845, 0.054440786160250616, 0.30434503952723996, -1.017410038992182, -2.0289135444050945, -1.017076254901686]
   1: 8284.877807447845
   2: 8161.708700091349
   3: 8156.208523569669
   4: 8156.116764659553
   5: 8156.116675394206
   6: 8156.116675393959
f_12: 8156.116675393959 [1.3236493916636756, 0.7142754704851885]
varyβ = true, obj₀ = 10404.356312407219, β = [0.20939665067030727, 0.054411231835824886, 0.3042078506972754, -1.0170635346399617, -2.0281059261519183, -1.016666961042768]
   1: 8284.15587979821
   2: 8161.5390145667825
   3: 8156.0921687832415
   4: 8156.002155450718
   5: 8156.002070048327
   6: 8156.002070048107
f_13: 8156.002070048107 [1.318465324298597, 0.7088555585343973]
varyβ = true, obj₀ = 10404.53296455974, β = [0.20956728950785716, 0.05436677355363324, 0.3040009712301877, -1.0165433887574258, -2.0269623799216827, -1.0160886153662267]
   1: 8284.027704597496
   2: 8161.308933576595
   3: 8155.84427374012
   4: 8155.753676771535
   5: 8155.753590029702
   6: 8155.753590029469
f_14: 8155.753590029469 [1.320716938857388, 0.7017015224792295]
varyβ = true, obj₀ = 10406.08059012708, β = [0.20947609033479664, 0.054381708709041736, 0.30407138125848593, -1.0167138677796181, -2.027198583515564, -1.0162057468818277]
   1: 8283.987555395504
   2: 8160.874305438914
   3: 8155.367283709278
   4: 8155.275310500365
   5: 8155.275220625106
   6: 8155.275220624852
f_15: 8155.275220624852 [1.3263558167622518, 0.6878017722665803]
varyβ = true, obj₀ = 10409.588006839367, β = [0.2092533197475376, 0.05442052407951385, 0.30425395783378195, -1.0171578113767867, -2.027878717803296, -1.0165449133816418]
   1: 8284.070875041212
   2: 8160.104155117698
   3: 8154.505123130821
   4: 8154.410093357733
   5: 8154.409996159085
   6: 8154.409996158775
f_16: 8154.409996158775 [1.3385857165571091, 0.660407803021366]
varyβ = true, obj₀ = 10422.176224407653, β = [0.2087679854479831, 0.054504010629259254, 0.30464664470271935, -1.0181072531047377, -2.029334902760042, -1.0172718754794445]
   1: 8286.148930569901
   2: 8159.408927883642
   3: 8153.501004776301
   4: 8153.395082917684
   5: 8153.394956151873
   6: 8153.394956151293
f_17: 8153.394956151293 [1.3758192930075792, 0.6133582463315629]
varyβ = true, obj₀ = 10427.671854015107, β = [0.20732305195162948, 0.05476609730564688, 0.3058757219214882, -1.0210628480448425, -2.034408391105722, -1.019819941218809]
   1: 8286.691066751828
   2: 8158.869556551096
   3: 8152.851268971666
   4: 8152.74108389247
   5: 8152.740942646492
   6: 8152.740942645721
f_18: 8152.740942645721 [1.3951508993559725, 0.563096398842232]
varyβ = true, obj₀ = 10412.134716826527, β = [0.20647531795319488, 0.0548743832097482, 0.30638967606898054, -1.0222230837232897, -2.0355303719223774, -1.0203755723551062]
   1: 8282.017684910908
   2: 8157.5220168154065
   3: 8151.862297991303
   4: 8151.764833400248
   5: 8151.764726543592
   6: 8151.764726543154
f_19: 8151.764726543154 [1.3676323395139596, 0.5091235621060084]
varyβ = true, obj₀ = 10364.111639276893, β = [0.20730625718507492, 0.05461186270246281, 0.3051748195857039, -1.0191804210341766, -2.027680778962547, -1.0164075726118529]
   1: 8271.970907130364
   2: 8157.5555577531
   3: 8152.875645706553
   4: 8152.808979106382
   5: 8152.808935981127
   6: 8152.808935981069
f_20: 8152.808935981069 [1.2677640577776739, 0.47512314632935615]
varyβ = true, obj₀ = 10420.049185016687, β = [0.2106775805664361, 0.053768043388408154, 0.30122955512263344, -1.0092843628515935, -2.006396181498495, -1.0056787436393213]
   1: 8285.752253358396
   2: 8159.0739780960885
   3: 8152.980811898105
   4: 8152.866577307153
   5: 8152.866415002163
   6: 8152.86641500103
f_21: 8152.86641500103 [1.4147967067023959, 0.4710990869768519]
varyβ = true, obj₀ = 10395.48938932329, β = [0.2054165953692425, 0.054920824377008, 0.30662779729746137, -1.0225711061446965, -2.033140785136617, -1.0191682418196022]
   1: 8278.256525895518
   2: 8157.154672639189
   3: 8151.854900402225
   4: 8151.769667985884
   5: 8151.769591419974
   6: 8151.769591419772
f_22: 8151.769591419772 [1.325886065404328, 0.5275226916888389]
varyβ = true, obj₀ = 10406.84702257384, β = [0.20886367217970872, 0.05430317854236076, 0.3037280282631728, -1.0156611090051952, -2.0211403904573264, -1.0131090903115638]
   1: 8281.195202177089
   2: 8157.492470479971
   3: 8151.835428349153
   4: 8151.737863139613
   5: 8151.737755041395
   6: 8151.737755040944
f_23: 8151.737755040944 [1.3668062704599326, 0.49860621312361686]
varyβ = true, obj₀ = 10396.62318319631, β = [0.20729825419661121, 0.054593630044578394, 0.3050918349444385, -1.018955366189819, -2.0268543124582865, -1.0159907029292399]
   1: 8278.460002109783
   2: 8157.053323937607
   3: 8151.673309975078
   4: 8151.585246119005
   5: 8151.585161942601
   6: 8151.58516194234
f_24: 8151.58516194234 [1.3397371137748282, 0.49349219821758816]
varyβ = true, obj₀ = 10393.60931624831, β = [0.2082549659763669, 0.054376602982506315, 0.30407785019276673, -1.016462598791784, -2.021611817000961, -1.0133460788416402]
   1: 8277.918179767254
   2: 8157.043669989343
   3: 8151.68944059344
   4: 8151.602147764987
   5: 8151.602064961424
   6: 8151.602064961172
f_25: 8151.602064961172 [1.3375752795348084, 0.48631052204007374]
varyβ = true, obj₀ = 10397.764486556898, β = [0.20830542908618235, 0.05435096553106824, 0.3039593106172651, -1.016155180970176, -2.0207219752066155, -1.0128978657446543]
   1: 8278.935236559155
   2: 8157.145230819007
   3: 8151.691194952189
   4: 8151.600591190683
   5: 8151.600500691267
   6: 8151.60050069096
f_26: 8151.60050069096 [1.3469240905642434, 0.4913480516208264]
varyβ = true, obj₀ = 10395.774178197511, β = [0.20798988251473483, 0.054430597969235724, 0.3043309127060538, -1.0170803419710324, -2.022785425820781, -1.013938218472325]
   1: 8278.406326249222
   2: 8157.064425497865
   3: 8151.6719513436765
   4: 8151.583465843254
   5: 8151.583380717896
   6: 8151.583380717632
f_27: 8151.583380717632 [1.3395817610360452, 0.49733727686890145]
varyβ = true, obj₀ = 10395.443943189504, β = [0.20827453095426368, 0.05437984247156536, 0.3040922464549829, -1.016506043864169, -2.0218392967293224, -1.013460616264679]
   1: 8278.362655736923
   2: 8157.064544983574
   3: 8151.672160260856
   4: 8151.58367501763
   5: 8151.583589892588
   6: 8151.583589892321
f_28: 8151.583589892321 [1.3392718767160512, 0.49802026421008855]
varyβ = true, obj₀ = 10395.535382999346, β = [0.20828805461930472, 0.05437818021174276, 0.3040843349543604, -1.0164879468260042, -2.021827115712798, -1.0134544452465297]
   1: 8278.380710340887
   2: 8157.066651045122
   3: 8151.671999799076
   4: 8151.583432800491
   5: 8151.583347442415
   6: 8151.583347442148
f_29: 8151.583347442148 [1.3397210818303735, 0.49695538534953754]
varyβ = true, obj₀ = 10395.855932030907, β = [0.20826817633963607, 0.05438050242433693, 0.30409540944507163, -1.0165130522337495, -2.0218397675906115, -1.0134608700296197]
   1: 8278.456305937294
   2: 8157.074106468333
   3: 8151.672400464071
   4: 8151.583600201322
   5: 8151.5835143054655
   6: 8151.583514305196
f_30: 8151.583514305196 [1.3404085881839323, 0.49725511024642816]
varyβ = true, obj₀ = 10395.450813218617, β = [0.20824473328012896, 0.05438627252763261, 0.30412234667904087, -1.0165800287330096, -2.021986355042864, -1.0135347946044966]
   1: 8278.353361600997
   2: 8157.063550036673
   3: 8151.671954135129
   4: 8151.583485621317
   5: 8151.58340047098
   6: 8151.583400470715
f_31: 8151.583400470715 [1.3395728772296231, 0.49622017422590264]
varyβ = true, obj₀ = 10395.726653540596, β = [0.2082707914372589, 0.05437848364454214, 0.3040861090956531, -1.0164886208435095, -2.021763716970758, -1.0134225446121197]
   1: 8278.426177247662
   2: 8157.071791450571
   3: 8151.6723002428325
   4: 8151.583568754747
   5: 8151.583482990038
   6: 8151.583482989765
f_32: 8151.583482989765 [1.3403128971075373, 0.49649467429289207]
varyβ = true, obj₀ = 10395.480743445069, β = [0.2082453864901393, 0.05438464054170947, 0.30411486152871514, -1.0165600218525852, -2.021918394381647, -1.0135005493945077]
   1: 8278.364600272997
   2: 8157.064659600624
   3: 8151.671928507909
   4: 8151.583425349218
   5: 8151.583340139288
   6: 8151.583340139021
f_33: 8151.583340139021 [1.3395583304721455, 0.4968332751125267]
varyβ = true, obj₀ = 10395.442267530794, β = [0.20827353831163917, 0.0543790769603831, 0.30408876413394115, -1.0164964269722956, -2.021801661225551, -1.0134416550289664]
   1: 8278.359702896334
   2: 8157.064704875616
   3: 8151.671931416344
   4: 8151.5834264301775
   5: 8151.58334121545
   6: 8151.583341215184
f_34: 8151.583341215184 [1.3395286260146129, 0.49690214197830906]
varyβ = true, obj₀ = 10395.479774749134, β = [0.20827484776294208, 0.05437892187903991, 0.30408802496957216, -1.0164947473254018, -2.0218007362745256, -1.0134411854625793]
   1: 8278.368369189171
   2: 8157.065584741343
   3: 8151.671959910017
   4: 8151.583426365367
   5: 8151.583341082105
   6: 8151.5833410818395
f_35: 8151.5833410818395 [1.3396255162939181, 0.496866607759477]
varyβ = true, obj₀ = 10395.478844626039, β = [0.20827126293062556, 0.05437964575020565, 0.3040914187934968, -1.0165030367587673, -2.0218162635246677, -1.0134490187641243]
   1: 8278.367429147567
   2: 8157.06544565772
   3: 8151.671955750738
   4: 8151.583426429475
   5: 8151.583341154113
   6: 8151.583341153844
f_36: 8151.583341153844 [1.339626756927264, 0.4968025697344976]
varyβ = true, obj₀ = 10395.453884480128, β = [0.20827098617862877, 0.054379581686990676, 0.3040911316853379, -1.0165021974731276, -2.021812254820609, -1.0134469996600357]
   1: 8278.361837064556
   2: 8157.064895678978
   3: 8151.671936855652
   4: 8151.583425364516
   5: 8151.583340132135
   6: 8151.583340131867
f_37: 8151.583340131867 [1.339563900001637, 0.49683278385267693]
varyβ = true, obj₀ = 10395.451247329112, β = [0.20827333787817207, 0.054379120359622626, 0.30408896728574886, -1.016496926541366, -2.0218026564175355, -1.0134421570121597]
   1: 8278.361538693684
   2: 8157.064912191759
   3: 8151.671937454685
   4: 8151.583425366179
   5: 8151.583340132134
   6: 8151.583340131869

````





The optimization process is summarized by

````julia
julia> mdl.LMM.optsum
Initial parameter vector: [1.0, 1.0]
Initial objective value:  8201.848559060623

Optimizer (from NLopt):   LN_BOBYQA
Lower bounds:             [0.0, 0.0]
ftol_rel:                 1.0e-12
ftol_abs:                 1.0e-8
xtol_rel:                 0.0
xtol_abs:                 [1.0e-10, 1.0e-10]
initial_step:             [0.75, 0.75]
maxfeval:                 -1

Function evaluations:     37
Final parameter vector:   [1.339563900001637, 0.49683278385267693]
Final objective value:    8151.583340131867
Return code:              FTOL_REACHED


````





As one would hope, given the name of the option, this fit is comparatively fast.
````julia
julia> @time(fit!(GeneralizedLinearMixedModel(@formula(r2 ~ 1 + a + g + b + s + (1 | id) + (1 | item)), 
        dat[:VerbAgg], Bernoulli()), fast=true))
  0.849437 seconds (2.87 M allocations: 64.163 MiB, 1.33% gc time)
Generalized Linear Mixed Model fit by maximum likelihood (nAGQ = 1)
  r2 ~ 1 + a + g + b + s + (1 | id) + (1 | item)
  Distribution: Distributions.Bernoulli{Float64}
  Link: GLM.LogitLink()

  Deviance: 8151.5833

Variance components:
          Column    Variance   Std.Dev. 
 id   (Intercept)  1.79443144 1.3395639
 item (Intercept)  0.24684282 0.4968328
 Number of obs: 7584; levels of grouping factors: 316, 24

Fixed-effects parameters:
──────────────────────────────────────────────────────
               Estimate  Std.Error    z value  P(>|z|)
──────────────────────────────────────────────────────
(Intercept)   0.208273   0.387547    0.537414   0.5910
a             0.0543791  0.0160145   3.39561    0.0007
g: M          0.304089   0.182791    1.66359    0.0962
b: scold     -1.0165     0.246175   -4.12917    <1e-4 
b: shout     -2.0218     0.247803   -8.15891    <1e-15
s: self      -1.01344    0.201588   -5.02728    <1e-6 
──────────────────────────────────────────────────────

````





The alternative algorithm is to use PIRLS to find the conditional mode of the random effects, given $\beta$ and $\theta$ and then use the general nonlinear optimizer to fit with respect to both $\beta$ and $\theta$.
Because it is slower to incorporate the $\beta$ parameters in the general nonlinear optimization, the fast fit is performed first and used to determine starting estimates for the more general optimization.

````julia
julia> @time mdl1 = fit!(GeneralizedLinearMixedModel(@formula(r2 ~ 1+a+g+b+s+(1|id)+(1|item)),
        dat[:VerbAgg], Bernoulli()))
  2.728348 seconds (12.98 M allocations: 164.780 MiB, 1.34% gc time)
Generalized Linear Mixed Model fit by maximum likelihood (nAGQ = 1)
  r2 ~ 1 + a + g + b + s + (1 | id) + (1 | item)
  Distribution: Distributions.Bernoulli{Float64}
  Link: GLM.LogitLink()

  Deviance: 8151.3997

Variance components:
          Column    Variance   Std.Dev. 
 id   (Intercept)  1.79477044 1.3396904
 item (Intercept)  0.24529882 0.4952765
 Number of obs: 7584; levels of grouping factors: 316, 24

Fixed-effects parameters:
──────────────────────────────────────────────────────
               Estimate  Std.Error    z value  P(>|z|)
──────────────────────────────────────────────────────
(Intercept)   0.199084   0.38773     0.513461   0.6076
a             0.0574292  0.0160358   3.58131    0.0003
g: M          0.320644   0.183024    1.75192    0.0798
b: scold     -1.05895    0.245739   -4.30925    <1e-4 
b: shout     -2.10546    0.247388   -8.51074    <1e-16
s: self      -1.05535    0.20124    -5.24423    <1e-6 
──────────────────────────────────────────────────────

````





This fit provided slightly better results (Laplace approximation to the deviance of 8151.400 versus 8151.583) but took 6 times as long.
That is not terribly important when the times involved are a few seconds but can be important when the fit requires many hours or days of computing time.

The comparison of the slow and fast fit is available in the optimization summary after the slow fit.

````julia
julia> mdl1.LMM.optsum
Initial parameter vector: [0.20827333787911337, 0.054379120359429704, 0.3040889672848241, -1.016496926539076, -2.021802656413764, -1.013442157010256, 1.3395638999730468, 0.49683278387583557]
Initial objective value:  8151.583340131868

Optimizer (from NLopt):   LN_BOBYQA
Lower bounds:             [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 0.0, 0.0]
ftol_rel:                 1.0e-12
ftol_abs:                 1.0e-8
xtol_rel:                 0.0
xtol_abs:                 [1.0e-10, 1.0e-10]
initial_step:             [0.1291823170508878, 0.005338175496666596, 0.06093022804397936, 0.08205826499805279, 0.08260097852350845, 0.06719613447778985, 0.05, 0.05]
maxfeval:                 -1

Function evaluations:     178
Final parameter vector:   [0.19908440905470953, 0.05742923240217957, 0.32064368389297343, -1.0589508896133553, -2.105456449464348, -1.0553501482027767, 1.339690426159839, 0.49527651021787783]
Final objective value:    8151.399720690761
Return code:              FTOL_REACHED


````

