var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "MixedModels.jl Documentation",
    "title": "MixedModels.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#MixedModels.jl-Documentation-1",
    "page": "MixedModels.jl Documentation",
    "title": "MixedModels.jl Documentation",
    "category": "section",
    "text": "CurrentModule = MixedModelsMixedModels.jl is a Julia package providing capabilities for fitting and examining linear and generalized linear mixed-effect models. It is similar in scope to the lme4 package for R.Pages = [\n        \"constructors.md\",\n        \"optimization.md\",\n        \"GaussHermite.md\",\n        \"bootstrap.md\",\n]\nDepth = 2"
},

{
    "location": "constructors/#",
    "page": "Model constructors",
    "title": "Model constructors",
    "category": "page",
    "text": ""
},

{
    "location": "constructors/#MixedModels.LinearMixedModel",
    "page": "Model constructors",
    "title": "MixedModels.LinearMixedModel",
    "category": "type",
    "text": "LinearMixedModel\n\nLinear mixed-effects model representation\n\nFields\n\nformula: the formula for the model\nallterms: a vector of random-effects terms, the fixed-effects terms and the response\nsqrtwts: vector of square roots of the case weights.  Can be empty.\nA: an nt × nt symmetric BlockMatrix of matrices representing hcat(Z,X,y)\'hcat(Z,X,y)\nL: a nt × nt BlockMatrix - the lower Cholesky factor of Λ\'AΛ+I\noptsum: an OptSummary object\n\nProperties\n\nθ or theta: the covariance parameter vector used to form λ\nβ or beta: the fixed-effects coefficient vector\nλ or lambda: a vector of lower triangular matrices repeated on the diagonal blocks of Λ\nσ or sigma: current value of the standard deviation of the per-observation noise\nb: random effects on the original scale, as a vector of matrices\nreterms: a Vector{ReMat{T}} of random-effects terms.\nfeterms: a Vector{FeMat{T}} of the fixed-effects model matrix and the response\nu: random effects on the orthogonal scale, as a vector of matrices\nlowerbd: lower bounds on the elements of θ\nX: the fixed-effects model matrix\ny: the response vector\n\n\n\n\n\n"
},

{
    "location": "constructors/#Model-constructors-1",
    "page": "Model constructors",
    "title": "Model constructors",
    "category": "section",
    "text": "The LinearMixedModel type represents a linear mixed-effects model. Typically it is constructed from a Formula and an appropriate data type, usually a DataFrame.LinearMixedModel"
},

{
    "location": "constructors/#Examples-of-linear-mixed-effects-model-fits-1",
    "page": "Model constructors",
    "title": "Examples of linear mixed-effects model fits",
    "category": "section",
    "text": "For illustration, several data sets from the lme4 package for R are made available in .rda format in this package. These include the Dyestuff and Dyestuff2 data sets.julia> using DataFrames, MixedModels, RData, StatsBase\n\njulia> datf = joinpath(dirname(pathof(MixedModels)),\"..\",\"test\",\"dat.rda\")\n\"/home/travis/build/JuliaStats/MixedModels.jl/src/../test/dat.rda\"\n\njulia> const dat = Dict(Symbol(k)=>v for (k,v) in load(datf))\nError: SystemError: opening file \"/home/travis/build/JuliaStats/MixedModels.jl/src/../test/dat.rda\": No such file or directory\n\njulia> describe(dat[:Dyestuff])\nError: UndefVarError: dat not defined\nThe columns in these data sets have been renamed for convenience. The response is always named Y. Potential grouping factors for random-effects terms are named G, H, etc. Numeric covariates are named starting with U. Categorical covariates not suitable as grouping factors are named starting with A."
},

{
    "location": "constructors/#Models-with-simple,-scalar-random-effects-1",
    "page": "Model constructors",
    "title": "Models with simple, scalar random effects",
    "category": "section",
    "text": "The formula language in Julia is similar to that in R except that the formula must be enclosed in a call to the @formula macro. A basic model with simple, scalar random effects for the levels of G (the batch of an intermediate product, in this case) is declared and fit asjulia> fm1 = fit!(LinearMixedModel(@formula(Y ~ 1 + (1|G)),\n    dat[:Dyestuff]))\nError: UndefVarError: dat not defined\nAn alternative expression isjulia> fm1 = fit(MixedModel, @formula(Y ~ 1 + (1|G)), dat[:Dyestuff])\nError: UndefVarError: dat not defined\n(If you are new to Julia you may find that this first fit takes an unexpectedly long time, due to Just-In-Time (JIT) compilation of the code. The second and subsequent calls to such functions are much faster.)julia> @time fit(MixedModel, @formula(Y ~ 1 + (1|G)), dat[:Dyestuff2]);\nError: UndefVarError: dat not defined\nBy default, the model fit is by maximum likelihood.  To use the REML criterion instead, add the optional named argument REML = true to the call to fit!julia> fm1R = fit(MixedModel, @formula(Y ~ 1 + (1|G)),\n    dat[:Dyestuff], REML=true)\nError: UndefVarError: dat not defined\n"
},

{
    "location": "constructors/#Simple,-scalar-random-effects-1",
    "page": "Model constructors",
    "title": "Simple, scalar random effects",
    "category": "section",
    "text": "A simple, scalar random effects term in a mixed-effects model formula is of the form (1|G). All random effects terms end with |G where G is the grouping factor for the random effect. The name or, more generally, the expression G should evaluate to a categorical array that has a distinct set of levels. The random effects are associated with the levels of the grouping factor.A scalar random effect is, as the name implies, one scalar value for each level of the grouping factor. A simple, scalar random effects term is of the form, (1|G). It corresponds to a shift in the intercept for each level of the grouping factor."
},

{
    "location": "constructors/#Models-with-vector-valued-random-effects-1",
    "page": "Model constructors",
    "title": "Models with vector-valued random effects",
    "category": "section",
    "text": "The sleepstudy data are observations of reaction time, Y, on several subjects, G, after 0 to 9 days of sleep deprivation, U. A model with random intercepts and random slopes for each subject, allowing for within-subject correlation of the slope and intercept, is fit asjulia> fm2 = fit(MixedModel, @formula(Y ~ 1+U + (1+U|G)), dat[:sleepstudy])\nError: UndefVarError: dat not defined\n"
},

{
    "location": "constructors/#Models-with-multiple,-scalar-random-effects-terms-1",
    "page": "Model constructors",
    "title": "Models with multiple, scalar random-effects terms",
    "category": "section",
    "text": "A model for the Penicillin data incorporates random effects for the plate, G, and for the sample, H. As every sample is used on every plate these two factors are crossed.julia> fm4 = fit(MixedModel,@formula(Y ~ 1+(1|G)+(1|H)),dat[:Penicillin])\nError: UndefVarError: dat not defined\nIn contrast the sample, G, grouping factor is nested within the batch, H, grouping factor in the Pastes data. That is, each level of G occurs in conjunction with only one level of H.julia> fm5 = fit(MixedModel, @formula(Y ~ 1 + (1|G) + (1|H)), dat[:Pastes])\nError: UndefVarError: dat not defined\nIn observational studies it is common to encounter partially crossed grouping factors. For example, the InstEval data are course evaluations by students, G, of instructors, H. Additional covariates include the academic department, I, in which the course was given and A, whether or not it was a service course.julia> fm6 = fit(MixedModel,@formula(Y ~ 1+ A*I +(1|G)+(1|H)),dat[:InstEval])\nError: UndefVarError: dat not defined\n"
},

{
    "location": "constructors/#MixedModels.zerocorr!",
    "page": "Model constructors",
    "title": "MixedModels.zerocorr!",
    "category": "function",
    "text": "zerocorr!(m::LinearMixedModel[, trmnms::Vector{Symbol}])\n\nRewrite the random effects specification for the grouping factors in trmnms to zero correlation parameter.\n\nThe default for trmnms is all the names of random-effects terms.\n\nA random effects term is in the zero correlation parameter configuration when the off-diagonal elements of λ are all zero - hence there are no correlation parameters in that term being estimated.\n\n\n\n\n\n"
},

{
    "location": "constructors/#MixedModels.zerocorr",
    "page": "Model constructors",
    "title": "MixedModels.zerocorr",
    "category": "function",
    "text": "zerocorr(term::RandomEffectsTerm)\n\nRemove correlations between random effects in term.\n\n\n\n\n\n"
},

{
    "location": "constructors/#Simplifying-the-random-effect-correlation-structure-1",
    "page": "Model constructors",
    "title": "Simplifying the random effect correlation structure",
    "category": "section",
    "text": "MixedEffects.jl estimates not only the variance of the effects for each random effect level, but also the correlation between the random effects for different predictors. So, for the model of the sleepstudy data above, one of the parameters that is estimated is the correlation between each subject\'s random intercept (i.e., their baseline reaction time) and slope (i.e., their particular change in reaction time over days of sleep deprivation). In some cases, you may wish to simplify the random effects structure by removing these correlation parameters. This often arises when there are many random effects you want to estimate (as is common in psychological experiments with many conditions and covariates), since the number of random effects parameters increases as the square of the number of predictors, making these models difficult to estimate from limited data.A model with uncorrelated random effects for the intercept and slope by subject is fit asjulia> fm3 = fit!(zerocorr!(LinearMixedModel(@formula(Y ~ 1+U+(1+U|G)),\n    dat[:sleepstudy])))\nError: UndefVarError: dat not defined\nNote that the use of zerocorr! requires the model to be constructed, then altered to eliminate the correlation of the random effects, then fit with a call to the mutating function, fit!.zerocorr!The special syntax zerocorr can be applied to individual random effects terms inside the @formula:zerocorrjulia> fit(MixedModel, @formula(Y ~ 1 + U + zerocorr(1+U|G)), dat[:sleepstudy])\nError: UndefVarError: dat not defined\nAlternatively, correlations between parameters can be removed by including them as separate random effects terms:julia> fit(MixedModel, @formula(Y ~ 1+U+(1|G)+(0+U|G)), dat[:sleepstudy])\nError: UndefVarError: dat not defined\nNote that it is necessary to explicitly block the inclusion of an intercept term by adding 0 in the random-effects term (0+U|G).Finally, for predictors that are categorical, MixedModels.jl will estimate correlations between each level. Notice the large number of correlation parameters if we treat U as a categorical variable by giving it contrasts:julia> fit(MixedModel, @formula(Y ~ 1+U+(1+U|G)), dat[:sleepstudy],\n    contrasts=Dict(:U=>DummyCoding()))\nError: UndefVarError: dat not defined\nSeparating the 1 and U random effects into separate terms removes the correlations between the intercept and the levels of U, but not between the levels themselves:julia> fit(MixedModel, @formula(Y ~ 1+U+(1|G)+(0+U|G)), dat[:sleepstudy],\n    contrasts=Dict(:U=>DummyCoding()))\nError: UndefVarError: dat not defined\nBut using zerocorr on the individual terms (or zerocorr! on the constructed model object as above) does remove the correlations between the levels:julia> fit(MixedModel, @formula(Y ~ 1+U+zerocorr(1+U|G)), dat[:sleepstudy],\n    contrasts=Dict(:U=>DummyCoding()))\nError: UndefVarError: dat not defined\n\njulia> fit(MixedModel, @formula(Y ~ 1+U+(1|G)+zerocorr(0+U|G)),dat[:sleepstudy],\n    contrasts=Dict(:U=>DummyCoding()))\nError: UndefVarError: dat not defined\n"
},

{
    "location": "constructors/#MixedModels.GeneralizedLinearMixedModel",
    "page": "Model constructors",
    "title": "MixedModels.GeneralizedLinearMixedModel",
    "category": "type",
    "text": "GeneralizedLinearMixedModel\n\nGeneralized linear mixed-effects model representation\n\nFields\n\nLMM: a LinearMixedModel - the local approximation to the GLMM.\nβ: the pivoted and possibly truncated fixed-effects vector\nβ₀: similar to β. Used in the PIRLS algorithm if step-halving is needed.\nθ: covariance parameter vector\nb: similar to u, equivalent to broadcast!(*, b, LMM.Λ, u)\nu: a vector of matrices of random effects\nu₀: similar to u.  Used in the PIRLS algorithm if step-halving is needed.\nresp: a GlmResp object\nη: the linear predictor\nwt: vector of prior case weights, a value of T[] indicates equal weights.\n\nThe following fields are used in adaptive Gauss-Hermite quadrature, which applies only to models with a single random-effects term, in which case their lengths are the number of levels in the grouping factor for that term.  Otherwise they are zero-length vectors.\n\ndevc: vector of deviance components\ndevc0: vector of deviance components at offset of zero\nsd: approximate standard deviation of the conditional density\nmult: multiplier\n\nProperties\n\nIn addition to the fieldnames, the following names are also accessible through the . extractor\n\ntheta: synonym for θ\nbeta: synonym for β\nσ or sigma: common scale parameter (value is NaN for distributions without a scale parameter)\nlowerbd: vector of lower bounds on the combined elements of β and θ\nformula, trms, A, L, and optsum: fields of the LMM field\nX: fixed-effects model matrix\ny: response vector\n\n\n\n\n\n"
},

{
    "location": "constructors/#Fitting-generalized-linear-mixed-models-1",
    "page": "Model constructors",
    "title": "Fitting generalized linear mixed models",
    "category": "section",
    "text": "To create a GLMM representationGeneralizedLinearMixedModelthe distribution family for the response, and possibly the link function, must be specified.julia> verbaggform = @formula(r2 ~ 1+a+g+b+s+m+(1|id)+(1|item));\n\njulia> gm1 = fit(MixedModel, verbaggform, dat[:VerbAgg], Bernoulli())\nError: UndefVarError: dat not defined\nThe canonical link, which is GLM.LogitLink for the Bernoulli distribution, is used if no explicit link is specified.Note that, in keeping with convention in the GLM package, the distribution family for a binary (i.e. 0/1) response is the Bernoulli distribution. The Binomial distribution is only used when the response is the fraction of trials returning a positive, in which case the number of trials must be specified as the case weights."
},

{
    "location": "constructors/#Optional-arguments-to-fit!-1",
    "page": "Model constructors",
    "title": "Optional arguments to fit!",
    "category": "section",
    "text": "An alternative approach is to create the GeneralizedLinearMixedModel object then call fit! on it. In this form optional arguments fast and/or nAGQ can be passed to the optimization process.As the name implies, fast=true, provides a faster but somewhat less accurate fit. These fits may suffice for model comparisons.julia> gm1a = fit!(GeneralizedLinearMixedModel(verbaggform, dat[:VerbAgg],\n    Bernoulli()), fast=true)\nError: UndefVarError: dat not defined\n\njulia> deviance(gm1a) - deviance(gm1)\nError: UndefVarError: gm1a not defined\n\njulia> @time fit!(GeneralizedLinearMixedModel(verbaggform, dat[:VerbAgg],\n    Bernoulli()));\nError: UndefVarError: dat not defined\n\njulia> @time fit!(GeneralizedLinearMixedModel(verbaggform, dat[:VerbAgg],\n    Bernoulli()), fast=true);\nError: UndefVarError: dat not defined\nThe optional argument nAGQ=k causes evaluation of the deviance function to use a k point adaptive Gauss-Hermite quadrature rule. This method only applies to models with a single, simple, scalar random-effects term, such asjulia> contraform = @formula(use ~ 1+a+abs2(a)+l+urb+(1|d));\n\njulia> @time gm2 = fit!(GeneralizedLinearMixedModel(contraform,\n    dat[:Contraception], Bernoulli()), nAGQ=9)\nError: UndefVarError: dat not defined\n\njulia> @time deviance(fit!(GeneralizedLinearMixedModel(contraform,\n    dat[:Contraception], Bernoulli()), nAGQ=9, fast=true))\nError: UndefVarError: dat not defined\n\njulia> @time deviance(fit!(GeneralizedLinearMixedModel(contraform,\n    dat[:Contraception], Bernoulli())))\nError: UndefVarError: dat not defined\n\njulia> @time deviance(fit!(GeneralizedLinearMixedModel(contraform,\n    dat[:Contraception], Bernoulli()), fast=true))\nError: UndefVarError: dat not defined\n"
},

{
    "location": "constructors/#Extractor-functions-1",
    "page": "Model constructors",
    "title": "Extractor functions",
    "category": "section",
    "text": "LinearMixedModel and GeneralizedLinearMixedModel are subtypes of StatsBase.RegressionModel which, in turn, is a subtype of StatsBase.StatisticalModel. Many of the generic extractors defined in the StatsBase package have methods for these models."
},

{
    "location": "constructors/#StatsBase.loglikelihood-Tuple{StatisticalModel}",
    "page": "Model constructors",
    "title": "StatsBase.loglikelihood",
    "category": "method",
    "text": "loglikelihood(obj::StatisticalModel)\n\nReturn the log-likelihood of the model.\n\n\n\n\n\n"
},

{
    "location": "constructors/#StatsBase.aic",
    "page": "Model constructors",
    "title": "StatsBase.aic",
    "category": "function",
    "text": "aic(obj::StatisticalModel)\n\nAkaike\'s Information Criterion, defined as -2 log L + 2k, with L the likelihood of the model, and k its number of consumed degrees of freedom (as returned by dof).\n\n\n\n\n\n"
},

{
    "location": "constructors/#StatsBase.bic",
    "page": "Model constructors",
    "title": "StatsBase.bic",
    "category": "function",
    "text": "bic(obj::StatisticalModel)\n\nBayesian Information Criterion, defined as -2 log L + k log n, with L the likelihood of the model,  k its number of consumed degrees of freedom (as returned by dof), and n the number of observations (as returned by nobs).\n\n\n\n\n\n"
},

{
    "location": "constructors/#StatsBase.dof-Tuple{StatisticalModel}",
    "page": "Model constructors",
    "title": "StatsBase.dof",
    "category": "method",
    "text": "dof(obj::StatisticalModel)\n\nReturn the number of degrees of freedom consumed in the model, including when applicable the intercept and the distribution\'s dispersion parameter.\n\n\n\n\n\n"
},

{
    "location": "constructors/#StatsBase.nobs-Tuple{StatisticalModel}",
    "page": "Model constructors",
    "title": "StatsBase.nobs",
    "category": "method",
    "text": "nobs(obj::StatisticalModel)\n\nReturn the number of independent observations on which the model was fitted. Be careful when using this information, as the definition of an independent observation may vary depending on the model, on the format used to pass the data, on the sampling plan (if specified), etc.\n\n\n\n\n\n"
},

{
    "location": "constructors/#StatsBase.deviance-Tuple{StatisticalModel}",
    "page": "Model constructors",
    "title": "StatsBase.deviance",
    "category": "method",
    "text": "deviance(obj::StatisticalModel)\n\nReturn the deviance of the model relative to a reference, which is usually when applicable the saturated model. It is equal, up to a constant, to -2 log L, with L the likelihood of the model.\n\n\n\n\n\n"
},

{
    "location": "constructors/#MixedModels.objective",
    "page": "Model constructors",
    "title": "MixedModels.objective",
    "category": "function",
    "text": "objective(m::LinearMixedModel)\n\nReturn negative twice the log-likelihood of model m\n\n\n\n\n\n"
},

{
    "location": "constructors/#MixedModels.deviance!",
    "page": "Model constructors",
    "title": "MixedModels.deviance!",
    "category": "function",
    "text": "deviance!(m::GeneralizedLinearMixedModel, nAGQ=1)\n\nUpdate m.η, m.μ, etc., install the working response and working weights in m.LMM, update m.LMM.A and m.LMM.R, then evaluate the deviance.\n\n\n\n\n\n"
},

{
    "location": "constructors/#Model-fit-statistics-1",
    "page": "Model constructors",
    "title": "Model-fit statistics",
    "category": "section",
    "text": "The statistics describing the quality of the model fit includeloglikelihood(::StatisticalModel)\naic\nbic\ndof(::StatisticalModel)\nnobs(::StatisticalModel)julia> loglikelihood(fm1)\nError: UndefVarError: fm1 not defined\n\njulia> aic(fm1)\nError: UndefVarError: fm1 not defined\n\njulia> bic(fm1)\nError: UndefVarError: fm1 not defined\n\njulia> dof(fm1)   # 1 fixed effect, 2 variances\nError: UndefVarError: fm1 not defined\n\njulia> nobs(fm1)  # 30 observations\nError: UndefVarError: fm1 not defined\n\njulia> loglikelihood(gm1)\nError: UndefVarError: gm1 not defined\nIn general the deviance of a statistical model fit is negative twice the log-likelihood adjusting for the saturated model.deviance(::StatisticalModel)Because it is not clear what the saturated model corresponding to a particular LinearMixedModel should be, negative twice the log-likelihood is called the objective.objectiveThis value is also accessible as the deviance but the user should bear in mind that this doesn\'t have all the properties of a deviance which is corrected for the saturated model. For example, it is not necessarily non-negative.julia> objective(fm1)\nError: UndefVarError: fm1 not defined\n\njulia> deviance(fm1)\nError: UndefVarError: fm1 not defined\nThe value optimized when fitting a GeneralizedLinearMixedModel is the Laplace approximation to the deviance or an adaptive Gauss-Hermite evaluation.MixedModels.deviance!julia> MixedModels.deviance!(gm1)\nError: UndefVarError: gm1 not defined\n"
},

{
    "location": "constructors/#StatsBase.coef",
    "page": "Model constructors",
    "title": "StatsBase.coef",
    "category": "function",
    "text": "coef(obj::StatisticalModel)\n\nReturn the coefficients of the model.\n\n\n\n\n\n"
},

{
    "location": "constructors/#MixedModels.fixef",
    "page": "Model constructors",
    "title": "MixedModels.fixef",
    "category": "function",
    "text": "fixef(m::MixedModel, permuted=true)\n\nReturn the fixed-effects parameter vector estimate of m.\n\nIf permuted is true the vector elements are permuted according to first(m.feterms).piv and truncated to the rank of that term.\n\n\n\n\n\n"
},

{
    "location": "constructors/#StatsBase.vcov",
    "page": "Model constructors",
    "title": "StatsBase.vcov",
    "category": "function",
    "text": "vcov(obj::StatisticalModel)\n\nReturn the variance-covariance matrix for the coefficients of the model.\n\n\n\n\n\n"
},

{
    "location": "constructors/#StatsBase.stderror",
    "page": "Model constructors",
    "title": "StatsBase.stderror",
    "category": "function",
    "text": "stderror(obj::StatisticalModel)\n\nReturn the standard errors for the coefficients of the model.\n\n\n\n\n\n"
},

{
    "location": "constructors/#StatsBase.coeftable",
    "page": "Model constructors",
    "title": "StatsBase.coeftable",
    "category": "function",
    "text": "coeftable(obj::StatisticalModel; level::Real=0.95)\n\nReturn a table of class CoefTable with coefficients and related statistics. level determines the level for confidence intervals (by default, 95%).\n\n\n\n\n\n"
},

{
    "location": "constructors/#Fixed-effects-parameter-estimates-1",
    "page": "Model constructors",
    "title": "Fixed-effects parameter estimates",
    "category": "section",
    "text": "The coef and fixef extractors both return the maximum likelihood estimates of the fixed-effects coefficients.coef\nfixefjulia> show(coef(fm1))\nError: UndefVarError: fm1 not defined\n\njulia> show(fixef(fm1))\nError: UndefVarError: fm1 not defined\n\njulia> show(fixef(gm1))\nError: UndefVarError: gm1 not defined\nAn alternative extractor for the fixed-effects coefficient is the β property. Properties whose names are Greek letters usually have an alternative spelling, which is the name of the Greek letter.julia> show(fm1.β)\nError: UndefVarError: fm1 not defined\n\njulia> show(fm1.beta)\nError: UndefVarError: fm1 not defined\n\njulia> show(gm1.β)\nError: UndefVarError: gm1 not defined\nA full list of property names is returned by propertynamesjulia> propertynames(fm1)\nError: UndefVarError: fm1 not defined\n\njulia> propertynames(gm1)\nError: UndefVarError: gm1 not defined\nThe variance-covariance matrix of the fixed-effects coefficients is returned byvcovjulia> vcov(fm2)\nError: UndefVarError: fm2 not defined\n\njulia> vcov(gm1)\nError: UndefVarError: gm1 not defined\nThe standard errors are the square roots of the diagonal elements of the estimated variance-covariance matrix of the fixed-effects coefficient estimators.stderrorjulia> show(StatsBase.stderror(fm2))\nError: UndefVarError: fm2 not defined\n\njulia> show(StatsBase.stderror(gm1))\nError: UndefVarError: gm1 not defined\nFinally, the coeftable generic produces a table of coefficient estimates, their standard errors, and their ratio. The p-values quoted here should be regarded as approximations.coeftablejulia> coeftable(fm2)\nError: UndefVarError: fm2 not defined\n"
},

{
    "location": "constructors/#MixedModels.VarCorr",
    "page": "Model constructors",
    "title": "MixedModels.VarCorr",
    "category": "type",
    "text": "VarCorr\n\nInformation from the fitted random-effects variance-covariance matrices.\n\nMembers\n\nσρ: a NamedTuple of NamedTuples as returned from σρs\ns: the estimate of the per-observation dispersion parameter\n\nThe main purpose of defining this type is to isolate the logic in the show method.\n\n\n\n\n\n"
},

{
    "location": "constructors/#MixedModels.varest",
    "page": "Model constructors",
    "title": "MixedModels.varest",
    "category": "function",
    "text": "varest(m::LinearMixedModel)\n\nReturns the estimate of σ², the variance of the conditional distribution of Y given B.\n\n\n\n\n\n"
},

{
    "location": "constructors/#MixedModels.sdest",
    "page": "Model constructors",
    "title": "MixedModels.sdest",
    "category": "function",
    "text": "sdest(m::LinearMixedModel)\n\nReturn the estimate of σ, the standard deviation of the per-observation noise.\n\n\n\n\n\n"
},

{
    "location": "constructors/#Covariance-parameter-estimates-1",
    "page": "Model constructors",
    "title": "Covariance parameter estimates",
    "category": "section",
    "text": "The covariance parameters estimates, in the form shown in the model summary, are a VarCorr objectVarCorrjulia> VarCorr(fm2)\nError: UndefVarError: fm2 not defined\n\njulia> VarCorr(gm1)\nError: UndefVarError: gm1 not defined\nIndividual components are returned by other extractorsvarest\nsdestjulia> varest(fm2)\nError: UndefVarError: fm2 not defined\n\njulia> sdest(fm2)\nError: UndefVarError: fm2 not defined\n\njulia> fm2.σ\nError: UndefVarError: fm2 not defined\n"
},

{
    "location": "constructors/#MixedModels.ranef",
    "page": "Model constructors",
    "title": "MixedModels.ranef",
    "category": "function",
    "text": "ranef(m::LinearMixedModel; uscale=false, named=true)\n\nReturn, as a Vector{Vector{T}} (Vector{NamedVector{T}} if named=true), the conditional modes of the random effects in model m.\n\nIf uscale is true the random effects are on the spherical (i.e. u) scale, otherwise on the original scale.\n\n\n\n\n\n"
},

{
    "location": "constructors/#MixedModels.condVar",
    "page": "Model constructors",
    "title": "MixedModels.condVar",
    "category": "function",
    "text": "condVar(m::LinearMixedModel)\n\nReturn the conditional variances matrices of the random effects.\n\nThe random effects are returned by ranef as a vector of length k, where k is the number of random effects terms.  The ith element is a matrix of size vᵢ × ℓᵢ  where vᵢ is the size of the vector-valued random effects for each of the ℓᵢ levels of the grouping factor.  Technically those values are the modes of the conditional distribution of the random effects given the observed data.\n\nThis function returns an array of k three dimensional arrays, where the ith array is of size vᵢ × vᵢ × ℓᵢ.  These are the diagonal blocks from the conditional variance-covariance matrix,\n\ns² Λ(Λ\'Z\'ZΛ + I)⁻¹Λ\'\n\n\n\n\n\n"
},

{
    "location": "constructors/#Conditional-modes-of-the-random-effects-1",
    "page": "Model constructors",
    "title": "Conditional modes of the random effects",
    "category": "section",
    "text": "The ranef extractorranefjulia> ranef(fm1)\nError: UndefVarError: fm1 not defined\n\njulia> fm1.b\nError: UndefVarError: fm1 not defined\nreturns the conditional modes of the random effects given the observed data. That is, these are the values that maximize the conditional density of the random effects given the observed data. For a LinearMixedModel these are also the conditional mean values.These are sometimes called the best linear unbiased predictors or BLUPs but that name is not particularly meaningful.At a superficial level these can be considered as the \"estimates\" of the random effects, with a bit of hand waving, but pursuing this analogy too far usually results in confusion.The corresponding conditional variances are returned bycondVarjulia> condVar(fm1)\nError: UndefVarError: fm1 not defined\n"
},

{
    "location": "optimization/#",
    "page": "Details of the parameter estimation",
    "title": "Details of the parameter estimation",
    "category": "page",
    "text": ""
},

{
    "location": "optimization/#Details-of-the-parameter-estimation-1",
    "page": "Details of the parameter estimation",
    "title": "Details of the parameter estimation",
    "category": "section",
    "text": ""
},

{
    "location": "optimization/#The-probability-model-1",
    "page": "Details of the parameter estimation",
    "title": "The probability model",
    "category": "section",
    "text": "Maximum likelihood estimates are based on the probability model for the observed responses. In the probability model the distribution of the responses is expressed as a function of one or more parameters.For a continuous distribution the probability density is a function of the responses, given the parameters. The likelihood function is the same expression as the probability density but regarding the observed values as fixed and the parameters as varying.In general a mixed-effects model incorporates two random variables: mathcalB, the q-dimensional vector of random effects, and mathcalY, the n-dimensional response vector. The value, bf y, of mathcalY is observed; the value, bf b, of mathcalB is not."
},

{
    "location": "optimization/#Linear-Mixed-Effects-Models-1",
    "page": "Details of the parameter estimation",
    "title": "Linear Mixed-Effects Models",
    "category": "section",
    "text": "In a linear mixed model the unconditional distribution of mathcalB and the conditional distribution, (mathcalY  mathcalB=bfb), are both multivariate Gaussian distributions,beginequation\nbeginaligned\n  (mathcalY  mathcalB=bfb) simmathcalN(bf Xbeta + Z bsigma^2bfI)\n  mathcalBsimmathcalN(bf0Sigma_theta) \nendaligned\nendequationThe conditional mean of mathcal Y, given mathcal B=bf b, is the linear predictor, bf Xbfbeta+bf Zbf b, which depends on the p-dimensional fixed-effects parameter, bf beta, and on bf b. The model matrices, bf X and bf Z, of dimension ntimes p and ntimes q, respectively, are determined from the formula for the model and the values of covariates. Although the matrix bf Z can be large (i.e. both n and q can be large), it is sparse (i.e. most of the elements in the matrix are zero).The relative covariance factor, Lambda_theta, is a qtimes q lower-triangular matrix, depending on the variance-component parameter, bftheta, and generating the symmetric qtimes q variance-covariance matrix, Sigma_theta, asbeginequation\nSigma_theta=sigma^2Lambda_thetaLambda_theta\nendequationThe spherical random effects, mathcalUsimmathcalN(bf0sigma^2bfI_q), determine mathcal B according tobeginequation\nmathcalB=Lambda_thetamathcalU\nendequationThe penalized residual sum of squares (PRSS),beginequation\nr^2(thetabetabfu)=bfy - bfXbeta -bfZLambda_thetabfu^2+bfu^2\nendequationis the sum of the residual sum of squares, measuring fidelity of the model to the data, and a penalty on the size of bf u, measuring the complexity of the model. Minimizing r^2 with respect to bf u,beginequation\nr^2_betatheta =min_bfuleft(bfy -bfXbeta -bfZLambda_thetabfu^2+bfu^2right)\nendequationis a direct (i.e. non-iterative) computation. The particular method used to solve this generates a blocked Choleksy factor, bfL_theta, which is an lower triangular qtimes q matrix satisfyingbeginequation\nbfL_thetabfL_theta=Lambda_thetabfZbfZLambda_theta+bfI_q \nendequationwhere bf I_q is the qtimes q identity matrix.Negative twice the log-likelihood of the parameters, given the data, bf y, isbeginequation\nd(bfthetabfbetasigmabf y)\n=nlog(2pisigma^2)+log(bf L_theta^2)+fracr^2_betathetasigma^2\nendequationwhere bf L_theta denotes the determinant of bf L_theta. Because bf L_theta is triangular, its determinant is the product of its diagonal elements.Because the conditional mean, bfmu_mathcal Ymathcal B=bf b=bf Xbfbeta+bf ZLambda_thetabf u, is a linear function of both bfbeta and bf u, minimization of the PRSS with respect to both bfbeta and bf u to producebeginequation\nr^2_theta =min_bfbetabf uleft(bf y -bf Xbfbeta -bf ZLambda_thetabf u^2+bf u^2right)\nendequationis also a direct calculation. The values of bf u and bfbeta that provide this minimum are called, respectively, the conditional mode, tildebf u_theta, of the spherical random effects and the conditional estimate, widehatbfbeta_theta, of the fixed effects. At the conditional estimate of the fixed effects the objective isbeginequation\nd(bfthetawidehatbeta_thetasigmabf y)\n=nlog(2pisigma^2)+log(bf L_theta^2)+fracr^2_thetasigma^2\nendequationMinimizing this expression with respect to sigma^2 produces the conditional estimatebeginequation\nwidehatsigma^2_theta=fracr^2_thetan\nendequationwhich provides the profiled log-likelihood on the deviance scale asbeginequation\ntilded(thetabf y)=d(thetawidehatbeta_thetawidehatsigma_thetabf y)\n=log(bf L_theta^2)+nleft1+logleft(frac2pi r^2_thetanright)right\nendequationa function of bftheta alone.The MLE of bftheta, written widehatbftheta, is the value that minimizes this profiled objective. We determine this value by numerical optimization. In the process of evaluating tilded(widehatthetabf y) we determine widehatbeta=widehatbeta_widehattheta, tildebf u_widehattheta and r^2_widehattheta, from which we can evaluate widehatsigma=sqrtr^2_widehatthetan.The elements of the conditional mode of mathcal B, evaluated at the parameter estimates,beginequation\ntildebf b_widehattheta=Lambda_widehatthetatildebf u_widehattheta\nendequationare sometimes called the best linear unbiased predictors or BLUPs of the random effects. Although BLUPs an appealing acronym, I don’t find the term particularly instructive (what is a “linear unbiased predictor” and in what sense are these the “best”?) and prefer the term “conditional modes”, because these are the values of bf b that maximize the density of the conditional distribution mathcalB  mathcalY = bf y. For a linear mixed model, where all the conditional and unconditional distributions are Gaussian, these values are also the conditional means."
},

{
    "location": "optimization/#MixedModels.ReMat",
    "page": "Details of the parameter estimation",
    "title": "MixedModels.ReMat",
    "category": "type",
    "text": "ReMat{T,S} <: AbstractMatrix{T}\n\nA section of a model matrix generated by a random-effects term.\n\nFields\n\ntrm: the grouping factor as a StatsModels.CategoricalTerm\nrefs: indices into the levels of the grouping factor as a Vector{Int32}\nz: transpose of the model matrix generated by the left-hand side of the term\nwtz: a weighted copy of z (z and wtz are the same object for unweighted cases)\nλ: a LowerTriangular matrix of size S×S\ninds: a Vector{Int} of linear indices of the potential nonzeros in λ\nadjA: the adjoint of the matrix as a SparseMatrixCSC{T}\n\n\n\n\n\n"
},

{
    "location": "optimization/#Internal-structure-of-\\Lambda_\\theta-and-\\bf-Z-1",
    "page": "Details of the parameter estimation",
    "title": "Internal structure of Lambda_theta and bf Z",
    "category": "section",
    "text": "In the types of LinearMixedModel available through the MixedModels package, groups of random effects and the corresponding columns of the model matrix, bf Z, are associated with random-effects terms in the model formula.For the simple examplejulia> fm1 = fit(MixedModel, @formula(Y ~ 1 + (1|G)), dat[:Dyestuff])\nError: UndefVarError: dat not defined\nthe only random effects term in the formula is (1|G), a simple, scalar random-effects term.julia> t1 = first(fm1.reterms);\nError: UndefVarError: fm1 not defined\n\njulia> Int.(t1)  # convert to integers for more compact display\nError: UndefVarError: t1 not defined\nReMatThis RandomEffectsTerm contributes a block of columns to the model matrix bf Z and a diagonal block to Lambda_theta. In this case the diagonal block of Lambda_theta (which is also the only block) is a multiple of the 6times6 identity matrix where the multiple isjulia> t1.λ\nError: UndefVarError: t1 not defined\nBecause there is only one random-effects term in the model, the matrix bf Z is the indicators matrix shown as the result of Matrix(t1), but stored in a special sparse format. Furthermore, there is only one block in Lambda_theta.For a vector-valued random-effects term, as injulia> fm2 = fit(MixedModel, @formula(Y ~ 1+U+(1+U|G)), dat[:sleepstudy])\nError: UndefVarError: dat not defined\nthe model matrix bf Z is of the formjulia> t21 = first(fm2.reterms);\nError: UndefVarError: fm2 not defined\n\njulia> Int.(t21) # convert to integers for more compact display\nError: UndefVarError: t21 not defined\nand Lambda_theta is a 36times36 block diagonal matrix with 18 diagonal blocks, all of the formjulia> t21.λ\nError: UndefVarError: t21 not defined\nThe theta vector isjulia> MixedModels.getθ(t21)\nError: UndefVarError: t21 not defined\nRandom-effects terms in the model formula that have the same grouping factor are amagamated into a single ReMat object.julia> fm3 = fit!(zerocorr!(LinearMixedModel(@formula(Y ~ 1+U+(1+U|G)),\n    dat[:sleepstudy])))\nError: UndefVarError: dat not defined\n\njulia> t31 = first(fm3.reterms);\nError: UndefVarError: fm3 not defined\n\njulia> Int.(t31)\nError: UndefVarError: t31 not defined\nNote that we could also have achieved this by re-fitting (a copy of) fm2.julia> fm3alt = fit!(zerocorr!(deepcopy(fm2)))\nError: UndefVarError: fm2 not defined\nFor this model the matrix bf Z is the same as that of model fm2 but the diagonal blocks of Lambda_theta are themselves diagonal.julia> t31.λ\nError: UndefVarError: t31 not defined\n\njulia> MixedModels.getθ(t31)\nError: UndefVarError: t31 not defined\nRandom-effects terms with distinct grouping factors generate distinct elements of the trms member of the LinearMixedModel object. Multiple ReMat objects are sorted by decreasing numbers of random effects.julia> fm4 = fit(MixedModel, @formula(Y ~ 1 + (1|H) + (1|G)),\n    dat[:Penicillin])\nError: UndefVarError: dat not defined\n\njulia> Int.(first(fm4.reterms))\nError: UndefVarError: fm4 not defined\n\njulia> Int.(last(fm4.reterms))\nError: UndefVarError: fm4 not defined\nNote that the first ReMat in fm4.terms corresponds to grouping factor G even though the term (1|G) occurs in the formula after (1|H)."
},

{
    "location": "optimization/#MixedModels.OptSummary",
    "page": "Details of the parameter estimation",
    "title": "MixedModels.OptSummary",
    "category": "type",
    "text": "OptSummary\n\nSummary of an NLopt optimization\n\nFields\n\ninitial: a copy of the initial parameter values in the optimization\nlowerbd: lower bounds on the parameter values\nftol_rel: as in NLopt\nftol_abs: as in NLopt\nxtol_rel: as in NLopt\nxtol_abs: as in NLopt\ninitial_step: as in NLopt\nmaxfeval: as in NLopt\nfinal: a copy of the final parameter values from the optimization\nfmin: the final value of the objective\nfeval: the number of function evaluations\noptimizer: the name of the optimizer used, as a Symbol\nreturnvalue: the return value, as a Symbol\nnAGQ: number of adaptive Gauss-Hermite quadrature points in deviance evaluation for GLMMs\nREML: use the REML criterion for LMM fits\n\nThe latter field doesn\'t really belong here but it has to be in a mutable struct in case it is changed.\n\n\n\n\n\n"
},

{
    "location": "optimization/#Progress-of-the-optimization-1",
    "page": "Details of the parameter estimation",
    "title": "Progress of the optimization",
    "category": "section",
    "text": "An optional named argument, verbose=true, in the call to fit for a LinearMixedModel causes printing of the objective and the theta parameter at each evaluation during the optimization.julia> fit(MixedModel, @formula(Y ~ 1 + (1|G)), dat[:Dyestuff],\n    verbose=true);\nError: UndefVarError: dat not defined\n\njulia> fit(MixedModel, @formula(Y ~ 1 + U + (1+U|G)), dat[:sleepstudy],\n    verbose=true);\nError: UndefVarError: dat not defined\nA shorter summary of the optimization process is always available as anOptSummaryobject, which is the optsum member of the LinearMixedModel.julia> fm2.optsum\nError: UndefVarError: fm2 not defined\n"
},

{
    "location": "optimization/#Modifying-the-optimization-process-1",
    "page": "Details of the parameter estimation",
    "title": "Modifying the optimization process",
    "category": "section",
    "text": "The OptSummary object contains both input and output fields for the optimizer. To modify the optimization process the input fields can be changed after constructing the model but before fitting it.Suppose, for example, that the user wishes to try a Nelder-Mead optimization method instead of the default BOBYQA (Bounded Optimization BY Quadratic Approximation) method.julia> fm2 = LinearMixedModel(@formula(Y ~ 1+U+(1+U|G)), dat[:sleepstudy]);\nError: UndefVarError: dat not defined\n\njulia> fm2.optsum.optimizer = :LN_NELDERMEAD;\nError: UndefVarError: fm2 not defined\n\njulia> fit!(fm2)\nError: UndefVarError: fm2 not defined\n\njulia> fm2.optsum\nError: UndefVarError: fm2 not defined\nThe parameter estimates are quite similar to those using :LN_BOBYQA but at the expense of 140 functions evaluations for :LN_NELDERMEAD versus 57 for :LN_BOBYQA.See the documentation for the NLopt package for details about the various settings."
},

{
    "location": "optimization/#MixedModels.issingular",
    "page": "Details of the parameter estimation",
    "title": "MixedModels.issingular",
    "category": "function",
    "text": "issingular(m::LinearMixedModel, θ=m.θ)\n\nTest whether the model m is singular if the parameter vector is θ.\n\n\n\n\n\n"
},

{
    "location": "optimization/#Convergence-to-singular-covariance-matrices-1",
    "page": "Details of the parameter estimation",
    "title": "Convergence to singular covariance matrices",
    "category": "section",
    "text": "To ensure identifiability of Sigma_theta=sigma^2Lambda_theta Lambda_theta, the elements of theta corresponding to diagonal elements of Lambda_theta are constrained to be non-negative. For example, in a trivial case of a single, simple, scalar, random-effects term as in fm1, the one-dimensional theta vector is the ratio of the standard deviation of the random effects to the standard deviation of the response. It happens that -theta produces the same log-likelihood but, by convention, we define the standard deviation to be the positive square root of the variance. Requiring the diagonal elements of Lambda_theta to be non-negative is a generalization of using this positive square root.If the optimization converges on the boundary of the feasible region, that is if one or more of the diagonal elements of Lambda_theta is zero at convergence, the covariance matrix Sigma_theta will be singular. This means that there will be linear combinations of random effects that are constant. Usually convergence to a singular covariance matrix is a sign of an over-specified model.Singularity can be checked with the issingular predicate function.issingularjulia> issingular(fm2)\nError: UndefVarError: fm2 not defined\n"
},

{
    "location": "optimization/#Distributions.Bernoulli",
    "page": "Details of the parameter estimation",
    "title": "Distributions.Bernoulli",
    "category": "type",
    "text": "Bernoulli(p)\n\nA Bernoulli distribution is parameterized by a success rate p, which takes value 1 with probability p and 0 with probability 1-p.\n\nP(X = k) = begincases\n1 - p  quad textfor  k = 0 \np  quad textfor  k = 1\nendcases\n\nBernoulli()    # Bernoulli distribution with p = 0.5\nBernoulli(p)   # Bernoulli distribution with success rate p\n\nparams(d)      # Get the parameters, i.e. (p,)\nsuccprob(d)    # Get the success rate, i.e. p\nfailprob(d)    # Get the failure rate, i.e. 1 - p\n\nExternal links:\n\nBernoulli distribution on Wikipedia\n\n\n\n\n\n"
},

{
    "location": "optimization/#Distributions.Poisson",
    "page": "Details of the parameter estimation",
    "title": "Distributions.Poisson",
    "category": "type",
    "text": "Poisson(λ)\n\nA Poisson distribution descibes the number of independent events occurring within a unit time interval, given the average rate of occurrence λ.\n\nP(X = k) = fraclambda^kk e^-lambda quad text for  k = 012ldots\n\nPoisson()        # Poisson distribution with rate parameter 1\nPoisson(lambda)       # Poisson distribution with rate parameter lambda\n\nparams(d)        # Get the parameters, i.e. (λ,)\nmean(d)          # Get the mean arrival rate, i.e. λ\n\nExternal links:\n\nPoisson distribution on Wikipedia\n\n\n\n\n\n"
},

{
    "location": "optimization/#Generalized-Linear-Mixed-Effects-Models-1",
    "page": "Details of the parameter estimation",
    "title": "Generalized Linear Mixed-Effects Models",
    "category": "section",
    "text": "In a generalized linear model the responses are modelled as coming from a particular distribution, such as Bernoulli for binary responses or Poisson for responses that represent counts. The scalar distributions of individual responses differ only in their means, which are determined by a linear predictor expression eta=bf Xbeta, where, as before, bf X is a model matrix derived from the values of covariates and beta is a vector of coefficients.The unconstrained components of eta are mapped to the, possiby constrained, components of the mean response, mu, via a scalar function, g^-1, applied to each component of eta. For historical reasons, the inverse of this function, taking components of mu to the corresponding component of eta is called the link function and more frequently used map from eta to mu is the inverse link.A generalized linear mixed-effects model (GLMM) is defined, for the purposes of this package, bybeginequation\nbeginaligned\n  (mathcalY  mathcalB=bfb) simmathcalD(bfg^-1(Xbeta + Z b)phi)\n  mathcalBsimmathcalN(bf0Sigma_theta) \nendaligned\nendequationwhere mathcalD indicates the distribution family parameterized by the mean and, when needed, a common scale parameter, phi. (There is no scale parameter for Bernoulli or for Poisson. Specifying the mean completely determines the distribution.)Bernoulli\nPoissonA GeneralizedLinearMixedModel object is generated from a formula, data frame and distribution family.julia> const vaform = @formula(r2 ~ 1 + a + g + b + s + (1|id) + (1|item));\n\njulia> mdl = GeneralizedLinearMixedModel(vaform, dat[:VerbAgg], Bernoulli());\nError: UndefVarError: dat not defined\n\njulia> typeof(mdl)\nError: UndefVarError: mdl not defined\nA separate call to fit! can be used to fit the model. This involves optimizing an objective function, the Laplace approximation to the deviance, with respect to the parameters, which are beta, the fixed-effects coefficients, and theta, the covariance parameters. The starting estimate for beta is determined by fitting a GLM to the fixed-effects part of the formulajulia> mdl.β\nError: UndefVarError: mdl not defined\nand the starting estimate for theta, which is a vector of the two standard deviations of the random effects, is chosen to bejulia> mdl.θ\nError: UndefVarError: mdl not defined\nThe Laplace approximation to the deviance requires determining the conditional modes of the random effects. These are the values that maximize the conditional density of the random effects, given the model parameters and the data. This is done using Penalized Iteratively Reweighted Least Squares (PIRLS). In most cases PIRLS is fast and stable. It is simply a penalized version of the IRLS algorithm used in fitting GLMs.The distinction between the \"fast\" and \"slow\" algorithms in the MixedModels package (nAGQ=0 or nAGQ=1 in lme4) is whether the fixed-effects parameters, beta, are optimized in PIRLS or in the nonlinear optimizer. In a call to the pirls! function the first argument is a GeneralizedLinearMixedModel, which is modified during the function call. (By convention, the names of such mutating functions end in ! as a warning to the user that they can modify an argument, usually the first argument.) The second and third arguments are optional logical values indicating if beta is to be varied and if verbose output is to be printed.julia> pirls!(mdl, true, true)\nError: UndefVarError: mdl not defined\njulia> deviance(mdl)\nError: UndefVarError: mdl not defined\njulia> mdl.β\nError: UndefVarError: mdl not defined\njulia> mdl.θ # current values of the standard deviations of the random effects\nError: UndefVarError: mdl not defined\nIf the optimization with respect to beta is performed within PIRLS then the nonlinear optimization of the Laplace approximation to the deviance requires optimization with respect to theta only. This is the \"fast\" algorithm. Given a value of theta, PIRLS is used to determine the conditional estimate of beta and the conditional mode of the random effects, b.julia> mdl.b # conditional modes of b\nError: UndefVarError: mdl not defined\njulia> fit!(mdl, fast=true, verbose=true);\nError: UndefVarError: mdl not defined\nThe optimization process is summarized byjulia> mdl.LMM.optsum\nError: UndefVarError: mdl not defined\nAs one would hope, given the name of the option, this fit is comparatively fast.julia> @time(fit!(GeneralizedLinearMixedModel(vaform,\n    dat[:VerbAgg], Bernoulli()), fast=true))\nError: UndefVarError: dat not defined\nThe alternative algorithm is to use PIRLS to find the conditional mode of the random effects, given beta and theta and then use the general nonlinear optimizer to fit with respect to both beta and theta. Because it is slower to incorporate the beta parameters in the general nonlinear optimization, the fast fit is performed first and used to determine starting estimates for the more general optimization.julia> @time mdl1 = fit(MixedModel, vaform, dat[:VerbAgg], Bernoulli())\nError: UndefVarError: dat not defined\nThis fit provided slightly better results (Laplace approximation to the deviance of 8151.400 versus 8151.583) but took 6 times as long. That is not terribly important when the times involved are a few seconds but can be important when the fit requires many hours or days of computing time.The comparison of the slow and fast fit is available in the optimization summary after the slow fit.julia> mdl1.LMM.optsum\nError: UndefVarError: mdl1 not defined\n"
},

{
    "location": "GaussHermite/#",
    "page": "Normalized Gauss-Hermite Quadrature",
    "title": "Normalized Gauss-Hermite Quadrature",
    "category": "page",
    "text": ""
},

{
    "location": "GaussHermite/#Normalized-Gauss-Hermite-Quadrature-1",
    "page": "Normalized Gauss-Hermite Quadrature",
    "title": "Normalized Gauss-Hermite Quadrature",
    "category": "section",
    "text": "Gaussian Quadrature rules provide sets of x values, called abscissae, and weights, w, to approximate an integral with respect to a weight function, g(x). For a kth order rule the approximation isint f(x)g(x)dx approx sum_i=1^k w_i f(x_i)For the Gauss-Hermite rule the weight function isg(x) = e^-x^2and the domain of integration is (-infty infty). A slight variation of this is the normalized Gauss-Hermite rule for which the weight function is the standard normal densityg(z) = phi(z) = frace^-z^22sqrt2piThus, the expected value of f(z), where mathcalZsimmathscrN(01), is approximated asmathbbEf=int_-infty^infty f(z) phi(z)dzapproxsum_i=1^k w_if(z_i) Naturally, there is a caveat. For the approximation to be accurate the function f(z) must behave like a low-order polynomial over the range of interest. More formally, a kth order rule is exact when f is a k-1 order polynomial."
},

{
    "location": "GaussHermite/#MixedModels.GHnorm",
    "page": "Normalized Gauss-Hermite Quadrature",
    "title": "MixedModels.GHnorm",
    "category": "function",
    "text": "GHnorm(k::Int)\n\nReturn the (unique) GaussHermiteNormalized{k} object.\n\nThe function values are stored (memoized) when first evaluated.  Subsequent evaluations for the same k have very low overhead.\n\n\n\n\n\n"
},

{
    "location": "GaussHermite/#Evaluating-the-weights-and-abscissae-1",
    "page": "Normalized Gauss-Hermite Quadrature",
    "title": "Evaluating the weights and abscissae",
    "category": "section",
    "text": "In the Golub-Welsch algorithm the abscissae for a particular Gaussian quadrature rule are determined as the eigenvalues of a symmetric tri-diagonal matrix and the weights are derived from the squares of the first row of the matrix of eigenvectors. For a kth order normalized Gauss-Hermite rule the tridiagonal matrix has zeros on the diagonal and the square roots of 1:k-1 on the super- and sub-diagonal, e.g.julia> using LinearAlgebra, Gadfly\n\njulia> sym3 = SymTridiagonal(zeros(3), sqrt.(1:2))\n3×3 LinearAlgebra.SymTridiagonal{Float64,Array{Float64,1}}:\n 0.0  1.0       ⋅     \n 1.0  0.0      1.41421\n  ⋅   1.41421  0.0    \n\njulia> ev = eigen(sym3);\n\njulia> show(ev.values)\n[-1.7320508075688739, 1.1102230246251565e-15, 1.7320508075688774]\njulia> show(abs2.(ev.vectors[1,:]))\n[0.16666666666666743, 0.6666666666666657, 0.16666666666666677]As a function of k this can be written asfunction gausshermitenorm(k)\n    ev = eigen(SymTridiagonal(zeros(k), sqrt.(1:k-1)))\n    ev.values, abs2.(ev.vectors[1,:])\nendprovidingjulia> gausshermitenorm(3)\n([-1.7320508075688739, 1.1102230246251565e-15, 1.7320508075688774], [0.16666666666666743, 0.6666666666666657, 0.16666666666666677])\nThe weights and positions are often shown as a lollipop plot. For the 9th order rule these are (Image: Lollipop plot of 9th order normalized Gauss-Hermite rule)Notice that the magnitudes of the weights drop quite dramatically away from zero, even on a logarithmic scale (Image: Lollipop plot of 9th order normalized Gauss-Hermite rule (logarithmic scale)The definition of MixedModels.GHnorm is similar to the gausshermitenorm function with some extra provisions for ensuring symmetry of the abscissae and the weights and for caching values once they have been calculated.GHnormjulia> using MixedModels\n\njulia> GHnorm(3)\nMixedModels.GaussHermiteNormalized{3}([-1.7320508075688772, 0.0, 1.7320508075688772], [0.16666666666666666, 0.6666666666666666, 0.16666666666666666])\nBy the properties of the normal distribution, when mathcalXsimmathscrN(mu sigma^2)mathbbEg(x) approx sum_i=1^k g(mu + sigma z_i)w_iFor example, mathbbEmathcalX^2 where mathcalXsimmathcalN(2 3^2) isjulia> μ = 2; σ = 3; ghn3 = GHnorm(3);\n\njulia> sum(@. ghn3.w * abs2(μ + σ * ghn3.z))  # should be μ² + σ² = 13\n13.0\n(In general a dot, \'.\', after the function name in a function call, as in abs2.(...), or before an operator creates a fused vectorized evaluation in Julia. The macro @. has the effect of vectorizing all operations in the subsequent expression.)"
},

{
    "location": "GaussHermite/#Application-to-a-model-for-contraception-use-1",
    "page": "Normalized Gauss-Hermite Quadrature",
    "title": "Application to a model for contraception use",
    "category": "section",
    "text": "A binary response is a \"Yes\"/\"No\" type of answer. For example, in a 1989 fertility survey of women in Bangladesh (reported in Huq, N. M. and Cleland, J., 1990) one response of interest was whether the woman used artificial contraception. Several covariates were recorded including the woman\'s age (centered at the mean), the number of live children the woman has had (in 4 categories: 0, 1, 2, and 3 or more), whether she lived in an urban setting, and the district in which she lived. The version of the data used here is that used in review of multilevel modeling software conducted by the Center for Multilevel Modelling, currently at University of Bristol (http://www.bristol.ac.uk/cmm/learning/mmsoftware/data-rev.html). These data are available as the Contraception data frame in the test data for the MixedModels package.julia> using DataFrames, PooledArrays, RData\n\njulia> const dat = Dict(Symbol(k)=>v for (k,v) in\n    load(joinpath(dirname(pathof(MixedModels)), \"..\", \"test\", \"dat.rda\")));\nError: SystemError: opening file \"/home/travis/build/JuliaStats/MixedModels.jl/src/../test/dat.rda\": No such file or directory\n\njulia> contra = dat[:Contraception];\nError: UndefVarError: dat not defined\n\njulia> contra.urbdist = PooledArray(string.(contra.urb, contra.d));\nError: UndefVarError: contra not defined\n\njulia> describe(contra)\nError: UndefVarError: contra not defined\nA smoothed scatterplot of contraception use versus ageError: UndefVarError: contra not definedshows that the proportion of women using artificial contraception is approximately quadratic in age.A model with fixed-effects for age, age squared, number of live children and urban location and with random effects for district, is fit asjulia> const form1 = @formula use ~ 1 + a + abs2(a) + l + urb + (1|d);\n\njulia> m1 = fit!(GeneralizedLinearMixedModel(form1, contra,\n    Bernoulli()), fast=true)\nError: UndefVarError: contra not defined\nFor a model such as m1, which has a single, scalar random-effects term, the unscaled conditional density of the spherical random effects variable, mathcalU, given the observed data, mathcalY=mathbfy_0, can be expressed as a product of scalar density functions, f_i(u_i) i=1dotsq. In the PIRLS algorithm, which determines the conditional mode vector, tildemathbfu, the optimization is performed on the deviance scale,D(mathbfu)=-2sum_i=1^q log(f_i(u_i))The objective, D, consists of two parts: the sum of the (squared) deviance residuals, measuring fidelity to the data, and the squared length of mathbfu, which is the penalty. In the PIRLS algorithm, only the sum of these components is needed. To use Gauss-Hermite quadrature the contributions of each of the u_ii=1dotsq should be separately evaluated.julia> const devc0 = map!(abs2, m1.devc0, m1.u[1]);  # start with uᵢ²\nError: UndefVarError: m1 not defined\n\njulia> const devresid = m1.resp.devresid;   # n-dimensional vector of deviance residuals\nError: UndefVarError: m1 not defined\n\njulia> const refs = first(m1.LMM.reterms).refs;  # n-dimensional vector of indices in 1:q\nError: UndefVarError: m1 not defined\n\njulia> for (dr, i) in zip(devresid, refs)\n    devc0[i] += dr\nend\nError: UndefVarError: devresid not defined\n\njulia> show(devc0)\nError: UndefVarError: devc0 not defined\nOne thing to notice is that, even on the deviance scale, the contributions of different districts can be of different magnitudes. This is primarily due to different sample sizes in the different districts.julia> using FreqTables\n\njulia> freqtable(contra, :d)\'\nError: UndefVarError: contra not defined\nBecause the first district has one of the largest sample sizes and the third district has the smallest sample size, these two will be used for illustration. For a range of u values, evaluate the individual components of the deviance and store them in a matrix.const devc = m1.devc;Error: UndefVarError: m1 not definedconst xvals = -5.0:2.0^(-4):5.0;\nconst uv = vec(m1.u[1]);Error: UndefVarError: m1 not definedconst u₀ = vec(m1.u₀[1]);Error: UndefVarError: m1 not definedresults = zeros(length(devc0), length(xvals))Error: UndefVarError: devc0 not definedfor (j, u) in enumerate(xvals)\n    fill!(devc, abs2(u))\n    fill!(uv, u)\n    MixedModels.updateη!(m1)\n    for (dr, i) in zip(devresid, refs)\n        devc[i] += dr\n    end\n    copyto!(view(results, :, j), devc)\nendError: UndefVarError: devc not definedA plot of the deviance contribution versus u_1Error: UndefVarError: results not definedshows that the deviance contribution is very close to a quadratic. This is also true for u_3Error: UndefVarError: results not definedThe PIRLS algorithm provides the locations of the minima of these scalar functions, stored asjulia> m1.u₀[1]\nError: UndefVarError: m1 not defined\nthe minima themselves, evaluated as devc0 above, and a horizontal scale, which is the inverse of diagonal of the Cholesky factor. As shown below, this is an estimate of the conditional standard deviations of the components of mathcalU.julia> const s = inv.(m1.LMM.L[Block(1,1)].diag);\nError: UndefVarError: m1 not defined\n\njulia> s\'\nError: UndefVarError: s not defined\nThe curves can be put on a common scale, corresponding to the standard normal, asjulia> for (j, z) in enumerate(xvals)\n    @. uv = u₀ + z * s\n    MixedModels.updateη!(m1)\n    @. devc = abs2(uv) - devc0\n    for (dr, i) in zip(devresid, refs)\n        devc[i] += dr\n    end\n    copyto!(view(results, :, j), devc)\nend\nError: UndefVarError: s not defined\nError: UndefVarError: results not definedError: UndefVarError: results not definedOn the original density scale these becomejulia> for (j, z) in enumerate(xvals)\n    @. uv = u₀ + z * s\n    MixedModels.updateη!(m1)\n    @. devc = abs2(uv) - devc0\n    for (dr, i) in zip(devresid, refs)\n        devc[i] += dr\n    end\n    copyto!(view(results, :, j), @. exp(-devc/2))\nend\nError: UndefVarError: s not defined\nError: UndefVarError: results not definedError: UndefVarError: results not definedand the function to be integrated with the normalized Gauss-Hermite rule isjulia> for (j, z) in enumerate(xvals)\n    @. uv = u₀ + z * s\n    MixedModels.updateη!(m1)\n    @. devc = abs2(uv) - devc0\n    for (dr, i) in zip(devresid, refs)\n        devc[i] += dr\n    end\n    copyto!(view(results, :, j), @. exp((abs2(z) - devc)/2))\nend\nError: UndefVarError: s not defined\nError: UndefVarError: results not definedError: UndefVarError: results not defined"
},

{
    "location": "bootstrap/#",
    "page": "Parametric bootstrap for linear mixed-effects models",
    "title": "Parametric bootstrap for linear mixed-effects models",
    "category": "page",
    "text": ""
},

{
    "location": "bootstrap/#MixedModels.parametricbootstrap",
    "page": "Parametric bootstrap for linear mixed-effects models",
    "title": "MixedModels.parametricbootstrap",
    "category": "function",
    "text": "parametricbootstrap(rng::AbstractRNG, nsamp::Integer, m::LinearMixedModel;\n    β = m.β, σ = m.σ, θ = m.θ)\nparametricbootstrap(nsamp::Integer, m::LinearMixedModel;\n    β = m.β, σ = m.σ, θ = m.θ)\n\nPerform nsamp parametric bootstrap replication fits of m, returning a MixedModelBootstrap.\n\nThe default random number generator is Random.GLOBAL_RNG.\n\nNamed Arguments\n\nβ, σ, and θ are the values of m\'s parameters for simulating the responses.\n\n\n\n\n\n"
},

{
    "location": "bootstrap/#Parametric-bootstrap-for-linear-mixed-effects-models-1",
    "page": "Parametric bootstrap for linear mixed-effects models",
    "title": "Parametric bootstrap for linear mixed-effects models",
    "category": "section",
    "text": "Julia is well-suited to implementing bootstrapping and other simulation-based methods for statistical models. The parametricbootstrap function in the MixedModels package provides an efficient parametric bootstrap for linear mixed-effects models.parametricbootstrap"
},

{
    "location": "bootstrap/#The-parametric-bootstrap-1",
    "page": "Parametric bootstrap for linear mixed-effects models",
    "title": "The parametric bootstrap",
    "category": "section",
    "text": "Bootstrapping is a family of procedures for generating sample values of a statistic, allowing for visualization of the distribution of the statistic or for inference from this sample of values.A parametric bootstrap is used with a parametric model, m, that has been fit to data. The procedure is to simulate n response vectors from m using the estimated parameter values and refit m to these responses in turn, accumulating the statistics of interest at each iteration.The parameters of a LinearMixedModel object are the fixed-effects parameters, β, the standard deviation, σ, of the per-observation noise, and the covariance parameter, θ, that defines the variance-covariance matrices of the random effects.For example, a simple linear mixed-effects model for the Dyestuff data in the lme4 package for R is fit byjulia> using DataFrames, Gadfly, MixedModels, Random, RData\n\njulia> datf = joinpath(dirname(pathof(MixedModels)), \"..\", \"test\", \"dat.rda\");\n\njulia> const dat = Dict(Symbol(k)=>v for (k,v) in load(datf));\nError: SystemError: opening file \"/home/travis/build/JuliaStats/MixedModels.jl/src/../test/dat.rda\": No such file or directory\njulia> ds = rename!(dat[:Dyestuff], [:Batch, :Yield])  # the Dyestuff data\nError: UndefVarError: dat not defined\n\njulia> m1 = fit(MixedModel, @formula(Yield ~ 1 + (1 | Batch)), ds)\nError: UndefVarError: ds not defined\nTo bootstrap the model parameters, first initialize a random number generatorjulia> const rng = MersenneTwister(1234321);\nthen create a bootstrap samplejulia> samp = parametricbootstrap(rng, 10_000, m1);\nError: UndefVarError: m1 not defined\n\njulia> propertynames(samp)\nError: UndefVarError: samp not defined\nAs shown above, the sample has several named properties, which allow for convenient extraction of information.  For example, a density plot of the estimates of σ, the residual standard deviation, can be created asplot(x=samp.σ, Geom.density, Guide.xlabel(\"Parametric bootstrap estimates of σ\"))Error: UndefVarError: samp not definedFor the estimates of the intercept parameter, the getproperty extractor must be usedplot(x = first.(samp.β), Geom.density, Guide.xlabel(\"Parametric bootstrap estimates of β₁\"))Error: UndefVarError: samp not definedThe σs property contains the estimates of the standard deviation of the random effects in a hierarchical format.julia> typeof(samp.σs)\nError: UndefVarError: samp not defined\nThis is to allow for random effects associated with more than one grouping factor. If we only have one grouping factor for random effects, which is the case here, we can use the first extractor, as injulia> first(samp.σs)\nError: UndefVarError: samp not defined\nor, to carry this one step further,plot(x=first.(first(samp.σs)), Geom.density,\n    Guide.xlabel(\"Parametric bootstrap estimates of σ₁\"))Error: UndefVarError: samp not definedNotice that this density plot has a spike, or mode, at zero. Although this mode appears to be diffuse, this is an artifact of the way that density plots are created. In fact, it is a pulse, as can be seen from a histogram.plot(x=first.(first(samp.σs)), Geom.histogram,\n    Guide.xlabel(\"Parametric bootstrap estimates of σ₁\"))Error: UndefVarError: samp not definedA value of zero for the standard deviation of the random effects is an example of a singular covariance. It is easy to detect the singularity in the case of a scalar random-effects term. However, it is not as straightforward to detect singularity in vector-valued random-effects terms.For example, if we bootstrap a model fit to the sleepstudy datajulia> m2 = fit(MixedModel, @formula(Y ~ 1+U+(1+U|G)), dat[:sleepstudy])\nError: UndefVarError: dat not defined\njulia> samp2 = parametricbootstrap(rng, 10_000, m2);\nError: UndefVarError: m2 not defined\nthe singularity can be exhibited as a standard deviation of zero or as a correlation of pm1. The σρs property of the sample is a vector of named tuplesjulia> σρ = first(samp2.σρs);\nError: UndefVarError: samp2 not defined\n\njulia> typeof(σρ)\nError: UndefVarError: σρ not defined\nwhere the first element of the tuple is itself a tuple of standard deviations and the second (also the last) element of the tuple is the correlation.A histogram of the estimated correlations from the bootstrap sample has a spike at +1.ρs = first.(last.(σρ))Error: UndefVarError: σρ not definedplot(x = ρs, Geom.histogram,\n    Guide.xlabel(\"Parametric bootstrap samples of correlation of random effects\"))Error: UndefVarError: ρs not definedor, as a count,julia> sum(ρs .≈ 1)\nError: UndefVarError: ρs not defined\nClose examination of the histogram shows a few values of -1.julia> sum(ρs .≈ -1)\nError: UndefVarError: ρs not defined\nFurthermore there are even a few cases where the estimate of the standard deviation of the random effect for the intercept is zero.julia> sum(first.(first.(first.(σρ))) .≈ 0)\nError: UndefVarError: σρ not defined\nThere is a general condition to check for singularity of an estimated covariance matrix or matrices in a bootstrap sample. The parameter optimized in the estimation is θ, the relative covariance parameter. Some of the elements of this parameter vector must be non-negative and, when one of these components is approximately zero, one of the covariance matrices will be singular.The boundary values are available asjulia> samp2.m.optsum.lowerbd\nError: UndefVarError: samp2 not defined\nso the check on singularity becomesjulia> sum(θ -> any(θ .≈ samp2.m.optsum.lowerbd), samp2.θ)\nError: UndefVarError: samp2 not defined\nThe issingular method for a LinearMixedModel object that tests if a parameter vector θ corresponds to a boundary or singular fit. The default value of θ is that from the model but another value can be given as the second argument.This operation is encapsulated in a method for the issingular function.julia> sum(issingular(samp2))\nError: UndefVarError: samp2 not defined\n"
},

]}
