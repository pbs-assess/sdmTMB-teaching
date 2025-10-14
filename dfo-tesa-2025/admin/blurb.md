## Spatial and spatiotemporal modelling with sdmTMB

This course will cover introductory to advanced uses of the sdmTMB R package. sdmTMB is an R package for fitting spatial and spatiotemporal generalized linear mixed effect models (GLMMs) using TMB (Template Model Builder) and the SPDE (Stochastic Partial Differential Equation) approach to approximating Gaussian random fields with Gaussian Markov random fields. A common application is spatially explicit species distribution modeling (SDM). Other fisheries applications include index standardization, combining survey data, tracking species redistribution, and simulating data to evaluate survey designs. This course will consist of 2.5 days of instruction and exercises, plus a half day featuring talks on recent applications of sdmTMB and opportunities for participants to present and discuss their own use cases or modeling challenges. 

Topics are expected to include:

- An introduction to Gaussian random fields, Gaussian Markov random fields, and the SPDE approach
- Spatial GLMMs
- Spatiotemporal GLMMs
- SPDE mesh design
- Standard delta or hurdle models for zero inflation
- Poisson-link delta models
- Penalized smoothers in sdmTMB
- Time-varying and spatially varying coefficient models
- Forecasting
- Index standardization and calculation of distribution metrics (e.g., center of gravity)
- Model checking and comparison (e.g., residuals, marginal and conditional AIC, deviance explained, self-simulation testing)
- Simulating new or conditioned datasets with sdmTMB (e.g., to test survey designs or model estimation performance)
- Anisotropy and physical barriers to correlation
- Priors and penalized complexity priors
- Integrating survey data with sdmTMB
- Approaches to using sdmTMB to model multivariate data (e.g., length or age composition or multiple species) and when it makes sense to use tinyVAST

Attendees should have an intermediate knowledge of R and some experience with GLMs and mixed effects models (e.g., using glm(), mgcv, lme4, or glmmTMB). Prior experience with sdmTMB is not needed, but participants with prior sdmTMB experience who are interested in the more advanced topics are also welcome.

