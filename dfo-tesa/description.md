# Spatial and spatiotemporal models with sdmTMB for estimating species distribution and survey index standardization

In this workshop, we will introduce fitting spatial and spatiotemporal
generalized linear mixed effects models (GLMMs) with Gaussian random fields.
Such models are increasingly used for species distribution modelling and survey
index standardization—particularly to control for survey protocol
inconsistencies or to incorporate climate variables. Such models are also
sometimes referred to by fancy terms such as "geostatistical", spatial "latent
variables", "spatially explicit", "predictive-process", and "SPDE". This
workshop will introduce this class of models and focus on fitting and
understanding them with the sdmTMB R package—a package that has a flexible
and user-friendly interface that should be familiar to users of packages such
as mgcv, lme4, or glmmTMB.

Topics will include:
* An introduction to Gaussian random fields
* An introduction to the sdmTMB R package
* Spatial random field GLMMs
* Spatiotemporal random field GLMMs
* Two-part delta or hurdle models for zero inflation
* Time-varying and spatially varying coefficients
* Forecasting
* Index standardization
* Model checking and comparison
* Depending on attendee interest as gauged before the workshop, we can touch on
  or substitute other advanced topics including simulation, anisotropy, center
  of gravity, physical barriers to spatial correlation, 'stitching' nearby
  comparable surveys, fitting the models with Stan for Bayesian inference, or
  modelling presence-only data.

Attendees should have an intermediate knowledge of R and some experience with
GLMs, GAMs, and mixed effects models (e.g., using glm(), mgcv, lme4, or
glmmTMB). While all examples will use sdmTMB, similarly structured models can
be fit with INLA or VAST and so the concepts will also be applicable to
participants using these packages. The workshop will consist of presentations,
tutorials, and group exercises with the support of instructors.

Preprint describing sdmTMB: https://doi.org/10.1101/2022.03.24.485545
sdmTMB documentation site: https://pbs-assess.github.io/sdmTMB
