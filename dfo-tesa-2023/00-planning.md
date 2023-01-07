9-11:30 - 1 hour lecture - 10 minute break - 50 minute break - 30 minutes

40 minute exercises? 20 minutes after

up to 1:15 lectures, 10 minutes,

# Day 1

### Lectures:

-   Intro to RFs
-   Intro to sdmTMB
-   Spatial example

First day 1:15...

10 minute break

exercises 40-45

### Exercises:

-   Matern correlation
-   Simulating spatial RFs
-   Simulating spatiotemporal RFs
-   Fitting a basic spatial model (plus print(), sanity(), tidy())

# Day 2

### Lectures:

-   Recap
-   Spatiotemporal example extended from before
-   Basic diagnostics
    -   reading print(), sanity(), QQ plots
    -   parameter estimates make sense - e.g., magnitude of spatial SD and spatiotemporal SD
    -   assumptions we make here (like constant spatiotemporal SD, shared or not shared range...)
-   Anistropy!!
-   Get into AR(1), RW, extra_time, ...
-   Comparing models - cut stacking... this time

### Exercises:

-   Fitting a spatiotemporal model
-   Trying different field structures
-   Inspecting output from that model in many ways
    -   sanity(), residuals(), simulate() [e.g., proportion zeros]

# Day 3

### Lectures:

-   Recap
-   Families (10 minutes)
-   Time-varying (15 min)
-   Spatial varying example... (20 min)
-   Calculating an index with that model: (20 min) TODO: expand a bit? COG? ADD A FAIR BIT

### Exercises:

-   Spatially varying trend and plot it...
-   Create an index with delta_gamma and once with Tweedie

# Day 4

-   Recap
-   Extrapolation, interpolation, and forecasting
-   Barrier models
-   Presence only?
-   Grab bag of other stuff you can do we haven't covered (add stacking)
-   Priors (especially PC prior); touch on other functionality like Bayesian?

### Exercises:

-   From beginning to end with new dataset... and minimal hand holding
    -   find something say where depth is too big

    -   use another synoptic trawl dataset WCVI

    -   prompt: try anisotropy here plot_anisotropy()

    -   try share_range = FALSE
-   A series of models with problems: find those problems (e.g., with sanity(), print(), tidy(), and/or residuals()), suggest how they might be solved.

# TODO

-   Add sanity() in lectures and lessons

-   Do we cover residuals() well? ...

-   Turn 05-extra (words of wisdom) into trouble shooting, diagnostics for day 4

-   Installation instructions before

-   Divide people into groups and decide on group leaders

-   Quarto work on DFO?

-   Try INLA install on DFO

-   Anisotropy! and the new plot_anistropy()

-   Google doc... open to anyone

-   Delta models! - add

    -   universe of families

    -   when use various ones

    -   delta models

need slides on families and delta models for day 3

add anisotropy() probably in day 2 slides

add residuals stuff probably day 2?

day 3: show how can interact with delta version of day 2?
