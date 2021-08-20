# **safir**

<!-- badges: start -->

[![R build
status](https://github.com/mrc-ide/safir/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/safir/actions)
[![CodeFactor](https://www.codefactor.io/repository/github/mrc-ide/safir/badge)](https://www.codefactor.io/repository/github/mrc-ide/safir)
[![codecov.io](https://codecov.io/github/mrc-ide/safir/coverage.svg?branch=main)](https://codecov.io/github/mrc-ide/safir?branch=main)
<!-- badges: end -->

**safir**: **s**quire **a**nd **f**riends **i**ndividual **r**ewrite
(and maintains the squire naming
[theme](https://en.wikipedia.org/wiki/Knights_of_the_Round_Table#Safir)).

It uses the [`{individual}`](https://github.com/mrc-ide/individual)
package to specify and run the simulation.

**safir** can run individual based stochastic versions of the **squire** and
**nimue** models, which recover trajectories from those aggregated population
models (see package vignettes). It also implements a model for vaccination
with an arbitrary number of doses and antibody titre dynamics which affect
protective efficacy as well as efficacy against severe disease outcomes. 

## Installation

``` r
install_github('mrc-ide/safir')
library(safir)
```

## Documentation

[`{safir}`](https://github.com/mrc-ide/safir) is documented on a
[dedicated website](https://mrc-ide.github.io/safir).

This includes the following vignettes:

-   **`Squire Validation Run`**: comparison of
    [`{safir}`](https://github.com/mrc-ide/safir) to
    [`{squire}`](https://github.com/mrc-ide/squire).
-   **`Nimue Validation Run`**: comparison of
    [`{safir}`](https://github.com/mrc-ide/safir) to
    [`{nimue}`](https://github.com/mrc-ide/nimue).
-   **`Vaccine model with multiple doses`**: runs the model for
    explicit tracking of antibody titre for a 2 dose vaccine.
-   **`Code Design`**: a diary of ideas, notes, etc. as
    [`{safir}`](https://github.com/mrc-ide/safir) is developed.

## License

MIT © Imperial College of Science, Technology and Medicine
