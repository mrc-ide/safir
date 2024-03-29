% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{make_rerlang}
\alias{make_rerlang}
\title{Make Erlang waiting time distribution}
\usage{
make_rerlang(mu, dt, shape = 2, shift = 0L)
}
\arguments{
\item{mu}{Mean duration}

\item{dt}{Size of time step}

\item{shape}{Shape parameter of Erlang distribution}

\item{shift}{number of time steps to add to sampled value}
}
\description{
Random draws from Erlang distribution
scaled by time step size.
}
