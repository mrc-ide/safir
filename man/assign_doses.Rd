% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution_vaccines.R
\name{assign_doses}
\alias{assign_doses}
\title{Assign doses to eligible persons based on available supply (multi-dose, no types)}
\usage{
assign_doses(doses_left, events, dose, eligible, parameters)
}
\arguments{
\item{doses_left}{the number of vaccine doses available}

\item{events}{a list of model events}

\item{dose}{what dose to assign for}

\item{eligible}{output of \code{\link{target_pop}}}

\item{parameters}{a list of model parameters}
}
\description{
Given a list of bitsets indicating who is eligible for this dose from \code{\link{target_pop}},
and the number of doses available today \code{doses_left}, schedule individuals for doses.
If fewer doses are available than eligible persons, the doses are apportioned
proportional to age group size.
}
