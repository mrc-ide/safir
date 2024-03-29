% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/events_vaccination.R
\name{create_vaccination_dose_listener}
\alias{create_vaccination_dose_listener}
\title{Create listener for vaccination dose (multi-dose, no types)}
\usage{
create_vaccination_dose_listener(variables, parameters, dose)
}
\arguments{
\item{variables}{a named list of variables}

\item{parameters}{the parameters}

\item{dose}{integer dose}
}
\description{
Updates state when a vaccine dose is given. It does not schedule future events.
}
