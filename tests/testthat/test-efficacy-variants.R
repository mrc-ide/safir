test_that("vaccine efficacy for identical variants returns same as default (one variant) model", {

  time_period <- 1000
  pop <- get_population("AFG")

  psq <- get_parameters(
    iso3c = "AFG",
    population = pop$n,
    R0 = 2,
    time_period = time_period,
    max_age = NULL,
    dt = 1
  )

  vaccine_parameters <- get_vaccine_ab_titre_parameters(vaccine = "Pfizer")
  voc_pars <- get_voc_parameters(safir_parameters = psq, voc_types = c("wt", "alpha"), voc_trajectory = NULL, vaccine_parameters = vaccine_parameters)

  n <- 50
  ab_titre <- log(10^rnorm(n = n, mean = log10(vaccine_parameters$mu_ab[1]),sd = vaccine_parameters$std10))
  ab_titre[sample.int(n = n,size = 20,replace = FALSE)] <- -Inf

  ef_inf <- vaccine_efficacy_infection(ab_titre = ab_titre,parameters = vaccine_parameters)
  ef_inf_voc <- vaccine_efficacy_infection_voc(ab_titre = ab_titre,parameters = voc_pars)
  ef_inf_voc_cpp <- vaccine_efficacy_infection_voc_cpp(ab_titre = ab_titre,parameters = voc_pars)

  expect_equal(ef_inf_voc[1, ], ef_inf)
  expect_equal(ef_inf_voc[2, ], ef_inf)
  expect_equal(ef_inf_voc_cpp[1, ], ef_inf)
  expect_equal(ef_inf_voc_cpp[2, ], ef_inf)

  ef_severe <- vaccine_efficacy_severe(ab_titre = ab_titre, ef_infection = ef_inf, parameters = vaccine_parameters)
  ef_severe_voc <- vaccine_efficacy_severe_voc(ab_titre = ab_titre, ef_infection = ef_inf_voc, parameters = voc_pars)
  ef_severe_voc_cpp <- vaccine_efficacy_severe_voc_cpp(ab_titre = ab_titre, ef_infection = ef_inf, parameters = voc_pars)

  expect_equal(ef_severe_voc[1, ], ef_severe)
  expect_equal(ef_severe_voc[2, ], ef_severe)
  expect_equal(ef_severe_voc_cpp[1, ], ef_severe)
  expect_equal(ef_severe_voc_cpp[2, ], ef_severe)

})
