# Our testing strategy

Testing stochastic code is difficult. There are loads of execution paths to check and recreating steering tests down them is not always easy.

To get robust tests which cover the critical execution paths and states, we have decided to combine unit-tests (with mocks) and integration tests.

## Unit-tests

These are the smallest, most common and easiest tests to write. Each unit-test is designed to test *one* aspect of *one* function.

We want to use them for:

* processes
* event listeners

### Why?

*Flagging regressions*. As safir gets bigger, it will become harder to change things without breaking something else. For example, it can feel scary to update the R version, package dependencies (like squire or individual), or simply re-write a function.

Good unit-tests should show you which functions break and how. So we can improve the codebase more confidently.

*Recreating bugs*. When we notice something is wrong, a unit-test can recreate the issue with one click. So we can spend more time developing a fix than recreating the issue with ad-hoc scripts.

Unit-tests also tell you if you've broken something else with your fix!

*Explaining usage*. When we come across a mysterious function, the unit-test will show the kinds of inputs and outputs we expect.

*Checking edge-cases*. Mocking gives us the flexibility to check weird execution paths that integration tests struggle to check. In a few lines we can check that a function will not explode with an empty vector, or with integers vs reals.

*Determinism*. We can make our functions deterministic by mocking trusted random number generators. e.g. `stats::runif`

### Example

Here's an example unit-test for the initialisation process (courtesy of OJ):

```r
test_that("test create_setup_process", {
   # create our dummy variables for create_setup_process
   # these are mock arguments to be passed to the function
   # Because create_steup_process is istelf a wrapper to a function that takes
   # an argument api, we will mock the arguments passed to the api
   individuals <- list(
      human = mockery::mock()
   )
   states <- list(
      E = mockery::mock()
   )
   variables <- list(discrete_age = mockery::mock())
   events <- list(mild_infection = mockery::mock(),
                  severe_infection = mockery::mock())
   process <- create_setup_process(
      individuals,
      states,
      events,
      variables
   )
   # here we know mock what our api will be doing
   # To test that our function processes both infection outcomes (severe and mild)
   # let's make sure that our mocked api will test this

   # We create our api list. This must have mocked named elements that match
   # each function required by api in the create_setup_process function:
   api <- list(
      # we need a get_state function to return the indices of the individuals in
      # state E. Here let's pretend that individuals 5:8 are in E
      get_state = mockery::mock(
         c(5, 6, 7, 8), # exposed states
      ),
      # we also need to get their discrete ages. Let's pretend that they have
      # discrete ages 1, 2, 3, 5
      get_variable = mockery::mock(
         c(1, 2, 3, 5), # age of our individuals
      ),
      # we also need get_parameters. This should return in create_setup_process
      # a list that can be indexed for both prob_hosp and dur_E. To make the outcome
      # of the bernoulli predictable let's have the probs for individuals 5 and 6
      # (so ages 1, 2) be 0 and then for individuals 7, 8 (so ages 3, 5) be 1
      get_parameters = mockery::mock(list(
         # note here 0.5 is for an individual age 4 but we don't have anyone age 4
         # so this wouldn't impact the predictability of the outcome
         prob_hosp = c(0, 0, 1, 0.5, 1),
         # And we aso need the duration dur_E that is used in the erlang draws.
         # We could mock it, or we have dur_E be equal to 0 to contol the outcome
         dur_E = 0 # set to zero to control behaviour of erlang without mocking it
      )),
      schedule = mockery::mock()
   )

   # now that we have our api let's call our process function
   process(api)

   # The process should have called api$schedule twice. The first call is to
   # schedule the severe infection for individuals 7,8 who schedule in 0 days
   mockery::expect_args(
      api$schedule,
      1,
      event = events$severe_infection,
      target = c(7, 8),
      delay = c(0,0)
   )

   # The process should have called api$schedule twice. The second call is to
   # schedule the severe infection for individuals 5,6 who schedule in 0 days
   mockery::expect_args(
      api$schedule,
      2,
      event = events$mild_infection,
      target = c(5, 6),
      delay = c(0, 0)
   )
})
```

## Integration tests

There should be relatively few of these. They check that everything is wired up and interacts as expected. Each test can involve executing several functions for one test.

We want to use them for high-level functions like:

 * create_processes (collects all the processes and wires them up to the state)
 * create_event_based_processes (collects events and listeners and wires them up)
 * run_simulation (collects all the model components, parameterises them and runs the model)

### Why?

*Large refactors*. If we want to re-wire large parts of the code base, these tests are going to tell us which components are affected.

*Low power*. It's tricky to steer an integration test down a particular execution path. You have to set up the whole state to do it - without mocks!

So we should only create them for simple, high-level integration functions. We should try to push logic down into smaller, unit-testable functions.

### Example

Here's an example integration test for create_event_based_processes

```r
test_that("create_event_based_processes assigns a listener to each event", {

  # Set up the whole simulation state
  pop <- get_population("ATG")
  pop$n <- as.integer(pop$n/100) # reduce population size for speed
  parameters <- get_parameters("ATG")
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(pop, parameters)
  human <- create_human(states, variables, events)

  # Integrate the events and the listeners
  create_event_based_processes(
    human,
    states,
    variables,
    events,
    parameters
  )

  # Check that all events have been assigned a listener
  for (event in events) {
    expect_gt(length(event$listeners), 0)
  }
})
```

As we introduce new testing strategies, we can document them here.
