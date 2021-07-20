# Notes

1. If we need to include onward infectiousness as something dependent upon ab titre
for breakthrough infections, it will need to interact with this: 
```
# Group infection by age
ages <- variables$discrete_age$get_values(infectious)
inf_ages <- tab_bins(a = ages, nbins = parameters$N_age)
```
If it can be turned into a vector of values of relative infectiousness the same
length as `ages` then we just need to make a modified `tab_bins` that handles that
relative multiplicative contribution of each individual.

2. One thing to note is that when a large amount of doses are available from the
start of the simulation there may be wasted doses during the lag between people
getting their first dose, and when persons can get (become eligible for) their second dose.
The prioritization step will only increment when both are fulfilled.

3. To test if 100% effective vaccines are having an effect, just set `ab_50 = 2e-12`
or some small value, which ensures that for almost any positive Ab titre the efficacy
will be nearly 1.
