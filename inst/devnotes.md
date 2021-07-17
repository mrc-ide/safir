# Notes

If we need to include onward infectiousness as something dependent upon ab titre
for breakthrough infections, it will need to interact with this: 
```
# Group infection by age
ages <- variables$discrete_age$get_values(infectious)
inf_ages <- tab_bins(a = ages, nbins = parameters$N_age)
```
If it can be turned into a vector of values of relative infectiousness the same
length as `ages` then we just need to make a modified `tab_bins` that handles that
relative multiplicative contribution of each individual.
