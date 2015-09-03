compLionMort
============

Comparison of age- and sex-specific mortality among lion populations

We implement a Bayesian hierachical framework to infer age-specific survival 
of both sexes using data from two populations of African Lions. One population
thrives almost undisturbed by humans, whereas the other population should
have high levels of anthropogenic mortality.

Age-specific survival of lion males has never been published because 
male lions disperse at maturity and are therefore commonly lost for data 
collection. In order to study age- and sex-specific survival, we therefore had to develop
a method to infer survival from lifespan data that have a strong sex-bias in right-censored 
records. Our model takes into account sex- and age-specific probabilities of dispersal when 
proposing Siler mortality parameters.

To run the model download runsModel.R, fcts01.R, and datSH.txt. The script runsModel.R sourses fcts01.R and reads in the data from datSH.txt.



