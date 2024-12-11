This repository contains Jolly Seber models fit using the sufficient statistics instead of using the complete data likelihood. This approach is much faster, but can only be used if there are no individual effects on any parameters.

This repository contains M0 and Mt for Jolly Seber with or without removals (same code handles all cases). 
The model is the model from Schofield et al. (2009), section 2.4 without density dependent effects used in their example with two differences:
1. I remove underflow in the marked individual observation model
2. I modified the binomial size parameter for the unmarked survivor distribution to properly account for unmarked removals (see model file comments).
