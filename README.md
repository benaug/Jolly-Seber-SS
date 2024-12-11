This repository contains Jolly Seber models fit using the sufficient statistics instead of using the complete data likelihood. This approach is much faster, but can only be used if there are no individual effects on any parameters.

This repository contains M0 and Mt for Jolly Seber with or without removals (same code handles all cases). 
The model is the model from Schofield et al. (2009), section 2.4 without density dependent effects used in their example with two differences:
1. I remove underflow in the marked individual observation model
2. I modified the binomial size parameter for the unmarked survivor distribution to properly account for unmarked removals (see model file comments).


https://link.springer.com/article/10.1007/s10651-007-0069-1
@article{schofield2009flexible,
  title={Flexible hierarchical mark-recapture modeling for open populations using WinBUGS},
  author={Schofield, Matthew R and Barker, Richard J and MacKenzie, Darryl I},
  journal={Environmental and Ecological Statistics},
  volume={16},
  pages={369--387},
  year={2009},
  publisher={Springer}
}