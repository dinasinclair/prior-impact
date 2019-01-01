# prior-impact
Using Bayesian modeling to better prioritize and access impact evaluations.

This set of files examines the scenario: imagine you have *I* cities in which you can implement a program, and you have
priors on how well each city will perform. If you can do pilot studies to update those priors on a limited subset *K* of the 
possible cities, which cities should you perform the pilot in?

To decide which cities to pilot, we use a bayesian hierarchical model coded in R and Stan. In the simplest (1D) case, the
cities are not grouped in any way, and we use a one-layer random effects model. In the more relevant (2D) case, the 
cities are grouped together through some predetermined underlying structure, and we use a two-layer random effects model.

The repo contains the following files/folders:
* CitiesMainCode.R files, which run the full analysis from data generation to the question "how frequently do we change our minds if
we run the pilot on this city subset"
* RandomEffectsModel.stan files, which contain the underlying random effects model bayesian hierarchical code called 
in CitiesMainCode.R files.
* EvaluatingConvergence.R files, which help visualize the results of the Stan REM models and determine if the models are converging 
properly.
* explanatory_code contains R markdown files and their corresponding pdfs walking through the code and math behind some of the main files.
* single_layer_model_results contains R markdown files with graphs of results using the simpler 1D model.
