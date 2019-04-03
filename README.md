# ppta-2019
Code supporting the Liao and Zigler (2019) paper "Posterior Predictive Treatment Assignment Methods for Causal Inference in the Context of Time-Varying Treatments"

## Simulation code

- runIPW_boot.R: Performs IPW, SW, OW and PPTA with bootstrapped SEs (user may toggle command on or off)on all replications of a single dataset (savedata), paralellized to run on multiple cores. 

### Data generation
- data_generation_homo.R: Contains code to simulate datasets with homogeneous treatment effects. Outputs a list with number of elements equal to the number of replications specified. 
- data_generation_hetero.R: Same as data_generation_homo.R, for heterogeneous treatment effects. 
- data_3.R and data_5.R: 200 replications from the same homogeneous treatment effect data-generating mechanism utilized in the simulation section of the paper. 
- data_overlap_3.R and data_overlap_5.R: Same as data_3.R and data_5.R, for heterogeneous treatment effects.

## Application code

- runIPW_boot_app.R: A modification of runIPW_boot.R to run on the application data, 
- simulated_app_data.R: application data frame in the same format as analyzed with runIPW_boot_app.R, but with simulated outcomes and census data. 
