Smooth Empirical Likelihood with Shares
=======================================

R code for A. Cosma, A. V. Kostyrka, G. Tripathi (2018) “Inference in Conditional Moment Restriction Models When There is Selection due to Stratification”.

University of Luxembourg, Faculty of Law, Economics and Finance (FDEF), Centre for Research in Economics and Management (CREA).

## File tree
* _functions.R_ — functions provided by the authors that are used in all simulations (data generation, simulation, estimation, testing).
* _OLS-GMM-strata-all-cases.R_ — complete R script to replicate **all** OLS- and GMM-based results (benchmark models). Outputs _runlogs.txt_.
* _SEL-strata-one-case.R_ — complete R script to replicate **one** selected SEL-based result. Outputs an _.RData_ file with name containing parameters used for simulation.
* _runlogs.txt_ — authors’ output from _OLS-GMM-strata-all-cases.R_ that was used to generate tables in the paper.
* _EL/scel.R_ — Art Owen’s (2015-03) original code from http://statweb.stanford.edu/~owen/empirical/.
* _EL/scelcount.R_ — Art Owen’s (2017-02) original code from http://statweb.stanford.edu/~owen/empirical/.
* _EL/scelcount2.R_ — authors’ modification of _EL/scelcount2.R_ allowing for an unconditional moment restriction thanks to an additional argument—a shift.

## How does one run the code?

It is simple! Just download the entire repository or clone it, and then used your favourite R IDE or console to launch the scripts!

By default, it will occupy all available cores on a machine (unless it is run on Windows); one might want to change the number of allocated cores manually!

For console invocation on one machine:

`Rscript OLS-GMM-strata-all-cases.R`

`Rscript SEL-strata-one-case.R`

Or within R:

`> source("OLS-GMM-strata-all-cases.R")`

`> source("SEL-strata-one-case.R")`

For console invocation on an MPI cluster (in this example, using the node file provided by the OAR scheduler):

1. Modify line 10 of _SEL-strata-one-case.R_ to `multimachine <- TRUE`;
2. `mpirun -np 1 -machinefile $OAR_NODE_FILE Rscript SEL-strata-one-case.R`

Warning: OLS and GMM perform extremely fast even on low-end machines. On the other hand, SEL is computationally more intensive, and even moreso with aggregate share estimation (due to an extra nested optimisation routine), so the authors used a cluster to run all SEL-based scripts. However, even if you do not have access to an MPI cluster, by default, it is assumed that you are running your script locally on your multi-core machine. Even if you are using Windows, infamous for its lack of parallel support in R (at the time of writing), it will correctly use only one core. By default, tracing information is on, so the output window might fill up with messages from the optimiser. This is done in order to allow one to get a ballpark estimate of necessary time (on average, with default initial values, the derivative-free Nelder—Mead optimiser converges after 200 evaluations or so); therefore, one can measure the time between two function evaluations (based on the output refresh rate) and multiply it by 200 to obtain the right order of magnitude for the job time.

Approximate times for one Monte-Carlo simulation (the authors used 1000 in their paper):

* Without aggrgate share estimation (`shares <- FALSE` on line 7 of _SEL-strata-one-case.R_): 0.036, 0.129, 0.523 minutes for N=50, 150, 500;
* With aggrgate share estimation (`shares <- TRUE` on line 7 of _SEL-strata-one-case.R_): 16.17, 45.70, 149.4 minutes for N=50, 150, 500.

## Improvements and suggestions?

You can find the authors’ email addresses on the website of the University of Luxembourg: https://uni.lu (on their personal pages; links valid as of 2018-05-15):

* https://wwwen.uni.lu/recherche/fdef/crea/people/antonio_cosma (corresponding author)
* https://wwwen.uni.lu/recherche/fdef/crea/people/andrei_kostyrka (code repository maintainer)
* https://wwwen.uni.lu/recherche/fdef/crea/people/gautam_tripathi
