# cbq: An R Package for Conditional Binary Quantile Models

[![Build Status](https://www.r-pkg.org/badges/version/cbq)](https://www.r-pkg.org/badges/version/cbq)
[![Build Status](https://www.r-pkg.org/badges/last-release/cbq)](https://www.r-pkg.org/badges/last-release/cbq)
[![Build Status](http://cranlogs.r-pkg.org/badges/grand-total/cbq)](http://cranlogs.r-pkg.org/badges/grand-total/cbq)

Introduction
------------

The conditional binary quantile (CBQ) models extend the simple version
of binary quantile models for analyzing discrete choices (including but
beyond binary choices) with varying choice alternatives. Each quantile
estimation represents a local inspection of effects at the specified
quantile of the choice probabilities. In the simple binary setting, the
features of the choice alternatives within each choice set are assumed
to be the same. However, in reality we oftentimes observe individuals
who are facing different sets and numbers of choice alterantives. The
CBQ models are developed to solve this problem by introducing a
conditional multinomial structure for modeling varying choice
alternatives. Even though the CBQ models are called "binary", they
actually belong to a more general family of dicrete choice models. I
refer the readers to [Lu (2020)](https://doi.org/10.1017/pan.2019.29)
for the details of the estimation.

This demo page is intended to introduce an R package named `cbq`
associated with the CBQ models. Compared to the original version in [Lu
(2020)](https://doi.org/10.1017/pan.2019.29), the package additionally
enables estimation of mixed effects (including both random and fixed
effects) in the CBQ models and allows for more flexible specifications
that can be determiend by the users of the package. Since the package
has been uploaded to the Comprehensive R Archive Network: [link to the
package in CRAN](https://CRAN.R-project.org/package=cbq), the package is
ready for use after download. To download the package, simply execute in
the R console the following line of code:

    install.packages("cbq")

Or alternatively, users can download a developing version of the package
from my GitHub repositories by using the following codes:

    # Make sure that the following packages have been installed in your local R environment
    if(!require(rstan)) install.packages("rstan")

    # Install cbq from github
    if(!require(devtools)) install.packages("devtools")
    devtools::install_github("xiao-lu-research/cbq")

You can always reach me by xiao.lu.research \[at\] gmail.com whenever
you have questions about the method and the package. In what follows, I
demonstrate step by step how we can use the package. It is important to
note that to reduce the estimation time the datasets used in the
following examples are only subsets of data used in the original studies
and only one chain has been specified to sample the model. Therefore,
the estimates do not necessarily reflect the results of the original
studies and I will refrain from interpreting the results (see [Lu
(2020)](https://doi.org/10.1017/pan.2019.29) for the full estimations).
This, however, should not matter so much for our pedagogical purpose of
showing how we can use the package. Now, let us get started.

[check [here](https://www.xiao-lu.com/files/cbq_demo.html) for the full page of the demo (including plots).]

Package loading and introduction of datasets
--------------------------------------------

To load the package and get help of the pacakge, we can use the
following codes.

    # Load the package
    library(cbq)

    # Get help from R documents
    ?cbq

The datasets used in this demo comes from three published papers in
political science:

\[1\] data\_legislature: Rasmussen, A. (2010). Early Conclusion in
Bicameral Bargaining: Evidence from the Co-decision Legislative
Procedure of the European Union. European Union Politics 12 (1), 41-64.

\[2\] data\_vote: Alvarez, R. M. and J. Nagler (1995). Economics, Issues
and the Perot Candidacy: Voter Choice in the 1992 Presidential Election.
American Journal of Political Science 39 (3), 714-744.

\[3\] data\_coalition: Martin, L. W. and R. T. Stevenson (2001).
Government Formation in Parliamentary Democracies. American Journal of
Political Science 45 (1), 33-50.

Let us first load the data and check the dimensions of the datasets and
the variables included.

    # Legislative data
    load("data_legislature.RData")
    # US vote choice data
    load("data_vote.RData")
    # Coalition formation data
    load("data_coalition.RData")

    # Dimensions of the datasets
    dim(data_legislature)

    ## [1] 385  13

    dim(data_vote)

    ## [1] 2727   41

    dim(data_coalition)

    ## [1] 33256    16

    # Variables included. Note that in all three datasets, the first variable is the dependent variable.
    names(data_legislature)

    ##  [1] "choice"      "var6"        "var3"        "var5"        "var10"      
    ##  [6] "var15"       "urgcy_yes"   "urgcy_noleg" "var4"        "var11"      
    ## [11] "var14"       "var7"        "var17"

    names(data_vote)

    ##  [1] "choice"           "dist"             "respfinp.bush"   
    ##  [4] "natlec.bush"      "respgjob.bush"    "resphlth.bush"   
    ##  [7] "respblk.bush"     "respab.bush"      "east.bush"       
    ## [10] "south.bush"       "west.bush"        "newvoter.bush"   
    ## [13] "termlim.bush"     "deficit.bush"     "dem.bush"        
    ## [16] "rep.bush"         "women.bush"       "educ.bush"       
    ## [19] "age1829.bush"     "age3044.bush"     "age4559.bush"    
    ## [22] "respfinp.clinton" "natlec.clinton"   "respgjob.clinton"
    ## [25] "resphlth.clinton" "respblk.clinton"  "respab.clinton"  
    ## [28] "east.clinton"     "south.clinton"    "west.clinton"    
    ## [31] "newvoter.clinton" "termlim.clinton"  "deficit.clinton" 
    ## [34] "dem.clinton"      "rep.clinton"      "women.clinton"   
    ## [37] "educ.clinton"     "age1829.clinton"  "age3044.clinton" 
    ## [40] "age4559.clinton"  "chid"

    names(data_coalition)

    ##  [1] "realg"    "minor"    "minwin"   "numpar"   "dompar"   "median"  
    ##  [7] "gdiv1"    "mgodiv1"  "prevpm"   "sq"       "mginvest" "anmax2"  
    ## [13] "pro"      "anti"     "country"  "case"

    # Have a look at the first 5 rows of the datasets.
    head(data_legislature,5)

    ##    choice var6 var3 var5 var10 var15 urgcy_yes urgcy_noleg var4 var11
    ## 1       0 0.44 0.44    1     0     1         0           0  130     3
    ## 2       0 0.44 0.44    1     0     1         0           0  130     2
    ## 3       0 0.44 0.44    1     0     1         0           0  130     2
    ## 9       0 0.40 0.40    1     0     1         0           0  130     1
    ## 10      0 0.19 0.00    0     0     1         0           0  130     0
    ##    var14 var7 var17
    ## 1      0    0     1
    ## 2      0    1     1
    ## 3      0    0     1
    ## 9      0    1     1
    ## 10     0    1     1

    head(data_vote,5)

    ##           choice   dist respfinp.bush natlec.bush respgjob.bush
    ## 1.Bush      TRUE 0.1024             1           3             6
    ## 1.Clinton  FALSE 4.0804             0           0             0
    ## 1.Perot    FALSE 0.2601             0           0             0
    ## 2.Bush      TRUE 0.1024             3           5             6
    ## 2.Clinton  FALSE 4.0804             0           0             0
    ##           resphlth.bush respblk.bush respab.bush east.bush south.bush
    ## 1.Bush                7            6           2         0          0
    ## 1.Clinton             0            0           0         0          0
    ## 1.Perot               0            0           0         0          0
    ## 2.Bush                6            6           3         0          1
    ## 2.Clinton             0            0           0         0          0
    ##           west.bush newvoter.bush termlim.bush deficit.bush dem.bush
    ## 1.Bush            0             0            1            0        0
    ## 1.Clinton         0             0            0            0        0
    ## 1.Perot           0             0            0            0        0
    ## 2.Bush            0             0            1            1        0
    ## 2.Clinton         0             0            0            0        0
    ##           rep.bush women.bush educ.bush age1829.bush age3044.bush
    ## 1.Bush           1          1         3            0            1
    ## 1.Clinton        0          0         0            0            0
    ## 1.Perot          0          0         0            0            0
    ## 2.Bush           1          1         4            0            0
    ## 2.Clinton        0          0         0            0            0
    ##           age4559.bush respfinp.clinton natlec.clinton respgjob.clinton
    ## 1.Bush               0                0              0                0
    ## 1.Clinton            0                1              3                6
    ## 1.Perot              0                0              0                0
    ## 2.Bush               0                0              0                0
    ## 2.Clinton            0                3              5                6
    ##           resphlth.clinton respblk.clinton respab.clinton east.clinton
    ## 1.Bush                   0               0              0            0
    ## 1.Clinton                7               6              2            0
    ## 1.Perot                  0               0              0            0
    ## 2.Bush                   0               0              0            0
    ## 2.Clinton                6               6              3            0
    ##           south.clinton west.clinton newvoter.clinton termlim.clinton
    ## 1.Bush                0            0                0               0
    ## 1.Clinton             0            0                0               1
    ## 1.Perot               0            0                0               0
    ## 2.Bush                0            0                0               0
    ## 2.Clinton             1            0                0               1
    ##           deficit.clinton dem.clinton rep.clinton women.clinton
    ## 1.Bush                  0           0           0             0
    ## 1.Clinton               0           0           1             1
    ## 1.Perot                 0           0           0             0
    ## 2.Bush                  0           0           0             0
    ## 2.Clinton               1           0           1             1
    ##           educ.clinton age1829.clinton age3044.clinton age4559.clinton
    ## 1.Bush               0               0               0               0
    ## 1.Clinton            3               0               1               0
    ## 1.Perot              0               0               0               0
    ## 2.Bush               0               0               0               0
    ## 2.Clinton            4               0               0               0
    ##           chid
    ## 1.Bush       1
    ## 1.Clinton    1
    ## 1.Perot      1
    ## 2.Bush       2
    ## 2.Clinton    2

    head(data_coalition,5)

    ##   realg minor minwin numpar dompar median    gdiv1   mgodiv1 prevpm sq
    ## 1     0     1      0      1      0      0 0.000000 0.2308820      0  0
    ## 2     0     1      0      1      0      0 0.000000 0.0485294      0  0
    ## 3     0     1      0      1      0      0 0.000000 0.0485294      0  0
    ## 4     0     1      0      2      0      0 0.822794 0.2308820      0  0
    ## 5     0     1      0      2      0      0 0.279412 0.4227940      0  0
    ##   mginvest anmax2 pro anti country case
    ## 1        0      0   0    0      11    3
    ## 2        0      0   0    0      11    3
    ## 3        0      0   0    0      11    3
    ## 4        0      0   0    0      11    3
    ## 5        0      0   0    0      11    3

The correspondence between the variable names and their meanings is
shown below. In `data_legislature` (unit of analysis: legislative
bills), the variables correspond to

1.  choice: Early agreement (dummy)
2.  var6: Policy distance rapporteur-EP median
3.  var3: Policy distance rapporteur-EP median \* Big party group
4.  var5: Big party group
5.  var10: Country coherence of key negotiators
6.  var15: Closeness of working relationship (calendar year)
7.  urgcy\_yes: Urgency of legislation: existing legislation expires;
8.  urgcy\_noleg: No existing legislation
9.  var4: Workload (pending legislation during first reading)
10. var11: Issue areas consulted (no. of consultative committees)
11. var14: New act
12. var7: Directive
13. var17: Inter-institutional disagreement

In `data_vote` (unit of analysis: individual vote choice):

1.  choice: 1 Bush 2 Clinton 3 Perot
2.  educ: Respondents education
3.  east: Region (East)
4.  south: Region (South)
5.  west: Region (West)
6.  women: Gender (Female)
7.  respfinp: Felt personal finance were worse
8.  natlec: Felt national economy was worse
9.  bclibdis: Ideological Distance
10. gblibdis: Ideological Distance
11. rplibdis: Ideological Distance
12. dem: Democrat
13. rep: Republican
14. respgjob: Oppose government jobs
15. resphlth: Oppose government health care
16. respblk: Oppose government minority assistance
17. respab: Abortion
18. termlim: Term limits
19. age1829: Age: 18-29
20. age3044: Age: 30-44
21. age4559: Age: 45-59
22. newvoter: New or returning voter
23. deficit: Felt deficit was a major problem
24. chid: Choice set indicator

In `data_coalition` (unit of analysis: potential coalitions):

1.  realg: Really formed government  
2.  minor: Minority Coalition
3.  Minwin: Minimal Winning Coalition
4.  numpar: Number of Parties in the Coalition
5.  dompar: Largest Party in the Coalition
6.  median: Median Party in the Coalition
7.  gdiv1: Ideological Divisions in the Coalition
8.  mgodiv1: Ideological Divisions within Majority Opposition
9.  prevpm: Previous Prime Minister in the Coalition
10. sq: Incumbent Coalition
11. mginvest: Minority Coalition where Investiture Vote Required
12. anmax2: Anti-System Presence in the Coalition
13. pro: Pre-Electoral Pact associated with the Coalition
14. anti: Anti-Pact associated with the Coalition
15. case: choice set indicator

Note that all the datasets are in a long-format with the binary-coded (0
or 1) dependent variable and an additional indicator variable showing
the choice sets each row of observation belongs to. Each row in the data
represents a choice alternative.

Estimation of the models
------------------------

### Binary choice: data\_legislature

In a simple binary choice situation, `cbq` only requres specifying the
discrete choice variable on the left and a bunch of independent
variables on the right, seperated by `~`. We can also specify a set of
controling parameters in the function as displayed below.

    # Estimate the model
    model1 <- cbq(choice ~ var6 + var3 + var5 + var10 + var15 + urgcy_yes + urgcy_noleg + var4 + var11 +var14 + var7 + var17, 
                 data = data_legislature, 
                 q = 0.5, # We first specify a 0.5 quantile (median)
                 vi = FALSE, # Indicating whether to use variational inference to approximate the results. The default is FALSE, meaning using MCMC sampling procedures.
                 nsim = 1000, # Running the model for 1000 iterations
                 grad_samples = 1, # Passed to vb (positive integer), the number of samples for Monte Carlo estimate of gradients, defaulting to 1.
                 elbo_samples = 100, # Passed to vb (positive integer), the number of samples for Monte Carlo estimate of ELBO (objective function), defaulting to 100. (ELBO stands for "the evidence lower bound".)
                 tol_rel_obj = 0.01, # Passed to vb (positive double), the convergence tolerance on the relative norm of the objective, defaulting to 0.01.
                 output_samples = 2000, # Passed to vb (positive integer), number of posterior samples to draw and save, defaults to 2000.
                 burnin = NULL, # The number of burnin iterations. If unspecified, hald of the iterations will be burnin iterations.
                 thin = 1, # Thinning parameter.
                 CIsize = 0.95, # The size of confidence interval.
                 nchain = 1, # The number of parallel chains.
                 seeds = 12345, # Random seeds to replicate the results.
                 inverse_distr = FALSE, # If FALSE, the ALD will not be reversed. The default is FALSE.
                 offset = 1e-20, # Offset values to enhance sampling stability. The default value is 1e-20.
                 mc_core = TRUE # Indicating whether the estimation will be run in multiple parallel chains. The default is TRUE.
                 )

    ## 
    ## SAMPLING FOR MODEL 'cbqbv' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.00024 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.4 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 19.4393 seconds (Warm-up)
    ## Chain 1:                5.61687 seconds (Sampling)
    ## Chain 1:                25.0562 seconds (Total)
    ## Chain 1:

We can check the estimated model by using the following three functions,
where the first generates the summary table, the second presents the
coefficient matrix and the third checks model convergence.

    # Print the results
    print(model1)

    ## Conditional binary quantile regression 
    ## 
    ## Call:
    ## cbq(formula = choice ~ var6 + var3 + var5 + var10 + var15 + urgcy_yes + 
    ##     urgcy_noleg + var4 + var11 + var14 + var7 + var17, data = data_legislature, 
    ##     q = 0.5, vi = FALSE, nsim = 1000, grad_samples = 1, elbo_samples = 100, 
    ##     tol_rel_obj = 0.01, output_samples = 2000, burnin = NULL, 
    ##     thin = 1, CIsize = 0.95, nchain = 1, seeds = 12345, inverse_distr = FALSE, 
    ##     offset = 1e-20, mc_core = TRUE)
    ## 
    ## MCMC run for 1000 iterations, with 500 used. 
    ## 
    ## Coefficients:
    ##             Estimate      LB     UB
    ## (Intercept)   -4.819  -7.789 -1.688
    ## var6          -0.107  -1.981  1.512
    ## var3          -5.301 -10.841 -1.132
    ## var5           3.076   0.964  5.789
    ## var10          1.323   0.151  2.560
    ## var15          0.821   0.545  1.151
    ## urgcy_yes      2.752   1.152  4.484
    ## urgcy_noleg    1.048  -0.105  2.353
    ## var4           0.015  -0.001  0.032
    ## var11          0.415   0.140  0.685
    ## var14         -1.803  -2.997 -0.721
    ## var7          -0.619  -1.376  0.191
    ## var17         -2.876  -4.089 -1.800

    print(model1, "mcmc") # Print MCMC

    ## Inference for Stan model: cbqbv.
    ## 1 chains, each with iter=1000; warmup=500; thin=1; 
    ## post-warmup draws per chain=500, total post-warmup draws=500.
    ## 
    ##             mean se_mean   sd    2.5%     25%     50%     75%   97.5%
    ## beta[1]    -4.82    0.10 1.55   -7.79   -5.77   -4.87   -3.84   -1.69
    ## beta[2]    -0.11    0.06 0.91   -1.98   -0.66   -0.08    0.51    1.51
    ## beta[3]    -5.30    0.18 2.45  -10.84   -6.62   -5.29   -3.49   -1.13
    ## beta[4]     3.08    0.08 1.19    0.96    2.21    3.02    3.80    5.79
    ## beta[5]     1.32    0.02 0.65    0.15    0.83    1.32    1.76    2.56
    ## beta[6]     0.82    0.01 0.15    0.54    0.72    0.81    0.93    1.15
    ## beta[7]     2.75    0.04 0.90    1.15    2.10    2.70    3.39    4.48
    ## beta[8]     1.05    0.03 0.62   -0.11    0.60    1.07    1.44    2.35
    ## beta[9]     0.02    0.00 0.01    0.00    0.01    0.02    0.02    0.03
    ## beta[10]    0.41    0.01 0.14    0.14    0.32    0.41    0.52    0.68
    ## beta[11]   -1.80    0.03 0.60   -3.00   -2.21   -1.79   -1.37   -0.72
    ## beta[12]   -0.62    0.02 0.43   -1.38   -0.93   -0.62   -0.33    0.19
    ## beta[13]   -2.88    0.02 0.61   -4.09   -3.26   -2.84   -2.45   -1.80
    ## lp__     -195.67    0.19 2.54 -201.24 -197.29 -195.45 -193.87 -191.54
    ##          n_eff Rhat
    ## beta[1]    254 1.00
    ## beta[2]    269 1.01
    ## beta[3]    177 1.02
    ## beta[4]    205 1.01
    ## beta[5]    930 1.00
    ## beta[6]    399 1.00
    ## beta[7]    556 1.00
    ## beta[8]    466 1.00
    ## beta[9]    602 1.00
    ## beta[10]   748 1.00
    ## beta[11]   450 1.00
    ## beta[12]   737 1.00
    ## beta[13]   752 1.00
    ## lp__       182 1.00
    ## 
    ## Samples were drawn using NUTS(diag_e) at Fri Feb 28 07:49:49 2020.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

    # Examine the coefficient estimates
    coef(model1)

    ##                Estimate            LB          UB
    ## (Intercept) -4.81901733  -7.789123864 -1.68751112
    ## var6        -0.10674397  -1.981043295  1.51233659
    ## var3        -5.30076900 -10.840607199 -1.13156396
    ## var5         3.07629825   0.963926360  5.78875404
    ## var10        1.32320303   0.150698098  2.55984875
    ## var15        0.82061505   0.544685268  1.15065942
    ## urgcy_yes    2.75164526   1.151922105  4.48366812
    ## urgcy_noleg  1.04807005  -0.105176133  2.35284643
    ## var4         0.01518161  -0.001493851  0.03206085
    ## var11        0.41465409   0.139894025  0.68455239
    ## var14       -1.80290166  -2.996690275 -0.72144895
    ## var7        -0.61902091  -1.375919925  0.19057199
    ## var17       -2.87639000  -4.088896423 -1.79992731

    # Check convergence (need to increase the number of chains for a better inspection)
    plot(model1)

    ## 'pars' not specified. Showing first 10 parameters by default.

![](cbq_demo_md_files/figure-markdown_strict/unnamed-chunk-4-1.png)

    plot(model1, type = "coef") # coefficient plots

    ## 'pars' not specified. Showing first 10 parameters by default.

    ## ci_level: 0.8 (80% intervals)

    ## outer_level: 0.95 (95% intervals)

![](cbq_demo_md_files/figure-markdown_strict/unnamed-chunk-4-2.png)
Additionally, we can also use the `predict` function to predict
probabilities based on the estimates:

    head(predict(model1),10)

    ##          mean       lower      upper
    ## 1  0.11380583 0.057162135 0.20669083
    ## 2  0.06787422 0.034680334 0.12051771
    ## 3  0.09249612 0.046550278 0.16751722
    ## 9  0.06146571 0.032320547 0.10970117
    ## 10 0.03132339 0.012996228 0.06401769
    ## 11 0.02739166 0.009872572 0.06895325
    ## 12 0.33752477 0.094446186 0.77365148
    ## 14 0.06694196 0.036189839 0.12299252
    ## 23 0.03782285 0.017075335 0.07475665
    ## 61 0.01280077 0.003689737 0.04289610

### Multiple choices with no variation in choice alternatives: data\_vote

When choice sets contain more than two choice alternatives, we need to
specify the conditionality of the estimation by using choice-set
indicators. Every choice set has a unique indicator so that the
estimator can tell the number of alternatives in a choice set and how
they are correlated with each other. In the package, the conditionality
of the estimation on the choice-set indicators is formulated by a
vertical bar immediately following the specification of the covariates:
`y ~ x|indicator variable`. In the data `data_vote`, all voters face
three choices: Bush, Clinton, or Perot. Thus, there is no variation in
the choice alternatives.

    # Take a subset of data_vote with only the first 30 choice sets.
    obs <- unique(data_vote$chid)[1:30]
    data <- subset(data_vote, data_vote$chid %in% obs)
    # Estimate the model
    model2 <- cbq(choice ~ dist + respfinp.bush + natlec.bush + respgjob.bush + resphlth.bush + 
      respblk.bush + respab.bush + east.bush + south.bush + west.bush + newvoter.bush + termlim.bush +
      deficit.bush + dem.bush + rep.bush + women.bush + educ.bush + age1829.bush + age3044.bush + age4559.bush
      + respfinp.clinton + natlec.clinton + respgjob.clinton + resphlth.clinton + 
      respblk.clinton + respab.clinton + east.clinton + south.clinton + west.clinton + newvoter.clinton + termlim.clinton +
      deficit.clinton + dem.clinton + rep.clinton + women.clinton + educ.clinton + age1829.clinton + age3044.clinton + age4559.clinton|chid, 
                 data = data, 
                 q = 0.5, # We first specify a 0.5 quantile (median)
                 vi = FALSE, # Indicating whether to use variational inference to approximate the results. The default is FALSE, meaning using MCMC sampling procedures.
                 nsim = 1000, # Running the model for 1000 iterations
                 grad_samples = 1, # Passed to vb (positive integer), the number of samples for Monte Carlo estimate of gradients, defaulting to 1.
                 elbo_samples = 100, # Passed to vb (positive integer), the number of samples for Monte Carlo estimate of ELBO (objective function), defaulting to 100. (ELBO stands for "the evidence lower bound".)
                 tol_rel_obj = 0.01, # Passed to vb (positive double), the convergence tolerance on the relative norm of the objective, defaulting to 0.01.
                 output_samples = 2000, # Passed to vb (positive integer), number of posterior samples to draw and save, defaults to 2000.
                 burnin = NULL, # The number of burnin iterations. If unspecified, hald of the iterations will be burnin iterations.
                 thin = 1, # Thinning parameter.
                 CIsize = 0.95, # The size of confidence interval.
                 nchain = 1, # The number of parallel chains.
                 seeds = 12345, # Random seeds to replicate the results.
                 inverse_distr = FALSE, # If FALSE, the ALD will not be reversed. The default is FALSE.
                 offset = 1e-20, # Offset values to enhance sampling stability. The default value is 1e-20.
                 mc_core = TRUE # Indicating whether the estimation will be run in multiple parallel chains. The default is TRUE.
                 )

    ## 
    ## SAMPLING FOR MODEL 'cbqdv' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000101 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.01 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 13.3507 seconds (Warm-up)
    ## Chain 1:                9.12629 seconds (Sampling)
    ## Chain 1:                22.477 seconds (Total)
    ## Chain 1:

Again, we can print and check the results:

    # Print the results
    print(model2)

    ## Conditional binary quantile regression 
    ## 
    ## Call:
    ## cbq(formula = choice ~ dist + respfinp.bush + natlec.bush + respgjob.bush + 
    ##     resphlth.bush + respblk.bush + respab.bush + east.bush + 
    ##     south.bush + west.bush + newvoter.bush + termlim.bush + deficit.bush + 
    ##     dem.bush + rep.bush + women.bush + educ.bush + age1829.bush + 
    ##     age3044.bush + age4559.bush + respfinp.clinton + natlec.clinton + 
    ##     respgjob.clinton + resphlth.clinton + respblk.clinton + respab.clinton + 
    ##     east.clinton + south.clinton + west.clinton + newvoter.clinton + 
    ##     termlim.clinton + deficit.clinton + dem.clinton + rep.clinton + 
    ##     women.clinton + educ.clinton + age1829.clinton + age3044.clinton + 
    ##     age4559.clinton | chid, data = data, q = 0.5, vi = FALSE, 
    ##     nsim = 1000, grad_samples = 1, elbo_samples = 100, tol_rel_obj = 0.01, 
    ##     output_samples = 2000, burnin = NULL, thin = 1, CIsize = 0.95, 
    ##     nchain = 1, seeds = 12345, inverse_distr = FALSE, offset = 1e-20, 
    ##     mc_core = TRUE)
    ## 
    ## MCMC run for 1000 iterations, with 500 used. 
    ## 
    ## Coefficients:
    ##                  Estimate      LB     UB
    ## (Intercept)        -0.479 -18.605 17.562
    ## dist               -4.370  -7.622 -1.200
    ## respfinp.bush       0.026  -8.334  8.159
    ## natlec.bush         1.943  -5.049 10.107
    ## respgjob.bush      -4.104 -13.777  5.876
    ## resphlth.bush       8.049   1.551 16.239
    ## respblk.bush        2.521  -4.066  9.157
    ## respab.bush        -8.437 -18.956  2.756
    ## east.bush           0.197 -20.620 21.052
    ## south.bush         -2.544 -16.326 13.411
    ## west.bush           0.568 -19.309 18.164
    ## newvoter.bush     -18.320 -35.316 -2.029
    ## termlim.bush       14.233  -0.062 30.434
    ## deficit.bush        2.701 -12.479 19.278
    ## dem.bush           -6.877 -24.636 10.360
    ## rep.bush            1.876 -15.349 18.358
    ## women.bush         -7.943 -24.863  8.251
    ## educ.bush           1.229  -4.714  7.214
    ## age1829.bush       -3.376 -18.269 10.491
    ## age3044.bush       -9.406 -24.070  5.298
    ## age4559.bush       -1.213 -18.227 16.455
    ## respfinp.clinton    7.272  -0.247 15.116
    ## natlec.clinton     15.014   4.177 26.429
    ## respgjob.clinton    0.742  -8.222  9.647
    ## resphlth.clinton   -8.309 -18.370  1.671
    ## respblk.clinton    -9.619 -17.307 -2.105
    ## respab.clinton    -12.711 -22.120 -3.063
    ## east.clinton        0.088 -17.483 19.294
    ## south.clinton       4.546 -11.185 19.891
    ## west.clinton        0.601 -19.954 20.072
    ## newvoter.clinton    3.711 -13.586 21.084
    ## termlim.clinton    11.461  -5.122 28.713
    ## deficit.clinton     1.392 -16.900 20.148
    ## dem.clinton         0.480 -16.593 17.192
    ## rep.clinton        -5.419 -24.619 11.895
    ## women.clinton      -2.445 -16.074 12.248
    ## educ.clinton       -1.084  -8.176  5.569
    ## age1829.clinton     2.084 -17.905 19.647
    ## age3044.clinton     8.759  -6.554 24.020
    ## age4559.clinton   -12.088 -32.012  7.926

    # Examine the coefficient estimates
    coef(model2)

    ##                      Estimate           LB        UB
    ## (Intercept)       -0.47941091 -18.60510943 17.561755
    ## dist              -4.37028757  -7.62249724 -1.200419
    ## respfinp.bush      0.02591841  -8.33411246  8.159270
    ## natlec.bush        1.94325237  -5.04920767 10.106515
    ## respgjob.bush     -4.10370732 -13.77735874  5.876419
    ## resphlth.bush      8.04873622   1.55133355 16.238642
    ## respblk.bush       2.52103798  -4.06559218  9.156750
    ## respab.bush       -8.43717837 -18.95640055  2.756086
    ## east.bush          0.19729292 -20.61968336 21.052040
    ## south.bush        -2.54426682 -16.32606814 13.410857
    ## west.bush          0.56801922 -19.30880849 18.163992
    ## newvoter.bush    -18.31982898 -35.31613178 -2.028722
    ## termlim.bush      14.23272959  -0.06247274 30.433986
    ## deficit.bush       2.70076293 -12.47899293 19.278460
    ## dem.bush          -6.87657078 -24.63594384 10.360383
    ## rep.bush           1.87552775 -15.34945743 18.358498
    ## women.bush        -7.94310923 -24.86301063  8.250544
    ## educ.bush          1.22889041  -4.71411577  7.214400
    ## age1829.bush      -3.37571072 -18.26943905 10.490509
    ## age3044.bush      -9.40626799 -24.07022801  5.298092
    ## age4559.bush      -1.21293745 -18.22703475 16.455109
    ## respfinp.clinton   7.27152339  -0.24654546 15.116165
    ## natlec.clinton    15.01412804   4.17748475 26.428846
    ## respgjob.clinton   0.74187366  -8.22223339  9.646937
    ## resphlth.clinton  -8.30939569 -18.37012475  1.671152
    ## respblk.clinton   -9.61906772 -17.30683492 -2.105019
    ## respab.clinton   -12.71097424 -22.11968637 -3.062707
    ## east.clinton       0.08820642 -17.48301848 19.294444
    ## south.clinton      4.54575408 -11.18500317 19.890875
    ## west.clinton       0.60066725 -19.95356252 20.072065
    ## newvoter.clinton   3.71060644 -13.58612762 21.084146
    ## termlim.clinton   11.46085601  -5.12158256 28.713361
    ## deficit.clinton    1.39208823 -16.89978817 20.147732
    ## dem.clinton        0.47973858 -16.59275549 17.192418
    ## rep.clinton       -5.41885430 -24.61937187 11.895420
    ## women.clinton     -2.44525615 -16.07398715 12.248482
    ## educ.clinton      -1.08367819  -8.17618982  5.569185
    ## age1829.clinton    2.08426498 -17.90539149 19.646744
    ## age3044.clinton    8.75886057  -6.55430448 24.019660
    ## age4559.clinton  -12.08822556 -32.01246932  7.925659

    # Check convergence (need to increase the number of chains for better inspection)
    plot(model2)

    ## 'pars' not specified. Showing first 10 parameters by default.

![](cbq_demo_md_files/figure-markdown_strict/unnamed-chunk-6-1.png)

### Multiple choices with varying choice alternatives: data\_coalition

When the choice alternatives vary in different choice sets, such as
government formation where different countries and different elections
have different potential coalitions, we need to additionally controling
for this variation in estimating the choice probabilities by quantiles.
While it can be quite complicated in math, in the package we only need
to add an indicating variable at the end of the variable specifation to
deal with varying choice alternatives. The following shows the example
of coalition formation.

    # Take a subset of data_coalition with only the first 10 choice sets.
    obs <- unique(data_coalition$case)[1:10]
    data <- subset(data_coalition, data_coalition$case %in% obs)
    # Estimate the model
    model3 <- cbq(realg ~ minor+ minwin + numpar + dompar + median + gdiv1 + mgodiv1 + prevpm + sq + mginvest + anmax2 + pro + anti|case, 
                 data = data, 
                 q = 0.5, # We first specify a 0.5 quantile (median)
                 vi = FALSE, # Indicating whether to use variational inference to approximate the results. The default is FALSE, meaning using MCMC sampling procedures.
                 nsim = 1000, # Running the model for 1000 iterations
                 grad_samples = 1, # Passed to vb (positive integer), the number of samples for Monte Carlo estimate of gradients, defaulting to 1.
                 elbo_samples = 100, # Passed to vb (positive integer), the number of samples for Monte Carlo estimate of ELBO (objective function), defaulting to 100. (ELBO stands for "the evidence lower bound".)
                 tol_rel_obj = 0.01, # Passed to vb (positive double), the convergence tolerance on the relative norm of the objective, defaulting to 0.01.
                 output_samples = 2000, # Passed to vb (positive integer), number of posterior samples to draw and save, defaults to 2000.
                 burnin = NULL, # The number of burnin iterations. If unspecified, hald of the iterations will be burnin iterations.
                 thin = 1, # Thinning parameter.
                 CIsize = 0.95, # The size of confidence interval.
                 nchain = 1, # The number of parallel chains.
                 seeds = 12345, # Random seeds to replicate the results.
                 inverse_distr = FALSE, # If FALSE, the ALD will not be reversed. The default is FALSE.
                 offset = 1e-20, # Offset values to enhance sampling stability. The default value is 1e-20.
                 mc_core = TRUE # Indicating whether the estimation will be run in multiple parallel chains. The default is TRUE.
                 )

    ## 
    ## SAMPLING FOR MODEL 'cbqdv' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.00192 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 19.2 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 24.7547 seconds (Warm-up)
    ## Chain 1:                20.5171 seconds (Sampling)
    ## Chain 1:                45.2718 seconds (Total)
    ## Chain 1:

Show the results:

    # Print the results
    print(model3)

    ## Conditional binary quantile regression 
    ## 
    ## Call:
    ## cbq(formula = realg ~ minor + minwin + numpar + dompar + median + 
    ##     gdiv1 + mgodiv1 + prevpm + sq + mginvest + anmax2 + pro + 
    ##     anti | case, data = data, q = 0.5, vi = FALSE, nsim = 1000, 
    ##     grad_samples = 1, elbo_samples = 100, tol_rel_obj = 0.01, 
    ##     output_samples = 2000, burnin = NULL, thin = 1, CIsize = 0.95, 
    ##     nchain = 1, seeds = 12345, inverse_distr = FALSE, offset = 1e-20, 
    ##     mc_core = TRUE)
    ## 
    ## MCMC run for 1000 iterations, with 500 used. 
    ## 
    ## Coefficients:
    ##             Estimate      LB     UB
    ## (Intercept)    0.207 -17.679 21.587
    ## minor         -0.454 -14.518 13.984
    ## minwin         5.133  -5.368 17.808
    ## numpar        -9.085 -20.658 -0.973
    ## dompar        10.043  -5.136 25.936
    ## median        -5.248 -15.850  2.234
    ## gdiv1         -6.204 -17.231  2.100
    ## mgodiv1       -8.517 -20.489  3.618
    ## prevpm        10.218  -3.715 25.715
    ## sq             6.432   2.950 11.340
    ## mginvest       0.108 -18.867 20.372
    ## anmax2        -0.237 -21.439 20.068
    ## pro            0.440 -19.952 22.420
    ## anti         -11.721 -26.300 -2.249

    # Examine the coefficient estimates
    coef(model3)

    ##                Estimate         LB         UB
    ## (Intercept)   0.2073859 -17.678637 21.5871038
    ## minor        -0.4544294 -14.518391 13.9837767
    ## minwin        5.1325820  -5.367509 17.8083680
    ## numpar       -9.0853487 -20.658318 -0.9726962
    ## dompar       10.0429065  -5.135943 25.9359808
    ## median       -5.2478665 -15.850199  2.2340944
    ## gdiv1        -6.2040626 -17.230546  2.1003287
    ## mgodiv1      -8.5169400 -20.489068  3.6180528
    ## prevpm       10.2175928  -3.715190 25.7148654
    ## sq            6.4318213   2.949932 11.3404063
    ## mginvest      0.1079683 -18.866675 20.3722624
    ## anmax2       -0.2370202 -21.439364 20.0683179
    ## pro           0.4400939 -19.952044 22.4199402
    ## anti        -11.7211375 -26.300406 -2.2485107

    # Check convergence (need to increase the number of chains for better inspection)
    plot(model3)

    ## 'pars' not specified. Showing first 10 parameters by default.

![](cbq_demo_md_files/figure-markdown_strict/unnamed-chunk-7-1.png)

### Incorporating mixed effects (fixed and/or random effects)

We are also able to incorporate mixed effects (fixed and/or random
effects) into the estimation of CBQ. To do so, we only need to modify
the syntax to
`y ~ x|indicator|variable of fixed effects|variable of random effects`.
In case we are only interested in either fixed or random effects, we can
replace the other with 1 so that it will not be estimated. For example,
if we are interested in the random effects, we can specify the formular
as `y ~ x|indicator|1|variable of random effects`. We consider again the
example of coalition formation by incorporating the country-fixed
effects. Because the estimation can be rather time consuming, I only
include the first two independent variables.

    # Take a subset of data_coalition
    obs <- unique(data_coalition$case)[c(1,40,50)]
    data <- subset(data_coalition, data_coalition$case %in% obs)
    # Estimate the model
    model4 <- cbq(realg ~ minor + minwin|case|country|1, 
                 data = data, 
                 q = 0.5, # We first specify a 0.5 quantile (median)
                 vi = FALSE, # Indicating whether to use variational inference to approximate the results. The default is FALSE, meaning using MCMC sampling procedures.
                 nsim = 1000, # Running the model for 1000 iterations
                 grad_samples = 1, # Passed to vb (positive integer), the number of samples for Monte Carlo estimate of gradients, defaulting to 1.
                 elbo_samples = 100, # Passed to vb (positive integer), the number of samples for Monte Carlo estimate of ELBO (objective function), defaulting to 100. (ELBO stands for "the evidence lower bound".)
                 tol_rel_obj = 0.01, # Passed to vb (positive double), the convergence tolerance on the relative norm of the objective, defaulting to 0.01.
                 output_samples = 2000, # Passed to vb (positive integer), number of posterior samples to draw and save, defaults to 2000.
                 burnin = NULL, # The number of burnin iterations. If unspecified, hald of the iterations will be burnin iterations.
                 thin = 1, # Thinning parameter.
                 CIsize = 0.95, # The size of confidence interval.
                 nchain = 1, # The number of parallel chains.
                 seeds = 12345, # Random seeds to replicate the results.
                 inverse_distr = FALSE, # If FALSE, the ALD will not be reversed. The default is FALSE.
                 offset = 1e-20, # Offset values to enhance sampling stability. The default value is 1e-20.
                 mc_core = TRUE # Indicating whether the estimation will be run in multiple parallel chains. The default is TRUE.
                 )

    ## 
    ## SAMPLING FOR MODEL 'cbqfixdv' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.079709 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 797.09 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 154.766 seconds (Warm-up)
    ## Chain 1:                130.542 seconds (Sampling)
    ## Chain 1:                285.308 seconds (Total)
    ## Chain 1:

Show the results:

    # Print the results
    print(model4)

    ## Conditional binary quantile regression 
    ## 
    ## Call:
    ## cbq(formula = realg ~ minor + minwin | case | country | 1, data = data, 
    ##     q = 0.5, vi = FALSE, nsim = 1000, grad_samples = 1, elbo_samples = 100, 
    ##     tol_rel_obj = 0.01, output_samples = 2000, burnin = NULL, 
    ##     thin = 1, CIsize = 0.95, nchain = 1, seeds = 12345, inverse_distr = FALSE, 
    ##     offset = 1e-20, mc_core = TRUE)
    ## 
    ## MCMC run for 1000 iterations, with 500 used. 
    ## 
    ## Coefficients:
    ##             Estimate      LB     UB
    ## (Intercept)   -0.531 -20.977 18.531
    ## minor         12.734   4.823 24.724
    ## minwin        -5.210 -24.359 10.200

    # Examine the coefficient estimates
    coef(model4)

    ##               Estimate         LB       UB
    ## (Intercept) -0.5312046 -20.977191 18.53107
    ## minor       12.7336211   4.823037 24.72359
    ## minwin      -5.2096248 -24.358964 10.19952

    # Check convergence (need to increase the number of chains for better inspection)
    plot(model4)

![](cbq_demo_md_files/figure-markdown_strict/unnamed-chunk-8-1.png)

In a nutshell, the model specification in `cbq` has the following
patterns:

1.  Simple binary choices: `y ~ x`.
2.  Multiple choices with or without varying choice alternatives:
    `y ~ x|choice-set indicator`.
3.  Incorporating mixed-effects:
    `y ~ x|choice-set indicator|fixed-effects|random-effects`.

Finally, as a concluding remark, the above examples show the basic
functions of the package, which can be further extended in the updated
versions. Some details of the functions, such as naming of the
variables, will be refined in the renewed version. For the
interpretation and additional visualization of the results, you can
check [Lu (2020)](https://doi.org/10.1017/pan.2019.29) and other related
works.
