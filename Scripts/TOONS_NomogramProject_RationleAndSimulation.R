# TOONS_NomogramProject_RationleAndSimulation.R
# ------------------------------------------------------------------------------
# Project examining the effect of information formats on decision making
# The Ottawa Operative Neuropsychology System (TOONS)
# Bryce P Mulligan, PhD, CPsych
# 19 April 2026
# ------------------------------------------------------------------------------
## ReadMe

# This project was meant as an experimental validation -- using undergraduate
# volunteers -- of the decision-support method proposed by Mulligan & Carniello
# (2023). In brief, the volunteers were asked to make decisions about whether or
# not they would consider changing majors based on information that was provided
# either (1) verbally, (2) using words and numbers, or (3) using words, numbers,
# and Bayesian nomograms. This study was meant as a proof of principle that
# could inform the development of future studies, ultimately eventuating in a
# study recruiting epilepsy neurosurgical candidates.

# The present R script walks through (rationale for) data simulation and model
# development for this project. Once the data are collected, the same script
# will be used to analyze the actual data using procedures validated on
# simulated data (i.e. wherein the input parameters are known).

# This R script is meant to document our rationale and plans for data
# collection, analysis, and visualisation. Page numbers, "R code" headers, and
# some text excerpts borrow heavily from -- and with sincere appreciation for --
# Richard McElreath's Statistical Rethinking (2nd Edition) and the accompanying
# code file. These are included to illustrate and provide context to our plans,
# expectations, and interpretations. In no way do we intend to portray excerpted
# content as our own original intellectual property.

# This script is intended as a record of our thinking at the time of creation
# and upload. It will further serve as a useful benchmark against which we can
# gauge the evolution of our understanding of and approach to these clinical
# issues.

# Finally, all data included in public repositories are simulated (i.e. no data
# from actual persons is included).
# ------------------------------------------------------------------------------
# remove all elements for a clean start (careful, now!)

# rm(list=ls(all=TRUE))  # clears the R environment
# dev.off()            # clears all plots
# cat("\014")            # clears R terminal
# ------------------------------------------------------------------------------
## Load packages from library
library(rethinking)
# ------------------------------------------------------------------------------
# Revamped April 2025 to include estimation of parameter covariance
# (see Statistical Rethinking, 2nd Ed., Chapter 14: Adventures in Covariance)
# ------------------------------------------------------------------------------
## Simulate data
# ------------------------------------------------------------------------------
### 1) Simulate Populuation
# ------------------------------------------------------------------------------
# A) Participants (intercepts)
# mBias <- 0      # average participant bias (participant intercept)
# sdBias <- 1   # sd of participant bias (sd of participant intercepts)

# B) Treatments (slopes)
ProbSet <- c(0.05, 0.10, 0.05, 0.40, 0.05, 0.80) # ideal treatment effect
# ProbSet <- c(0.50, 0.50, 0.50, 0.50, 0.50, 0.50)   # null treatment effect
 # convert to log-odds using logit function
# see: https://www.statisticshowto.com/log-odds/
loProbSet <- logit(ProbSet) # average treatment effect (slope)
sdProbSet <- 0.5 # std dev in slopes (treatments)

# C) Variance-Covariance Matrix
# These values define the entire population of participants. To use these values
# to simulate a sample of participants, we’ll need to build them into a
# 7-dimensional multivariate Gaussian distribution. This means we need a vector
# of 7 means and 7-by-7 matrix of variances and covariances. The means are
# easiest. The vector we need is just:
Mu <- loProbSet

# That’s it. The values are the the mean slopes, the difference in response bias
# associated with each of the 6 treatments.

# The matrix of variances and covariances is arranged like this:

#      variance of intercepts       | covariance of intercepts & slopes
# _____________________________________________________________________
# covariance of intercepts & slopes |       variance of slopes

# To build this matrix with R code, there are several options. One common way to
# build the covariance matrix is conceptually very useful, because it treats the
# standard deviations and correlations separately. Then it matrix multiplies
# them to produce the covariance matrix. We’re going to use this approach later
# on, to define priors, so it’s worth seeing it now. Here’s how it’s done:
## R code 14.5
sigmas <- c(rep(sdProbSet,6)) # standard deviations

# correlation matrix (manually created to reflect plausible/ideal treatment
# effects)
Rho <- matrix( c(1, -0.1, 0.2, -0.2, 0.4, -0.4,
                 -0.1, 1, -0.2, 0.2, -0.4, 0.4, 
                 0.2, -0.2, 1, -0.4, 0.6, -0.6, 
                 -0.2, 0.2, -0.4, 1, -0.6, 0.6, 
                 0.4, -0.4, 0.6, -0.6, 1, -0.8, 
                 -0.4, 0.4, -0.6, 0.6, -0.8, 1), nrow=6 ) 

# now matrix multiply to get covariance matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)
#If you are not sure what diag(sigmas) accomplishes, then try typing just
#diag(sigmas) at the R prompt.
# ------------------------------------------------------------------------------
### 1) Simulate Participant Characteristics (intercepts and slopes)
# ------------------------------------------------------------------------------
#Now we’re ready to simulate some participant, each with its own slope.
#Let’s define the number of participants:
N <- 100        # number of participants

# And to simulate their properties, we just sample randomly from the
# multivariate Gaussian distribution defined by Mu and Sigma:
## R code 14.7
library(MASS)
set.seed(5) # used to replicate example
vary_effects <- mvrnorm( N , Mu , Sigma )

# Note the set.seed(5) line above. That’s there so you can replicate the precise
# results in the example figures. The particular number, 5, produces a
# particular sequence of random numbers. Each unique number generates a unique
# sequence. Including a set.seed line like this in your code allows others to
# exactly replicate your analyses. Later you’ll want to repeat the example
# without repeating the set.seed call, or with a different number, so you can
# appreciate the variation across simulations.

# Look at the contents of vary_effects now. It should be a matrix with 100 rows
# and 6 columns. Each row is a participant. Each column is a treatment.
vary_effects
plot(as.data.frame(vary_effects))

# For transparency, let’s split these columns apart into nicely named vectors:
## R code 14.8
Tx_1 <- vary_effects[,1]
Tx_2 <- vary_effects[,2]
Tx_3 <- vary_effects[,3]
Tx_4 <- vary_effects[,4]
Tx_5 <- vary_effects[,5]
Tx_6 <- vary_effects[,6]

# To visualize these slopes, go ahead and plot them against one
# another. This code will also show the distribution’s contours:
## R code 14.9
plot( Tx_4 , Tx_6 , col=rangi2 ,
      xlab="Effect of Num/Inc Treatment (Tx_4)" , 
      ylab="Effect of Nom/Inc Treatment (Tx_6)" )
# overlay population distribution
library(ellipse)
for ( l in c(0.1,0.3,0.5,0.8,0.99) )
  lines(ellipse(c(Sigma[4,6], Sigma[6,4]),centre=c(Mu[4], Mu[6]),level=l),col=col.alpha("black",0.2))
# 100 participants sampled from a statistical population. The horizontal
# axis is the Effect of Numbers/Inc Treatment for each participant. The vertical
# axis is the Effect of Nomogram/Inc Treatment for each participant. The gray
# ellipses illustrate the multivariate Gaussian population of intercepts and
# slopes.

# In any particular simulation, the correlation may not be as obvious. But on
# average, the correlations between the effects will reflect the input
# correlation matrix (Rho), and you’ll be able to see this in the scatterplot.
# The contour lines in the plot, produced by the ellipse package (make sure you
# install it), display the multivariate Gaussian population of effects that the
# 100 participants were sampled from.

# 14.1.2. Simulate observations. We’re almost done simulating. What we did above
# was simulate individual participants and their average properties (i.e.
# response to different treatments). Now all that remains is to simulate the
# individual experimental trials and resultant response data. The code below
# simulates 10 trials for each participant. The script generates a binary
# decision (no/yes, 0/1) for each trial. Then it they are combined  into a
# common data frame.
## R code 14.10
set.seed(22)
TrialN <- 12    # number of trials per participant
ID <- rep(1:N, each=TrialN)  # participant ID numbers

Trial <- rep(1:TrialN,N) # 10 trials for each of 100 participants
Cond <- rep(1:3, length=TrialN*N) # format of information (word, num, nomo)
Dir <- rep(0:1, length=TrialN*N)  # direction of information nudge (decrease, increase)
Treatment <- rep(1:6, length=TrialN*N)  # all combinations of Cond and Dir

# Simulate individual decisions
DecisionProbM <- cbind(vary_effects,vary_effects) # each participant gets 2 trials of each treatment
DecisionProb <- as.vector(t(DecisionProbM))
Decision <- rbinom(TrialN*N, 1, prob=(inv_logit(DecisionProb)))

# D) bind into dataframe
d <- as.data.frame(cbind(ID, #Bias, 
                         Trial, Cond, Dir, Treatment, 
                         DecisionProb,
                         Decision))
precis(d)


# E) plot raw data
labs <- c("Wor/Dec","Wor/Inc","Num/Dec","Num/Inc","Nomo/Dec","Nomo/Inc")
# This is a plot of the raw simulated participant trajectories. Each line is a 
# single simulated participant.
# First, create a dataset of mean values for each participant.
plot_dat <- aggregate(x = d$Decision,               # Specify data column
                      by = list(d$ID, d$Treatment),  # Specify group indicator
                      FUN = mean)                     # Specify function (i.e. mean)
# Reshape from long to wide format
plot_dat_w <- reshape(plot_dat, idvar = "Group.1", 
                      timevar = "Group.2", 
                      direction = "wide") 
# Drop first column and convert to matrix
plot_dat_mat <- as.matrix(subset(plot_dat_w, select = -c(1)))
dens( plot_dat_mat[1,], ylim=c(0,7)) # plot first row
for ( i in 2:nrow(plot_dat_mat) ) dens( plot_dat_mat[i,] , add=TRUE) # overlay other rows

# Plot the trajectories
plot( NULL , xlab="Treatment" , ylab="proportion 'change major' decisions" ,
      ylim=c(0,1) , xaxt="n" , xlim=c(1,6) )
axis( 1 , at=levels(as.factor(Treatment)) , labels=labs )
# grau makes the lines transparent so can see more/less common trajectories
for ( i in 1:nrow(plot_dat_mat) ) lines( 1:6 , plot_dat_mat[i,] , col=grau(0.25) , lwd=2 )

# ------------------------------------------------------------------------------
## Modeling

## R code 11.10
# trimmed data list
dat_list <- list(
  Decision = d$Decision,
  ID = d$ID,
  Treatment = as.integer(d$Treatment) )

# Now we can start the Markov chain. I’ll add log_lik=TRUE to the call, so that
# ulam computes the values necessary for PSIS and WAIC.
# ## R code 11.11
mario1.0 <- ulam(
  alist(
    Decision ~ dbinom( 1 , p ) ,
    logit(p) <- a[ID] + b[Treatment] ,
    a[ID] ~ dnorm( 0 , 1.5 ),
    b[Treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list , chains=4 , cores=4, iter=1e4, log_lik=TRUE )

precis( mario1.0 , depth=2 )
trankplot(mario1.0)
model <- mario1.0


########## Try multilevel model allowing for partial pooling of the
# participant-specific intercepts in order to allow for varying response bias.
# (model 13.4, p.415)
mario2.0 <- ulam(
  alist(
    Decision ~ dbinom( 1 , p ) ,
    logit(p) <- a[ID] + b[Treatment] ,
    b[Treatment] ~ dnorm( 0 , 0.5 ),
    a[ID] ~ dnorm( a_bar , sigma_a ),
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1)
  ) , data=dat_list , chains=4 , cores=4, iter=1e4, log_lik=TRUE )

precis( mario2.0 , depth=2 )
trankplot(mario2.0)
model <- mario2.0


# coeftab(mario1.0, mario2.0)
# ------------------------------------------------------------------------------
# ########## Try another multilevel model, this time including partial pooling
# of the treatment effects in an attempt to minimize the confound introduced by
# the participant-specific variation in Treatment on Decision outcomes.
# (model 13.6, p.419)
set.seed(4387510) # for reproducibility
mario3.0 <- ulam(
  alist(
    Decision ~ dbinom( 1 , p ) ,
    logit(p) <- a[ID] + b[Treatment] ,
    b[Treatment] ~ dnorm( 0 , sigma_b ),
    a[ID] ~ dnorm( a_bar , sigma_a ),
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1),
    sigma_b ~ dexp(1)
  ) , data=dat_list , chains=2 , cores=4, iter=1e4, log_lik=TRUE,
      control=list(adapt_delta=0.99) )

precis( mario3.0 , depth=2 )
trankplot(mario3.0)
WAIC(mario3.0)
model <- mario3.0

# coeftab(mario2.0, mario3.0)


# ####### (p.420) Estimating this model results in a number of divergent
# transitions. This is common in multilevel models and usually happens when the
# posterior distribution is very steep in some region of parameter space. First,
# we can try to tune the simulation so that it doesn't over shoot and valley
# wall. This means doing warmup with an increase in Stan’s
# target acceptance rate. This is controlled by the adapt_delta control
# parameter. The ulam default is 0.95, which means that it aims to attain a 95%
# acceptance rate. It tries this during the warmup phase, adjusting the step
# size of each leapfrog step (go back to Chapter 9 if these terms aren’t
# familiar). When adapt_delta is set high, it results in a smaller step size,
# which means a more accurate approximation of the curved surface. It can also
# mean slower exploration of the distribution.

# # Increasing adapt_delta will often, but not always, help with divergent
# # transitions. For example, model mario3.0 in the previous section presented
# divergent(mario3.0)
# # divergent transitions. We can re-run the model, using a higher target
# # acceptance rate, with:
# ## R code 13.28
# mario3.1 <- ulam( mario3.0 , chains=4 , cores=4 , iter=1e4 ,
#                   control=list(adapt_delta=0.99) )
# divergent(mario3.1)
# # So that did help. But sometimes this won’t be enough. And
# # while the divergences are reduced, the chain still isn’t very efficient—look at
# # the precis output and notice that many of the n_eff values are still far below
# # the true number of samples (10000 in this case: 4 chains, 2500 from each).
# precis(mario3.1, depth=2)
# ------------------------------------------------------------------------------
# We can do much better with the non-centered version of the model. What we want
# is a version of mario3.0 in which we get the parameters out of the adaptive
# priors and instead into the linear model.

mario3.2 <- ulam(
  alist(
    Decision ~ dbinom( 1 , p ) ,
    logit(p) <- a_bar + z[ID]*sigma_a +  # participant intercepts
                x[Treatment]*sigma_b ,   # Treatment intercepts
    z[ID] ~ dnorm( 0 , 1 ),
    x[Treatment] ~ dnorm( 0 , 1 ),
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1),
    sigma_b ~ dexp(1),
    gq> vector[ID]:a <<- a_bar + z*sigma_a,
    gq> vector[Treatment]:b <<- x*sigma_b
  ) , data=dat_list , chains=4 , cores=4, iter=1e4, warmup=1e3, log_lik=TRUE )

# divergent(mario3.2) # that did 'er!
# precis(mario3.2, depth=2)
# and it is much more efficient, as shown in the the ess_bulk (effective sample
# size) approaching the true number of samples.
model <- mario3.2

# ------------------------------------------------------------------------------
# (Based on model m14.2, described beginning on p.447)
# Now we'll model both types of clusters (ID and Treatment) and place varying
# effects on the intercepts and slopes.

# ## R code 14.18
# mario4.0 <- ulam(
#   alist(
#     Decision ~ dbinom( 1 , p ) ,
#     logit(p) <- b[Treatment] + a[ID,Treatment],
#     
#     # adaptive priors
#     # mulitvariate Gaussian prior is 6-dimensional because there are 6 treatments
#     vector[6]:a[ID] ~ multi_normal(0,Rho_ID,sigma_ID),
#     
#     # fixed priors
#     b[Treatment] ~ dnorm(0,1),
#     sigma_ID ~ dexp(1),
#     Rho_ID ~ dlkjcorr(4)
#   ) , data=dat_list , chains=4 , cores=2, iter=1e4, warmup=1e3, log_lik=TRUE,
#       control=list(adapt_delta=0.99)) # to deal with divergent transitions
# 
# precis( mario4.0 , depth=2 )
# # precis( mario4.0 , depth=3 )
# # trankplot(mario4.0)
# 
# model <- mario4.0

# ------------------------------------------------------------------------------
# Use non-centered parameterization to address divergent transitions
## R code 14.19
set.seed(4387510)
mario4.1 <- ulam(
  alist(
    Decision ~ dbinom( 1 , p ) ,
    logit(p) <- b[Treatment] + a[ID,Treatment],
    
    # adaptive priors - non-centered
    transpars> matrix[ID,6]:a <-
      compose_noncentered( sigma_ID , L_Rho_ID , z_ID ),
    matrix[6,ID]:z_ID ~ normal( 0 , 1 ),
    
    # fixed priors
    b[Treatment] ~ dnorm(0,1),
    vector[6]:sigma_ID ~ dexp(1),
    cholesky_factor_corr[6]:L_Rho_ID ~ lkj_corr_cholesky( 2 ),
    
    # compute ordinary correlation matrices from Cholesky factors
    gq> matrix[6,6]:Rho_ID <<- Chol_to_Corr(L_Rho_ID)
  ) , data=dat_list , chains=4 , cores=2, iter=1e4, warmup=1e3, log_lik=TRUE)

trankplot(mario4.1)
traceplot(mario4.1)
WAIC(mario4.1)

compare( mario3.2 , mario4.1)

# The difference between the models is about 8.2 and the standard error is about
# 4.95. If we imagine the 99% interval of the difference (corresponding to a
# z-score of about 2.6), it'll be about:
8.2 + c(-1,1)*4.95*2.6
# So these models are not necessarily easy to distinguish by expected
# out-of-sample accuracy. Model mario4.1 is marginally better in this respect.
# However, we know that mario4.1 is theoretically correct because we simulated
# the data. Both theoretical and empirical considerations like these need to
# factor in any decision about selecting/developing a model to analyze a
# particular dataset. In short, WAIC cannot be used to infer causation.

# You might be able to see all of this better, if we plot the compare table:
## R code 7.29
plot(compare( mario3.2 , mario4.1))

# The filled points are the in-sample deviance values. The open points are the
# WAIC values. Notice that naturally each model does better in-sample than it is
# expected to do out-ofsample. The line segments show the standard error of each
# WAIC. These are the values in the column labeled dSE in the table above. So
# you can probably see that mario4.1 is marginally better than mario3.2. What we
# really want however is the standard error of the difference in WAIC between
# the two models. That is shown by the lighter line segment with the triangle on
# it, between mario4.1 and mario3.0.

model <- mario4.1
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
## Extract samples and visualizize posterior
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# This is the guts of the tide prediction engine. We’ll need to do a little work
# to interpret it. The first 100 "a" parameters are the intercepts unique to
# each participant. Each of these expresses the tendency of each individual to
# make a "change majors" response. Let’s look at these on the outcome scale:
## R code 11.12
post <- extract.samples(model)

# inverse-logit to transform back to outcome (probability) scale
p_Decision_a <- inv_logit( post$a ) 
plot( precis( as.data.frame(p_Decision_a) ) , xlim=c(0,1), xlab="Probability" )

# Each row is a participant, the numbers corresponding to the values in ID.
# There are substantial differences among the participants in their baseline
# tendencies, as well as subtantial variation within participants. This is
# exactly the kind of effect that makes pure experiments difficult in the
# behavioral sciences. Having repeat measurements, like in this experiment, and
# measuring them is very useful.

# Now let’s consider the treatment effects, hopefully estimated more precisely
# because the model could subtract out the variation among participants. 
# (inverse-logit to transform back to outcome (probability) scale)
p_Decision_b <- inv_logit( post$b )
plot( precis( as.data.frame(p_Decision_b) ) , 
      xlim=c(0,1) , xlab="Probability", labels=labs)
# I’ve added treatment labels in place of the parameter names.
# (i.e. words/decrease; words/increase;
#       numbers/decrease; numbers/increase;
#       nomogram/decrease; nomogram/increase)

# Here are the treatment effects on the logit (log-odds) scale:
## R code 11.13
plot( precis( model , depth=2 , pars="b" ) , xlab="Log-Odds", labels=labs )

# To understand these distributions, it’ll help to consider our expectations.
# What we are looking for is evidence that the participants choose to the
# "change majors" option more when they are presented with information that (1)
# suggests an increased chance of career satisfaction AND (2) is presented in an
# informative/comprehensible/convincing way. In particular, we are expecting
# that providing numbers and/or nomograms suggesting the benefits of changing
# majors will result in more "change majors" decisions than providing words
# alone. This implies comparing the fourth and sixth row to the second row. By
# way of a manipulation check, we will also verify that the information format
# (words, numbers, nomograms) did not have any influence on decisions to "change
# majors" when the information suggested that staying in the current major was
# advantageous (i.e. compare third and fifth row with first row). Let’s
# calculate the differences between information formats on the log-odds scale.
## R code 11.14
diffs <- list(
  db31 = post$b[,3] - post$b[,1],
  db51 = post$b[,5] - post$b[,1],
  db42 = post$b[,4] - post$b[,2],
  db62 = post$b[,6] - post$b[,2] )
plot( precis(diffs) , xlab="Log-Odds Difference Score")

# These are the contrasts between the information format treatments. The scale
# is log-odds of making a "change major" decision. Remember the tide engine!
# db31 is the difference between number/words treatments when the information
# was unfavourable to a major-switch. So if there is evidence of more "change
# major" choice when numbers are presented, this will show up here as a larger
# difference, consistent with deciding to switch more when numbers are present.
# Clearly that is not the case, as the estimates are near zero (and the 89%
# compatability interval overlaps with zero). db42 is the same difference
# (numbers/words), but for when the information is favourable to a "major
# change". There is indeed evidence that participants made more "change major"
# decisions when favourable information that included numbers was presented.

# Now let’s consider a posterior prediction check. Let’s summarize the
# proportions of "change major" decisions for each treatment and then plot
# against the posterior predictions. First, to calculate the proportion in
# treatment:
prM <- by( d$Decision , d$Treatment , mean )
prI <- by( d$Decision , d$Treatment , HPDI )

# Each column is a treatment. The cells contain the mean proportions of
# decisions that were to change major. The model will make predictions for these
# values, so we can see how the posterior predictions look against the raw data.
# Remember that we don’t want an exact match—that would mean overfitting. But we
# would like to understand how the model sees the data and learn from any
# anomalies.

# We can compute posterior predictions using link:
## R code 11.17
dat <- list( Treatment=rep(1:6) )
p_post <- link( model , data=dat )
p_mu <- apply( p_post , 2 , mean )
p_ci <- apply( p_post , 2 , HPDI )

# This is a plot of the raw data (closed circles), the set probabilities for
# data simulation, which are the parameters we are trying to recover (red
# asterisks), and posterior predictions (open circles, with 89% HPDI).
plot(prM , ylim=c(0,1), pch=16, xlab="Treatment" , xaxt = "n", # no x-axis
     ylab="proportion 'change major' decisions" )
axis( 1 , at=levels(as.factor(Treatment)) , labels=labs )
abline( h=0.5 , lty=2 )  # chance responding
abline( h=mean(d$Decision), col="blue", lwd=2)  # overall mean Decision
# abline( h=mean(d$Bias), col="green", lwd=2)  # overall mean Bias
points(ProbSet, pch=8, col="red")
points(p_mu, pch=1)
for ( i in 1:6) {
  arrows(i, p_ci[1,i], i, p_ci[2,i], 
                       length=0.05, angle=90, code=3)
}
# ------------------------------------------------------------------------------

# Eventually, add multilevel posterior prediction for new clusters, as
# described/illustrated beginning on p.428