Ch 11 part I: GLMs and binomial regression
================

``` r
knitr::opts_chunk$set(
  #tidy = TRUE,
  #tidy.opts = list(width.cutoff = 95),
  message = FALSE,
  warning = TRUE, 
  eval = TRUE
  # root.dir = here()
)
```

Generalized linear models for counts.

Logistic and aggregated binomial regression are generalizations of
binomial regression. If we want to model some observed counts as a
function of some predictor variables we can start with: y(count) ~
Binomial(n,p) p is the proabbility of a success on some “trial” and n is
the number of trials. Logistic regression is binomial when outcomes can
take values 0, 1. When individual trials with the same covariate are
aggregated together it’s “aggregated binomial regression” – outcome can
take any value 0, to n, the number of trials. Both use the logit link.
As in ch 10, a link function maps the linear space of the predictor
linear model betaX + alpha onto the non linear space of the parameter in
a probability distribution.

see p 325 for experiment description.

first quickly plot the actual outcome data per actor bc i dont
understand what i am acually trying to model here. ok understand now;
each actor was tested multiple times under different treatment
conditions; the outcome variable is the percentage of time the actor
pulled left but pulling left is a different thing depending on the
treatment. the point of hte treatment vars for no food present is to
determine the handedness of the animal / prefernce if there is one.

1 = 2 food R, no partner  
2 = 2 food L, no partner  
3 = 2 food R, partner  
4 = 2 food L, partner

Key - there is always food for the animal on both sides they can pull
left or right and get food. Only by pulling right do they help out the
other guy.

``` r
## R code 11.1
library(rethinking)
r = rangi2
data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition
```

``` r
library(ggplot2); theme_set(theme_bw())


ggplot(d, aes(x = as.factor(treatment), y = pulled_left, fill = as.factor(treatment))) + 
  geom_col() + 
  facet_wrap(~actor)
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
unique(d$treatment)
```

    ## [1] 1 2 3 4

``` r
table(d$treatment, d$prosoc_left)
```

    ##    
    ##       0   1
    ##   1 126   0
    ##   2   0 126
    ##   3 126   0
    ##   4   0 126

The model:

L\[i\] ~ Binomial( 1, p\[i\] )  
logit(p\[i\]) = a\[ actor\[i\] \] + B\[ treatment\[i\] \]  
a\[j\] ~ tbd  
B\[k\] ~ tbd

L is a 0 , 1 indicator outcome for “pulled\_left” so it could also be
writeen as bernoulli(pi) == Binomial(1, pi).  
There are 7 alpha parameters one for each chimp.  
4 treatment parameters , one for each of the combined factor indicating
presence of partner and choise of prosocial option.

Specifying priors - consider small model with single a variable.

**remmber: logit means log odds** so we will have to go back and forth
between the log odds and hte outcome probability scale with an inverse
logit transform.

``` r
## R code 11.4
# very wide prior 
m11.1 <- quap(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a ,
        a ~ dnorm( 0 , 10 )
    ) , data=d )

## R code 11.5
set.seed(1999)
prior <- extract.prior( m11.1 , n=1e4 )
```

Now need to convert the prior to the outcome scale by calculating the
inverse function. This shows a flat prior in logit space is not a flat
prior in the outcome probability space–the model thinks that the chimps
either never or always pull left lever.

``` r
## R code 11.6
p <- inv_logit( prior$a )
dens( p , adj=0.1 )
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Adding treatment effects–show what a flat prior for beta implies:  
L\[i\] ~ Binomial( 1, p\[i\] )  
logit(p\[i\]) = a\[ actor\[i\] \] + B\[ treatment\[i\] \]  
a\[j\] ~ N(0, - )  
B\[k\] ~ N(0, -)

``` r
## R code 11.7
# make sure i understand: 
# f11.2 = alist(
#   pulled_left ~ dbinom(1, p)
#   logit(p) <- a + B[treatment],
#   a ~ dnorm(0,1.5), 
#   B ~ dnorm(0,10)
# )


m11.2 <- quap(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a + b[treatment] ,
        a ~ dnorm( 0 , 1.5 ),
        b[treatment] ~ dnorm( 0 , 10 )
    ) , data=d )
set.seed(1999)
prior <- extract.prior( m11.2 , n=1e4 )
names(prior)
```

    ## [1] "a" "b"

``` r
hist(prior$a, col = r)
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
par(mfrow = c(2,2))
sapply(1:4, function(x) hist(prior$b[ ,x], col=r))
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

    ##          [,1]           [,2]           [,3]           [,4]          
    ## breaks   Numeric,16     Numeric,18     Numeric,16     Numeric,17    
    ## counts   Integer,15     Integer,17     Integer,15     Integer,16    
    ## density  Numeric,15     Numeric,17     Numeric,15     Numeric,16    
    ## mids     Numeric,15     Numeric,17     Numeric,15     Numeric,16    
    ## xname    "prior$b[, x]" "prior$b[, x]" "prior$b[, x]" "prior$b[, x]"
    ## equidist TRUE           TRUE           TRUE           TRUE

Probability of pulling left for each treatment

``` r
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )
```

prior differences among treatments – The logit scale piles all
probability at 0 and 1. With a flat prior model thinks before seeing
dataa treatments are completely alike or completely different.

``` r
## R code 11.8
dens( abs( p[,1] - p[,2] ) , adj=0.1 )
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

With a N(0,0.5) prior the avg difference expected is around 10%.

``` r
## R code 11.9
m11.3 <- quap(
    alist(
        pulled_left ~ dbinom( 1 , p ) ,
        logit(p) <- a + b[treatment] ,
        a ~ dnorm( 0 , 1.5 ),
        b[treatment] ~ dnorm( 0 , 0.5 )
    ) , data=d )
set.seed(1999)
prior <- extract.prior( m11.3 , n=1e4 )
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )
dens( abs( p[,1] - p[,2] ) , adj=0.1 )
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Strong effects still possible. but these regularizing priors will help
reduce overfitting.

### fit HMC stan model.

If this was our central research question its reasonable but this takes
ridiculously long to fit. this would not scale well to 20 thousand
separately fit models.

``` r
## R code 11.10
# trimmed data list
dat_list <- list(
    pulled_left = d$pulled_left,
    actor = d$actor,
    treatment = as.integer(d$treatment) 
    )

## R code 11.11
f11.4 = alist(
  pulled_left ~ dbinom(1, p), 
  logit(p) <- a[actor] + b[treatment], 
  a[actor] ~ dnorm(0, 1.5), 
  b[treatment] ~ dnorm(0, 0.5)
)


m11.4 <- ulam(f11.4, data=dat_list , chains=4 , log_lik=TRUE)
```

    ## Running /Library/Frameworks/R.framework/Resources/bin/R CMD SHLIB foo.c
    ## clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/Rcpp/include/"  -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/"  -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/unsupported"  -I"/Users/matthewmule/Library/R/3.5/library/BH/include" -I"/Users/matthewmule/Library/R/3.5/library/StanHeaders/include/src/"  -I"/Users/matthewmule/Library/R/3.5/library/StanHeaders/include/"  -I"/Users/matthewmule/Library/R/3.5/library/rstan/include" -DEIGEN_NO_DEBUG  -D_REENTRANT  -DBOOST_DISABLE_ASSERTS -DBOOST_PENDING_INTEGER_LOG2_HPP -include stan/math/prim/mat/fun/Eigen.hpp   -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -I/usr/local/include   -fPIC  -Wall -g -O2  -c foo.c -o foo.o
    ## In file included from <built-in>:1:
    ## In file included from /Users/matthewmule/Library/R/3.5/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:13:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Dense:1:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Core:88:
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:613:1: error: unknown type name 'namespace'
    ## namespace Eigen {
    ## ^
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:613:16: error: expected ';' after top level declarator
    ## namespace Eigen {
    ##                ^
    ##                ;
    ## In file included from <built-in>:1:
    ## In file included from /Users/matthewmule/Library/R/3.5/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:13:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Dense:1:
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Core:96:10: fatal error: 'complex' file not found
    ## #include <complex>
    ##          ^~~~~~~~~
    ## 3 errors generated.
    ## make: *** [foo.o] Error 1
    ## 
    ## SAMPLING FOR MODEL '82480dff1a626a42c2ca9de938d65b9d' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 8.6e-05 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.86 seconds.
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
    ## Chain 1:  Elapsed Time: 0.470824 seconds (Warm-up)
    ## Chain 1:                0.386418 seconds (Sampling)
    ## Chain 1:                0.857242 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL '82480dff1a626a42c2ca9de938d65b9d' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 5.6e-05 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.56 seconds.
    ## Chain 2: Adjust your expectations accordingly!
    ## Chain 2: 
    ## Chain 2: 
    ## Chain 2: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 2: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 2: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 2: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 2: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 2: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 2: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 2: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 2: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 2: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 2: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 2: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 2: 
    ## Chain 2:  Elapsed Time: 0.519243 seconds (Warm-up)
    ## Chain 2:                0.321154 seconds (Sampling)
    ## Chain 2:                0.840397 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL '82480dff1a626a42c2ca9de938d65b9d' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 6e-05 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.6 seconds.
    ## Chain 3: Adjust your expectations accordingly!
    ## Chain 3: 
    ## Chain 3: 
    ## Chain 3: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 3: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 3: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 3: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 3: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 3: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 3: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 3: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 3: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 3: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 3: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 3: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 3: 
    ## Chain 3:  Elapsed Time: 0.491864 seconds (Warm-up)
    ## Chain 3:                0.37279 seconds (Sampling)
    ## Chain 3:                0.864654 seconds (Total)
    ## Chain 3: 
    ## 
    ## SAMPLING FOR MODEL '82480dff1a626a42c2ca9de938d65b9d' NOW (CHAIN 4).
    ## Chain 4: 
    ## Chain 4: Gradient evaluation took 5.9e-05 seconds
    ## Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.59 seconds.
    ## Chain 4: Adjust your expectations accordingly!
    ## Chain 4: 
    ## Chain 4: 
    ## Chain 4: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 4: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 4: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 4: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 4: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 4: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 4: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 4: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 4: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 4: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 4: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 4: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 4: 
    ## Chain 4:  Elapsed Time: 0.42815 seconds (Warm-up)
    ## Chain 4:                0.340795 seconds (Sampling)
    ## Chain 4:                0.768945 seconds (Total)
    ## Chain 4:

``` r
precis( m11.4 , depth=2 )
```

    ##             mean        sd        5.5%       94.5%     n_eff     Rhat4
    ## a[1] -0.45642629 0.3477435 -1.01318764  0.09523512  614.8645 1.0000542
    ## a[2]  3.88682837 0.7227322  2.80040209  5.11843419 1471.2153 1.0009275
    ## a[3] -0.75063656 0.3363033 -1.31493135 -0.23030060  648.5380 1.0015048
    ## a[4] -0.76124169 0.3421155 -1.30706452 -0.20818540  566.1462 1.0025290
    ## a[5] -0.46091281 0.3369121 -1.00762028  0.06616788  631.0298 0.9992719
    ## a[6]  0.46639103 0.3419080 -0.08814385  1.00046011  543.3119 1.0016696
    ## a[7]  1.94304138 0.4348036  1.25204631  2.68025204  821.7180 1.0014386
    ## b[1] -0.03074828 0.2949858 -0.49803371  0.45067799  565.9269 1.0013012
    ## b[2]  0.48927316 0.2919955  0.01573074  0.95312804  516.1493 1.0015648
    ## b[3] -0.37039709 0.2929954 -0.83303525  0.09380429  542.2040 0.9999995
    ## b[4]  0.37756874 0.2932367 -0.10293106  0.84357819  557.5963 1.0000941

### Examining model posterior distributions

We need to extract samples from the model fit and look at them on the
outcome scale with the inverse logit transform. Remember: model was
logit() so the betas are in terms of log odds; want to go back to
probblity with the inverse logit transform (the exp).

``` r
## R code 11.12
post <- extract.samples(m11.4)
p_left <- inv_logit( post$a )
plot(post$a, inv_logit(post$a), col = r , cex = 0.1)
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Now we extract alpha- the intercept term one for eah animal tells us if
there is a handedness perference. At 0.5 on the outcome
scale

``` r
plot( precis( as.data.frame(p_left) ) , xlim=c(0,1) )
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Now at the effects: **Note since we adjusted coefficients for the
variation in handedness across actors with the random intercept, these
estimates adjust for handedness**

LN means prosocial left no partner – meaning the “prosocial option” so
the chimp slid food down the table but there wasnt another chimp on the
other side.

``` r
## R code 11.13
labs <- c("R/N","L/N","R/P","L/P")
plot( precis( m11.4 , depth=2 , pars="b" ) , labels=labs )
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Now at the contrasts db13 is difference between beta 1 and beta 3. – do
the chimps choose pro social option more when there is another chimp on
the other side of the table?

contrasts between partner no partner treatment.

db13 says difference between no artner pertner treatments when prosocial
option is on the right. If there is evidence of more prosocial choice
when artner on the right will show up as a larger difference consistent
with pulling right more when partner is prestne. this hows pulled left
more when partner was absent…

``` r
## R code 11.14
diffs <- list(
    db13 = post$b[,1] - post$b[,3],
    db24 = post$b[,2] - post$b[,4] )
plot( precis(diffs) )
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

Below filled points means partner present.

``` r
## R code 11.15
pl <- by( d$pulled_left , list( d$actor , d$treatment ) , mean )
pl[1,]
```

    ##         1         2         3         4 
    ## 0.3333333 0.5000000 0.2777778 0.5555556

``` r
## R code 11.16
plot( NULL , xlim=c(1,28) , ylim=c(0,1) , xlab="" ,
    ylab="proportion left lever" , xaxt="n" , yaxt="n" )
axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
abline( h=0.5 , lty=2 )
for ( j in 1:7 ) abline( v=(j-1)*4+4.5 , lwd=0.5 )
for ( j in 1:7 ) text( (j-1)*4+2.5 , 1.1 , concat("actor ",j) , xpd=TRUE )
for ( j in (1:7)[-2] ) {
    lines( (j-1)*4+c(1,3) , pl[j,c(1,3)] , lwd=2 , col=rangi2 )
    lines( (j-1)*4+c(2,4) , pl[j,c(2,4)] , lwd=2 , col=rangi2 )
}
points( 1:28 , t(pl) , pch=16 , col="white" , cex=1.7 )
points( 1:28 , t(pl) , pch=c(1,1,16,16) , col=rangi2 , lwd=2 )
yoff <- 0.01
text( 1 , pl[1,1]-yoff , "R/N" , pos=1 , cex=0.8 )
text( 2 , pl[1,2]+yoff , "L/N" , pos=3 , cex=0.8 )
text( 3 , pl[1,3]-yoff , "R/P" , pos=1 , cex=0.8 )
text( 4 , pl[1,4]+yoff , "L/P" , pos=3 , cex=0.8 )
mtext( "observed proportions\n" )
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
## R code 11.17
dat <- list( actor=rep(1:7,each=4) , treatment=rep(1:4,times=7) )
p_post <- link( m11.4 , data=dat )
p_mu <- apply( p_post , 2 , mean )
p_ci <- apply( p_post , 2 , PI )
```

### Understanding relative and absolute effect size estimates from count glm

absolute effects – difference in probabiltiy of an event – this is the
outcome scale (probability pulled left).

relative effects – the changes in **odds** as a proportion. e.g. the
odds of an event double. These can be hard to interpret if the base rate
is low (doubling the odds from 1/ 1000 to 2/ 1000) but help in
understanding conditional effects.

To calculate relative odds, we exponentiate the estimate of the effect
size. Below we can calculate the contrast effect size of treatment 4 vs
2 in terms of relative odds.

The change in probability puled left is the contrast on the outcome
scale.

``` r
## R code 11.23
post <- extract.samples(m11.4)
mean( post$b[ ,4] - post$b[ ,2] )
```

    ## [1] -0.1117044

The outcome contrast on the relative odds scale: effect of switching
from treatent 2 to 4 is 0.92; can also interpret as a 100-x = 0.92 = 8%
reduction in odds:

``` r
mean( exp(post$b[ ,4] - post$b[ ,2]) )
```

    ## [1] 0.9275074

Note because the 4 vs 2 contrast is our effect of interest we are
calculating that on the outcome scale, then exponentiating that; need to
remember log rules–this is not equivalent to:

``` r
mean( exp(post$b[ ,4]) - exp(post$b[ ,2]))
```

    ## [1] -0.179353

### Aggregated binomial regression

Below we aggregate counts (left pull) per individual. We don’t care
about the order of individual pulls of the lever. The same information
is contained in a count of how many times each individual pulled the
left hand lever for each combination of predictor variables.

``` r
## R code 11.24
suppressMessages(library(tidyverse))
data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition
d$side <- d$prosoc_left + 1 # right 1, left 2
d$cond <- d$condition + 1 # no partner 1, partner 2

d_aggregated = d %>%
  group_by(treatment, actor, side, condition) %>% 
  summarize(pulled_left = sum(pulled_left))

d_aggregated
```

    ## # A tibble: 28 x 5
    ## # Groups:   treatment, actor, side [28]
    ##    treatment actor  side condition pulled_left
    ##        <dbl> <int> <dbl>     <int>       <int>
    ##  1         1     1     1         0           6
    ##  2         1     2     1         0          18
    ##  3         1     3     1         0           5
    ##  4         1     4     1         0           6
    ##  5         1     5     1         0           6
    ##  6         1     6     1         0          14
    ##  7         1     7     1         0          14
    ##  8         2     1     2         0           9
    ##  9         2     2     2         0          18
    ## 10         2     3     2         0          11
    ## # … with 18 more rows

### Fit and compare the aggregated binomial to the individual count binomial.

``` r
## R code 11.25
dat <- with( d_aggregated , list(
    left_pulls = pulled_left,
    treatment = treatment,
    actor = actor,
    side = side,
    cond = condition ) )

dat
```

    ## $left_pulls
    ##  [1]  6 18  5  6  6 14 14  9 18 11  9 10 11 15  5 18  3  2  5 10 17 10 18
    ## [24]  6  8  9 11 18
    ## 
    ## $treatment
    ##  [1] 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4
    ## 
    ## $actor
    ##  [1] 1 2 3 4 5 6 7 1 2 3 4 5 6 7 1 2 3 4 5 6 7 1 2 3 4 5 6 7
    ## 
    ## $side
    ##  [1] 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2
    ## 
    ## $cond
    ##  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1

here the 18 is the sum for each individual; doesnt have to be hard
coded.

``` r
f11.6 = alist(
  left_pulls ~ dbinom( 18 , p ) ,
  logit(p) <- a[actor] + b[treatment] ,
  a[actor] ~ dnorm( 0 , 1.5 ) ,
  b[treatment] ~ dnorm( 0 , 0.5 )
  )
m11.6 <- ulam(f11.6,data=dat , chains=4 , log_lik=TRUE )
```

    ## Running /Library/Frameworks/R.framework/Resources/bin/R CMD SHLIB foo.c
    ## clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/Rcpp/include/"  -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/"  -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/unsupported"  -I"/Users/matthewmule/Library/R/3.5/library/BH/include" -I"/Users/matthewmule/Library/R/3.5/library/StanHeaders/include/src/"  -I"/Users/matthewmule/Library/R/3.5/library/StanHeaders/include/"  -I"/Users/matthewmule/Library/R/3.5/library/rstan/include" -DEIGEN_NO_DEBUG  -D_REENTRANT  -DBOOST_DISABLE_ASSERTS -DBOOST_PENDING_INTEGER_LOG2_HPP -include stan/math/prim/mat/fun/Eigen.hpp   -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -I/usr/local/include   -fPIC  -Wall -g -O2  -c foo.c -o foo.o
    ## In file included from <built-in>:1:
    ## In file included from /Users/matthewmule/Library/R/3.5/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:13:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Dense:1:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Core:88:
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:613:1: error: unknown type name 'namespace'
    ## namespace Eigen {
    ## ^
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:613:16: error: expected ';' after top level declarator
    ## namespace Eigen {
    ##                ^
    ##                ;
    ## In file included from <built-in>:1:
    ## In file included from /Users/matthewmule/Library/R/3.5/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:13:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Dense:1:
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Core:96:10: fatal error: 'complex' file not found
    ## #include <complex>
    ##          ^~~~~~~~~
    ## 3 errors generated.
    ## make: *** [foo.o] Error 1
    ## 
    ## SAMPLING FOR MODEL '79cee3368c3337e807538171c896ad26' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 2.2e-05 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
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
    ## Chain 1:  Elapsed Time: 0.05911 seconds (Warm-up)
    ## Chain 1:                0.058062 seconds (Sampling)
    ## Chain 1:                0.117172 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL '79cee3368c3337e807538171c896ad26' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 1e-05 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.1 seconds.
    ## Chain 2: Adjust your expectations accordingly!
    ## Chain 2: 
    ## Chain 2: 
    ## Chain 2: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 2: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 2: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 2: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 2: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 2: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 2: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 2: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 2: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 2: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 2: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 2: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 2: 
    ## Chain 2:  Elapsed Time: 0.061813 seconds (Warm-up)
    ## Chain 2:                0.050896 seconds (Sampling)
    ## Chain 2:                0.112709 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL '79cee3368c3337e807538171c896ad26' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 9e-06 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
    ## Chain 3: Adjust your expectations accordingly!
    ## Chain 3: 
    ## Chain 3: 
    ## Chain 3: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 3: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 3: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 3: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 3: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 3: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 3: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 3: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 3: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 3: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 3: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 3: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 3: 
    ## Chain 3:  Elapsed Time: 0.058207 seconds (Warm-up)
    ## Chain 3:                0.061188 seconds (Sampling)
    ## Chain 3:                0.119395 seconds (Total)
    ## Chain 3: 
    ## 
    ## SAMPLING FOR MODEL '79cee3368c3337e807538171c896ad26' NOW (CHAIN 4).
    ## Chain 4: 
    ## Chain 4: Gradient evaluation took 9e-06 seconds
    ## Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
    ## Chain 4: Adjust your expectations accordingly!
    ## Chain 4: 
    ## Chain 4: 
    ## Chain 4: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 4: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 4: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 4: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 4: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 4: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 4: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 4: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 4: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 4: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 4: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 4: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 4: 
    ## Chain 4:  Elapsed Time: 0.052763 seconds (Warm-up)
    ## Chain 4:                0.052984 seconds (Sampling)
    ## Chain 4:                0.105747 seconds (Total)
    ## Chain 4:

``` r
precis(m11.6, depth = 2)
```

    ##             mean        sd        5.5%       94.5%     n_eff     Rhat4
    ## a[1] -0.45566843 0.3285378 -0.99894864  0.05226855  611.0715 1.0009858
    ## a[2]  3.87384548 0.7423892  2.76072113  5.19184124 1426.8074 0.9996103
    ## a[3] -0.76472427 0.3380879 -1.31119996 -0.23881974  772.9595 1.0022930
    ## a[4] -0.75616625 0.3487291 -1.33549890 -0.23290084  650.6202 1.0066766
    ## a[5] -0.46983253 0.3349363 -1.02586959  0.06349543  657.3546 1.0034095
    ## a[6]  0.45585377 0.3353734 -0.08085424  0.99417729  657.9981 1.0050892
    ## a[7]  1.93502783 0.4086033  1.27074640  2.60738870 1000.3562 1.0019003
    ## b[1] -0.02534619 0.2803133 -0.46975028  0.42789439  561.0507 1.0039193
    ## b[2]  0.49344655 0.2907439  0.03240391  0.96437497  556.9953 1.0043014
    ## b[3] -0.36540258 0.2886930 -0.83242753  0.09660509  634.4884 1.0035907
    ## b[4]  0.39238967 0.2830651 -0.04822704  0.85151421  652.5963 1.0035876

``` r
precis(m11.4, depth = 2)
```

    ##             mean        sd        5.5%       94.5%     n_eff     Rhat4
    ## a[1] -0.45642629 0.3477435 -1.01318764  0.09523512  614.8645 1.0000542
    ## a[2]  3.88682837 0.7227322  2.80040209  5.11843419 1471.2153 1.0009275
    ## a[3] -0.75063656 0.3363033 -1.31493135 -0.23030060  648.5380 1.0015048
    ## a[4] -0.76124169 0.3421155 -1.30706452 -0.20818540  566.1462 1.0025290
    ## a[5] -0.46091281 0.3369121 -1.00762028  0.06616788  631.0298 0.9992719
    ## a[6]  0.46639103 0.3419080 -0.08814385  1.00046011  543.3119 1.0016696
    ## a[7]  1.94304138 0.4348036  1.25204631  2.68025204  821.7180 1.0014386
    ## b[1] -0.03074828 0.2949858 -0.49803371  0.45067799  565.9269 1.0013012
    ## b[2]  0.48927316 0.2919955  0.01573074  0.95312804  516.1493 1.0015648
    ## b[3] -0.37039709 0.2929954 -0.83303525  0.09380429  542.2040 0.9999995
    ## b[4]  0.37756874 0.2932367 -0.10293106  0.84357819  557.5963 1.0000941

Above we get the same estimates, but since the aggregated model uses a
term for the different ways the counts could be realized, information
criteria is not comparable between the aggregated binomial binomial
models

``` r
## R code 11.26
compare( m11.6 , m11.4 , func=PSIS )
```

    ## Warning in compare(m11.6, m11.4, func = PSIS): Different numbers of observations found for at least two models.
    ## Model comparison is valid only for models fit to exactly the same observations.
    ## Number of observations for each model:
    ## m11.6 28 
    ## m11.4 504

    ##           PSIS        SE    dPSIS      dSE    pPSIS       weight
    ## m11.6 114.0425  8.483187   0.0000       NA 8.311594 1.000000e+00
    ## m11.4 532.2816 18.917739 418.2391 9.488178 8.511503 1.515433e-91

this is because there are more possible ways to see the aggregated data
(see below). This difference doesn’t change inference.

``` r
## R code 11.27
# deviance of aggregated 6-in-9
-2*dbinom(6,9,0.2,log=TRUE)
```

    ## [1] 11.79048

``` r
# deviance of dis-aggregated
-2*sum(dbern(c(1,1,1,1,1,1,0,0,0),0.2,log=TRUE))
```

    ## [1] 20.65212

The WAIC and psis penaltys have influential vlaues in the aggregated
models because instead of calculating something akin to leave one out
cross validation, they are calculating leave 18 out in a sense.

The aggregated format can be better for multilevel
models.

<!-- ```{r} -->

<!-- ## R code 11.18 -->

<!-- d$side <- d$prosoc_left + 1 # right 1, left 2 -->

<!-- d$cond <- d$condition + 1 # no partner 1, partner 2 -->

<!-- ## R code 11.19 -->

<!-- dat_list2 <- list( -->

<!--     pulled_left = d$pulled_left, -->

<!--     actor = d$actor, -->

<!--     side = d$side, -->

<!--     cond = d$cond ) -->

<!-- m11.5 <- ulam( -->

<!--     alist( -->

<!--         pulled_left ~ dbinom( 1 , p ) , -->

<!--         logit(p) <- a[actor] + bs[side] + bc[cond] , -->

<!--         a[actor] ~ dnorm( 0 , 1.5 ), -->

<!--         bs[side] ~ dnorm( 0 , 0.5 ), -->

<!--         bc[cond] ~ dnorm( 0 , 0.5 ) -->

<!--     ) , data=dat_list2 , chains=4 , log_lik=TRUE ) -->

<!-- ## R code 11.20 -->

<!-- compare( m11.5 , m11.4 , func=PSIS ) -->

<!-- ## R code 11.21 -->

<!-- post <- extract.samples( m11.4 , clean=FALSE ) -->

<!-- str(post) -->

<!-- ## R code 11.22 -->

<!-- m11.4_stan_code <- stancode(m11.4) -->

<!-- m11.4_stan <- stan( model_code=m11.4_stan_code , data=dat_list , chains=4 ) -->

<!-- compare( m11.4_stan , m11.4 ) -->

<!-- ``` -->

### Aggregated binmoial model - Berkley admissions.

This is a small dataset with department, acceptance rate, applicants and
sex.

``` r
## R code 11.28
suppressMessages(library(rethinking))
data(UCBadmit)
d <- UCBadmit
d
```

    ##    dept applicant.gender admit reject applications
    ## 1     A             male   512    313          825
    ## 2     A           female    89     19          108
    ## 3     B             male   353    207          560
    ## 4     B           female    17      8           25
    ## 5     C             male   120    205          325
    ## 6     C           female   202    391          593
    ## 7     D             male   138    279          417
    ## 8     D           female   131    244          375
    ## 9     E             male    53    138          191
    ## 10    E           female    94    299          393
    ## 11    F             male    22    351          373
    ## 12    F           female    24    317          341

These 12 rows are a summary of ~4500 applications

``` r
sum(d$applications)
```

    ## [1] 4526

### A model of admit ~ sex.

we make a index variable for sex, and model admit status as a function
of sex; here we are saying a indicator term indexed over levels of sex
predicts admit status. Admit is binomially distributed with probability
p with total count “applications”.

``` r
## R code 11.29
dat_list <- list(
    admit = d$admit,
    applications = d$applications,
    gid = ifelse( d$applicant.gender=="male" , 1 , 2 )
)
# formula 
f11.7 = alist(
  admit ~ dbinom( applications , p ) ,
  logit(p) <- a[gid],
  a[gid] ~ dnorm( 0 , 1.5 )
)
# fit model 
m11.7 <- ulam(f11.7 , data=dat_list , chains=4 )
```

    ## Running /Library/Frameworks/R.framework/Resources/bin/R CMD SHLIB foo.c
    ## clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/Rcpp/include/"  -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/"  -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/unsupported"  -I"/Users/matthewmule/Library/R/3.5/library/BH/include" -I"/Users/matthewmule/Library/R/3.5/library/StanHeaders/include/src/"  -I"/Users/matthewmule/Library/R/3.5/library/StanHeaders/include/"  -I"/Users/matthewmule/Library/R/3.5/library/rstan/include" -DEIGEN_NO_DEBUG  -D_REENTRANT  -DBOOST_DISABLE_ASSERTS -DBOOST_PENDING_INTEGER_LOG2_HPP -include stan/math/prim/mat/fun/Eigen.hpp   -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -I/usr/local/include   -fPIC  -Wall -g -O2  -c foo.c -o foo.o
    ## In file included from <built-in>:1:
    ## In file included from /Users/matthewmule/Library/R/3.5/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:13:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Dense:1:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Core:88:
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:613:1: error: unknown type name 'namespace'
    ## namespace Eigen {
    ## ^
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:613:16: error: expected ';' after top level declarator
    ## namespace Eigen {
    ##                ^
    ##                ;
    ## In file included from <built-in>:1:
    ## In file included from /Users/matthewmule/Library/R/3.5/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:13:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Dense:1:
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Core:96:10: fatal error: 'complex' file not found
    ## #include <complex>
    ##          ^~~~~~~~~
    ## 3 errors generated.
    ## make: *** [foo.o] Error 1
    ## 
    ## SAMPLING FOR MODEL '92bf1800293ab05bf80648b931ab8346' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 1.8e-05 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.18 seconds.
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
    ## Chain 1:  Elapsed Time: 0.016544 seconds (Warm-up)
    ## Chain 1:                0.016084 seconds (Sampling)
    ## Chain 1:                0.032628 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL '92bf1800293ab05bf80648b931ab8346' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 7e-06 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.07 seconds.
    ## Chain 2: Adjust your expectations accordingly!
    ## Chain 2: 
    ## Chain 2: 
    ## Chain 2: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 2: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 2: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 2: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 2: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 2: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 2: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 2: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 2: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 2: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 2: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 2: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 2: 
    ## Chain 2:  Elapsed Time: 0.016289 seconds (Warm-up)
    ## Chain 2:                0.016817 seconds (Sampling)
    ## Chain 2:                0.033106 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL '92bf1800293ab05bf80648b931ab8346' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 6e-06 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.
    ## Chain 3: Adjust your expectations accordingly!
    ## Chain 3: 
    ## Chain 3: 
    ## Chain 3: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 3: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 3: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 3: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 3: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 3: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 3: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 3: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 3: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 3: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 3: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 3: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 3: 
    ## Chain 3:  Elapsed Time: 0.018473 seconds (Warm-up)
    ## Chain 3:                0.015829 seconds (Sampling)
    ## Chain 3:                0.034302 seconds (Total)
    ## Chain 3: 
    ## 
    ## SAMPLING FOR MODEL '92bf1800293ab05bf80648b931ab8346' NOW (CHAIN 4).
    ## Chain 4: 
    ## Chain 4: Gradient evaluation took 7e-06 seconds
    ## Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.07 seconds.
    ## Chain 4: Adjust your expectations accordingly!
    ## Chain 4: 
    ## Chain 4: 
    ## Chain 4: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 4: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 4: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 4: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 4: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 4: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 4: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 4: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 4: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 4: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 4: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 4: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 4: 
    ## Chain 4:  Elapsed Time: 0.016636 seconds (Warm-up)
    ## Chain 4:                0.016723 seconds (Sampling)
    ## Chain 4:                0.033359 seconds (Total)
    ## Chain 4:

``` r
precis( m11.7 , depth=2 )
```

    ##            mean         sd       5.5%      94.5%    n_eff     Rhat4
    ## a[1] -0.2211643 0.04028955 -0.2868443 -0.1590310 1187.649 0.9993454
    ## a[2] -0.8304286 0.04919629 -0.9085865 -0.7509627 1438.134 1.0010663

The mean for probability of admisssion on the log odds scale is higher
for males. Compute the contrast effect sizes on the logit scale and
outcome scale.

The model was fit on the log odds (logit) scale so the contrast on the
posterior tells the log odds difference between male and female, below:
`diff_a`. If we want to *back transform that to a probability scale* we
take the inverse logit transform to get the difference in probability of
admission `diff_p`

``` r
## R code 11.30
post <- extract.samples(m11.7)
diff_a <- post$a[,1] - post$a[,2]
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis( list( diff_a=diff_a , diff_p=diff_p ) )
```

    ##             mean         sd      5.5%     94.5%   histogram
    ## diff_a 0.6092642 0.06296702 0.5045635 0.7073150  ▁▁▁▂▇▇▅▂▁▁
    ## diff_p 0.1413006 0.01425302 0.1178644 0.1638095 ▁▁▁▂▃▇▇▅▂▁▁

``` r
par(mfrow = c(1,2))
hist(diff_a , col = r, breaks = 50, xlab = "M-F difference in log odds");
hist(diff_p, col = r, breaks = 50, xlab = "M-F difference in prob admit (inv_logit)")
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

``` r
## R code 11.31
postcheck( m11.7 )
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

### Fitting the model with another index variable for department.

``` r
## R code 11.32
dat_list$dept_id <- rep(1:6,each=2)
dat_list
```

    ## $admit
    ##  [1] 512  89 353  17 120 202 138 131  53  94  22  24
    ## 
    ## $applications
    ##  [1] 825 108 560  25 325 593 417 375 191 393 373 341
    ## 
    ## $gid
    ##  [1] 1 2 1 2 1 2 1 2 1 2 1 2
    ## 
    ## $dept_id
    ##  [1] 1 1 2 2 3 3 4 4 5 5 6 6

``` r
f11.8 = alist(
  admit ~ dbinom(applications, p ),
  # likelihood 
  logit(p) <- a[gid] + delta[dept_id], 
  # priors 
  a[gid] ~ dnorm( 0 , 1.5 ),
  delta[dept_id] ~ dnorm( 0 , 1.5 )
)

m11.8 <- ulam(flist = f11.8, data=dat_list , chains=4 , iter=4000 )
```

    ## Running /Library/Frameworks/R.framework/Resources/bin/R CMD SHLIB foo.c
    ## clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/Rcpp/include/"  -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/"  -I"/Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/unsupported"  -I"/Users/matthewmule/Library/R/3.5/library/BH/include" -I"/Users/matthewmule/Library/R/3.5/library/StanHeaders/include/src/"  -I"/Users/matthewmule/Library/R/3.5/library/StanHeaders/include/"  -I"/Users/matthewmule/Library/R/3.5/library/rstan/include" -DEIGEN_NO_DEBUG  -D_REENTRANT  -DBOOST_DISABLE_ASSERTS -DBOOST_PENDING_INTEGER_LOG2_HPP -include stan/math/prim/mat/fun/Eigen.hpp   -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -I/usr/local/include   -fPIC  -Wall -g -O2  -c foo.c -o foo.o
    ## In file included from <built-in>:1:
    ## In file included from /Users/matthewmule/Library/R/3.5/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:13:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Dense:1:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Core:88:
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:613:1: error: unknown type name 'namespace'
    ## namespace Eigen {
    ## ^
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:613:16: error: expected ';' after top level declarator
    ## namespace Eigen {
    ##                ^
    ##                ;
    ## In file included from <built-in>:1:
    ## In file included from /Users/matthewmule/Library/R/3.5/library/StanHeaders/include/stan/math/prim/mat/fun/Eigen.hpp:13:
    ## In file included from /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Dense:1:
    ## /Library/Frameworks/R.framework/Versions/3.5/Resources/library/RcppEigen/include/Eigen/Core:96:10: fatal error: 'complex' file not found
    ## #include <complex>
    ##          ^~~~~~~~~
    ## 3 errors generated.
    ## make: *** [foo.o] Error 1
    ## 
    ## SAMPLING FOR MODEL 'eefef2dd38cd530185967a7d0d74fe47' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 1.8e-05 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.18 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4000 [  0%]  (Warmup)
    ## Chain 1: Iteration:  400 / 4000 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  800 / 4000 [ 20%]  (Warmup)
    ## Chain 1: Iteration: 1200 / 4000 [ 30%]  (Warmup)
    ## Chain 1: Iteration: 1600 / 4000 [ 40%]  (Warmup)
    ## Chain 1: Iteration: 2000 / 4000 [ 50%]  (Warmup)
    ## Chain 1: Iteration: 2001 / 4000 [ 50%]  (Sampling)
    ## Chain 1: Iteration: 2400 / 4000 [ 60%]  (Sampling)
    ## Chain 1: Iteration: 2800 / 4000 [ 70%]  (Sampling)
    ## Chain 1: Iteration: 3200 / 4000 [ 80%]  (Sampling)
    ## Chain 1: Iteration: 3600 / 4000 [ 90%]  (Sampling)
    ## Chain 1: Iteration: 4000 / 4000 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 0.400316 seconds (Warm-up)
    ## Chain 1:                0.396502 seconds (Sampling)
    ## Chain 1:                0.796818 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'eefef2dd38cd530185967a7d0d74fe47' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 9e-06 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
    ## Chain 2: Adjust your expectations accordingly!
    ## Chain 2: 
    ## Chain 2: 
    ## Chain 2: Iteration:    1 / 4000 [  0%]  (Warmup)
    ## Chain 2: Iteration:  400 / 4000 [ 10%]  (Warmup)
    ## Chain 2: Iteration:  800 / 4000 [ 20%]  (Warmup)
    ## Chain 2: Iteration: 1200 / 4000 [ 30%]  (Warmup)
    ## Chain 2: Iteration: 1600 / 4000 [ 40%]  (Warmup)
    ## Chain 2: Iteration: 2000 / 4000 [ 50%]  (Warmup)
    ## Chain 2: Iteration: 2001 / 4000 [ 50%]  (Sampling)
    ## Chain 2: Iteration: 2400 / 4000 [ 60%]  (Sampling)
    ## Chain 2: Iteration: 2800 / 4000 [ 70%]  (Sampling)
    ## Chain 2: Iteration: 3200 / 4000 [ 80%]  (Sampling)
    ## Chain 2: Iteration: 3600 / 4000 [ 90%]  (Sampling)
    ## Chain 2: Iteration: 4000 / 4000 [100%]  (Sampling)
    ## Chain 2: 
    ## Chain 2:  Elapsed Time: 0.397734 seconds (Warm-up)
    ## Chain 2:                0.409652 seconds (Sampling)
    ## Chain 2:                0.807386 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL 'eefef2dd38cd530185967a7d0d74fe47' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 7e-06 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.07 seconds.
    ## Chain 3: Adjust your expectations accordingly!
    ## Chain 3: 
    ## Chain 3: 
    ## Chain 3: Iteration:    1 / 4000 [  0%]  (Warmup)
    ## Chain 3: Iteration:  400 / 4000 [ 10%]  (Warmup)
    ## Chain 3: Iteration:  800 / 4000 [ 20%]  (Warmup)
    ## Chain 3: Iteration: 1200 / 4000 [ 30%]  (Warmup)
    ## Chain 3: Iteration: 1600 / 4000 [ 40%]  (Warmup)
    ## Chain 3: Iteration: 2000 / 4000 [ 50%]  (Warmup)
    ## Chain 3: Iteration: 2001 / 4000 [ 50%]  (Sampling)
    ## Chain 3: Iteration: 2400 / 4000 [ 60%]  (Sampling)
    ## Chain 3: Iteration: 2800 / 4000 [ 70%]  (Sampling)
    ## Chain 3: Iteration: 3200 / 4000 [ 80%]  (Sampling)
    ## Chain 3: Iteration: 3600 / 4000 [ 90%]  (Sampling)
    ## Chain 3: Iteration: 4000 / 4000 [100%]  (Sampling)
    ## Chain 3: 
    ## Chain 3:  Elapsed Time: 0.403313 seconds (Warm-up)
    ## Chain 3:                0.400877 seconds (Sampling)
    ## Chain 3:                0.80419 seconds (Total)
    ## Chain 3: 
    ## 
    ## SAMPLING FOR MODEL 'eefef2dd38cd530185967a7d0d74fe47' NOW (CHAIN 4).
    ## Chain 4: 
    ## Chain 4: Gradient evaluation took 8e-06 seconds
    ## Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.08 seconds.
    ## Chain 4: Adjust your expectations accordingly!
    ## Chain 4: 
    ## Chain 4: 
    ## Chain 4: Iteration:    1 / 4000 [  0%]  (Warmup)
    ## Chain 4: Iteration:  400 / 4000 [ 10%]  (Warmup)
    ## Chain 4: Iteration:  800 / 4000 [ 20%]  (Warmup)
    ## Chain 4: Iteration: 1200 / 4000 [ 30%]  (Warmup)
    ## Chain 4: Iteration: 1600 / 4000 [ 40%]  (Warmup)
    ## Chain 4: Iteration: 2000 / 4000 [ 50%]  (Warmup)
    ## Chain 4: Iteration: 2001 / 4000 [ 50%]  (Sampling)
    ## Chain 4: Iteration: 2400 / 4000 [ 60%]  (Sampling)
    ## Chain 4: Iteration: 2800 / 4000 [ 70%]  (Sampling)
    ## Chain 4: Iteration: 3200 / 4000 [ 80%]  (Sampling)
    ## Chain 4: Iteration: 3600 / 4000 [ 90%]  (Sampling)
    ## Chain 4: Iteration: 4000 / 4000 [100%]  (Sampling)
    ## Chain 4: 
    ## Chain 4:  Elapsed Time: 0.408898 seconds (Warm-up)
    ## Chain 4:                0.40813 seconds (Sampling)
    ## Chain 4:                0.817028 seconds (Total)
    ## Chain 4:

``` r
precis( m11.8 , depth=2 )
```

    ##                mean        sd       5.5%      94.5%    n_eff    Rhat4
    ## a[1]     -0.5343631 0.5133460 -1.3508307  0.3056072 638.1812 1.007069
    ## a[2]     -0.4389693 0.5129825 -1.2447279  0.3997600 631.4575 1.007397
    ## delta[1]  1.1146057 0.5158696  0.2629082  1.9328841 651.8555 1.007085
    ## delta[2]  1.0696739 0.5182031  0.2250171  1.8937450 640.3082 1.006979
    ## delta[3] -0.1435763 0.5139504 -0.9830603  0.6668315 643.8968 1.006990
    ## delta[4] -0.1759430 0.5161819 -1.0234409  0.6374747 638.7620 1.007313
    ## delta[5] -0.6181643 0.5172998 -1.4637785  0.1979002 645.1350 1.006962
    ## delta[6] -2.1763005 0.5249216 -3.0344782 -1.3475191 687.2707 1.006555

``` r
postcheck(m11.8)
```

![](ch11_part_I_GLMs_and_binomial_regression_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

``` r
## R code 11.33
post <- extract.samples(m11.8)
diff_a <- post$a[,1] - post$a[,2]
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis( list( diff_a=diff_a , diff_p=diff_p ) )
```

    ##               mean         sd        5.5%       94.5%      histogram
    ## diff_a -0.09539386 0.08198346 -0.22651691 0.036964009  ▁▁▁▂▅▇▇▅▂▁▁▁▁
    ## diff_p -0.02133179 0.01858574 -0.05204162 0.008223685 ▁▁▁▂▃▅▇▇▅▂▁▁▁▁

``` r
sessionInfo()
```

    ## R version 3.5.3 Patched (2019-03-11 r77192)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Mojave 10.14.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] forcats_0.4.0        stringr_1.4.0        dplyr_0.8.5         
    ##  [4] purrr_0.3.3          readr_1.3.1          tidyr_1.0.2         
    ##  [7] tibble_2.1.1         tidyverse_1.2.1      rethinking_2.12     
    ## [10] rstan_2.19.3         ggplot2_3.1.1        StanHeaders_2.21.0-1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.1         lubridate_1.7.4    mvtnorm_1.0-10    
    ##  [4] lattice_0.20-38    prettyunits_1.0.2  ps_1.3.0          
    ##  [7] utf8_1.1.4         assertthat_0.2.1   digest_0.6.25     
    ## [10] cellranger_1.1.0   R6_2.4.0           plyr_1.8.4        
    ## [13] backports_1.1.4    stats4_3.5.3       evaluate_0.14     
    ## [16] coda_0.19-2        httr_1.4.0         pillar_1.4.1      
    ## [19] rlang_0.4.5        readxl_1.3.1       lazyeval_0.2.2    
    ## [22] rstudioapi_0.10    callr_3.2.0        rmarkdown_1.13    
    ## [25] labeling_0.3       loo_2.3.1          munsell_0.5.0     
    ## [28] broom_0.5.2        compiler_3.5.3     modelr_0.1.4      
    ## [31] xfun_0.7           pkgconfig_2.0.2    pkgbuild_1.0.3    
    ## [34] shape_1.4.4        htmltools_0.3.6    tidyselect_0.2.5  
    ## [37] gridExtra_2.3      codetools_0.2-16   matrixStats_0.54.0
    ## [40] fansi_0.4.0        crayon_1.3.4       withr_2.1.2       
    ## [43] MASS_7.3-51.1      grid_3.5.3         nlme_3.1-137      
    ## [46] jsonlite_1.6       gtable_0.3.0       lifecycle_0.1.0   
    ## [49] magrittr_2.0.1     scales_1.0.0       cli_1.1.0         
    ## [52] stringi_1.4.3      xml2_1.2.0         generics_0.0.2    
    ## [55] vctrs_0.2.4        tools_3.5.3        glue_1.3.1        
    ## [58] hms_0.4.2          processx_3.3.1     yaml_2.2.0        
    ## [61] inline_0.3.15      colorspace_1.4-1   rvest_0.3.4       
    ## [64] knitr_1.23         haven_2.1.0
