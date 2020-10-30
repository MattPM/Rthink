Chapter 2
================

<!--   output: -->

<!--   html_document: -->

<!--     df_print: paged -->

<!-- editor_options: -->

<!--   chunk_output_type: inline -->

``` r
suppressMessages(library(tidyverse))
suppressMessages(library(rethinking))
```

A likelihood is a distribution function assigned to an observed
variable. A non-Bayesian likelihood is a function of the parameters
conditional only on the observed data l(param|data) . In bayesian stats
it is reasonable to write p(data | param).

Binomial distibution of 6 globe tosses the d in dbinom means density.
the partner functions include dbinom, rbinom - r for random distribution
of a given size, and pbinom for cumulative probability.

``` r
dbinom(x = 5, size = 10, prob = 0.7)
```

    ## [1] 0.1029193

``` r
dbinom(x = 9, size = 10, prob = 0.3)
```

    ## [1] 0.000137781

``` r
#  a random binomial sample 
rbinom(n = 50,size = 10, 0.7)
```

    ##  [1] 4 8 8 8 6 6 6 6 8 8 6 6 6 9 8 7 9 6 6 7 7 9 4 7 7 6 6 9 8 6 9 8 6 5 8
    ## [36] 6 8 7 6 5 6 7 9 6 4 5 7 8 6 5

The binomial distribution is the maximum entropy (Ch 10) way to count
binary events.

The distributions assigned to observed variables here are the parameters
- in the case of the binomial the parameter p was the probability of
sampling a head on a coin toss. Generalize to things like variation
among group average treatment effect from a contrast etc.

For every parameter in the bayesian “golem” I am constructing, I must
provide a distribution of prior plausibility - an initial plausibility
assigmnent for each posible value of the parameter.

Bayes theorm summarized is  
**Posterior = probability of data x prior / average probability of
data**

For the water/land example where W and L represent landing on water r
land from the globe toss:

Pr(p|W,L) = Pr(W,L | p ) \* Pr(p) / Pr(W,L)

the denominator is the confusing part. It is the probabiliity of data
averaged over the prior. The format for the denominator is an
expectation **E()** means take the **expectation** - this average is
commonly called a **marginal likelihood** the whole point of this is to
standardize the posterior so that it sums to 1. In mathematical stats
format:  
E(Pr(W,L|p) -\> this means computing the average over a continuous
distribution of values by calculating the integral: ∫PR(W,L|p)Pr(p)dp  
So the overall equation would be  
Pr(p|W,L) = Pr(W,L | p ) \* Pr(p) / ∫PR(W,L|p)Pr(p)dp

This is all unimportant formal stuff the key is posterior is
proportional to the prior \* p(data) because for each specific value of
p, the number of paths through the garden of forking data is the prior
number of paths times the new number of paths. Average on the bottom
turns it into a probability.

Bayesian analysis is not about Baayes theorm it is about quantifying
uncertainty about parameters and models.

### Conditioning engines

The action/ verb of the motor of the bayesian machine can be thought of
as “conditioning”. Conditioning engines are used for approximating the
conditioning that would be done with integral calculus which is not
posible for most models. **Conditioning engines**  
1\) Grid approximation  
2\) Quadratic Approximation  
3\) Markov Chain Monte Carlo (MCMC)

Grid approximation is kind of like fundamental theorm of calculus in
that we can approximate continuous parameters by considering a finite
grid of parameter values. This does not scale to many parameter models.

### Grid approximation

``` r
# define grid 
p_grid <- seq( from=0 , to=1 , length.out=20 )

# define prior 
prior <- rep( 1 , 20 )

# compute likelihood (6 waters, 9 tosses calculated at each p above in grid)
likelihood <- dbinom(x = 6, size = 9, prob = p_grid)

# compute product of likelihood and prior 
unstd.posterior <- likelihood * prior

# standardize the posterior, so it sums to 1 
posterior <- unstd.posterior / sum(unstd.posterior)

plot( x = p_grid , y = posterior , type="b" , xlab="probability of water" , ylab="posterior probability" )
```

![](ch2_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
plot( x = p_grid , y = likelihood , type="b" , xlab="probability of water" , ylab="Likelihood " )
```

![](ch2_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

### Quadratic approximation

Approximate the posterior distribution with a Gaussian distribution
which can conveniently be described by computing only the mean and
variance. The name Quadratic comes from parabolic, the shape of the
logarithm of a gaussian distribution…

This just uses a optomization algorithm to find the peak of
posterior(the mode) and then estimate the curvature near the peak -
algos do this by knowing the slope at each point. use the **quap**
function.

``` r
globe.qa <- quap(
  alist( 
    # binomial likelihood 
    W ~ dbinom(W+L, p), 
    
    # uniform prior
    p ~ dunif(0,1)
    ),
  # specify input data 
  data=list(W=6,L=3) 
)

# display summary of quadratic approximation
precis( globe.qa )
```

    ##        mean        sd      5.5%     94.5%
    ## p 0.6666669 0.1571337 0.4155369 0.9177969

Assuming a Gaussian posterior, it is maximized at 0.67 with a SD of
0.16.

With more data this gets much better. It can do things like assign
positive probability to p = 1 which we know is false because we observed
both W and L.

Note this is often equivalent to MLE with a uniform prior / lots of data
- see Huber chapter on MLE.

An aside about **Hessians** which appear in error messages: *A Hessian
is a square matrix of second derivatives. In the quadratic approximation
it is second derivatives of the log of posterior probability with
respect to the parameters. It turns out that these derivatives are
sufficient to describe a Gaussian distribution, because the logarithm of
a Gaussian distribution is just a parabola. **Parabolas have no
derivatives beyond the second**, so once we know the center of the
parabola (the posterior mode) and its second derivative, we know
everything about it. And indeed the second derivative (with respect to
the outcome) of the logarithm of a Gaussian distribution is proportional
to its inverse squared standard deviation (its “precision”: page 79). So
knowing the standard deviation tells us everything about its shape. The
standard deviation is typically computed from the Hessian, so computing
the Hessian is nearly always a necessary step. But sometimes the
computation goes wrong, and your golem will choke while trying to
compute the Hessian. In those cases, you have several options. Not all
hope is lost. But for now it’s enough to recognize the term and
associate it with an attempt to find the standard deviation for a
quadratic approximation.*

### Markov Chain Monte Carlo

Instead of computing or infering the posterior distribution directly,
draw samples from the posterior (in proportion to their relative
frequency).

multilevel models do not always allow us to write down a single unified
function forr the posterior distribution. The function to maximize when
computing the maximum a posterori estimation is not known and must be
computed in peices.

``` r
n_samples = 1000 
p = rep( NA , n_samples ) 
p[1] = 0.5
W = 6 
L = 3 

for ( i in 2:n_samples ) {
  p_new <- rnorm( n = 1 ,mean =  p[i-1] , sd = 0.1 )
  if ( p_new < 0 ) p_new <- abs( p_new )
  if ( p_new > 1 ) p_new <- 2 - p_new
  q0 <- dbinom(x =  W , size =  W+L , prob =  p[i-1] )
  q1 <- dbinom( x = W , size = W+L ,prob =  p_new )
  p[i] <- ifelse( runif(1) < q1/q0 , p_new , p[i-1] )
}
```

*The values in p are samples from the posterior distribution. To compare
to the analytical posterior:*

``` r
dens( p , xlim=c(0,1) )
curve( dbeta( x , W+1 , L+1 ) , lty=2 , add=TRUE )
```

![](ch2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

*It’s weird. But it works. I’ll explain this algorithm, the Metropolis
algorithm, in Chapter 9.*

## Summary

The target of Bayesian inference is a posterior probability distribution
= the relative number of ways each conjectured cause of the data could
have produced the data. These relative numbers of ways represent
plausabilities of the conjectures. The plausibilities are updated in
light of new observations by bayesian updating which is multiplying the
prior \* likelihood.

A Bayesian model is a composite of variables and distributional
definitions of the variables. The prior provides plausibility of each
plausible value for the parameters before accounting for the data. The
plausibliities after accounting for the data are computed by Bayes
theorm, resulting in a posterior distribution. These are fit with the 3
options above.