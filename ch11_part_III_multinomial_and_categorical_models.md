Ch 11 part III: Multinomial and categorical models
================

### Multinomial models

Imagine blue / white marble example had 3 colors. WE need a way to
represent more than 2 unordered events with constant probability across
uneven trials. The maximum entropy distribution for this is the
**Multinomial distribution**.

When each event is isolated on a single row this can be called
**categorical regression** and in machine learning it is typically
called **the maximum entropy classifier**.

2 approaches:  
1\. model multinomial likelihood using a logit link  
2\. transform into a series of poisson likelihoods

### The multinomial link function

**multinomial logit** also known as the **softmax function**  
it takes a vector of scores for each K type and computes the probability
of a particular event type:  
p(k | score 1 score 2 … ) = exp( s\[k\] ) / exp ( sum( s\[i\] ) )

In multinomial model we need k -1 linear models for K types of events
and we can also have different predictors for each k

### multinomial – predictors matched to outcomes

model career choice for individuals using expected income. b\*income
appears in each model, but its a different income value for each
career.  
example - 3 careers

``` r
suppressMessages(library(rethinking))
## R code 11.55
# simulate career choices among 500 individuals

N <- 500             # number of individuals
income <- c(1,2,5)   # expected income of each career
score <- 0.5*income  # scores for each career, based on income

# next line converts scores to probabilities
p <- softmax(score[1],score[2],score[3])
p
```

    ## [1] 0.09962365 0.16425163 0.73612472

now simulate choice: outcome career holds event type values, not counts

``` r
# make empty vector of choices for each individual
career <- rep(NA,N)  

# sample chosen career for each individual
set.seed(34302)
for ( i in 1:N ) career[i] <- sample( 1:3 , size=1 , prob=p )
plot(career, cex = 0.5)
```

![](ch11_part_III_multinomial_and_categorical_models_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

### Fit model (stan code)

We use `dcategorical()` for the likelihood function, the multinomial
logistic regression distribution.

the “pivot” is holding one score constant in this case the score for
career 3.

``` r
## R code 11.56
code_m11.13 <- "data{
    int N; // number of individuals
    int K; // number of possible careers
    int career[N]; // outcome
    vector[K] career_income;
}
parameters{
    vector[K-1] a; // intercepts
    real<lower=0> b; // association of income with choice
}
model{
    vector[K] p;
    vector[K] s;
    a ~ normal( 0 , 1 );
    b ~ normal( 0 , 0.5 );
    s[1] = a[1] + b*career_income[1];
    s[2] = a[2] + b*career_income[2];
    s[3] = 0; // pivot
    p = softmax( s );
    career ~ categorical( p );
}
"

# set up data list and invoke stan 
## R code 11.57
dat_list <- list( N=N , K=3 , career=career , career_income=income )
m11.13 <- stan( model_code=code_m11.13 , data=dat_list , chains=4 )
```

    ## Trying to compile a simple C file

    ## Warning: There were 51 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ## Warning: Examine the pairs() plot to diagnose sampling problems

``` r
precis( m11.13 , depth = 2 )
```

Really hard to interpret coefs so below we do a counterfactual
simulation. Compare to a counterfactual career where the income is
changed (by doubling the income of career 2).

``` r
## R code 11.58
post <- extract.samples( m11.13 )

# set up logit scores
s1 <- with( post , a[,1] + b*income[1] )
s2_orig <- with( post , a[,2] + b*income[2] )
s2_new <- with( post , a[,2] + b*income[2]*2 )

par(mfrow = c(1,3))
hist(s1);hist(s2_orig);hist(s2_new)
```

![](ch11_part_III_multinomial_and_categorical_models_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

How much does the probability of career 2 selection increase when we
double its expected income

``` r
# compute probabilities for original and counterfactual
p_orig <- sapply( 1:length(post$b) , function(i)
    softmax( c(s1[i],s2_orig[i],0) ) )
p_new <- sapply( 1:length(post$b) , function(i)
    softmax( c(s1[i],s2_new[i],0) ) )

# summarize
p_diff <- p_new[2,] - p_orig[2,]
precis( p_diff )
```

    ##              mean         sd        5.5%      94.5% histogram
    ## p_diff 0.03655682 0.03167913 0.001975706 0.09680895 ▇▅▂▂▁▁▁▁▁

I get something different from the printed book here which says 0.13.

### multinomial regression with parameters matched to observations

model career choice ~ that persons family income.

``` r
## R code 11.59
N <- 500
# simulate family incomes for each individual
family_income <- runif(N)
# assign a unique coefficient for each type of event
b <- c(-2,0,2)
career <- rep(NA,N)  # empty vector of choices for each individual
for ( i in 1:N ) {
    score <- 0.5*(1:3) + b*family_income[i]
    p <- softmax(score[1],score[2],score[3])
    career[i] <- sample( 1:3 , size=1 , prob=p )
}
```

so this is what we’re
modeling:

``` r
plot(career ~ family_income)
```

![](ch11_part_III_multinomial_and_categorical_models_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

This is how we would build that model. note the loop needed in stan to
make a model for each observation.

``` r
code_m11.14 <- "data{
    int N; // number of observations
    int K; // number of outcome values
    int career[N]; // outcome
    real family_income[N];
}
parameters{
    vector[K-1] a; // intercepts
    vector[K-1] b; // coefficients on family income
}
model{
    vector[K] p;
    vector[K] s;
    a ~ normal(0,1.5);
    b ~ normal(0,1);
    for ( i in 1:N ) {
        for ( j in 1:(K-1) ) s[j] = a[j] + b[j]*family_income[i];
        s[K] = 0; // the pivot
        p = softmax( s );
        career[i] ~ categorical( p );
    }
}
"

dat_list <- list( N=N , K=3 , career=career , family_income=family_income )
m11.14 <- stan( model_code=code_m11.14 , data=dat_list , chains=4 )
# precis( m11.14 , 2 )
```

### multinomial regression as Poisson regression

``` r
## R code 11.60
library(rethinking)
data(UCBadmit)
d <- UCBadmit

## R code 11.61
# binomial model of overall admission probability
fb = alist(
  admit ~ dbinom(applications,p), 
  logit(p) <- a, 
  a ~ dnorm( 0 , 1.5 )
)

m_binom <- quap(flist = fb, data = d)

# Poisson model of overall admission rate and rejection rate
# 'reject' is a reserved word in Stan, cannot use as variable name
dat <- list( 
  admit=d$admit , 
  rej=d$reject 
  )

fl = alist(
  admit ~ dpois(lambda1),  
  rej ~ dpois(lambda2),  
  log(lambda1) <- a1,
  log(lambda2) <- a2,
  c(a1,a2) ~ dnorm(0,1.5)
)
m_pois <- ulam(fl, data=dat , chains=3 , cores=3 )
```

    ## Trying to compile a simple C file

the inferred binomial probability of admission across the dataset is the
inverse function of the log odds.

``` r
## R code 11.62
inv_logit(coef(m_binom))
```

    ##         a 
    ## 0.3878045

in the poisson model padmis = lambda1 / lambda 1 + lambda 2

we exponentiate the coefficients to get the probabiliites

``` r
## R code 11.63
k <- coef(m_pois)
a1 <- k['a1']
a2 <- k['a2']
exp(a1)/(exp(a1)+exp(a2))
```

    ##        a1 
    ## 0.3874246

cool, it’s the same. See overthinking box p 365
