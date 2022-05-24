Ch 11 p set
================

``` r
suppressMessages(library(rethinking))
suppressMessages(library(magrittr))
```

11E1 - log odds of pr = 0.35

``` r
log(0.35  / (1- 0.35))
```

    ## [1] -0.6190392

11E2 - log odds of event = 3.2 what is p

``` r
inv_logit(3.2)
```

    ## [1] 0.9608343

to get a feel for log odds vs p:  
at 50% p, log odds is 0 at 0.1 its around -2 and at 0.9 its around 2.

``` r
pvec = seq(0, 1, by = 0.1)
plot(pvec, log(pvec / (1-pvec)))
```

![](ch11_pset_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

11E3 logistic regression beta = 1.7 what does it say about proportional
change in odds of the event?

logistic regression coefs are in units of log odds, so 1.7 is log odds.
The odds is the exponentiated coefficient (transforms to the outcome
scale)

``` r
exp(1.7)
```

    ## [1] 5.473947

11E4  
Poisson regression requires an offset because the totla counts per class
(exposure) can be different, i.e. sampling one hospital for events for a
year and another for 2 years. Or in scRNAseq the rate parameter for a
gene tells you something but to conrast, the total genes per cell (log
UMI) needs to be accounted for.

11M1  
the likelihood is the probability of the data given the model so in the
case of aggregated binomial globe toss, 3 water one land, that is 3 W in
4 trials, by embedding this into the binomial probability distribution
we are aggregating the the data, the multiplicity term wiht n choose k
\* 1-k is multiplying the number of ways to observe, ie. W L W W | W W W
L | W W L W etc. we then plug in the probability. For individual counts,
we calculate p \* (1- p) for each separately. Basically p is independent
of the likeihood function.

11M2

Poisson regression coefficient 1.7:

y\[i\] = Poisson(lambda\[i\])  
log(lambda\[i\]) = a + Bx\[i\]

the a + bX term is equal to log(lambda)

``` r
exp(1.7)
```

    ## [1] 5.473947

Fr every 1 unit change in the predictor there is a 5.47 unit change in
the outcome variable. *note, in looking up the correct answer to this,
there was a slightly more nuanced answerusing ratio of means exp( a + Bx
+ 1) / exp(a + Bx); the top part of the fraction is the 1 unit of change
in predictor, solving you get exp(B), as i calculated above. but
important to note that the interpretation of this is the proportional
change in the expected count very similar to proportional change in log
odds from logistic regression*

11M3 explain why logit link appropriate for binomial GLM.

Because binomial is mapping the linar model output which is 0 , inf onto
2 possible outcomes, yes, no, heads, tails, 1, 0 etc. The odds makes
sense to use because we can write probability functions in terms of one
outcome like heads then the odds (logit) captures p ( 1-p ).

11M4 explain why log link appropriate for Poisson.

The log link constrains predictions from the linear model onto a space
that will not predict negative counts. When we get a prediction from the
linear model:

log(lambda) = a + Bx

then lambda = exp(a + Bx)  
exp is raising e to the power of the linear model, so say a+ Bx
evaluated to -5, the log link will make lambda at i =

``` r
exp(-5)
```

    ## [1] 0.006737947

meaning basically

``` r
2.71^-5
```

    ## [1] 0.006841535

*It is conventional to use a logit link for a binomial GLM because we
need to map the continuous linear model value to a probability parameter
that is bounded between zero and one.*

*It is conventional to use a log link for a Poisson GLM because we need
to map the continuous linear model value to a mean that must be
positive.*

11M5 could it make sense to use a logit link with a poisson GLM?  
using a logit link makes the outcome bounded between 0 and 1.  
basically if you want to put a upper bound on the estimates. see manual
ed 1 p73. One could potentially do this – a cell does have an upper
bound as to what the total counts could be, so we would expect the
counts of individual genes to be some fraction of that total; we’d have
to define this for eahc cell type and its more a function of the
technology used.

M116 State constraints under which poisson and binomial have maximum
entropy

The binomial distribution has maximum entropy for a binary outcome with
a constant expected rate, for example a globe toss with expected
probability of water 0.7 what is the likelihood of observing 5 waters in
a row – the likeliood function with the least assumptions is the one
that can happen the most ways which is the binomial.

A Poisson distribution has maximum entropy under the constraint that
there is a constant expected rate, as the binomial with low probability
of observing a “success” (count) in each trial. its a generalization of
the binomial.

11M8 Kline data without hawaii- fit an interaction model

likelihood for interaction model with a category and continuous var:  
y = a\[group\] + B\[group\] \* predictor

embed that in the glm using the link function on outcome  
log( y ) = …

``` r
data("Kline")
d = Kline
d = d[!d$culture == 'Hawaii', ]
d$pop = log(d$population)


# fit mcmc binomial model 
f1.1 = alist(
  # estimand
  total_tools ~ dpois(lambda), 
  # likelihood
  log(lambda) <- a[cid] + B[cid]*pop, 
  # priors 
  a[cid] ~ dnorm(3, 0.3), 
  B[cid] ~ dnorm(0, 0.1)
)

# specify data for stan
dat = list(
  total_tools = d$total_tools, 
  pop = d$pop, 
  cid = d$contact
)

# fit model 
m1.1 = ulam(flist = f1.1, data = dat,cores = 4)
```

    ## Trying to compile a simple C file

``` r
# coefficients
precis(m1.1, depth = 2)
```

Important – the model converts the ‘high’, ‘low’ factor to indicators
for contact ID; this screwed up the ability to use ‘high low’ in the
step idx = which(<m1.1@data$cid> ==ct). So probably best to use the data
that is stored in the stan model fit when using link.

``` r
suppressMessages(library(tidyverse))

# extract posterior 
post = extract.samples(m1.1)
m1.1@data
```

    ## $total_tools
    ## [1] 13 22 24 43 33 19 40 28 55
    ## 
    ## $pop
    ## [1] 7.003065 7.313220 8.188689 8.474494 8.909235 8.987197 9.126959 9.472705
    ## [9] 9.769956
    ## 
    ## $cid
    ## [1] 2 2 2 1 1 1 1 2 1

``` r
pop.seq = seq(min(d$pop), max(d$pop), length.out = 45)
pop.lim = c(min(d$pop), max(d$pop))
tool.seq = c(min(d$total_tools), max(d$total_tools))
## 

par(mfrow=c(1,2))
for (ct in c(1,2)) {
  
  idx = which(m1.1@data$cid ==ct)
  
  # compute the posterior for the levels of contact ID 
  mu = link(m1.1, data = data.frame(cid = ct, pop = pop.seq))
  plot( d$pop[idx], d$total_tools[idx], xlim = pop.lim, ylim = tool.seq)
  for (i in 1:20) {
    lines(x = pop.seq, y = mu[i, ], col = col.alpha('grey', 0.4))
  }
}
```

![](ch11_pset_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

another way to visualize the interaction

``` r
pop.seq = c(min(d$pop), max(d$pop))
tool.seq = c(min(d$total_tools), max(d$total_tools))

mu1 = link(m1.1,data = data.frame(cid = 1 , pop = pop.seq))
mu1.mean = apply(mu1, 2, mean)
mu1.pi = apply(mu1, 2, PI)

mu2 = link(m1.1,data = data.frame(cid = 2 , pop = pop.seq))
mu2.mean = apply(mu2, 2, mean)
mu2.pi = apply(mu2, 2, PI)

plot(dat$pop, dat$total_tools, col = rangi2, pch = 16)
lines(pop.seq, mu1.mean, col = 'grey')
shade(mu1.pi, lim = pop.lim)
lines(pop.seq, mu2.mean, col = 'black')
shade(mu2.pi, lim = pop.lim)
```

![](ch11_pset_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

compare to frequentist glm fit

``` r
# fit frequentist binomial glm with log link

d$con.id = ifelse(d$contact == 'low', 0, 1) 
f1 = total_tools ~  pop + con.id + pop:con.id
m2 = glm(formula = f1, family = poisson(link = 'log'),data = d)
precis(m2)
```

    ##                   mean        sd        5.5%     94.5%
    ## (Intercept)  1.3188919 0.8885545 -0.10118977 2.7389735
    ## pop          0.2174400 0.1075242  0.04559558 0.3892844
    ## con.id      -0.6938106 1.7735513 -3.52828812 2.1406670
    ## pop:con.id   0.1142124 0.1996576 -0.20487898 0.4333037

Interpretation of coefficients:  
(intercept) - Y intercept for the baseline reference group ‘low
contact’  
1.32  
b1 = pop = slope or the baseline group  
0.22

regression equation for low contact reference group  
y = 1.32 + (0.22 \* xi)

b2 = conid = offset for the y intercepy for the alternate group ‘high
contact’  
\-0.69  
b3 = offset to the slope of the alternative group compared to the
reference group high vs low  
0.11

regression equation for the low contact group  
y = (1.32 + -0.69 ) + (0.11 + 0.22 \*
    xi)

``` r
precis(m1.1, depth = 2)
```

    ##            mean         sd        5.5%      94.5%    n_eff     Rhat4
    ## a[1] 2.96182176 0.28072949  2.50262723 3.39491841 227.4898 1.0061881
    ## a[2] 2.86326689 0.27369350  2.45207797 3.31436592 196.1047 0.9999053
    ## B[1] 0.07427456 0.03150094  0.02381760 0.12431253 240.7092 1.0056578
    ## B[2] 0.02829464 0.03588428 -0.03456683 0.08395578 184.4798 0.9998972

``` r
par(mfrow = c(2,2))

plot(
  exp(predict(m2, newdata = data.frame(pop =  d$pop, con.id = 1))), 
  pch = 16, col = 'blue')
plot(
  exp(predict(m2, newdata = data.frame(pop =  d$pop, con.id = 0))), 
  pch = 16, col = 'blue')

plot(
  exp(mean(post$a[,1]) + mean(post$B[ ,1])*dat$pop), 
  pch = 16, col = rangi2
  )

plot(
  exp(mean(post$a[,2]) + mean(post$B[ ,2])*dat$pop), 
  pch = 16, col = rangi2
  )
```

![](ch11_pset_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

11H3

was struggling to get ulam to fit this prob bc of bad priors.

``` r
suppressMessages(library(rethinking))
suppressMessages(library(tidyverse))

data(salamanders)
d = salamanders

d$s = d$SALAMAN
d$p = (d$PCTCOVER - mean(d$PCTCOVER) ) / sd(d$PCTCOVER)
d$a = (d$FORESTAGE - mean(d$FORESTAGE)) / sd(d$FORESTAGE)

dlist = list( 
  sal = d$s, 
  pct = d$p,
  age = d$a
  )

pairs(dlist, col = rangi2, pch = 16)
```

![](ch11_pset_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

it looks like old growth forests have increased ground cover based on
pct ~ age. Once the forest hits a certain percentage the salmon
population can take off (it becomes more like gaussian dist).

There is a direct causal effect maybe of forest age on salamander, but
also an effect mediated through percent coverage.

``` r
f1 <- alist(
  sal ~ dpois(lambda), 
  log(lambda) <- a + b*pct,
  b ~ dnorm(0, 1),
  a ~ dnorm(0, sd = 1)
)

m1 <- ulam(flist = f1, data = dlist, chains = 4, cores = 4)
```

    ## Trying to compile a simple C file

``` r
trankplot(m1)
```

![](ch11_pset_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
precis(m1)
```

    ##       mean        sd      5.5%     94.5%    n_eff    Rhat4
    ## b 1.137885 0.1724942 0.8739675 1.4233835 601.7771 1.003465
    ## a 0.431753 0.1417133 0.2031508 0.6534903 547.0473 1.002759

``` r
# extract parameters 
post <- extract.samples(m1, n = 50)
head(post)
```

    ## $b
    ##  [1] 0.9555929 1.2825958 1.5468328 1.0634902 1.6107719 1.2803444 0.9223248
    ##  [8] 1.2085722 1.1313128 1.2252885 1.1511575 1.1768188 0.8740236 1.4068096
    ## [15] 1.3513733 1.4732654 1.0242242 1.2719923 1.1255519 1.2106576 0.9422423
    ## [22] 1.1671464 0.9763960 1.4933063 1.2223957 1.4298885 1.0049393 1.1589405
    ## [29] 1.2088814 0.9623257 1.1260404 1.2393074 1.3496300 0.9291431 1.2206181
    ## [36] 0.9812968 1.0723348 0.8271347 1.1728695 1.2393730 0.8634603 1.4524782
    ## [43] 1.0545215 0.9906052 1.0457508 0.9260136 1.2658268 1.2735204 1.3237081
    ## [50] 1.2750220
    ## 
    ## $a
    ##  [1] 0.50615166 0.42899471 0.10853533 0.42332302 0.32180233 0.39483531
    ##  [7] 0.69577035 0.38051789 0.38173482 0.53841695 0.53263463 0.22184067
    ## [13] 0.55313144 0.17025498 0.38720382 0.20839549 0.58550593 0.28768349
    ## [19] 0.27809062 0.31335406 0.60119497 0.22354595 0.54150297 0.03655037
    ## [25] 0.36761096 0.11095769 0.53365764 0.40977841 0.27563187 0.56524795
    ## [31] 0.43712756 0.38442131 0.27857972 0.63319300 0.34125209 0.53957136
    ## [37] 0.45018524 0.57289272 0.48164109 0.21482200 0.58337124 0.14387720
    ## [43] 0.53739571 0.57992584 0.42846915 0.40159840 0.43333539 0.43488343
    ## [49] 0.33202454 0.18127904

``` r
xseq = seq(-2,2, length.out = 30)
# summarize estimand based on new predictor values
mu = link(m1, n = 100, data = data.frame(pct = xseq))
mu.mean = apply(mu, 2, mean)
mu.PI = apply(mu, 2, PI)

# simulate uncertainty in estimates based on full model 
dsim = sim(fit = m1, data = data.frame(pct = xseq))
sim.pi = apply(dsim, 2, PI)

# plot summary
plot(sal ~ pct, dlist, col = rangi2, pch = 16)
lines(x = xseq,y = mu.mean, lwd = 3) # mean
shade(object = mu.PI, lim = xseq) # mean uncertainty
shade(object = sim.pi,lim = xseq)
```

![](ch11_pset_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Given the pairs plot the model `log(lambda) <- a + b*pct + b2*age` woudl
prob predict better but would be harder to understand the estimates.
Same prediction vs interpretability trade off but here we can see this
is because the 2 predictors are highly **collinear**. There is a
fundamental different scientific process for undertanding the links
between variables and doing predictive modeling.
