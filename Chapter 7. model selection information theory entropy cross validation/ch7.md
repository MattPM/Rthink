Chapter 7
================

<!--   output: -->

<!--   html_document: -->

<!--     df_print: paged -->

<!-- editor_options: -->

<!--   chunk_output_type: inline -->

``` r
suppressMessages(library(rethinking))
suppressMessages(library(magrittr))
```

This chapter connects the topics below to explain theory behind methods
that enable principled selection of the optimal model from a group of
candidate models.  
Cross validation vs information criteria, information theory,
regularization, overfitting, underfitting, data compression, minimum
description length, bias variance tradeoff, intro to information theory
to motivate maximum entropy probability distributions, Kullback-leiber
divergence, cross entropy, divergence vs. deviance, log pointwise
predictive density, Leave one out cross validation, pareto smoothed
importance sampling (PSIS), pareto distribution, Akaike information
criteria, ceviance information criteria, widely applicable information
criteria (WAIC), effective number of parameters, Bayesian information
criteria, bayes factor, model comparison, model averaging, robust
regression, student t. distribution.

We reframe the bias-variance tradeoff as simply underfitting(bias)
overfitting(variance). Bias variance tradeoff is just models bias in the
relationship it thinks could exist between the predictor and ouctomes, a
linear model is more ‘biased’ in that assumption; variance is the TEST
SET error from the model trained on training data.

A key insight: confounded models often produce better predictions than
non-confounded models. Trying to understand how parameters are wired
together mechanistically and cause an outcome (e.g. through the use of
DAGS or structural equation modeling) is a fundamentally different task
than making predictions. A simple DAG + GLM that we create based on
domain knowledge might help us understand how a system works. This model
will make worse predictions than a highly parameterized deep learning
model. We can use both tools in different
ways.

### Overview of information entropy, divergence vs deviance, predicting out of sample deviance without ‘out of sample’ data:

information theory gives us the shannon entropy, a measure of the
uncertainty in a target distribution; a model predicting something
difficult / hard to predict has higher uncertainty and shannon entropy
measures this precisely as  
H = -∑ pi \* log(pi)  
Accuracy is how far our model is from the true target distribution - we
measure this with something related to entropy , H. If the true
distribution of events occurs with P1 = 0.3 P2 = 0.7 and our model
predicts they happen with Q1 = 0.25 Q2 = 0.75, how much additional
uncertainty have we introduced? This is called the KullbackLeiber
divergence:

Dkl(p,q) = -∑ Pi \* log(Pi / Qi)

derivation:  
Dkl(p,q) = H(p,q) - H(p)  
the term H(p,q) is called the cross entropy; we are using a probability
distribution q t predict events that occur with a different probability
distribution p: the events arise with probability p but are expected
with probability q:  
H(p,q) = -∑ Pi \* log(Qi) \* notice it is Q in the log.

To evaluate a model we need a scoring system (a measure of accuracy)
that is maximized when we know the **true** model that generated the
data. This measure of accuracy turns out to be the joint probability of
the model estimates (multiply each estimated probability), usually
reported in log space. Example from weather prediction– compare 2
measures of accuracy 1) hit rate (percent of time correct over x days),
2) joint probability (the probability of predicting the exact sequence
of weather across the x days). In the book example p 203, a weather
person that always predicts 100% it will be sunny had a higher ‘hit’
rate but ignored the *cost* of predicting rain (this same idea discussed
in nate silvers book). The weather person who makes a decent guess each
day has a higher joint probability of getting the exact sequence.

How do we apply this to model fitting? We need a measure of distance
from the target that acccounts for the fact that some targets are harder
to hit than others. We measure this in units of uncertainty–these units
are in the form of information. **Information is the amount that
uncertainty is reduced when we learn an outcome.**

Now we need a measure of uncertainty that fits 3 criteria:  
1\) it’s continuous 2) it increases as the number of possible events /
was something could happen increases (a moving target as opposed to a
stationary target - 3 dimensions; predicting rain or sun vs rain, sun,
snow, hail)  
3\) It’s additive - it is the sum of the separate uncertainties
reflected in the model e.g. u(rain/shine) + u(hit/cold).  
This function is *Inforrmation Entropy* or Shannon Entropy.

\-∑ pi \* log(pi)

The true probability for rain and shine on mount washington p(rain) =
0.4, p(shine) = 0.6 The true probability for rain and shine in abu dabi
p(rain) = 0.01, p(shine) = 0.99

entropy for each:

``` r
pm = c(0.4, 0.6)
pa = c(0.01, 0.99)

# calculate Shannon entropy for mount washington 
-sum(pm * log(pm))
```

    ## [1] 0.6730117

``` r
# calculate Shannon entropy for Abu dabi 
-sum(pa * log(pa))
```

    ## [1] 0.05600153

Lower entropy means lower uncertainty in abu dabi because it rarely
rains, there is much less uncertainty about whether it will rain on any
given day.

### Why we need information theory

We need a “target” criteria that we want our models to be good at to
select the one that optimizes over / underfitting. Information theory
gives us a measurement scale for the distance between 2 probability
distributions accounting for the fact that some things are harder to
predict from others, some targets are harder to hit than others, some
estimates have higher uncertainty than others.

*Information is the amount of uncertainty that decreases when we learn
an outcome*

We measure this on the scale of
