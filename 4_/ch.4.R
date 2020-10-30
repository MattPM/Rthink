#


library(here)
library(rethinking)
library(tidyverse)


pos <- replicate( 1000 , sum( runif(16,-1,1) ) )


## 1 

curve( exp( -x^2 ) , from=-3 , to=3 )




# read in day 1 fold changes of cell frequency 
df = read_delim("../Flu_CITEseq_analysis/data/CHI_H1N1_data/flow/day1-day0.log10.txt", delim = "\t") %>%
  column_to_rownames("ID") %>%
  t %>%
  as.data.frame()

# cool histogram view 
precis(df)
#          mean   sd  5.5% 94.5%       histogram
# ID67     0.15 0.21 -0.14  0.49     ▁▁▅▇▃▂▅▂▁▁▁
# % CD86+ of total monocytes (ID67


# model of  ID 67
#hi ∼ Normal( µ , σ)
#µ ∼ Normal(1, 0.1)
#σ ∼ Uniform(0, 0.05)


# Prior predictive simulation 
# prior for the random variable of normally distributed means 
sample_mu <- rnorm( n = 1e4 ,mean =  0 ,sd =  0.1 ) 
# prior for the random variable of uniformly distributed SD 
sample_sigma <- runif( n = 1e4, min = 0, max = 0.05)
# 
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma ) 
dens( prior_h )

