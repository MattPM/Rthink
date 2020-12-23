#


library(here)
library(rethinking)
library(tidyverse)

## 1 
curve( exp( -x^2 ) , from=-3 , to=3 )

# read in day 1 fold changes of cell frequency 
df = read_delim("../Flu_CITEseq_analysis/data/CHI_H1N1_data/flow/day1-day0.log10.txt", delim = "\t") %>%
  column_to_rownames("ID") %>%
  t %>%
  as.data.frame()

# cool histogram view 
precis(df)
# ID67     0.15 0.21 -0.14  0.49     ▁▁▅▇▃▂▅▂▁▁▁
# % CD86+ of total monocytes (ID67


# model of  ID 67
#hi ∼ Normal( µ , σ)
#µ ∼ Normal(0, 0.15)
#σ ∼ Uniform(0, 0.4)


# Prior predictive simulation 
sample_mu <- rnorm( 1e4 , 178 , 20 ) 
sample_sigma <- runif( 1e4 , 0 , 50 )
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma ) 
dens( prior_h )

