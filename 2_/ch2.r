#


library(here)
library(rethinking)
library(tidyverse)

## 1 
ways <- c( 0 , 3 , 8 , 9 , 0 )
ways /  sum(ways)

# globe toss 
dbinom( 6 , size=9 , prob=0.5 )
dbinom(x = c(),  )