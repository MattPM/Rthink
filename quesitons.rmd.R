


predicting osteoclast frequency as a function of age and conditioning on sex when the actual relationship is that sex causes the age associated changes, e.g. osteolast percent O age and sex as a function of estrogen 

```{r}
dag.ex <- dagitty( "dag {
    S -> E
    A -> E -> O
    S -> E <- A
    A -> O
}")
drawdag(dag.ex)
adjustmentSets( dag.ex , exposure= "A" , outcome="O" )

```