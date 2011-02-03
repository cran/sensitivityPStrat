library(sensitivityPStrat)

data(vaccine.trial)

vaccine.trial.withNA <- vaccine.trial

set.seed(12345)
for(i in seq_len(20)) 
  vaccine.trial.withNA[sample(nrow(vaccine.trial), size=1, replace=TRUE),
                          sample(ncol(vaccine.trial), size=1)] <- NA

set.seed(12345)
est.bounds<-with(vaccine.trial,
                 sensitivityHHS(z=treatment, s=hiv.outcome, y=logVL,
                     selection="infected", groupings=c("placebo","vaccine"),
                     empty.principal.stratum=c("not infected","infected"),
                     N.boot=1000)
                )
est.bounds


set.seed(12345)
est.bounds<-with(vaccine.trial.withNA,
                 sensitivityHHS(z=treatment, s=hiv.outcome, y=logVL,
                     selection="infected", groupings=c("placebo","vaccine"),
                     empty.principal.stratum=c("not infected","infected"),     
                     na.rm=TRUE, N.boot=100)
                )
est.bounds
