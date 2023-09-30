library('tidyverse')
dat <- read_table(file="CA.fst.weir.fst", col_names=T)
mean(dat$WEIR_AND_COCKERHAM_FST)
# 0.1458919

dat <- read_table(file="WI.fst.weir.fst", col_names=T)
mean(dat$WEIR_AND_COCKERHAM_FST)
# 0.14
