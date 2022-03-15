# ## An example of a fitted brms model
# 
# data(mcompd)
#  
# sbp <- matrix(c(
#   -1, -1, -1,-1, 1,
#   1, -1, -1, -1, 0,
#   0, 1, -1, -1, 0,
#   0, 0, 1, -1, 0),
#   ncol = 5, byrow = TRUE)
# 
# compiltest <- compilr(data = mcompd[, 1:6], sbp = sbp, idvar = "ID")
# 
# bilr <- compiltest[[4]]
# wilr <- compiltest[[5]]
# 
# mcompd$bilr1 <- bilr[,1]
# mcompd$bilr2 <- bilr[,2]
# mcompd$bilr3 <- bilr[,3]
# mcompd$bilr4 <- bilr[,4]
# 
# mcompd$wilr1 <- wilr[,1]
# mcompd$wilr2 <- wilr[,2]
# mcompd$wilr3 <- wilr[,3]
# mcompd$wilr4 <- wilr[,4]
# 
# m <- brm(STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID), 
#          data = mcompd,
#          chain = 4, core = 8)
# 
# summary(m)
