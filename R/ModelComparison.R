
#######################################################################
######################### Model Comparison ############################
#######################################################################


set.seed(123)
imps_SMData <- TwoFCS(df, smtype="lm", smformula="AKT1~IL16+APC+ID3+SLA+CDK4+EGR1+TCF12+MCL1",
                      method = c("norm","norm","norm","norm","norm","norm","norm","norm",
                                 "","","","","","","","","","",
                                 "","","","","","","","","","",
                                 "","","","","","","","","","",
                                 "","","","","","","","","","",
                                 "","","","","","","","","","")
                      ,m=5)
summary(imps_SMData)

MIce_RF <- mice(df,m=5,maxit=50,meth='rf',seed=500)
summary(MIce_RF)
#get complete data (3rd out of 5)
MICE_RF_completeData <- complete(MIce_RF,3)
View(MICE_RF_completeData)

lm.mice.out <- with(MIce_RF, lm(AKT1~IL16+APC+ID3+SLA+CDK4+EGR1+TCF12+MCL1))
summary(lm.mice.out)

pool.mice <- pool(lm.mice.out)

summary(pool.mice)  # multiply imputed results


lm <- lm(AKT1~IL16+APC+ID3+SLA+CDK4+EGR1+TCF12+MCL1,data = dataLong)
summary(lm)
coef(summary(lm))[, 1:2]  # original results

s.real <- summary(lm(AKT1~IL16+APC+ID3+SLA+CDK4+EGR1+TCF12+MCL1,data = dataLong))$coef[, 1:2]
s.mice <- summary(pool.mice)[, 1:2]  # multiply imputed results

dt <- imps_SMData[["impDatasets"]][[3]]


SM_lm.mice.out <- with(dt, lm(AKT1~IL16+APC+ID3+SLA+CDK4+EGR1+TCF12+MCL1))
summary(SM_lm.mice.out)

SM_pool.mice <- pool(SM_lm.mice.out)

summary(SM_lm.mice.out)  # multiply imputed results

sum_miceOut[["coefficients"]]
sum_miceOut[["coefficients"]][,1]


allout <- cbind(s.real[, 1], s.mice[, 2], sum_miceOut[["coefficients"]][,1])
colnames(allout) <- c("Real Relationship", "mice", "FCS-Twofold")
allout
#######################################################################
######################### Model Comparison ############################
#######################################################################
