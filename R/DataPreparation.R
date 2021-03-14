#######################################################################
######################### Data Preparation ############################
#######################################################################

#install.packages("longitudinal")
library("longitudinal")

# load data sets
data(tcell)

tcell.44 <- combine.longitudinal(tcell.34, tcell.10)

dim(tcell.44)
typeof(tcell.44)

dataLong <- data.frame(rbind(tcell.44, tcell.44))
dim(dataLong)
typeof(dataLong)

dataLong <- dataLong[1:440,]
dim(dataLong)
colnames(dataLong)

df1 = subset(dataLong, select = c(IL16,APC,ID3,SLA,CDK4,EGR1,TCF12,MCL1) )
dim(df1)
df2 = subset(dataLong, select = -c(IL16,APC,ID3,SLA,CDK4,EGR1,TCF12,MCL1) )
dim(df2)

df1 <- prodNA(df1, noNA = 0.3)
colSums(is.na(df1))
df <- data.frame(cbind(df1, df2))
dim(df)
View(df)
#######################################################################
######################### Data Preparation ############################
#######################################################################