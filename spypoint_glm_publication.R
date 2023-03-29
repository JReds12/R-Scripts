# Spypoint GLM model testing in relation to Burmese python abundance in Key Largo

# 

## load libraries
library(AICcmodavg) # AIC testing for GLM models
library(pastecs)  # calculating mammal detection summary statistics
library(unmarked) # running pcount()

# load spypoint data
spypt.dat = read.csv("spypoint_mammal_occupancy_combined.csv") # for occupancy summary stats# was spypoint.occ.combined.corrected.csv
mam.occ.2016 = read.csv("spypoint_mammal_occupancy_2016.csv") # for 2016 GLMs# mam.occ.2016.corrected.3.csv
mam.occ.2021 = read.csv("spypoint_mammal_occupancy_2021.csv") # for 2021 GLMs # mam.occ.2021.corrected.3.csv

# calculate summaries of occupancy data ----
row.sum.ra = rowSums(spypt.dat[1:15, 3:33], na.rm = T)
stat.desc(row.sum.ra)

row.sum.op = rowSums(spypt.dat[16:30, 34:64], na.rm = T)
stat.desc(row.sum.op)

# create python presence variable for 2016 and 2021 data sets ----
mam.occ.2016$mam.sums = rowSums(mam.occ.2016[6:36])
mam.occ.2021$mam.sums = rowSums(mam.occ.2021[6:36])

#mam.occ.2016$python.presence = mam.occ.2016$python.yearsum 
#mam.occ.2016$python.presence[mam.occ.2016$python.presence > 0] = 1 

#mam.occ.2021$python.presence = mam.occ.2021$python.yearsum
#mam.occ.2021$python.presence[mam.occ.2021$python.presence > 0] = 1

# GLM Models ----

# 2016
glm.1 = glm(mam.sums ~ species + python.presence, data = mam.occ.2016, family = "poisson")
glm.2 = glm(mam.sums ~ species + python.sum, data = mam.occ.2016, family = "poisson")
glm.3 = glm(mam.sums ~ species, data = mam.occ.2016, family = "poisson")
glm.4 = glm(mam.sums ~ python.presence, data = mam.occ.2016, family = "poisson")
glm.5 = glm(mam.sums ~ python.sum, data = mam.occ.2016, family = "poisson")

models = list(glm.1, glm.2, glm.3, glm.4, glm.5)
aictab(models)

summary(glm.1)
summary(glm.2)
summary(glm.3)
summary(glm.4) # used in review 1
summary(glm.5) # used in review 1

boxplot(mam.occ.2016$mam.sums ~ mam.occ.2016$species)


# 2021 
glm.1 = glm(mam.sums ~ species +  python.presence, data = mam.occ.2021, family = "poisson")
glm.2 = glm(mam.sums ~ species +  python.sum, data = mam.occ.2021, family = "poisson")
glm.3 = glm(mam.sums ~ species, data = mam.occ.2021, family = "poisson")
glm.4 = glm(mam.sums ~ python.presence, data = mam.occ.2021, family = "poisson")
glm.5 = glm(mam.sums ~ python.sum, data = mam.occ.2021, family = "poisson")

models = list(glm.1, glm.2, glm.3, glm.4, glm.5)
aictab(models)

summary(glm.1)
summary(glm.2)
summary(glm.3)
summary(glm.4) # used in review 1
summary(glm.5) # used in review 1

boxplot(mam.occ.2021$mam.sums ~ mam.occ.2021$species)

######################################################################################################
y16 = matrix(c("2016"), 15, 1)
y21 = matrix(c("2021"), 15, 1)
years = as.data.frame(rbind(y16,y21))



#######################################################################################################
# pcount models ----
occ.total = as.data.frame(cbind(mam.occ.2016[6:36],mam.occ.2021[6:36]))
y16 = matrix(c("2016"), 30, 31)
y21 = matrix(c("2021"), 30, 31)
years = as.data.frame(cbind(y16,y21))

# mammals combined
spy.umf = unmarkedFramePCount(y = mam.occ.total,
                              siteCovs = mam.occ.2021[1:5],
                              obsCovs = list(year = years))

mod.0 = pcount(~ 1 ~ 1, data = spy.umf, mixture = "P", K = 30)
mod.1 = pcount(~ 1 ~ python.sum, data = spy.umf, mixture = "P", K = 30)
mod.2 = pcount(~ 1 ~ python.presence, data = spy.umf, mixture = "P", K = 30)
mod.3 = pcount(~ 1 ~ species, data = spy.umf, mixture = "P", K = 30)
mod.4 = pcount(~ 1 ~ python.sum + species, data = spy.umf, mixture = "P", K = 30)
mod.5 = pcount(~ 1 ~ python.presence + species, data = spy.umf, mixture = "P", K = 30)
mod.6 = pcount(~ year ~ 1 , data = spy.umf, mixture = "P", K = 30)
mod.7 = pcount(~ year-1 ~ python.sum, data = spy.umf, mixture = "P", K = 30)
mod.8 = pcount(~ year-1 ~ python.presence, data = spy.umf, mixture = "P", K= 30)
mod.9 = pcount(~ year-1 ~ species, data = spy.umf, mixture = "P", K = 30)

fl = fitList(mod.0, mod.1, mod.2, mod.3, mod.4, mod.5, mod.6, mod.7, mod.8, mod.9)
ms_occ = modSel(fl)  
ms_occ 

summary(mod.0)
summary(mod.1)
summary(mod.2)
summary(mod.3)
summary(mod.4)
summary(mod.5)
summary(mod.6)
summary(mod.7)
summary(mod.8)
summary(mod.9)

# raccoon ----
y16 = matrix(c("2016"), 15, 31)
y21 = matrix(c("2021"), 15, 31)
years = as.data.frame(cbind(y16,y21))

ra.spy.umf = unmarkedFramePCount(y = mam.occ.total[1:15,], siteCovs = mam.occ.2021[1:15, 1:5], obsCovs = list(year = years))

mod.0 = pcount(~ 1 ~ 1, data = ra.spy.umf, mixture = "P", K = 30)
mod.1 = pcount(~ 1 ~ python.sum, data = ra.spy.umf, mixture = "P", K = 30)
mod.2 = pcount(~ 1 ~ python.presence, data = ra.spy.umf, mixture = "P", K = 30)
mod.3 = pcount(~ year ~ 1 , data = ra.spy.umf, mixture = "P", K = 30)
mod.4 = pcount(~ year-1 ~ python.sum, data = ra.spy.umf, mixture = "P", K = 30)
mod.5 = pcount(~ year-1 ~ python.presence, data = ra.spy.umf, mixture = "P", K= 30)

fl = fitList(mod.0, mod.1, mod.2, mod.3, mod.4, mod.5)
ms_occ = modSel(fl)  
ms_occ 

summary(mod.0)
summary(mod.1)
summary(mod.2)
summary(mod.3)
summary(mod.4)
summary(mod.5)


# opossums
y16 = matrix(c("2016"), 15, 31)
y21 = matrix(c("2021"), 15, 31)
years = as.data.frame(cbind(y16,y21))

op.spy.umf = unmarkedFramePCount(y = mam.occ.total[16:30,], siteCovs = mam.occ.2021[16:30, 1:5], obsCovs = list(year = years))

mod.0 = pcount(~ 1 ~ 1, data = op.spy.umf, mixture = "P", K = 30)
mod.1 = pcount(~ 1 ~ python.sum, data = op.spy.umf, mixture = "P", K = 30)
mod.2 = pcount(~ 1 ~ python.presence, data = op.spy.umf, mixture = "P", K = 30)
mod.3 = pcount(~ year-1 ~ 1 , data = op.spy.umf, mixture = "P", K = 30)
mod.4 = pcount(~ year-1 ~ python.sum, data = op.spy.umf, mixture = "P", K = 30)
mod.5 = pcount(~ year-1 ~ python.presence, data = op.spy.umf, mixture = "P", K= 30)


fl = fitList(mod.0, mod.1, mod.2, mod.3, mod.4, mod.5)
ms_occ = modSel(fl)  
ms_occ 

summary(mod.0)
summary(mod.1)
summary(mod.2)
summary(mod.3)
summary(mod.4)
summary(mod.5)





####
spy.umf.2 = unmarkedFramePCount(y = mam.occ.total, siteCovs = mam.occ.2021[1:5])

mod.0 = pcount(~1 ~ 1, data = spy.umf.2, mixture = "P", K = 30)
mod.1 = pcount(~1 ~ python.sum, data = spy.umf.2, mixture = "P", K = 30)
mod.2 = pcount(~1 ~ species, data = spy.umf.2, mixture = "P", K = 30)
mod.3 = pcount(~1 ~ python.sum + species, data = spy.umf.2, mixture = "P", K = 30)

fl = fitList(mod.0, mod.1, mod.2, mod.3)
ms_occ = modSel(fl)  
ms_occ 


