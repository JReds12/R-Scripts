setwd("C:/Users/jredi/Documents/Mammal Python Occupancy/final data")

# Mammal Occupancy models in relation to Burmese python abundance in Key Largo
# 

# load libraries
library(pastecs)  # calculating mammal detection summary statistics
library(unmarked)  # creating occupancy models
#library(tidyverse) 

# import data  ----
env.covs = read.csv("CM_occ_dat_covs.csv", header = T)# load environmental variables 
  names(env.covs)  # only dist.to.development used in models


python.covs = read.csv("CLNWR_eMammal_python_variables.csv", header = T) # load python variables## originally CLNWR_eMammal_python_variables_3.csv
  names(python.covs)
  cor(python.covs[c('python.totalsum','python.dif','python.year')])
  
raccoon = read.csv("Northern_Raccoon_detection_history.csv", header = T) # load mammal detection histories
opossum = read.csv("Virginia_Opossum_detection_history.csv", header = T)  
squirrel = read.csv("Eastern_Gray_Squirrel_detection_history.csv", header = T)
cotton.mouse = read.csv("Cotton_Mouse_detection_history.csv", header = T)
woodrat = read.csv("Eastern_Woodrat_detection_history.csv", header = T)
black.rat = read.csv("House_rat_detection_history.csv", header = T)


# create function that formats data  ----
format.occ.dat = function(x) {
  mam.2016 = x[which(x$year == 2016), 3:10]
  mam.2018 = x[which(x$year == 2018), 3:10]
  mam.2020 = x[which(x$year == 2020), 3:10]
  mam.dat = merge(env.covs, python.covs, by = "Camera.ID")
  mam.dat = merge(mam.dat, mam.2016, by = "Camera.ID")
  mam.dat = merge(mam.dat, mam.2018, by = "Camera.ID")
  mam.dat = merge(mam.dat, mam.2020, by = "Camera.ID")
  return(mam.dat)
}  # combines variables and detection histories for models

# calculate summary statistics for mammal detection histories ----
stat.desc(rowSums(raccoon[,4:10], na.rm = T))
stat.desc(rowSums(opossum[,4:10], na.rm = T))
stat.desc(rowSums(squirrel[,4:10], na.rm = T))
stat.desc(rowSums(cotton.mouse[,4:10], na.rm = T))
stat.desc(rowSums(woodrat[,4:10], na.rm = T))
stat.desc(rowSums(black.rat[,4:10], na.rm = T))




# Raccoon Model Analysis ----
raccoon.dat = format.occ.dat(raccoon)
raccoon.umf = unmarkedMultFrame(y = raccoon.dat[28:48], siteCovs = raccoon.dat[1:27], numPrimary = 3)
occ.1 <- colext(~1, ~1, ~1 , ~1, raccoon.umf)
occ.2 <- colext(~1, ~1, ~1 , ~ scale(dist.to.development), raccoon.umf)
occ.3 <- colext(~1, ~1, ~1 , ~ scale(python.dif),  raccoon.umf)
occ.4 <- colext(~1, ~1, ~1 , ~ scale(python.dif) + scale(dist.to.development), raccoon.umf)
occ.5 <- colext(~1, ~1, ~ scale(python.totalsum), ~ 1, raccoon.umf)
occ.6 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(dist.to.development), raccoon.umf)
occ.7 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(python.dif), raccoon.umf)
occ.8 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(python.dif) + scale(dist.to.development), raccoon.umf)
occ.9 <- colext(~1, ~1, ~ 1, ~ scale(python.year), raccoon.umf)
occ.10 <- colext(~1, ~1,~ scale(python.year), ~ 1, raccoon.umf)

fl <- fitList(occ.1, occ.2, occ.3, occ.4, occ.5, occ.6, occ.7, occ.8, occ.9, occ.10)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.1)
summary(occ.2)
summary(occ.3)
summary(occ.4)
summary(occ.5)
summary(occ.6)
summary(occ.7)
summary(occ.8)
summary(occ.9)
summary(occ.10)

## Full time-dependent analysis ----
survey.years = matrix(c('2016','2018','2020'), 428, 3, byrow = T)

raccoon.time.umf = unmarkedMultFrame(y = raccoon.dat[28:48], siteCovs = raccoon.dat[1:27], yearlySiteCovs = list(year = survey.years), numPrimary = 3)

occ.0 = colext(~1, ~1, ~1, ~1, raccoon.time.umf)
occ.1 = colext(~1, ~1, ~1, ~year-1 + python.totalsum, raccoon.time.umf)
occ.2 = colext(~1, ~1, ~year-1, ~1, raccoon.time.umf)
occ.3 = colext(~1, ~1, ~year-1, ~year-1, raccoon.time.umf)
occ.4 = colext(~1, ~1, ~1, ~year, raccoon.time.umf)
occ.5 = colext(~1, ~1, ~year, ~1, raccoon.time.umf)

fl <- fitList(occ.0, occ.1, occ.2, occ.3)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.0)
summary(occ.1)
summary(occ.2)
  plogis(0.111)
  plogis(-0.89)

summary(occ.3)

## plotting
nd <- data.frame(year = c('2016','2018'))
E.ext <- predict(occ.2, type = 'ext', newdata = nd)

nd <- data.frame(year = c('2016','2018','2020'))
E.det <- predict(occ.1, type = 'det', newdata = nd)

op <- par(mfrow=c(1,2), mai=c(0.6, 0.6, 0.1, 0.1))
with(E.ext, { # Plot for extinction probability
  plot(1:2, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Extinction probability ( ', epsilon, ' )')),
       ylim=c(0,1))
  axis(1, at=1:2, labels=nd$year[1:2])
  arrows(1:2, lower, 1:2, upper, code=3, angle=90, length=0.03)
  #points((1:2)-0.1, 1-phi, col=1, lwd = 1, pch=16)
  legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

with(E.det, { # Plot for detection probability: note 10 years
  plot(1:3, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Detection probability ( ', p, ' )')),
       ylim=c(0,1))
  axis(1, at=1:3, labels=nd$year)
  arrows(1:3, lower, 1:3, upper, code=3, angle=90, length=0.03)
  #points((1:3)-0.1, p, col=1, lwd = 1, pch=16)
  legend(7.5, 1, c('Parameter','Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

par(op)

# Opossum Model Selection ----
#mam.dat = mam.dat.prep(opossum)
opossum.dat = format.occ.dat(opossum)
opossum.umf = unmarkedMultFrame(y = opossum.dat[28:48], siteCovs = opossum.dat[1:27], numPrimary = 3)
occ.1 <- colext(~1, ~1, ~1 , ~1, opossum.umf)
occ.2 <- colext(~1, ~1, ~1 , ~ scale(dist.to.development), opossum.umf)
occ.3 <- colext(~1, ~1, ~1 , ~ scale(python.dif),  opossum.umf)
occ.4 <- colext(~1, ~1, ~1 , ~ scale(python.dif) + scale(dist.to.development), opossum.umf)
occ.5 <- colext(~1, ~1, ~ scale(python.totalsum) , ~ 1, opossum.umf)
occ.6 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(dist.to.development), opossum.umf)
occ.7 <- colext(~1, ~1, ~ scale(python.totalsum) , ~ scale(python.dif), opossum.umf)
occ.8 <- colext(~1, ~1, ~ scale(python.totalsum) , ~ scale(python.dif) + scale(dist.to.development), opossum.umf)
occ.9 <- colext(~1, ~1, ~ 1, ~ scale(python.year), opossum.umf)
occ.10 <- colext(~1, ~1,~ scale(python.year), ~ 1, opossum.umf)

fl <- fitList(occ.1, occ.2, occ.3, occ.4, occ.5, occ.6, occ.7, occ.8, occ.9, occ.10)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.1)
summary(occ.2)
summary(occ.3)
summary(occ.4)
summary(occ.5)
summary(occ.6)
summary(occ.7)
summary(occ.8)
summary(occ.9)
summary(occ.10)

## Full time-dependent analysis
survey.years = matrix(c('2016','2018','2020'), 428, 3, byrow = T)

opossum.time.umf = unmarkedMultFrame(y = opossum.dat[28:48], yearlySiteCovs = list(year = survey.years), numPrimary = 3)

occ.0 = colext(~1, ~1, ~1, ~1, opossum.time.umf)
occ.1 = colext(~1, ~1, ~1, ~year-1, opossum.time.umf)
occ.2 = colext(~1, ~1, ~year-1, ~1, opossum.time.umf)
occ.3 = colext(~1, ~1, ~year-1, ~year-1, opossum.time.umf)
#occ.4 = colext(~year-1, ~1, ~1, ~1, raccoon.time.umf)

fl <- fitList(occ.0, occ.1, occ.2, occ.3)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.0)
summary(occ.1)
summary(occ.2)
summary(occ.3)

## plotting
nd <- data.frame(year = c('2016','2018'))
E.ext <- predict(occ.2, type = 'ext', newdata = nd)

nd <- data.frame(year = c('2016','2018','2020'))
E.det <- predict(occ.1, type = 'det', newdata = nd)

op <- par(mfrow=c(1,2), mai=c(0.6, 0.6, 0.1, 0.1))
with(E.ext, { # Plot for extinction probability
  plot(1:2, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Extinction probability ( ', epsilon, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:2, labels=nd$year[1:2])
  arrows(1:2, lower, 1:2, upper, code=3, angle=90, length=0.03, col=4)
  legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

with(E.det, { # Plot for detection probability: note 10 years
  plot(1:3, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Detection probability ( ', p, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:3, labels=nd$year)
  arrows(1:3, lower, 1:3, upper, code=3, angle=90, length=0.03, col=4)
  points((1:3)-0.1, p, col=1, lwd = 1, pch=16)
  legend(7.5, 1, c('Parameter','Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

par(op)

# Squirrel Model Selection ----
#mam.dat = mam.dat.prep(squirrel)
squirrel.dat = format.occ.dat(squirrel)
squirrel.umf = unmarkedMultFrame(y = squirrel.dat[28:48], siteCovs = squirrel.dat[1:27], numPrimary = 3)
occ.1 <- colext(~1, ~1, ~1 , ~1, squirrel.umf)
occ.2 <- colext(~1, ~1, ~1 , ~ scale(dist.to.development), squirrel.umf)
occ.3 <- colext(~1, ~1, ~1 , ~ scale(python.dif),  squirrel.umf)
occ.4 <- colext(~1, ~1, ~1 , ~ scale(python.dif) + scale(dist.to.development), squirrel.umf)
occ.5 <- colext(~1, ~1, ~ scale(python.totalsum), ~ 1, squirrel.umf)
occ.6 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(dist.to.development), squirrel.umf)
occ.7 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(python.dif), squirrel.umf)
occ.8 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(python.dif) + scale(dist.to.development), squirrel.umf)
occ.9 <- colext(~1, ~1, ~ 1, ~ scale(python.year), squirrel.umf)
occ.10 <- colext(~1, ~1,~ scale(python.year), ~ 1, squirrel.umf)

fl <- fitList(occ.1, occ.2, occ.3, occ.4, occ.5, occ.6, occ.7, occ.8, occ.9, occ.10)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.1)
summary(occ.2)
summary(occ.3)
summary(occ.4)
summary(occ.5)
summary(occ.6)
summary(occ.7)
summary(occ.8)
summary(occ.9)
summary(occ.10)

## Full time-dependent analysis
survey.years = matrix(c('2016','2018','2020'), 428, 3, byrow = T)

squirrel.time.umf = unmarkedMultFrame(y = squirrel.dat[28:48], yearlySiteCovs = list(year = survey.years), numPrimary = 3)

occ.0 = colext(~1, ~1, ~1, ~1, squirrel.time.umf)
occ.1 = colext(~1, ~1, ~1, ~year-1, squirrel.time.umf)
occ.2 = colext(~1, ~1, ~year-1, ~1, squirrel.time.umf)
occ.3 = colext(~1, ~1, ~year-1, ~year-1, squirrel.time.umf)
#occ.4 = colext(~year-1, ~1, ~1, ~1, raccoon.time.umf)

fl <- fitList(occ.0, occ.1, occ.2, occ.3)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.0)
summary(occ.1)
summary(occ.2)
summary(occ.3)

## plotting
nd <- data.frame(year = c('2016','2018'))
E.ext <- predict(occ.2, type = 'ext', newdata = nd)

nd <- data.frame(year = c('2016','2018','2020'))
E.det <- predict(occ.1, type = 'det', newdata = nd)

op <- par(mfrow=c(1,2), mai=c(0.6, 0.6, 0.1, 0.1))
with(E.ext, { # Plot for extinction probability
  plot(1:2, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Extinction probability ( ', epsilon, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:9, labels=nd$year[1:9])
  arrows(1:9, lower, 1:9, upper, code=3, angle=90, length=0.03, col=4)
  points((1:9)-0.1, 1-phi, col=1, lwd = 1, pch=16)
  legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

with(E.det, { # Plot for detection probability: note 10 years
  plot(1:3, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Detection probability ( ', p, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:10, labels=nd$year)
  arrows(1:10, lower, 1:10, upper, code=3, angle=90, length=0.03, col=4)
  points((1:10)-0.1, p, col=1, lwd = 1, pch=16)
  legend(7.5, 1, c('Parameter','Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

par(op)


# Cotton Mouse Model Selection ----
#mam.dat = mam.dat.prep(cotton.mouse)
cotton.mouse.dat = format.occ.dat(cotton.mouse)
cotton.mouse.umf = unmarkedMultFrame(y = cotton.mouse.dat[28:48], siteCovs = cotton.mouse.dat[1:27], numPrimary = 3)
occ.1 <- colext(~1, ~1, ~1 , ~1, cotton.mouse.umf)
occ.2 <- colext(~1, ~1, ~1 , ~ scale(dist.to.development), cotton.mouse.umf)
occ.3 <- colext(~1, ~1, ~1 , ~ scale(python.dif),  cotton.mouse.umf)
occ.4 <- colext(~1, ~1, ~1 , ~ scale(python.dif) + scale(dist.to.development), cotton.mouse.umf)
occ.5 <- colext(~1, ~1, ~ scale(python.totalsum) , ~ 1, cotton.mouse.umf)
occ.6 <- colext(~1, ~1, ~ scale(python.totalsum) , ~ scale(dist.to.development), cotton.mouse.umf)
occ.7 <- colext(~1, ~1, ~ scale(python.totalsum) , ~ scale(python.dif), cotton.mouse.umf)
occ.8 <- colext(~1, ~1, ~ scale(python.totalsum) , ~ scale(python.dif) + scale(dist.to.development), cotton.mouse.umf)
occ.9 <- colext(~1, ~1, ~ 1, ~ scale(python.year), cotton.mouse.umf)
occ.10 <- colext(~1, ~1,~ scale(python.year), ~ 1, cotton.mouse.umf)

fl <- fitList(occ.1, occ.2, occ.3, occ.4, occ.5, occ.6, occ.7, occ.8, occ.9, occ.10)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.1)
summary(occ.2)
summary(occ.3)
summary(occ.4)
summary(occ.5)
summary(occ.6)
summary(occ.7)
summary(occ.8)
summary(occ.9)
summary(occ.10)


## Full time-dependent analysis
survey.years = matrix(c('2016','2018','2020'), 428, 3, byrow = T)

cotton.mouse.time.umf = unmarkedMultFrame(y = cotton.mouse.dat[28:48], yearlySiteCovs = list(year = survey.years), numPrimary = 3)

occ.0 = colext(~1, ~1, ~1, ~1, cotton.mouse.time.umf)
occ.1 = colext(~1, ~1, ~1, ~year-1, cotton.mouse.time.umf)
occ.2 = colext(~1, ~1, ~year-1, ~1, cotton.mouse.time.umf)
occ.3 = colext(~1, ~1, ~year-1, ~year-1, cotton.mouse.time.umf)
#occ.4 = colext(~year-1, ~1, ~1, ~1, raccoon.time.umf)

fl <- fitList(occ.0, occ.1, occ.2, occ.3)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.0)
summary(occ.1)
summary(occ.2)
summary(occ.3)

## plotting
nd <- data.frame(year = c('2016','2018'))
E.ext <- predict(occ.2, type = 'ext', newdata = nd)

nd <- data.frame(year = c('2016','2018','2020'))
E.det <- predict(occ.1, type = 'det', newdata = nd)

op <- par(mfrow=c(1,2), mai=c(0.6, 0.6, 0.1, 0.1))
with(E.ext, { # Plot for extinction probability
  plot(1:2, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Extinction probability ( ', epsilon, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:9, labels=nd$year[1:9])
  arrows(1:9, lower, 1:9, upper, code=3, angle=90, length=0.03, col=4)
  points((1:9)-0.1, 1-phi, col=1, lwd = 1, pch=16)
  legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

with(E.det, { # Plot for detection probability: note 10 years
  plot(1:3, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Detection probability ( ', p, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:10, labels=nd$year)
  arrows(1:10, lower, 1:10, upper, code=3, angle=90, length=0.03, col=4)
  points((1:10)-0.1, p, col=1, lwd = 1, pch=16)
  legend(7.5, 1, c('Parameter','Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

par(op)

# Woodrat Model Selection ----
#mam.dat = mam.dat.prep(woodrat)
woodrat.dat = format.occ.dat(woodrat)
woodrat.umf = unmarkedMultFrame(y = woodrat.dat[28:48], siteCovs = woodrat.dat[1:27], numPrimary = 3)
occ.1 <- colext(~1, ~1, ~1 , ~1, woodrat.umf)
occ.2 <- colext(~1, ~1, ~1 , ~ scale(dist.to.development), woodrat.umf)
occ.3 <- colext(~1, ~1, ~1 , ~ scale(python.dif),  woodrat.umf)
occ.4 <- colext(~1, ~1, ~1 , ~ scale(python.dif) + scale(dist.to.development), woodrat.umf)
occ.5 <- colext(~1, ~1, ~ scale(python.totalsum), ~ 1, woodrat.umf)
occ.6 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(dist.to.development), woodrat.umf)
occ.7 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(python.dif), woodrat.umf)
occ.8 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(python.dif) + scale(dist.to.development), woodrat.umf)
occ.9 <- colext(~1, ~1, ~ 1, ~ scale(python.year), woodrat.umf)
occ.10 <- colext(~1, ~1,~ scale(python.year), ~ 1, woodrat.umf)

fl <- fitList(occ.1, occ.2, occ.3, occ.4, occ.5, occ.6, occ.7, occ.8, occ.9, occ.10)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.1)
summary(occ.2)
summary(occ.3)
summary(occ.4)
summary(occ.5)
summary(occ.6)
summary(occ.7)
summary(occ.8)
summary(occ.9)
summary(occ.10)


## Full time-dependent analysis
survey.years = matrix(c('2016','2018','2020'), 428, 3, byrow = T)

woodrat.time.umf = unmarkedMultFrame(y = woodrat.dat[28:48], yearlySiteCovs = list(year = survey.years), numPrimary = 3)

occ.0 = colext(~1, ~1, ~1, ~1, woodrat.time.umf)
occ.1 = colext(~1, ~1, ~1, ~year-1, woodrat.time.umf)
occ.2 = colext(~1, ~1, ~year-1, ~1, woodrat.time.umf)
occ.3 = colext(~1, ~1, ~year-1, ~year-1, woodrat.time.umf)

fl <- fitList(occ.0, occ.1, occ.2, occ.3)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.0)
summary(occ.1)
summary(occ.2)
summary(occ.3)

## plotting
nd <- data.frame(year = c('2016','2018'))
E.ext <- predict(occ.2, type = 'ext', newdata = nd)

nd <- data.frame(year = c('2016','2018','2020'))
E.det <- predict(occ.1, type = 'det', newdata = nd)

op <- par(mfrow=c(1,2), mai=c(0.6, 0.6, 0.1, 0.1))
with(E.ext, { # Plot for extinction probability
  plot(1:2, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Extinction probability ( ', epsilon, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:2, labels=nd$year[1:2])
  arrows(1:2, lower, 1:2, upper, code=3, angle=90, length=0.03, col=4)
  legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

with(E.det, { # Plot for detection probability: note 10 years
  plot(1:3, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Detection probability ( ', p, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:3, labels=nd$year)
  arrows(1:3, lower, 1:3, upper, code=3, angle=90, length=0.03, col=4)
  legend(7.5, 1, c('Parameter','Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

par(op)


# Black Rat Model Selection ----
#mam.dat = mam.dat.prep(black.rat)
black.rat.dat = format.occ.dat(black.rat)
black.rat.umf = unmarkedMultFrame(y = black.rat.dat[28:48], siteCovs = black.rat.dat[1:27], numPrimary = 3)
occ.1 <- colext(~1, ~1, ~1 , ~1, black.rat.umf)
occ.2 <- colext(~1, ~1, ~1 , ~ scale(dist.to.development), black.rat.umf)
occ.3 <- colext(~1, ~1, ~1 , ~ scale(python.dif),  black.rat.umf)
occ.4 <- colext(~1, ~1, ~1 , ~ scale(python.dif) + scale(dist.to.development), black.rat.umf)
occ.5 <- colext(~1, ~1, ~ scale(python.totalsum), ~ 1, black.rat.umf)
occ.6 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(dist.to.development), black.rat.umf)
occ.7 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(python.dif), black.rat.umf)
occ.8 <- colext(~1, ~1, ~ scale(python.totalsum), ~ scale(python.dif) + scale(dist.to.development), black.rat.umf)
occ.9 <- colext(~1, ~1, ~ 1, ~ scale(python.year), black.rat.umf)
occ.10 <- colext(~1, ~1,~ scale(python.year), ~ 1, black.rat.umf)

fl <- fitList(occ.1, occ.2, occ.3, occ.4, occ.5, occ.6, occ.7, occ.8, occ.9, occ.10)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.1)
summary(occ.2)
summary(occ.3)
summary(occ.4)
summary(occ.5)
summary(occ.6)
summary(occ.7)
summary(occ.8)
summary(occ.9)
summary(occ.10)

## Full time-dependent analysis
survey.years = matrix(c('2016','2018','2020'), 428, 3, byrow = T)

black.rat.time.umf = unmarkedMultFrame(y = black.rat.dat[28:48], yearlySiteCovs = list(year = survey.years), numPrimary = 3)

occ.0 = colext(~1, ~1, ~1, ~1, black.rat.time.umf)
occ.1 = colext(~1, ~1, ~1, ~year-1, black.rat.time.umf)
occ.2 = colext(~1, ~1, ~year-1, ~1, black.rat.time.umf)
occ.3 = colext(~1, ~1, ~year-1, ~year-1, black.rat.time.umf)
#occ.4 = colext(~year-1, ~1, ~1, ~1, raccoon.time.umf)

fl <- fitList(occ.0, occ.1, occ.2, occ.3)
ms_occ <- modSel(fl)  
ms_occ 

summary(occ.0)
summary(occ.1)
summary(occ.2)
summary(occ.3)

## plotting
nd <- data.frame(year = c('2016','2018'))
E.ext <- predict(occ.2, type = 'ext', newdata = nd)

nd <- data.frame(year = c('2016','2018','2020'))
E.det <- predict(occ.1, type = 'det', newdata = nd)

op <- par(mfrow=c(1,2), mai=c(0.6, 0.6, 0.1, 0.1))
with(E.ext, { # Plot for extinction probability
  plot(1:2, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Extinction probability ( ', epsilon, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:9, labels=nd$year[1:9])
  arrows(1:9, lower, 1:9, upper, code=3, angle=90, length=0.03, col=4)
  points((1:9)-0.1, 1-phi, col=1, lwd = 1, pch=16)
  legend(7, 1, c('Parameter', 'Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

with(E.det, { # Plot for detection probability: note 10 years
  plot(1:3, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Detection probability ( ', p, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:10, labels=nd$year)
  arrows(1:10, lower, 1:10, upper, code=3, angle=90, length=0.03, col=4)
  points((1:10)-0.1, p, col=1, lwd = 1, pch=16)
  legend(7.5, 1, c('Parameter','Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})

par(op)

