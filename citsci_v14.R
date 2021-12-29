

library(dplyr); library(car); library(glmmTMB); library(DHARMa); library(MuMIn); library(bbmle)

# All subsets model selection ::
options(na.action = 'na.fail')

#####################################################################
## ANALYSIS :: STATES
#####################################################################

#####################################################################
## MODEL PREP :: CORRELATIONS, TRANSFORMATIONS

state.dat = read.csv('Input/analysisDFs/stateAnalysis.csv',
                     header = T, stringsAsFactors = F)

# Format GEOID with leading zeros
state.dat$state_code = stringr::str_pad(state.dat$state_code, 2, side = 'left', pad = '0')

# Transformations check
hist(state.dat$popDens) ## log10
hist(state.dat$eduDens) ## log1p
hist(state.dat$med_income)
hist(state.dat$bachDeg25p)
hist(state.dat$vacRecHs)
hist(state.dat$openWater) ## log1p
hist(state.dat$forest) ## log1p
hist(state.dat$plantCult) ## log1p
hist(state.dat$wetland) ## log1p
hist(state.dat$developed) ## log1p

# Center, scale, transform predictors
state.mod.df = state.dat

state.mod.df$area_km2 = log(state.mod.df$area_km2)  ## natural log transform area
state.mod.df$popDens = scale(log10(state.mod.df$popDens))
state.mod.df$eduDens = scale(log1p(state.mod.df$eduDens))
state.mod.df$med_income = scale(state.mod.df$med_income)
state.mod.df$bachDeg25p = scale(state.mod.df$bachDeg25p)
state.mod.df$vacRecHs = scale(state.mod.df$vacRecHs)
state.mod.df$openWater = scale(log1p(state.mod.df$openWater))
state.mod.df$developed = scale(log1p(state.mod.df$developed))
state.mod.df$forest = scale(log1p(state.mod.df$forest))
state.mod.df$plantCult = scale(log1p(state.mod.df$plantCult))
state.mod.df$wetland = scale(log1p(state.mod.df$wetland))

# CORRELATIONS
corMat = cor(state.mod.df[ , c('popDens', 'eduDens', 'med_income', 
                            'bachDeg25p', 'vacRecHs', 
                            'openWater', 'developed', 'forest',
                            'plantCult', 'wetland')])
corMat[upper.tri(corMat, diag = T)] = NA
corMat = as.data.frame(corMat)

# write.csv(corMat, 'Output/corMat/stateCorr.csv',
#           row.names = T)

rm(corMat)

#####################################################################
## MODELLING SITE COUNT :: SELECTING BEST ERROR STRUCTURE

# 1. First, determine what is the best error structure (from poisson, nb1, nb2, and the zero-inflated versions of these 3)

# poisson model	
state.fullmod.pois <- glmmTMB(nSites ~
                                      plantCult +
                                      bachDeg25p +
                                      forest +
                                      openWater +
                                      eduDens +
                                      popDens +
                                      vacRecHs +
                                      wetland +
                                offset(area_km2),
                              data = state.mod.df,
                              family = poisson(link = 'log'))

# neg bin 1 model									
state.fullmod.nb1 <- glmmTMB(nSites ~
                                     plantCult +
                                     bachDeg25p +
                                     forest +
                                     openWater +
                                     eduDens +
                                     popDens +
                                     vacRecHs +
                                     wetland +
                               offset(area_km2),
                             data = state.mod.df, 
                             family = nbinom1(link = 'log'))

# neg bin 2 model									
state.fullmod.nb2 <- glmmTMB(nSites ~
                                     plantCult +
                                     bachDeg25p +
                                     forest +
                                     openWater +
                                     eduDens +
                                     popDens +
                                     vacRecHs +
                                     wetland +
                               offset(area_km2),
                             data = state.mod.df,
                             family = nbinom2(link = 'log'))

# poisson model	with zeroinfl
state.fullmod.zipois <- glmmTMB(nSites ~
                                        plantCult +
                                        bachDeg25p +
                                        forest +
                                        openWater +
                                        eduDens +
                                        popDens +
                                        vacRecHs +
                                        wetland +
                                  offset(area_km2),
                                ziformula = ~1,
                                data = state.mod.df,
                                family = poisson(link = 'log'))

# neg bin 1 model with zeroinfl									
state.fullmod.zinb1 <- glmmTMB(nSites ~
                                       plantCult +
                                       bachDeg25p +
                                       forest +
                                       openWater +
                                       eduDens +
                                       popDens +
                                       vacRecHs +
                                       wetland +
                                 offset(area_km2),
                               ziformula = ~1,
                               data = state.mod.df, 
                               family = nbinom1(link = 'log'))

# neg bin 2 model with zeroinfl									
state.fullmod.zinb2 <- glmmTMB(nSites ~
                                       plantCult +
                                       bachDeg25p +
                                       forest +
                                       openWater +
                                       eduDens +
                                       popDens +
                                       vacRecHs +
                                       wetland +
                                 offset(area_km2),
                               ziformula = ~1,
                               data = state.mod.df, 
                               family = nbinom2(link = 'log'))


## use AICc to check with error structure is the best
AICctab(state.fullmod.pois, 
        state.fullmod.nb1, 
        state.fullmod.nb2, 
        state.fullmod.zipois, 
        state.fullmod.zinb1, 
        state.fullmod.zinb2, 
        base = T, 
        weights = T)

#                        AICc     dAICc    df weight
# state.fullmod.nb2       853.8      0.0 10 0.84  
# state.fullmod.zinb2     857.1      3.3 11 0.16  
# state.fullmod.nb1       873.6     19.8 10 <0.001
# state.fullmod.zinb1     876.9     23.1 11 <0.001
# state.fullmod.pois   182438.8 181585.0 9  <0.001
# state.fullmod.zipois 182442.0 181588.2 10 <0.001 

rm(list = setdiff(ls(), c('state.mod.df', 'state.dat', 'state.fullmod.nb2')))

#####################################################################
## MODELLING SITE COUNT :: MODEL SELECTION ON BEST-FIT MODEL

# All-subsets model selection - keep offset term in all models
mod.sel = dredge(state.fullmod.nb2,
				rank = 'AICc', fixed='cond(offset(area_km2))') # offset term as always included in every model				 

# Rename model selection object (ws=watershed level)	
state.mod.sel <- mod.sel

# Top 10 models
state.mod.sel[1:10]
				
# Top model
state.mod = get.models(state.mod.sel, subset = 1)[[1]]
summary(state.mod)
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.9667     0.1644 -24.125  < 2e-16 ***
#  forest        0.4724     0.1704   2.773  0.00555 ** 
#  popDens       1.2508     0.2024   6.179 6.46e-10 ***

# Top-ranked AICc model
# Null model to get R2
state.null = glmmTMB(nSites ~ 1,
                    data = state.mod.df,
                    family = nbinom2(link = 'log'))

MuMIn::r.squaredLR(state.mod, state.null)
# [1] 0.5437035

# Residual checks :: https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#installing-loading-and-citing-the-package
simulationOutput.st = simulateResiduals(fittedModel = state.mod)
windows(9,5)
plot(simulationOutput.st, quantreg=T) # for WebFigure 4

# Test for dispersion
testDispersion(simulationOutput.st)
# data:  simulationOutput
# dispersion = 2.4093, p-value = 0.136
# alternative hypothesis: two.sided

# Test zero inflation - probably not needed here (no zeros)
testZeroInflation(simulationOutput.st)
# data:  simulationOutput
# ratioObsSim = 0, p-value = 1
# alternative hypothesis: two.sided

# Get coefficients, se, and R2marginal for top 10 models (WebPanel 3)
# e.g., for 2nd ranked model, use subset=2; for 3rd ranked model, use subset=3, etc.
state.mod.rank2 = get.models(state.mod.sel, subset = 2)[[1]]
summary(state.mod.rank2)

MuMIn::r.squaredLR(state.mod.rank2, state.null)

rm(list = ls())

#####################################################################
## ANALYSIS :: COUNTY
#####################################################################

#####################################################################
## MODEL PREP :: CORRELATIONS, TRANSFORMATIONS

count.dat = read.csv('Input/analysisDFs/countyAnalysis.csv',
                     header = T, stringsAsFactors = F)

# Format GEOID with leading zeros
count.dat$GEOID = stringr::str_pad(count.dat$GEOID, 5, side = 'left', pad = '0')

# Transformations check
hist(count.dat$siteDens)
hist(count.dat$popDens) ## log10
hist(count.dat$eduDens) ## log1p
hist(count.dat$med_income)
hist(count.dat$bachDeg25p)
hist(count.dat$vacRecHs) ## log1p
hist(count.dat$openWater) ## log1p
hist(count.dat$developed) ## log1p
hist(count.dat$forest) ## log1p
hist(count.dat$plantCult) ## log1p
hist(count.dat$wetland) ## log1p
hist(count.dat$area_km2) ## log10

# Center, scale, transform predictors
count.mod.df = count.dat

count.mod.df$popDens = scale(log10(count.mod.df$popDens))
count.mod.df$eduDens = scale(log1p(count.mod.df$eduDens))
count.mod.df$med_income = scale(count.mod.df$med_income)
count.mod.df$bachDeg25p = scale(count.mod.df$bachDeg25p)
count.mod.df$vacRecHs = scale(log1p(count.mod.df$vacRecHs))
count.mod.df$openWater = scale(log1p(count.mod.df$openWater))
count.mod.df$developed = scale(log1p(count.mod.df$developed))
count.mod.df$forest = scale(log1p(count.mod.df$forest))
count.mod.df$plantCult = scale(log1p(count.mod.df$plantCult))
count.mod.df$wetland = scale(log1p(count.mod.df$wetland))
count.mod.df$area_km2 = log(count.mod.df$area_km2) ## natural log transform area


# CORRELATIONS
corMat = cor(count.mod.df[ , c('popDens', 'eduDens', 'med_income', 
                               'bachDeg25p', 'vacRecHs', 
                               'openWater', 'developed', 'forest',
                               'plantCult', 'wetland')])
corMat[upper.tri(corMat, diag = T)] = NA
corMat = as.data.frame(corMat)

# write.csv(corMat, 'Output/corMat/countyCorr.csv',
#           row.names = T)


#####################################################################
## SELECTING OPTIMAL MODEL STRUCTURE

# poisson model	
county.fullmod.pois <- glmmTMB(nSites ~
                                 plantCult +
                                 eduDens +
                                 bachDeg25p +
                                 forest +
                                 med_income +
                                 openWater +
                                 popDens +
                                 vacRecHs +
                                 wetland +
								 #developed +
                                 (1 | state) + 
                                 offset(area_km2),
                               data = count.mod.df,
                               family = poisson(link = 'log'))

# neg bin 1 model									
county.fullmod.nb1 <- glmmTMB(nSites ~
                                plantCult +
                                eduDens +
                                bachDeg25p +
                                forest +
                                med_income +
                                openWater +
                                popDens +
                                vacRecHs +
                                wetland +
								#developed +
                                (1 | state) +
                                offset(area_km2),
                              data = count.mod.df, 
                              family = nbinom1(link = 'log'))

# neg bin 2 model									
county.fullmod.nb2 <- glmmTMB(nSites ~
                                plantCult +
                                eduDens +
                                bachDeg25p +
                                forest +
                                med_income +
                                openWater +
                                popDens +
                                vacRecHs +
                                wetland +
								#developed +
                                (1 | state) +
                                offset(area_km2),
                              data = count.mod.df, 
                              family = nbinom2(link = 'log'))

# poisson model	with zeroinfl
county.fullmod.zipois <- glmmTMB(nSites ~
                                   plantCult +
                                   eduDens +
                                   bachDeg25p +
                                   forest +
                                   med_income +
                                   openWater +
                                   popDens +
                                   vacRecHs +
                                   wetland +
								   #developed +
                                   (1 | state) + 
                                   offset(area_km2),
                                 ziformula = ~ 1,
                                 data = count.mod.df,
                                 family = poisson(link = 'log'))

# neg bin 1 model with zeroinfl									
county.fullmod.zinb1 <- glmmTMB(nSites ~
                                  plantCult +
                                  eduDens +
                                  bachDeg25p +
                                  forest +
                                  med_income +
                                  openWater +
                                  popDens +
                                  vacRecHs +
                                  wetland +
								  #developed +
                                  (1 | state) + 
                                  offset(area_km2),
                                ziformula = ~ 1,
                                data = count.mod.df, 
                                family = nbinom1(link = 'log'))

# neg bin 2 model with zeroinfl									
county.fullmod.zinb2 <- glmmTMB(nSites ~
                                  plantCult +
                                  eduDens +
                                  bachDeg25p +
                                  forest +
                                  med_income +
                                  openWater +
                                  popDens +
                                  vacRecHs +
                                  wetland +
								  #developed +
                                  (1 | state) + 
                                  offset(area_km2),
                                ziformula = ~ 1,
                                data = count.mod.df, 
                                family = nbinom2(link = 'log'))

## use AICc to check with error structure is the best
AICctab(county.fullmod.pois, 
        county.fullmod.nb1, 
        county.fullmod.nb2, 
        county.fullmod.zipois, 
        county.fullmod.zinb1, 
        county.fullmod.zinb2, 
        base = T, 
        weights = T)

#                       AICc     dAICc    df weight
# county.fullmod.nb2     22230.5      0.0 12 0.73  
# county.fullmod.zinb2   22232.5      2.0 13 0.27  
# county.fullmod.nb1     23702.4   1471.9 12 <0.001
# county.fullmod.zinb1   23704.4   1473.9 13 <0.001
# county.fullmod.zipois 141173.7 118943.2 12 <0.001
# county.fullmod.pois   143584.6 121354.1 11 <0.001

rm(list = setdiff(ls(), c('count.mod.df', 'count.dat', 'county.fullmod.nb2')))

#####################################################################
## MODELLING SITE COUNT :: MODEL SELECTION ON BEST-FIT MODEL

# All-subsets model selection - keep offset term in all models
mod.sel = dredge(county.fullmod.nb2,
                 rank = 'AICc', fixed='cond(offset(area_km2))')

# Rename model selection object (county=county level)				 
county.mod.sel <- mod.sel
			
# Top 10 models
county.mod.sel[1:10]
			
# Top model
county.mod = get.models(county.mod.sel, subset = 1)[[1]]

summary(county.mod)
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -5.33091    0.18449 -28.896  < 2e-16 ***
# bachDeg25p   0.33182    0.02733  12.140  < 2e-16 ***
# eduDens      0.07235    0.03407   2.123  0.03371 *  
# forest       0.28727    0.03989   7.202 5.93e-13 ***
# plantCult   -0.09437    0.03755  -2.513  0.01197 *  
# popDens      1.09338    0.03798  28.791  < 2e-16 ***
# vacRecHs     0.34386    0.03054  11.261  < 2e-16 ***
# wetland      0.07830    0.03000   2.610  0.00907 **

r.squaredGLMM(county.mod)
#                R2m       R2c
# delta     0.4153530 0.8002972
# lognormal 0.4391956 0.8462370
# trigamma  0.3764012 0.7252453

# Checking residuals
simulationOutput.ct = simulateResiduals(fittedModel = county.mod)
windows(9,5)
plot(simulationOutput.ct, quantreg=T) # for WebFigure 4

# Test for overdispersion
testDispersion(simulationOutput.ct)
# data:  simulationOutput
# dispersion = 0.58973, p-value = 0.648
# alternative hypothesis: two.sided

# Test zero inflation 
testZeroInflation(simulationOutput.ct)
# data:  simulationOutput
# ratioObsSim = 1.016, p-value = 0.888
# alternative hypothesis: two.sided

# Get coefficients, se, and R2marginal for top 10 models (WebPanel 3)
# e.g., for 2nd ranked model, use subset=2; for 3rd ranked model, use subset=3, etc.
county.mod.rank2 = get.models(county.mod.sel, subset = 2)[[1]]
summary(county.mod.rank2)

r.squaredGLMM(county.mod.rank2)


rm(list = ls())

#####################################################################
## ANALYSIS :: WATERSHED
#####################################################################

#####################################################################
## MODEL PREP :: CORRELATIONS, TRANSFORMATIONS

hu8.dat = read.csv('Input/analysisDFs/watershedAnalysis.csv',
                   header = T, stringsAsFactors = F)

# Format HUC8 with leading zeros
hu8.dat$huc8 = stringr::str_pad(hu8.dat$huc8 , 8, side = 'left', pad = '0')
hu8.dat$huc2 = stringr::str_sub(hu8.dat$huc8 , 1, 2)


# Transformations check
hist(hu8.dat$popDens)  # log10
hist(hu8.dat$eduDens) ## log1p
hist(hu8.dat$med_income)
hist(hu8.dat$bachDeg25p)
hist(hu8.dat$vacRecHs)
hist(hu8.dat$openWater) ## log1p
hist(hu8.dat$developed) ## log1p
hist(hu8.dat$forest) ## log1p
hist(hu8.dat$plantCult) ## log1p
hist(hu8.dat$wetland) ## log1p

# Center, scale, transform predictors
hu8.mod.df = hu8.dat

hu8.mod.df$popDens = scale(log10(hu8.mod.df$popDens))
hu8.mod.df$eduDens = scale(log1p(hu8.mod.df$eduDens))
hu8.mod.df$med_income = scale(hu8.mod.df$med_income)
hu8.mod.df$bachDeg25p = scale(hu8.mod.df$bachDeg25p)
hu8.mod.df$vacRecHs = scale(hu8.mod.df$vacRecHs)
hu8.mod.df$openWater = scale(log1p(hu8.mod.df$openWater))
hu8.mod.df$developed = scale(log1p(hu8.mod.df$developed))
hu8.mod.df$forest = scale(log1p(hu8.mod.df$forest))
hu8.mod.df$plantCult = scale(log1p(hu8.mod.df$plantCult))
hu8.mod.df$wetland = scale(log1p(hu8.mod.df$wetland))
hu8.mod.df$area_km2 = log(hu8.mod.df$area_km2)

corMat = cor(hu8.mod.df[ , c('popDens', 'eduDens', 'med_income', 
                               'bachDeg25p', 'vacRecHs', 
                               'openWater', 'developed', 'forest',
                               'plantCult', 'wetland')])
corMat[upper.tri(corMat, diag = T)] = NA
corMat = as.data.frame(corMat)

# write.csv(corMat, 'Output/corMat/huc8Corr.csv',
#           row.names = T)


#####################################################################
## SELECTING OPTIMAL MODEL STRUCTURE

# poisson model	
ws.fullmod.pois <- glmmTMB(nSites ~
                             plantCult +
                             eduDens +
                             bachDeg25p +
                             forest +
                             med_income +
                             openWater +
                             popDens +
                             vacRecHs +
                             wetland +
                             (1 | huc2) + 
                             offset(area_km2),
                           data = hu8.mod.df, 
                           family = poisson(link = 'log'))

# neg bin 1 model									
ws.fullmod.nb1 <- glmmTMB(nSites ~
                                  plantCult +
                                  eduDens +
                                  bachDeg25p +
                                  forest +
                                  med_income +
                                  openWater +
                                  popDens +
                                  vacRecHs +
                                  wetland +
                            (1 | huc2) + 
                            offset(area_km2),
                          data = hu8.mod.df, 
                          family = nbinom1(link = 'log'))

# neg bin 2 model									
ws.fullmod.nb2 <- glmmTMB(nSites ~
                                  plantCult +
                                  eduDens +
                                  bachDeg25p +
                                  forest +
                                  med_income +
                                  openWater +
                                  popDens +
                                  vacRecHs +
                                  wetland +
                            (1 | huc2) + 
                            offset(area_km2),
                          data = hu8.mod.df, 
                          family = nbinom2(link = 'log'))

# poisson model	with zeroinfl
ws.fullmod.zipois <- glmmTMB(nSites ~
                                     plantCult +
                                     eduDens +
                                     bachDeg25p +
                                     forest +
                                     med_income +
                                     openWater +
                                     popDens +
                                     vacRecHs +
                                     wetland +
                               (1 | huc2) + 
                               offset(area_km2),
                             ziformula = ~ 1,
                             data = hu8.mod.df, 
                             family = poisson(link = 'log'))

# neg bin 1 model with zeroinfl									
ws.fullmod.zinb1 <- glmmTMB(nSites ~
                                    plantCult +
                                    eduDens +
                                    bachDeg25p +
                                    forest +
                                    med_income +
                                    openWater +
                                    popDens +
                                    vacRecHs +
                                    wetland +
                              (1 | huc2) + 
                              offset(area_km2),
                            ziformula = ~ 1,
                            data = hu8.mod.df, 
                            family = nbinom1(link = 'log'))

# neg bin 2 model with zeroinfl									
ws.fullmod.zinb2 <- glmmTMB(nSites ~
                                    plantCult +
                                    eduDens +
                                    bachDeg25p +
                                    forest +
                                    med_income +
                                    openWater +
                                    popDens +
                                    vacRecHs +
                                    wetland +
                              (1 | huc2) + 
                              offset(area_km2),
                            ziformula = ~ 1,
                            data = hu8.mod.df,
                            family = nbinom2(link = 'log'))


## use AICc to check with error structure is the best
AICctab(ws.fullmod.pois, 
        ws.fullmod.nb1, 
        ws.fullmod.nb2, 
        ws.fullmod.zipois, 
        ws.fullmod.zinb1, 
        ws.fullmod.zinb2, 
        base = T, 
        weights = T)
#                 AICc     dAICc    df weight
# ws.fullmod.nb2     17233.5      0.0 12 0.73  
# ws.fullmod.zinb2   17235.5      2.0 13 0.27  
# ws.fullmod.nb1     17681.6    448.1 12 <0.001
# ws.fullmod.zinb1   17683.6    450.1 13 <0.001
# ws.fullmod.zipois 202909.8 185676.3 12 <0.001
# ws.fullmod.pois   210194.8 192961.3 11 <0.001

rm(list = setdiff(ls(), c('hu8.mod.df', 'hu8.dat', 'ws.fullmod.nb2')))

#####################################################################
## MODELLING SITE COUNT :: MODEL SELECTION ON BEST-FIT MODEL

# All-subsets model selection - keep offset term in all models
mod.sel = dredge(ws.fullmod.nb2,
                 rank = 'AICc', fixed='cond(offset(area_km2))')

# Rename model selection object (ws=watershed level)				 
ws.mod.sel <- mod.sel

# Top 10 models
ws.mod.sel[1:10]
				 
# Top model
ws.mod = get.models(ws.mod.sel, subset = 1)[[1]]
summary(ws.mod)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -5.00458    0.23935 -20.909  < 2e-16 ***
# bachDeg25p   0.23940    0.03620   6.613 3.76e-11 ***
# forest       0.57539    0.04563  12.611  < 2e-16 ***
# openWater    0.32989    0.03727   8.851  < 2e-16 ***
# popDens      1.08819    0.05274  20.635  < 2e-16 ***
# vacRecHs     0.07868    0.04957   1.587  0.11242    
# wetland      0.11784    0.04432   2.659  0.00785 **

r.squaredGLMM(ws.mod)
#               R2m       R2c
# delta     0.6445247 0.9195283
# lognormal 0.6509990 0.9287650
# trigamma  0.6362515 0.9077252

simulationOutput.ws = simulateResiduals(fittedModel = ws.mod)
windows(9,5)
plot(simulationOutput.ws, quantreg=T) # for WebFigure 4

# Test for overdispersion
testDispersion(simulationOutput.ws)
# data:  simulationOutput
# dispersion = 0.25298, p-value = 0.808
# alternative hypothesis: two.sided

# Test zero inflation 
testZeroInflation(simulationOutput.ws)
# data:  simulationOutput
# ratioObsSim = 1.1516, p-value = 0.288
# alternative hypothesis: two.sided

# Get coefficients, se, and R2marginal for top 10 models (WebPanel 3)
# e.g., for 2nd ranked model, use subset=2; for 3rd ranked model, use subset=3, etc.
ws.mod.rank2 = get.models(ws.mod.sel, subset = 2)[[1]]
summary(ws.mod.rank2)

r.squaredGLMM(ws.mod.rank2)
