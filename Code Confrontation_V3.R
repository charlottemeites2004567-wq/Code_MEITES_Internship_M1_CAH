################################################################################
#                                                                              #
# EFFECTS OF PRENATAL ODOR LEARNING ON SOCIAL BEHAVIOURS IN PIGLETS              #
# Data: Happy Smeling Confrontation Test                                       #
# ID: Charlotte Meites                                                         #
################################################################################


citation()

# Packages ------------------------------------------------------------------------------------------------------------------------

install.packages(c("lmerTest", "car", "RVAideMemoire", "effects", "emmeans", "DHARMa", "openxlsx"), dependencies = TRUE)

# --- DATA MANIPULATION ---
library(tidyverse)     # Groups dplyr, tidyr, ggplot2, etc. (loads the ecosystem)
library(tidyr)         # For data manipulation / restructuring
library(dplyr)         # To manipulate, filter, and arrange data

# --- GRAPHICS AND VISUALIZATION ---
library(ggplot2)       # To create graphics/plots
library(ggsci)         # For themes and color palettes (scientific journals)
library(scales)        # To format axes and display color scales
library(RColorBrewer)  # To obtain predefined color palettes
library(viridis)       # For color palettes adapted to color-blindness and black & white printing
library(ggpubr)        # To create publication-ready graphics (and add stats)

# --- STATISTICS: LINEAR AND MIXED MODELS ---
library(lme4)          # For mixed models (linear and non-linear)
library(lmerTest)      # To improve visualization of mixed models and get p-values
library(glmmTMB)       # For generalized mixed models 
library(car)           # To run type II or III ANOVAs
library(lmtest)        # To test linear model assumptions 

# --- DIAGNOSTICS AND POST-HOC ---
library(RVAideMemoire) # For the "plotresid" function 
library(DHARMa)        # For residual diagnostics of mixed models (GLMM) via simulations
library(emmeans)       # To estimate marginal means and run post-hoc tests (pairwise comparisons)
library(effects)       # To calculate and display prediction effects of models (lmer, lm)



# Choose working directory ----------------------------------------------------------------------------------------------------

setwd("~/M1/Stage/Données")

# I. Importation and formatting of datasets ------------------------------------------------------------------------------------------

data_confrontation <- read.csv("data_confrontation.csv", header=T, stringsAsFactors = TRUE, dec = ",")

str(data_confrontation)

# transform variable as factor

data_confrontation$replicate  <- as.factor(data_confrontation$replicate)         

str(data_confrontation)
summary(data_confrontation)
names(data_confrontation)

# Create frequencies per min on total obs (10 min):

tot_obs = 10

data_confrontation$f_nbr_scream = data_confrontation$nbr_scream / tot_obs
data_confrontation$f_nbr_squeak = data_confrontation$nbr_squeak / tot_obs
data_confrontation$f_nbr_grunt = data_confrontation$nbr_grunt / tot_obs
data_confrontation$f_nbr_nose_disc_to_snout  = data_confrontation$nbr_nose_disc_to_snout / tot_obs
data_confrontation$f_nbr_nose_disc_to_head = data_confrontation$nbr_nose_disc_to_head/ tot_obs
data_confrontation$f_nbr_nose_disc_to_body = data_confrontation$nbr_nose_disc_to_body / tot_obs
data_confrontation$f_nbr_exploring_together  = data_confrontation$nbr_exploring_together/ tot_obs
data_confrontation$f_nbr_nudging = data_confrontation$nbr_nudging/ tot_obs
data_confrontation$f_nbr_agression = data_confrontation$nbr_aggression / tot_obs
data_confrontation$f_nbr_mounting  = data_confrontation$nbr_mounting/ tot_obs

summary(data_confrontation)
names(data_confrontation)


cols_numeriques <- names(data_confrontation)[sapply(data_confrontation, is.numeric)]

#create means 

data_confrontation <- data_confrontation %>%
  filter(!is.na(treatment), treatment != "")

cols_numeriques <- names(data_confrontation)[sapply(data_confrontation, is.numeric)]

# Mean + sd + n + SEM per treatment
data_confrontation_moyen <- data_confrontation %>%
  group_by(odor_presence, familiarity, treatment) %>%
  summarize(
    across(
      all_of(cols_numeriques),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd   = ~ sd(.x, na.rm = TRUE),
        n    = ~ sum(!is.na(.x)),
        sem  = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

str(data_confrontation_moyen)

# II. Descriptive analysis of data ------------------------------------------------------------------------------

# III. Statistical modeling and model validation -----------------------------------------------------------------------

# 5 models were tested : 
# treatment + odor_presence + familiarity + sexe)^4 + running_order + (1 | replicate)
# treatment + odor_presence + familiarity + sexe)^3 + running_order +(1 | replicate),data = data_confrontation)
# treatment + odor_presence + familiarity + sexe)^2 + running_order + (1 | replicate)
# treatment + odor_presence + familiarity)^2 + sexe + running_order + (1 | replicate)
# treatment*odor_presence + familiarity + sexe + running_order + (1 | replicate)

#final model : treatment + odor_presence + familiarity)^2 + sexe + running_order +(1 | replicate)
# other models didn't show any significant interactions

### scream frequency ---------------------------------------------------------------------------------------------------------------

# LM Models 

data_confrontation$f_nbr_scream_log = log(data_confrontation$f_nbr_scream + 1)
data_confrontation$f_nbr_scream_sqrt = sqrt(data_confrontation$f_nbr_scream)

mod_lm_scream_f =lm(f_nbr_scream~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed NO
# sqrt NO
# log NO

shapiro.test(residuals(mod_lm_scream_f))

par(mfrow=c(1,2))
hist (residuals(mod_lm_scream_f)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_scream_f)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_scream_f))

par(mfrow=c(1,1))
plot(residuals(mod_lm_scream_f)~fitted(mod_lm_scream_f)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 


# GLM model

mod_glmm_scream_f =glmmTMB(f_nbr_scream~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation, family = tweedie(link = "log"))

# GLM model validity

res_scream <- simulateResiduals(mod_glmm_scream_f, n = 1000)
plot(res_scream) # QQplot ok but issue with residuals vs predicted: decreasing trend --> inherent limitation of vocalization data (30% zeros, highly variable frequencies)
testDispersion(res_scream) #dispersion = 0.88 , p-value = 0.98  (dispersion around 1 and p>0.05) ok

par(mfrow=c(1,3))
plotResiduals(res_scream, data_confrontation$treatment) #ok
plotResiduals(res_scream, data_confrontation$odor_presence)#ok
plotResiduals(res_scream, data_confrontation$familiarity) #ok

par(mfrow=c(1,3))
plotResiduals(res_scream, data_confrontation$sexe)#ok
plotResiduals(res_scream, data_confrontation$running_order)
plotResiduals(res_scream, data_confrontation$replicate)#ok

#Final model 

mod_final_glm_scream_f = glmmTMB(f_nbr_scream~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation, family = tweedie(link = "log"))

summary(mod_final_glm_scream_f)
Anova(mod_final_glm_scream_f, type = 3)

"                            Chisq Df Pr(>Chisq)    
(Intercept)               25.4666  1  4.501e-07 ***
treatment                  0.8540  1  0.3554226    
odor_presence              9.2007  1  0.0024192 ** familiarity                6.6546  1  0.0098899 ** sexe                       0.0018  1  0.9664847    
running_order              5.8745  1  0.0153614 * replicate                  0.0216  1  0.8831134    
treatment:odor_presence    0.6905  1  0.4060057    
treatment:familiarity      0.9179  1  0.3380345    
odor_presence:familiarity 13.6815  1  0.0002166 ***"

mod_final_glm_scream_f_means <- emmeans(mod_final_glm_scream_f, pairwise ~ familiarity | odor_presence, adjust = "tukey")
mod_final_glm_scream_f_means$emmeans
mod_final_glm_scream_f_means$contrasts

"odor_presence = no_odour:
 contrast              estimate    SE  df z.ratio p.value
 familiar - unfamiliar    -1.27 0.512 Inf  -2.488  0.0129

odor_presence = odour:
 contrast              estimate    SE  df z.ratio p.value
 familiar - unfamiliar     1.36 0.508 Inf   2.671  0.0076"

### squeak frequency ---------------------------------------------------------------------------------------------------------------

# LM Models 

mod_lm_squeak_f =lm(f_nbr_squeak~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed Yes

shapiro.test(residuals(mod_lm_squeak_f))

par(mfrow=c(1,2))
hist (residuals(mod_lm_squeak_f)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_squeak_f)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_squeak_f))

par(mfrow=c(1,1))
plot(residuals(mod_lm_squeak_f)~fitted(mod_lm_squeak_f)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 

#Final model 

mod_final_lm_squeak_f =lm(f_nbr_squeak~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

summary(mod_final_lm_squeak_f)
Anova(mod_final_lm_squeak_f, type = 3)

"                           Sum Sq Df F value    Pr(>F)    
(Intercept)               1.5187  1 27.2594 1.221e-06 ***
treatment                  0.0170  1  0.3055    0.5819    
odor_presence              0.0000  1  0.0000    0.9958    
familiarity                0.1143  1  2.0519    0.1556    
sexe                       0.0547  1  0.9825    0.3244    
running_order              0.0007  1  0.0133    0.9086    
replicate                  0.0182  1  0.3267    0.5691    
treatment:odor_presence    0.0299  1  0.5360    0.4661    
treatment:familiarity      0.0029  1  0.0517    0.8206    
odor_presence:familiarity 0.0016  1  0.0294    0.8643    
Residuals                  4.7913 86  "

### grunt frequency ---------------------------------------------------------------------------------------------------------------

# LM Models 

mod_lm_grunt_f =lm(f_nbr_grunt~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed YES

shapiro.test(residuals(mod_lm_grunt_f))

par(mfrow=c(1,2))
hist (residuals(mod_lm_grunt_f)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_grunt_f)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_grunt_f))

par(mfrow=c(1,1))
plot(residuals(mod_lm_grunt_f)~fitted(mod_lm_grunt_f)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 

#Final model 

mod_final_lm_grunt_f =lm(f_nbr_grunt~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

summary(mod_final_lm_grunt_f)
Anova(mod_final_lm_grunt_f, type = 3)

"                           Sum Sq Df  F value    Pr(>F)    
(Intercept)               16.3179  1 102.2847 2.681e-16 ***
treatment                  0.4280  1   2.6830    0.1051    
odor_presence              0.0558  1   0.3497    0.5559    
familiarity                0.3554  1   2.2279    0.1392    
sexe                       0.1475  1   0.9245    0.3390    
running_order              0.0506  1   0.3171    0.5748    
replicate                  0.4030  1   2.5259    0.1157    
treatment:odor_presence    0.0448  1   0.2809    0.5975    
treatment:familiarity      0.0030  1   0.0189    0.8909    
odor_presence:familiarity  0.0012  1   0.0078    0.9298    
Residuals                 13.7199 86  "


### nose disc to snout frequency ---------------------------------------------------------------------------------------------------------------

# LM Models 

data_confrontation$f_nbr_nose_disc_to_snout_log = log(data_confrontation$f_nbr_nose_disc_to_snout + 1)
data_confrontation$f_nbr_nose_disc_to_snout_sqrt = sqrt(data_confrontation$f_nbr_nose_disc_to_snout)

mod_lm_nose_disc_to_snout_f =lm(f_nbr_nose_disc_to_snout_sqrt~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed NO
# sqrt YES
# log NO

shapiro.test(residuals(mod_lm_nose_disc_to_snout_f))

par(mfrow=c(1,2))
hist (residuals(mod_lm_nose_disc_to_snout_f)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_nose_disc_to_snout_f)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_nose_disc_to_snout_f))

par(mfrow=c(1,1))
plot(residuals(mod_lm_nose_disc_to_snout_f)~fitted(mod_lm_nose_disc_to_snout_f)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 


#Final model 

mod_final_lm_nose_disc_to_snout_f =lm(f_nbr_nose_disc_to_snout_sqrt~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

summary(mod_final_lm_nose_disc_to_snout_f)
Anova(mod_final_lm_nose_disc_to_snout_f, type = 3)

"                             Sum Sq Df F value    Pr(>F)    
(Intercept)               0.067345  1 32.6000 1.583e-07 ***
treatment                 0.000867  1  0.4197  0.518826    
odor_presence             0.000362  1  0.1755  0.676345    
familiarity               0.009116  1  4.4127  0.038599 * sexe                      0.000084  1  0.0405  0.841024    
running_order             0.013505  1  6.5376  0.012317 * replicate                 0.016057  1  7.7726  0.006528 ** treatment:odor_presence   0.000642  1  0.3109  0.578551    
treatment:familiarity     0.001933  1  0.9359  0.336049    
odor_presence:familiarity 0.000129  1  0.0626  0.803048    
Residuals                 0.177658 86                       "

### nose disc to snout latency ---------------------------------------------------------------------------------------------------------------

# LM Models 

data_confrontation$latency_nose_disc_to_snout_log = log(data_confrontation$latency_nose_disc_to_snout + 1)
data_confrontation$latency_nose_disc_to_snout_sqrt = sqrt(data_confrontation$latency_nose_disc_to_snout)

mod_lm_nose_disc_to_snout_l =lm(latency_nose_disc_to_snout_log~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed NO
# sqrt NO
# log OK

shapiro.test(residuals(mod_lm_nose_disc_to_snout_l))

par(mfrow=c(1,2))
hist (residuals(mod_lm_nose_disc_to_snout_l)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_nose_disc_to_snout_l)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_nose_disc_to_snout_l))

par(mfrow=c(1,1))
plot(residuals(mod_lm_nose_disc_to_snout_l)~fitted(mod_lm_nose_disc_to_snout_l)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 

#Final model 

mod_final_lm_nose_disc_to_snout_l =lm(latency_nose_disc_to_snout_log~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

summary(mod_final_lm_nose_disc_to_snout_l)
Anova(mod_final_lm_nose_disc_to_snout_l, type = 3)

"                           Sum Sq Df  F value  Pr(>F)    
(Intercept)               142.110  1 145.1530 < 2e-16 ***
treatment                   0.875  1   0.8935 0.34719    
odor_presence               0.444  1   0.4532 0.50264    
familiarity                 0.000  1   0.0000 0.99664    
sexe                        0.833  1   0.8506 0.35895    
running_order               1.424  1   1.4546 0.23110    
replicate                   0.011  1   0.0110 0.91686    
treatment:odor_presence     3.363  1   3.4354 0.06724 .  
treatment:familiarity       0.002  1   0.0025 0.96016    
odor_presence:familiarity   0.032  1   0.0332 0.85594    
Residuals                  84.197 86  "



### nose disc to body frequency ---------------------------------------------------------------------------------------------------------------

# LM Models 

data_confrontation$f_nbr_nose_disc_to_body_log = log(data_confrontation$f_nbr_nose_disc_to_body + 1)
data_confrontation$f_nbr_nose_disc_to_body_sqrt = sqrt(data_confrontation$f_nbr_nose_disc_to_body)

mod_lm_nose_disc_to_body_f =lm(f_nbr_nose_disc_to_body_sqrt~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed NO
# sqrt YES
# log NO

shapiro.test(residuals(mod_lm_nose_disc_to_body_f))

par(mfrow=c(1,2))
hist (residuals(mod_lm_nose_disc_to_body_f)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_nose_disc_to_body_f)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_nose_disc_to_body_f))

par(mfrow=c(1,1))
plot(residuals(mod_lm_nose_disc_to_body_f)~fitted(mod_lm_nose_disc_to_body_f)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 


#Final model 

mod_final_lm_nose_disc_to_body_f =lm(f_nbr_nose_disc_to_body_sqrt~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

summary(mod_final_lm_nose_disc_to_body_f)
Anova(mod_final_lm_nose_disc_to_body_f, type = 3)

"                Sum Sq Df F value    Pr(>F)    
(Intercept)               0.35642  1 76.5449 1.596e-13 ***
treatment                 0.00700  1  1.5028   0.22359    
odor_presence             0.00798  1  1.7134   0.19404    
familiarity               0.00473  1  1.0165   0.31618    
sexe                      0.01110  1  2.3832   0.12632    
running_order             0.01328  1  2.8529   0.09483 .  
replicate                 0.01413  1  3.0345   0.08509 .  
treatment:odor_presence   0.00668  1  1.4347   0.23430    
treatment:familiarity     0.01110  1  2.3832   0.12632    
odor_presence:familiarity 0.00258  1  0.5548   0.45841    
Residuals                 0.40045 86   
Residuals                 0.177658 86                       "

### nose disc to body latency ---------------------------------------------------------------------------------------------------------------

# LM Models 

data_confrontation$latency_nose_disc_to_body_log = log(data_confrontation$latency_nose_disc_to_body + 1)
data_confrontation$latency_nose_disc_to_body_sqrt = sqrt(data_confrontation$latency_nose_disc_to_body)

mod_lm_nose_disc_to_body_l =lm(latency_nose_disc_to_body_log~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed NO
# sqrt NO
# log OK

shapiro.test(residuals(mod_lm_nose_disc_to_body_l))

par(mfrow=c(1,2))
hist (residuals(mod_lm_nose_disc_to_snout_l)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_nose_disc_to_snout_l)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_nose_disc_to_snout_l))

par(mfrow=c(1,1))
plot(residuals(mod_lm_nose_disc_to_snout_l)~fitted(mod_lm_nose_disc_to_snout_l)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 

#Final model 

mod_final_lm_nose_disc_to_body_l =lm(latency_nose_disc_to_body_log~ (treatment + odor_presence + familiarity)^2  + sexe + running_order + replicate, data = data_confrontation)

summary(mod_final_lm_nose_disc_to_body_l)
Anova(mod_final_lm_nose_disc_to_body_l, type = 3)

"                    Sum Sq Df F value    Pr(>F)    
(Intercept)               61.834  1 74.6102 2.689e-13 ***
treatment                  0.017  1  0.0208    0.8857    
odor_presence              1.197  1  1.4446    0.2327    
familiarity                0.000  1  0.0000    0.9967    
sexe                       0.713  1  0.8604    0.3562    
running_order              1.806  1  2.1795    0.1435    
replicate                  0.219  1  0.2647    0.6082    
treatment:odor_presence    0.065  1  0.0783    0.7803    
treatment:familiarity      0.323  1  0.3901    0.5339    
odor_presence:familiarity  0.025  1  0.0306    0.8615    
Residuals                 71.273 86 "

### nose disc to head frequency ---------------------------------------------------------------------------------------------------------------

# LM Models 

data_confrontation$f_nbr_nose_disc_to_head_log = log(data_confrontation$f_nbr_nose_disc_to_head + 1)
data_confrontation$f_nbr_nose_disc_to_head_sqrt = sqrt(data_confrontation$f_nbr_nose_disc_to_head)

mod_lm_nose_disc_to_head_f =lm(f_nbr_nose_disc_to_head_sqrt~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed NO
# sqrt YES
# log NO

shapiro.test(residuals(mod_lm_nose_disc_to_head_f))

par(mfrow=c(1,2))
hist (residuals(mod_lm_nose_disc_to_head_f)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_nose_disc_to_head_f)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_nose_disc_to_head_f))

par(mfrow=c(1,1))
plot(residuals(mod_lm_nose_disc_to_head_f)~fitted(mod_lm_nose_disc_to_head_f)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 

#Final model 

mod_final_lm_nose_disc_to_head_f =lm(f_nbr_nose_disc_to_head_sqrt~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

summary(mod_final_lm_nose_disc_to_head_f)
Anova(mod_final_lm_nose_disc_to_head_f, type = 3)

"                             Sum Sq Df F value    Pr(>F)    
(Intercept)               0.163699  1 46.7645 1.102e-09 ***
treatment                 0.017276  1  4.9352  0.028942 * odor_presence             0.009203  1  2.6292  0.108577    
familiarity               0.017191  1  4.9110  0.029327 * sexe                      0.004091  1  1.1686  0.282714    
running_order             0.007612  1  2.1746  0.143961    
replicate                 0.040135  1 11.4655  0.001071 ** treatment:odor_presence   0.003851  1  1.1002  0.297166    
treatment:familiarity     0.007668  1  2.1907  0.142504    
odor_presence:familiarity 0.001217  1  0.3478  0.556908    
Residuals                 0.301043 86  "



### nose disc to head latency ---------------------------------------------------------------------------------------------------------------

# LM Models 
data_confrontation$latency_nose_disc_to_head_log = log(data_confrontation$latency_nose_disc_to_head + 1)
data_confrontation$latency_nose_disc_to_head_sqrt = sqrt(data_confrontation$latency_nose_disc_to_head)

mod_lm_nose_disc_to_head_l =lm(latency_nose_disc_to_head_sqrt~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed NO
#sqrt YES
#log NO

shapiro.test(residuals(mod_lm_nose_disc_to_head_l))

par(mfrow=c(1,2))
hist (residuals(mod_lm_nose_disc_to_head_l)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_nose_disc_to_head_l)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_nose_disc_to_head_l))

par(mfrow=c(1,1))
plot(residuals(mod_lm_nose_disc_to_head_l)~fitted(mod_lm_nose_disc_to_head_l)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 

#Final model 

mod_final_lm_nose_disc_to_head_l =lm(latency_nose_disc_to_head_sqrt~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

summary(mod_final_lm_nose_disc_to_head_l)
Anova(mod_final_lm_nose_disc_to_head_l, type = 3)

"               Sum Sq Df F value    Pr(>F)    
(Intercept)               158.55  1 24.6457 3.451e-06 ***
treatment                   4.16  1  0.6463   0.42366    
odor_presence              24.59  1  3.8226   0.05381 .  
familiarity                14.83  1  2.3051   0.13262    
sexe                        2.74  1  0.4258   0.51579    
running_order              15.30  1  2.3788   0.12667    
replicate                  19.52  1  3.0339   0.08512 .  
treatment:odor_presence     0.00  1  0.0002   0.98965    
treatment:familiarity      15.94  1  2.4783   0.11910    
odor_presence:familiarity   7.74  1  1.2027   0.27584    
Residuals                 553.25 86  "

### exploring together frequency ---------------------------------------------------------------------------------------------------------------

# LM Models 

mod_lm_exploring_together_f =lm(f_nbr_exploring_together~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed YES

shapiro.test(residuals(mod_lm_exploring_together_f))

par(mfrow=c(1,2))
hist (residuals(mod_lm_exploring_together_f)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_exploring_together_f)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_exploring_together_f))

par(mfrow=c(1,1))
plot(residuals(mod_lm_exploring_together_f)~fitted(mod_lm_exploring_together_f)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 

#Final model 

mod_final_lm_exploring_together_f =lm(f_nbr_exploring_together~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

summary(mod_final_lm_exploring_together_f)
Anova(mod_final_lm_exploring_together_f, type = 3)

"                  Sum Sq Df F value   Pr(>F)    
(Intercept)               0.013558  1 35.2981 5.86e-08 ***
treatment                 0.000337  1  0.8778  0.35144    
odor_presence             0.001132  1  2.9484  0.08956 .  
familiarity               0.000250  1  0.6503  0.42225    
sexe                      0.000154  1  0.4015  0.52799    
running_order             0.000848  1  2.2065  0.14109    
replicate                 0.000229  1  0.5973  0.44173    
treatment:odor_presence   0.000133  1  0.3463  0.55775    
treatment:familiarity     0.000017  1  0.0455  0.83167    
odor_presence:familiarity 0.000169  1  0.4398  0.50900    
Residuals                 0.033033 86     "



### exploring together duration ---------------------------------------------------------------------------------------------------------------

# LM Models 

data_confrontation$duration_exploring_together_log = log(data_confrontation$duration_exploring_together + 1)
data_confrontation$duration_exploring_together_sqrt = sqrt(data_confrontation$duration_exploring_together)


mod_lm_exploring_together_d =lm(duration_exploring_together_log~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed NO
# sqrt NO
# log NO

shapiro.test(residuals(mod_lm_exploring_together_d))

par(mfrow=c(1,2))
hist (residuals(mod_lm_exploring_together_d)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_exploring_together_d)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_exploring_together_d))

par(mfrow=c(1,1))
plot(residuals(mod_lm_exploring_together_d)~fitted(mod_lm_exploring_together_d)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 

# GLM model

mod_glm_exploring_together_d =glmmTMB(duration_exploring_together~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation,family = tweedie(link = "log") )

# GLM model validity

res_scream <- simulateResiduals(mod_glm_exploring_together_d, n = 1000)
plot(res_scream) # QQplot and DHARMa ok
testDispersion(res_scream) #dispersion = 0.85 , p-value = 0.47  (dispersion around 1 and p>0.05) ok

par(mfrow=c(1,3))
plotResiduals(res_scream, data_confrontation$treatment) #ok
plotResiduals(res_scream, data_confrontation$odor_presence)#ok
plotResiduals(res_scream, data_confrontation$familiarity) #ok

par(mfrow=c(1,3))
plotResiduals(res_scream, data_confrontation$sexe)#ok
plotResiduals(res_scream, data_confrontation$running_order)
plotResiduals(res_scream, data_confrontation$replicate)#ok

#Final model 

mod_final_glm_exploring_together_d =glmmTMB(duration_exploring_together~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation,family = tweedie(link = "log") )

summary(mod_final_glm_exploring_together_d)
Anova(mod_final_glm_exploring_together_d, type = 3)

"                        Chisq Df Pr(>Chisq)    
(Intercept)               422.0883  1    < 2e-16 ***
treatment                   0.0261  1     0.87162    
odor_presence               0.7425  1     0.38888    
familiarity                 0.9561  1     0.32817    
sexe                        0.8643  1     0.35253    
running_order               3.0643  1     0.08003 .  
replicate                   3.3879  1     0.06568 .  
treatment:odor_presence     0.1080  1     0.74239    
treatment:familiarity       0.1007  1     0.75099    
odor_presence:familiarity   0.3424  1     0.55846   "



### nudging frequency  ---------------------------------------------------------------------------------------------------------------

# LM Models 

data_confrontation$f_nbr_nudging_log = log(data_confrontation$f_nbr_nudging + 1)
data_confrontation$f_nbr_nudging_sqrt = sqrt(data_confrontation$f_nbr_nudging)


mod_lm_nudging_f=lm(f_nbr_nudging~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation)

# Lm Models Validity 
# non-transformed NO
# sqrt NO
# log NO

shapiro.test(residuals(mod_lm_nudging_f))

par(mfrow=c(1,2))
hist (residuals(mod_lm_nudging_f)
      , col='red'
      , xlab='Residual values'
      , ylab='Counts')
qqnorm(residuals(mod_lm_nudging_f)
       , col='red'
       ,pch=16)
qqline(residuals(mod_lm_nudging_f))

par(mfrow=c(1,1))
plot(residuals(mod_lm_nudging_f)~fitted(mod_lm_nudging_f)
     , col='red'
     , pch=16
     , xlab = "Fitted values",
     ylab = "Residuals", 
     main = "Homogeneity?")
abline(h = 0, v = 0, lty = 2) 

# GLM model

mod_glm_nudging_f =glmmTMB(f_nbr_nudging~ (treatment + odor_presence + familiarity)^2 + sexe + running_order + replicate, data = data_confrontation,family = tweedie(link = "log") )

# GLM model validity YES

res_scream <- simulateResiduals(mod_glm_nudging_f, n = 1000)
plot(res_scream) # QQplot and DHARMa ok
testDispersion(res_scream) #dispersion = 0.82 , p-value = 0.37  (dispersion around 1 and p>0.05) ok

par(mfrow=c(1,3))
plotResiduals(res_scream, data_confrontation$treatment) #ok
plotResiduals(res_scream, data_confrontation$odor_presence)#ok
plotResiduals(res_scream, data_confrontation$familiarity) #ok

par(mfrow=c(1,3))
plotResiduals(res_scream, data_confrontation$sexe)#ok
plotResiduals(res_scream, data_confrontation$running_order)