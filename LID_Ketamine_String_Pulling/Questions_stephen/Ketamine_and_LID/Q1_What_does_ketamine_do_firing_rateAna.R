##################################
# Analyze da cells
##################################
library(tidyverse)
library(officer)
library(rvg)
library(GGally)
require(ggpubr)
rm(list = ls()) # this clears variables in case any are already loaded. 
old <- theme_set(theme_minimal(base_size = 18)) # No grid lines: I like this. Clean.
clrs = c('#0072BD','#D95319','#EDB120','#7E2F8E') # Use rgb2hex in matlab. This is the matlab lines colormap.
clr_inc_dec = c('#0072BD','#D95319') # Use rgb2hex in matlab. This is the matlab lines colormap.
#clr_DA_nonDA = c('#7E2F8E','#77AC30')


To = read_csv('C:/Temp/Q1Ket.csv')

# To = read_csv('C:/Users/Stephen Cowen/Dropbox/Foldershare/Analysis_Results_Dropbox/DANA_VTA/Q1_DA_cell_response_transient/Q1_DA_cell.csv')
To$NeuronID = factor(To$NeuronID)
To$Session = factor(To$Session)
To$SilentBase = To$Frate_base < .01
To$SilentPost1 = To$Frate_post1 < .01
To$SilentPost2 = To$Frate_post2 < .01


# Summarize firing rates...
ggplot(To, aes(x = Frate_base)) + geom_histogram()

# Do firing rates change?
ggplot(To, aes(x = NeuronID, y = Frate_post1mbase)) + geom_point()
ggplot(To, aes(x = Frate_post1mbase)) + geom_histogram()
median(To$Frate_post1mbase,na.rm = T)
wilcox.test(To$Frate_post1mbase)
ggplot(To, aes(x = Frate_post2mbase)) + geom_histogram()
median(To$Frate_post2mbase,na.rm = T)
wilcox.test(To$Frate_post2mbase)
# Do number of silent cells change?
ns = c(B = sum(To$Frate_base < .01), P1 = sum(To$Frate_post1 < .01), P2 = sum(To$Frate_post2 < .01))
barplot(ns)
chisq.test(ns)

ggplot(To, aes(x = Frate_post2mbase)) + geom_histogram()


# Does Local Variance change?
# Doing the following eliminates all of the neurons involved - almost. MANY LV very small. 
To$LocVar_base[To$LocVar_base<0.1] = NA
To$LocVar_post1mbase[To$LocVar_post1<0.1] = NA
To$LocVar_post2mbase[To$LocVar_post2<0.1] = NA

ggplot(To, aes(x = LocVar_base)) + geom_histogram()

ggplot(To, aes(x = NeuronID, y = LocVar_post1mbase)) + geom_point()
ggplot(To, aes(x = LocVar_post1mbase)) + geom_histogram()
median(To$LocVar_post1mbase,na.rm = T)
wilcox.test(To$LocVar_post1mbase)
ggplot(To, aes(x = LocVar_post2mbase)) + geom_histogram()
median(To$LocVar_post2mbase,na.rm = T)
wilcox.test(To$LocVar_post2mbase)

# Is there a relationship between the effect on movement and the effect on firing rates?
ggplot(To, aes(x = Frate_post1mbase, y = Speed_change)) + geom_point() + geom_smooth(method = 'lm')
mod = lm(Speed_change~Frate_post1mbase + Frate_post2mbase, data = To)
summary(mod)
cor.test(To$Speed_change, To$Frate_post1mbase, method = "spearman", continuity = FALSE, conf.level = 0.95)

ggplot(To, aes(x = Frate_post1mbase, y = Jerk_change)) + geom_point() + geom_smooth(method = 'lm')
mod = lm(Jerk_change~Frate_post1mbase + Frate_post2mbase, data = To)
summary(mod)
cor.test(To$Jerk_change, To$Frate_post1mbase, method = "spearman", continuity = FALSE, conf.level = 0.95)

