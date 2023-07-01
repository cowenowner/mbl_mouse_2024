# Determine if some drugs, and especially prevastatin, reduced LID
# Data is from Mitch Bartlett and Torsten Falk 2022
# Mitch just updated the data on 3/20/2023
#
# Cowen 2023
#
#
# Goals: Do all 9 (0 7 10 for V K KP, V K KL, V K P ) combos but use 
# kruskal-Wallis followed by Wilcoxon for post-hoc with Holm correction.
# 
#
#
# As of 10/11/2022: These data are preliminary so hold off on analysis as
# Torsten just informed me that they are collecting more data. As a result,
# hypothesis testing is a bit premature.
# 
#
# from Torsten: So, the main prediction for this study was that P will block the
# long-term anti-LID effect of K, but not the immediate anti-LID effect, based
# on literature from K and depression.  After we saw the clear effect of P to
# block the long-term K effect, we wanted to do a second study to test a
# different type of statins (hydrophilic vs hydrophobic) to see if this effect
# is specific to the drug (P) or the class of statin drugs.  Here the prediction
# is less clear.  The second prediction was that L might behave differently and
# not block K, based on its ability to be anti-LID itself in established LID.
# Therefore, you are correct we could focus on day 14, as the main long-term
# time point, which might allow to get significance, that is not lost in a
# million comparisons.  And we could have a separate analysis of just day 0 to
# look at that unexpected sensitization effect of P.
# P = Pravastatine L = lovastatin - so will need 9 graphs.
# Do this for days 0, 7, 10 ONLY

library(tidyverse)
library(readxl)
library(ggpubr)
#install.packages("ggsignif")                      # Install ggsignif package
library("ggsignif")  
theme_set(theme_classic())


# Load data 
path = 'C:/Users/cowen/Documents/GitHub/LID_Ketamine_String_Pulling/Questions_stephen/Drugs_and_LID_2022/'
fname = 'KETStat_LAO Summary_mjb032023.xlsx'
full_fname = paste0(path,fname)
sheets = c('Day 0', 'Day 7', 'Day 14')
D = tibble()
for (iS in sheets)
{
  TMP <- read_excel(full_fname, sheet = iS, range = cell_cols("A:F"))
  TMP$Day = iS
  colnames(TMP)[1] <- colnames(TMP)[1][1]
  colnames(TMP)[2] <- colnames(TMP)[2][1]
  colnames(TMP)[3] <- colnames(TMP)[3][1]
  colnames(TMP)[4] <- 'KP'
  colnames(TMP)[5] <- colnames(TMP)[5][1]
  colnames(TMP)[6] <- 'KL'
  TMP$Subject = 1:length(TMP$V)
  D = rbind(D,TMP[,1:8])
}
D$Day = factor(D$Day)
# Great, now I have to make it long form and add a subject ID.
# Blech.
D_long <- gather(D, condition, AIMS, V:KL)
# Reorder the categories
D_long$condition <-factor(D_long$condition,levels = c('V','K','KP','KL','P','L'))
factor(D_long$condition)
# Pull out the days of interest.
D_long14_K_KP <- D_long %>% filter(Day == 'Day 14', condition %in% c('K', 'KP') )
D_long14_K_KL <- D_long %>% filter(Day == 'Day 14', condition %in% c('K', 'KL') )
D_long14_K_KP_KL <- D_long %>% filter(Day == 'Day 14', condition %in% c('K', 'KP', 'KL') )
D_long14_K_KP_V <- D_long %>% filter(Day == 'Day 14', condition %in% c('K', 'KP', 'V') )
# Visualize the entire dataset
ggplot(D_long,aes(Day, AIMS, color = condition)) + geom_boxplot() + geom_jitter(position = position_jitterdodge(jitter.width = .08), size = 3,  aes(color = condition))

# Compare and plot the K and KP groups...
ggplot(D_long14_K_KP,aes(condition, AIMS)) + geom_boxplot() + geom_jitter(width = .14, size = 4)
# Do the stats...
x <- D_long14_K_KP %>% filter(condition == 'K')
shapiro.test(x$AIMS) # test for normality. If p > 0.05 it's normal enough.
# As of 3/2023, it is.
t.test(AIMS~condition,data = D_long14_K_KP) # was considering Wilcox, but these data do not look bounded (the other groups have multiple points clustered at 0 or 40 - but not these.)

# look at K, KP, and V (I think that this is the full comparison - really should include vehicle.)
mod = aov(AIMS ~ condition, data = D_long14_K_KP_V)
summary(mod)
TukeyHSD(mod, conf.level=.95)

ggplot(data = D_long14_K_KP_V, aes(condition, AIMS)) + geom_boxplot() + geom_jitter(width = .14, size = 4)

mod = aov(AIMS ~ condition, data = D_long %>% filter(Day == day_txt, condition %in% c('V', 'K', 'KP')))
summary(mod)
TukeyHSD(mod, conf.level=.95)

D <-D_long %>% filter(Day == 'Day 7', condition %in% c('V')) 



# Now generate the figures for torsten.

uplim = 40
day_txt = 'Day 0'
txt_size = 16
ylab_txt = 'L-DOPA-Induced AIMs'

my_comparisons <- list( c('V', 'K'), c('K', 'KP'), c("V", "KP") )
D <-D_long %>% filter(Day == day_txt, condition %in% c('V', 'K', 'KP')) 
a0 <- ggplot(data = D, aes(condition, AIMS, fill = condition, line = 'black')) + 
  stat_summary(fun = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar",  width=.3, size = 1) + 
  #stat_compare_means(method = "t.test",comparisons = my_comparisons, label = "p.signif",label.y = c(33, 35, 37)) +
  scale_fill_manual(values=c('#7B7B7B','#0C65FE','#E10004')) + coord_cartesian(ylim=c(0, uplim)) + 
  labs(title=" ",x="", y = "",size = txt_size) + theme(text = element_text(size = 16),legend.position="none",plot.title = element_text(hjust = 0.5),axis.text.x=element_blank())

my_comparisons <- list( c('V', 'K'), c('K', 'KL'), c("V", "KL") )
b0 <- ggplot(data = D_long %>% filter(Day == day_txt, condition %in% c('V', 'K', 'KL')) , aes(condition, AIMS, fill = condition, line = 'black')) + 
  stat_summary(fun = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar",  width=.3, size = 1) + 
  #stat_compare_means(method = "t.test",comparisons = my_comparisons, label = "p.signif",label.y = c(33, 35, 37)) +
  scale_fill_manual(values=c('#7B7B7B','#0C65FE','#FA6B0C')) + coord_cartesian(ylim=c(0, uplim)) + 
  labs(title=day_txt, x="", y = "",size = txt_size) + theme(text = element_text(size = 16),legend.position="none",plot.title = element_text(hjust = 0.5),axis.text.x=element_blank())

my_comparisons <- list( c('V', 'P'), c('P', 'L'), c("V", "L") )
c0 <- ggplot(data = D_long %>% filter(Day == day_txt, condition %in% c('V', 'P', 'L')) , aes(condition, AIMS, fill = condition, line = 'black')) + 
  stat_summary(fun = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar",  width=.3, size = 1) + 
  #stat_compare_means(method = "t.test",comparisons = my_comparisons, label = "p.signif",label.y = c(33, 35, 37)) +
  scale_fill_manual(values=c('#7B7B7B','#688035','#69006B')) + coord_cartesian(ylim=c(0, uplim)) + 
  labs(title=" ",x="", y = "",size = txt_size) + theme(text = element_text(size = 16),legend.position="none",plot.title = element_text(hjust = 0.5),axis.text.x=element_blank())


day_txt = 'Day 7'

my_comparisons <- list( c('V', 'K'), c('K', 'KP'), c("V", "KP") )
D <- D_long %>% filter(Day == day_txt, condition %in% c('V', 'K', 'KP'))
a7 <- ggplot(data = D , aes(condition, AIMS, fill = condition, line = 'black')) + 
  stat_summary(fun = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar",  width=.3, size = 1) + 
  #stat_compare_means(method = "t.test",comparisons = my_comparisons, label = "p.signif",label.y = c(33, 35, 37)) +
  scale_fill_manual(values=c('#7B7B7B','#0C65FE','#E10004')) +  coord_cartesian(ylim=c(0, uplim)) + 
  labs(title=" ",x="", y = ylab_txt,size = txt_size) + theme(text = element_text(size = 16),legend.position="none",plot.title = element_text(hjust = 0.5),axis.text.x=element_blank(),axis.title.y=element_text(size = 8))

my_comparisons <- list( c('V', 'K'), c('K', 'KL'), c("V", "KL") )
b7 <- ggplot(data = D_long %>% filter(Day == day_txt, condition %in% c('V', 'K', 'KL')) , aes(condition, AIMS, fill = condition, line = 'black')) + 
  stat_summary(fun = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar",  width=.3, size = 1) + 
  #stat_compare_means(method = "t.test",comparisons = my_comparisons, label = "p.signif",label.y = c(33, 35, 37)) +
  scale_fill_manual(values=c('#7B7B7B','#0C65FE','#FA6B0C')) + coord_cartesian(ylim=c(0, uplim)) + 
  labs(title=day_txt, x="", y = "",size = txt_size) + theme(text = element_text(size = 16),legend.position="none",plot.title = element_text(hjust = 0.5),axis.text.x=element_blank())

my_comparisons <- list( c('V', 'P'), c('P', 'L'), c("V", "L") )
c7 <- ggplot(data = D_long %>% filter(Day == day_txt, condition %in% c('V', 'P', 'L')) , aes(condition, AIMS, fill = condition, line = 'black')) + 
  stat_summary(fun = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar",  width=.3, size = 1) + 
  #stat_compare_means(method = "t.test",comparisons = my_comparisons, label = "p.signif",label.y = c(33, 35, 37)) +
  scale_fill_manual(values=c('#7B7B7B','#688035','#69006B'))+ coord_cartesian(ylim=c(0, uplim)) + 
  labs(title=" ",x="", y = "",size = txt_size) + theme(text = element_text(size = 16),legend.position="none",plot.title = element_text(hjust = 0.5),axis.text.x=element_blank())

day_txt = 'Day 14'

my_comparisons <- list( c('V', 'K'), c('K', 'KP'), c("V", "KP") )
a14 <- ggplot(data = D_long %>% filter(Day == day_txt, condition %in% c('V', 'K', 'KP')) , aes(condition, AIMS, fill = condition, line = 'black')) + 
  stat_summary(fun = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar",  width=.3, size = 1) + 
  #stat_compare_means(method = "t.test",comparisons = my_comparisons, label = "p.signif",label.y = c(33, 35, 37)) +
  #geom_signif(comparisons = my_comparisons, map_signif_level = TRUE,test = 't.test') +
  scale_fill_manual(values=c('#7B7B7B','#0C65FE','#E10004')) + coord_cartesian(ylim=c(0, uplim)) + 
  labs(title=" ",x="", y = "",size = txt_size) + theme(text = element_text(size = 16),legend.position="none",plot.title = element_text(hjust = 0.5))


my_comparisons <- list( c('V', 'K'), c('K', 'KL'), c("V", "KL") )
b14 <- ggplot(data = D_long %>% filter(Day == day_txt, condition %in% c('V', 'K', 'KL')) , aes(condition, AIMS, fill = condition, line = 'black')) + 
  stat_summary(fun = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar",  width=.3, size = 1) + 
  #stat_compare_means(method = "t.test",comparisons = my_comparisons, label = "p.signif",label.y = c(33, 35, 37),tip_length = 0.20,margin_top = 0.00,size =4) +
  #geom_signif(comparisons = my_comparisons, map_signif_level = TRUE,test = 't.test') +
  scale_fill_manual(values=c('#7B7B7B','#0C65FE','#FA6B0C')) + coord_cartesian(ylim=c(0, uplim)) + 
  labs(title=day_txt, x="", y = "",size = txt_size) + theme(text = element_text(size = 16),legend.position="none",plot.title = element_text(hjust = 0.5),axis.text.x=element_blank())


my_comparisons <- list( c('V', 'P'), c('P', 'L'), c("V", "L") )
c14 <- ggplot(data = D_long %>% filter(Day == day_txt, condition %in% c('V', 'P', 'L')) , aes(condition, AIMS, fill = condition, line = 'black'),axis.text.x=element_blank()) + 
  stat_summary(fun = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar",  width=.3, size = 1) + 
  #stat_compare_means(method = "t.test",comparisons = my_comparisons, label = "p.signif",label.y = c(33, 35, 37)) +
  #geom_signif(comparisons = my_comparisons),map_signif_level = TRUE)+
  scale_fill_manual(values=c('#7B7B7B','#688035','#69006B')) + coord_cartesian(ylim=c(0, uplim))+ 
  labs(title=" ",x="", y = "",size = txt_size) + theme(text = element_text(size = 16),legend.position="none",plot.title = element_text(hjust = 0.5),axis.text.x=element_blank())


figure <- ggarrange(a0, b0, c0,a7, b7, c7, a14, b14,c14,
                    labels = c("A","B","C", "D", "E", "F","G", "H", "I" ),
                    ncol = 3, nrow = 3)
figure
ggsave('C:\\Temp\\test.pdf', width=6, height=5.27, device=cairo_pdf)
#install.packages('svglite')
library('svglite')
ggsave('C:\\Temp\\test.svg', width=6, height=5.27)

# P values

print("All Comparative Statistics")

D <-D_long %>% filter(Day == 'Day 0' & condition %in% c('V', 'K', 'KP')) 
print(TukeyHSD(aov(AIMS ~ condition, data = D), conf.level=.95))

D <-D_long %>% filter(Day == 'Day 0' & condition %in% c('V', 'K', 'KL')) 
print(TukeyHSD(aov(AIMS ~ condition, data = D), conf.level=.95))

D <- D_long %>% filter(Day == 'Day 0' & condition %in% c('V', 'P', 'L'))
print(TukeyHSD(aov(AIMS ~ condition, data = D), conf.level=.95))



D <-D_long %>% filter(Day == 'Day 7' & condition %in% c('V', 'K', 'KP')) 
print(TukeyHSD(aov(AIMS ~ condition, data = D), conf.level=.95))

D <-D_long %>% filter(Day == 'Day 7' & condition %in% c('V', 'K', 'KL')) 
print(TukeyHSD(aov(AIMS ~ condition, data = D), conf.level=.95))

D <- D_long %>% filter(Day == 'Day 7' & condition %in% c('V', 'P', 'L'))
print(TukeyHSD(aov(AIMS ~ condition, data = D), conf.level=.95))



D <-D_long %>% filter(Day == 'Day 14' & condition %in% c('V', 'K', 'KP')) 
print(TukeyHSD(aov(AIMS ~ condition, data = D), conf.level=.95))

D <-D_long %>% filter(Day == 'Day 14' & condition %in% c('V', 'K', 'KL')) 
print(TukeyHSD(aov(AIMS ~ condition, data = D), conf.level=.95))

D <- D_long %>% filter(Day == 'Day 14' & condition %in% c('V', 'P', 'L'))
print(TukeyHSD(aov(AIMS ~ condition, data = D), conf.level=.95))




# for (day_text in c('Day 0', 'Day 7', 'Day 14')){
#   # NOTE _ DID THIS WITH A FOR LOOP AND IT"S ALL FUCKED UP!!!! Only seems to get Day 14 despite the loop
# not using a loop as above wored fine! Is this another nasty ass bug in R?
#   print(day_text)
#   print("V K KP")
#   D <-D_long %>% filter(Day == day_txt & condition %in% c('V', 'K', 'KP')) 
#   mod = aov(AIMS ~ condition, data = D)
#   print(summary(mod))
#   print(TukeyHSD(mod, conf.level=.95))
# 
#   print(day_text)
#   print("V K KL")
#   D <-D_long %>% filter(Day == day_txt & condition %in% c('V', 'K', 'KL')) 
#   mod = aov(AIMS ~ condition, data = D)
#   print(summary(mod))
#   print(TukeyHSD(mod, conf.level=.95))
#   
#   
#   print(day_text)
#   print("V P L")
#   D <-D_long %>% filter(Day == day_txt & condition %in% c('V', 'P', 'L')) 
#   mod = aov(AIMS ~ condition, data = D)
#   print(summary(mod))
#   print(TukeyHSD(mod, conf.level=.95))
#   
#   
# }

