################################################################################
#-#   o~~     o~~     Tadpole Familiarity & Disturbance Cues   ~~o     ~~o   #-#
#-#        o~~               Kevin Bairos-Novak, 2020             ~~o        #-#
#-#     o~~     o~~      Run on a Mac with R version 3.5.2   ~~o     ~~o     #-#
################################################################################

### Set-up ###

setwd("~/Documents/Analytics/R/Disturbance Cues/Familiarity & Kinship Donors")
rm(list=ls())
options(contrasts = c("contr.sum", "contr.poly")) # get Type III Sum of Squares

require(readxl)    # for loading excel data
require(tidyverse) # for manipulating data
require(magrittr)  # for manipulating data
require(car)       # for conducting the Levene's F-test of homoscedasticity
require(lme4)      # for quasipoisson modelling
require(lmerTest)  # for quasipoisson modelling
require(grid)      # for plotting figures
require(gridExtra) # for plotting figures
require(ggplotify) # for plotting figures

# define a custom function to reorder factors
refact <- function(factor.var, new.order.vec = seq(length(levels(factor(factor.var))),1) ){
	x <- factor(factor.var, levels = levels(factor(factor.var))[new.order.vec])
	return(x)
}

# load the data and assign factor variables
data <- read_excel("~/Documents/Analytics/R/Disturbance Cues/Familiarity & Kinship Donors/Dryad/data.xlsx")
data <- data %>% mutate(trial = factor(trial), 
			kinship = factor(kinship), 
			familiarity = factor(familiarity),
			cue = factor(cue) %>% fct_relevel("W", "UC", "DC"))

################################################################################
# Main analysis 1: three-way design
################################################################################

# examine only undisturbed or disturbance cues (exclude control water, "W" cue)
# use familiarity and kinship distinctions

### pre-exposure analysis:

# remove water cues, as they cannot be grouped by the factor levels:
(m1.pre <- data %>% filter(cue!="W") %>%
	glm(data=., pre ~ cue * familiarity * kinship, family=quasipoisson)) %>%
	summary # significant familiarity x kinship iteraction, no other significant factors
m1.pre %>% drop1("familiarity:kinship", test="F") 
	# F(1,239) = 5.28, P = 0.023
m1.pre %>% drop1("kinship", test="F")
	# F(1,239) ≤ 0.52, P ≥ 0.47)

## split analysis for examining familiarity x kinship interaction:
# by kinship:
data %>% filter(cue!="W") %>%
	filter(kinship == "kin") %>%
	glm(data=., pre  ~ cue * familiarity, family=quasipoisson) %>% 
	# summary
	drop1("familiarity", test="F") # F(1,118) ≤ 2.00, P ≥ 0.16

data %>% filter(cue!="W") %>%
	filter(kinship == "non-kin") %>%
	glm(data=., pre  ~ cue * familiarity, family=quasipoisson) %>% 
	# summary
	drop1("familiarity", test="F") # F(1,121) ≤ 3.34, P ≥ 0.070

# by familiarity:
data %>% filter(cue!="W") %>%
	filter(familiarity == "familiar") %>%
	glm(data=., pre  ~ cue * kinship, family=quasipoisson) %>% 
	# summary
	drop1("kinship", test="F") # F(1,123) ≤ 1.26, P ≥ 0.26

data %>% filter(cue!="W") %>%
	filter(familiarity == "unfamiliar") %>%
	glm(data=., pre  ~ cue * kinship, family=quasipoisson) %>% 
	#summary
	drop1("kinship", test="F") # F(1,116) = 4.48, P = 0.036
data %>% filter(cue!="W" & familiarity == "unfamiliar") %>% 
	group_by(kinship) %>% 
	summarise(mean(pre), sd(pre)/n())

# overall, no significant pre-exposure differences after correcting for multiple splits


### proportional change in line crosses:

(m1.prop <- data %>% filter(cue!="W") %>%
		lm(data=., prop ~ cue * familiarity * kinship)) %>%
	summary # no factors significant except cue type
m1.prop %>% drop1("cue", test="F") # F(1,239) = 11.89, P = 0.0007
m1.prop %>% drop1("familiarity", test="F") # F(1,239) ≤ 1.99, P ≥ 0.16

### difference in line crosses:
(m1.diff <- data %>% filter(cue!="W") %>%
		lm(data=., diff ~ cue * familiarity * kinship)) %>%
	summary # same outcome
m1.diff %>% drop1("cue", test="F") # F(1,239) = 11.35, P = 0.0009
m1.diff %>% drop1("familiarity", test="F") # F(1,239) ≤ 2.77, P ≥ 0.097


### parametric assumptions:
# normality:
hist(resid(m1.prop))
qqnorm(resid(m1.prop)); qqline(resid(m1.prop))
ks.test(resid(m1.prop), "pnorm", mean = mean(resid(m1.prop)), sd = sd(resid(m1.prop)))
# homoscedasticity:
data %>% car::leveneTest(prop  ~ cue * familiarity * kinship, data=.)

# extra analysis: using quasibinomial to predict proportion of total line 
		# crosses occurring in post vs. pre-exposure period
d <- data %>% filter(cue!="W") %>%
	mutate(p = post/(pre+post))
glm(data=d, p ~ cue * familiarity * kinship, family=quasibinomial) %>%
	summary # cue still only significant parameter, p = 0.0012

### effect sizes:
data %>% filter(cue!="W") %>%
	group_by(familiarity, kinship, cue) %>%
	summarise(mean = mean(pre, na.rm=T), 
		  se = sd(pre, na.rm=T)/sqrt(n()),
		  n = n())
data %>% filter(cue!="W") %>%
	group_by(familiarity, kinship, cue) %>%
	summarise(mean = mean(post, na.rm=T), 
		  se = sd(post, na.rm=T)/sqrt(n()),
		  n = n())
data %>% filter(cue!="W") %>%
	group_by(familiarity, kinship, cue) %>%
	summarise(mean = mean(diff, na.rm=T), 
		  se = sd(diff, na.rm=T)/sqrt(n()),
		  n = n())
data %>% filter(cue=="W") %>%
	summarise(mean = mean(pre, na.rm=T), 
		  se = sd(pre, na.rm=T)/sqrt(n()),
		  n = n())
data %>% filter(cue=="W") %>%
	summarise(mean = mean(post, na.rm=T), 
		  se = sd(post, na.rm=T)/sqrt(n()),
		  n = n())
data %>% filter(cue=="W") %>%
	summarise(mean = mean(diff, na.rm=T), 
		  se = sd(diff, na.rm=T)/sqrt(n()),
		  n = n())

################################################################################
# Main analysis 2: one-way cue type
################################################################################

# examine all three cue types, including control water
# avoid familiarity and kinship distinctions (water cannot be familiar/kin)

### pre-exposure analysis:
(m2.pre <- data %>% glm(data=., pre ~ cue, family=quasipoisson)) %>%
	summary
m2.pre %>% drop1("cue", test="F") # F(2, 370) ≤ 0.07, p ≥ 0.93


### proportional change in line crosses:
(m2.prop <- data %>% #within({cue <- relevel(cue, ref="DC")}) %>%
	aov(data=., prop ~ cue)) %>%
	summary # F(2,270) = 13.0, P ≤ 0.0001

### difference in line crossses:
(m2.diff <- data %>% #within({cue <- relevel(cue, ref="DC")}) %>%
		aov(data=., diff ~ cue)) %>%
	summary # F(2,270) = 8.45, P = 0.0003

### parametric assumptions:
# normality:
hist(resid(m2.prop))
qqnorm(resid(m2.prop)); qqline(resid(m2.prop))
ks.test(resid(m2.prop), "pnorm", mean = mean(resid(m2.prop)), sd = sd(resid(m2.prop)))
# homoscedasticity:
data %>% car::leveneTest(prop  ~ cue * familiarity * kinship, data=.)

### post-hoc Tukey HSD tests:
TukeyHSD(m2.prop)
# W  vs. UC: P = 0.13
# W  vs. DC: P ≤ 0.0001
# UC vs. DC: P = 0.0053
TukeyHSD(m2.diff)
# W  vs. UC: P = 0.84
# W  vs. DC: P = 0.0005
# UC vs. DC: P = 0.0036

### effect sizes by cue type:
(sum.data <- data %>% group_by(cue) %>%
		summarise(mean = mean(prop, na.rm=T), 
			  # se = sd(prop, na.rm=T)/sqrt(n()), 
			  n = n()) )
sum.data$mean[3]/sum.data$mean[1] # ~23 fold greater response
sum.data$mean[3]/sum.data$mean[2] # ~2.5 fold greater response
sum.data$mean[3]*100 # 37% mean reduction due to DCs
data %>% filter(cue != "DC") %>% summarise(mean = mean(prop, na.rm=T))*100 
# 8% mean reduction due to UC/W
-37.30232/-8.355641 # proportional difference in response to DC vs. UC/W


################################################################################
# Figure 2: b&w and colour
################################################################################

# Create a custom ggtheme object for plotting barplots
ggtheme <-
	theme(	plot.margin = unit(c(0.5,0.5,0.5,0),"cm"),
	       axis.line = element_line(),
	       axis.text = element_text(size = 17, color = "black"),
	       axis.text.x=element_text(vjust=-0.1),
	       axis.title.x = element_text(margin = margin(t=20), size=20, angle=0),
	       axis.title.y = element_text(margin = margin(r=20), size=20, vjust=-5),
	       legend.position = c(0.5, 0.95), 
	       legend.direction = "horizontal",
	       legend.title = element_text(size=17),
	       legend.text  = element_text(size=17),
	       legend.justification = "center",
	       legend.key = element_rect(colour = "white", fill = NA),
	       panel.background = element_blank(),
	       panel.border = element_blank(),
	       panel.grid.major = element_blank(), # major panel grid lines
	       panel.grid.minor = element_blank(), # minor grid lines
	       panel.spacing = unit(-0.09, "lines"), # space between facet panels
	       strip.text.y = element_text(size = 15, angle = 0), # y-axis facet variable title
	       strip.text.x = element_text(size = 15), # x-axis facet variable title
	       strip.background = element_blank() # facet panel background
	)

## Figure 1:
p1 <- data %>% filter(cue!="W") %>%
	within({cue <- factor(cue, levels=c("UC","DC"), 
			      labels=c("Undisturbed cue", "Disturbance cue"))}) %>%
	mutate(kinfam = interaction(kinship, familiarity)) %>%
	group_by(kinfam, cue) %>%
	summarise(mean = mean(prop, na.rm=T), se = sd(prop, na.rm=T)/sqrt(n()), n = n()) %>%
	ggplot(aes(x = kinfam, y = mean, fill = cue)) + 
	ggtheme +
	geom_bar(position=position_dodge(), stat = "identity", col = "black") +
	geom_errorbar(position=position_dodge(width = 0.9), aes(ymin = mean - se, ymax = mean + se),
		      width = 0.3) +
	labs(y="Percent change in receiver line crosses\n", x=NULL) +
	scale_fill_manual(name = "Cue type: ", labels = c("Undisturbed cue ", "Disturbance cue"),
			  # values = c("bisque", "dodgerblue1", "light green")) +
			  values = c("grey", "gray27")) + # also this one:
	scale_y_continuous(limits = c(-0.58, 0.15), # -0.55 also good without connecting letters
			   breaks = seq(-0.5,0,by=0.1),
			   labels = scales::percent_format(accuracy = 1)) +
	scale_x_discrete(name = "Audience composition", labels = c("Familiar\nkin", "Familiar\nnon-kin", "Unfamiliar\nkin", "Unfamiliar\nnon-kin")) +
	geom_abline(aes(slope=0, intercept=0)) +
	geom_text(aes(y = ifelse(mean > 0, mean-se-0.04, 
				 ifelse(mean+se > 0, mean+se+0.04, 0.04)), 
				 label=paste0("n=",n,"")), 
		  position=position_dodge(width=0.9), size=5) +
	geom_text(aes(y=mean-se-0.04,
		      label = rep(c("A", "B"),4)), 
		  position=position_dodge(width=0.9), size=5)
p1

# Figure 2:
p2 <- data %>%
	group_by(cue) %>%
	summarise(mean = mean(prop, na.rm=T), se = sd(prop, na.rm=T)/sqrt(n()), n = n()) %>%
	ggplot(aes(x = cue, y = mean, fill = cue)) +
	ggtheme +
	geom_bar(position=position_dodge(width = 0.9), width=1, stat = "identity", col = "black") +
	geom_errorbar(position=position_dodge(width = 0.9), aes(ymin = mean - se, ymax = mean + se),
		      width = 0.3) +
	labs(y=NULL, x="Cue type") +
	scale_fill_manual(guide = FALSE, values = c("white", "grey", "gray27")) + # also this one:
	scale_y_continuous(limits = c(-0.58, 0.15), 
			   breaks = seq(-0.5,0,by=0.1),
			   labels = scales::percent_format(accuracy = 1)) +
	scale_x_discrete(labels = c("Control\nwater", "Undisturbed\ncue", "Disturbance\ncue")) +
	geom_abline(aes(slope=0, intercept=0)) +
	theme(axis.text.y=element_blank()) + 
	geom_text(aes(y = ifelse(mean > 0, mean-se-0.04, 
				 ifelse(mean+se > 0, mean+se+0.04, 0.04)), 
				 label=paste0("n=",n,"")), 
		  position=position_dodge(width=0.9), size=5) +
	geom_text(aes(y=mean-se-0.04,
		      label = c("a", "a","b")), 
		  position=position_dodge(width=0.9), size=5)
# p2

grid.arrange(as.grob(p1), as.grob(p2), layout_matrix=t(as.matrix(c(1,1,2))))

# save b & w
setEPS()
postscript("Fig2.eps", width = 14, height = 6)
grid.arrange(as.grob(p1), as.grob(p2), layout_matrix=t(as.matrix(c(1,1,2))))
dev.off()


# save colour
setEPS()
postscript("Fig2_colour.eps", width = 14, height = 6)
grid.arrange(as.grob(p1 +scale_fill_manual(name = "Cue type: ", labels = c("Undisturbed cue ", "Disturbance cue"),
					   values = c("bisque", "light green"))), 
	     as.grob(p2 +scale_fill_manual(guide=FALSE,
	     			      values = c("skyblue","bisque", "light green"))), 
	     layout_matrix=t(as.matrix(c(1,1,2))))
dev.off()

################################################################################
# Supplementary Figure 1: pre-exposure line crosses across treatments
################################################################################

# pre-exposure line crosses:
p1.pre <- data %>% filter(cue!="W") %>%
	within({cue <- factor(cue, levels=c("UC","DC"), 
			      labels=c("Undisturbed cue", "Disturbance cue"))}) %>%
	mutate(kinfam = interaction(kinship, familiarity)) %>%
	group_by(kinfam, cue) %>%
	summarise(mean = mean(pre, na.rm=T), se = sd(pre, na.rm=T)/sqrt(n()), n = n()) %>%
	ggplot(aes(x = kinfam, y = mean, fill = cue)) + 
	ggtheme +
	geom_bar(position=position_dodge(), stat = "identity", col = "black") +
	geom_errorbar(position=position_dodge(width = 0.9), aes(ymin = mean - se, ymax = mean + se),
		      width = 0.3) +
	labs(y="Pre-exposure line crosses\n", x=NULL) +
	scale_fill_manual(name = "Cue type: ", labels = c("Undisturbed cue ", "Disturbance cue"),
			  # values = c("bisque", "dodgerblue1", "light green")) +
			  values = c("bisque", "light green")) + # also this one:
	scale_y_continuous(limits = c(0, 29),
			   breaks = seq(0,25,by=5),
			   expand = c(0,0)) +
	scale_x_discrete(name = "Audience composition", labels = c("Familiar\nkin", "Familiar\nnon-kin", "Unfamiliar\nkin", "Unfamiliar\nnon-kin")) +
	geom_abline(aes(slope=0, intercept=0)) +
	geom_text(aes(y = mean+se+1,
		      label=paste0("n=",n,"")), 
		  position=position_dodge(width=0.9), size=5)
# p1.pre

p2.pre <- data %>%
	group_by(cue) %>%
	summarise(mean = mean(pre, na.rm=T), se = sd(pre, na.rm=T)/sqrt(n()), n = n()) %>%
	ggplot(aes(x = cue, y = mean, fill = cue)) +
	ggtheme +
	geom_bar(position=position_dodge(width = 0.9), width=1, stat = "identity", col = "black") +
	geom_errorbar(position=position_dodge(width = 0.9), aes(ymin = mean - se, ymax = mean + se),
		      width = 0.3) +
	labs(y=NULL, x="Cue type") +
	scale_fill_manual(guide = FALSE, values = c("skyblue","bisque", "light green")) + # also this one:
	scale_y_continuous(limits = c(0, 29),
			   breaks = seq(0,25,by=5),
			   expand = c(0,0)) +
	scale_x_discrete(labels = c("Control\nwater", "Undisturbed\ncue", "Disturbance\ncue")) +
	geom_abline(aes(slope=0, intercept=0)) +
	theme(axis.text.y=element_blank()) + 
	geom_text(aes(y = mean+se+1,
		      label=paste0("n=",n,"")), 
		  position=position_dodge(width=0.9), size=5)
# p2.pre

setEPS()
postscript("FigS1.eps", width = 14, height = 6)
grid.arrange(as.grob(p1.pre), as.grob(p2.pre), layout_matrix=t(as.matrix(c(1,1,2))))
dev.off()

################################################################################

sessionInfo()
# R version 3.5.2 (2018-12-20)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.2
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
# 	[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8
# 
# attached base packages:
# 	[1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# 	[1] ggplotify_0.0.3 gridExtra_2.3   lmerTest_3.1-0  lme4_1.1-19     Matrix_1.2-15  
# [6] car_3.0-2       carData_3.0-2   magrittr_1.5    forcats_0.3.0   stringr_1.4.0  
# [11] dplyr_0.8.3     purrr_0.3.2     readr_1.3.1     tidyr_1.0.0     tibble_2.1.3   
# [16] ggplot2_3.2.0   tidyverse_1.2.1 readxl_1.2.0   
# 
# loaded via a namespace (and not attached):
# 	[1] tidyselect_0.2.5   splines_3.5.2      haven_2.0.0        lattice_0.20-38   
# [5] colorspace_1.4-1   vctrs_0.2.0        generics_0.0.2     yaml_2.2.0        
# [9] gridGraphics_0.3-0 rlang_0.4.0        nloptr_1.2.1       pillar_1.4.2      
# [13] foreign_0.8-71     glue_1.3.1         withr_2.1.2        modelr_0.1.2      
# [17] rvcheck_0.1.3      lifecycle_0.1.0    munsell_0.5.0      gtable_0.3.0      
# [21] cellranger_1.1.0   rvest_0.3.2        zip_2.0.4          rio_0.5.16        
# [25] curl_3.3           broom_0.5.1        Rcpp_1.0.2         backports_1.1.4   
# [29] scales_1.0.0       jsonlite_1.6       abind_1.4-5        hms_0.4.2         
# [33] stringi_1.4.3      openxlsx_4.1.0     numDeriv_2016.8-1  cli_1.1.0         
# [37] tools_3.5.2        lazyeval_0.2.2     crayon_1.3.4       pkgconfig_2.0.2   
# [41] zeallot_0.1.0      MASS_7.3-51.1      data.table_1.12.0  xml2_1.2.0        
# [45] lubridate_1.7.4    minqa_1.2.4        assertthat_0.2.1   httr_1.4.0        
# [49] rstudioapi_0.8     R6_2.4.0           nlme_3.1-137       compiler_3.5.2 




