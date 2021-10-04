#Magdalena Slawinska
#BdPOX gsmax comparison with WT
#comparison of WT: gsmax physiological and anatomical



#libraries to load
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(ggpubr)
library(RColorBrewer)
library(dichromat)
library(nortest)
library(moments)
library(FSA)
library(rstatix)


##
#Functions
##

#pull the data out of a linear regression, and return important values (R-squares, slope, intercept and P value) at the top of a ggplot graph with the regression line
#modified after https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    geom_abline(linetype = "dashed")+
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
                       "Intercept =",signif(fit$coef[[1]],2 ),
                       " Slope =",signif(fit$coef[[2]], 2),
                       " P =",signif(summary(fit)$coef[2,4], 2)))
}

#getting equation of the linear model
lm_eqn <- function(m){
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~R2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        R2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

##
#reading in confocal data and calculations
setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/BdPOX_anatomy")

bdpox_confocal <- read.csv("bdpox_project_anatomical_gsmax.csv")

#constants for 25 *C
d <- 24.9*10^(-6)
v <- 24.4*10^(-3)

#coefficients for calculating gsmax from light microscopy data
a <- 0.9 #pore area
b <- 0.44 #pore length
c <- 0.13 #pore width

#anatomical gsmax calculation from confocal data using the exact pore area and Dow formula

bdpox_confocal <- bdpox_confocal %>%
  rename(Sample_id = ï..Sample_id) %>%
  mutate(GC_mean_width = (GC_left_width_middle + GC_right_width_middle)/2)

bdpox_confocal <- bdpox_confocal %>%
  group_by(Genotype, Individual) %>%
  summarize(GC_length_ind = mean(GC_length),
            GC_width_ind = mean(GC_mean_width),
            Pore_length_ind = mean(Pore_length),
            Pore_width_ind = mean(Pore_width), 
            Pore_depth_ind = mean(Pore_depth),
            Pore_area_ind = mean(Pore_area_polygon_selection))


#read in physiological and DIC data
bdpox_physi_DIC <- read.csv("bdpox_project_physi_DIC.csv")

bdpox_physi_DIC <- bdpox_physi_DIC %>%
  rename(Individual = ï..Individual)

#combine confocal, physiological and DIC data
bdpox_data <- merge(bdpox_confocal, bdpox_physi_DIC)

#calculate mean values for the individuals
bdpox_data <- bdpox_data %>%
  mutate(Pore_length_DIC = b*GC_length_DIC,
         Pore_width_DIC = c*GC_length_DIC,
         Pore_area_DIC = a*Pore_length_DIC*Pore_width_DIC,
         gsmax_confocal =(d*SD*Pore_area_ind)/(v*(Pore_depth_ind+(pi/2)*(sqrt(Pore_area_ind/pi)))),
         gsmax_DIC_center = (d*SD*Pore_area_DIC)/(v*((GC_width_center_DIC)+(pi/2)*(sqrt(Pore_area_DIC/pi)))),
         gsmax_DIC_poles = (d*SD*Pore_area_DIC)/(v*((GC_width_poles_DIC)+(pi/2)*(sqrt(Pore_area_DIC/pi)))))



#save output files
setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/BdPOX_anatomy/output")

#combined data with calculated values
write.csv(bdpox_data, file = "bdpox_project_all_data.csv", quote = FALSE, 
          row.names = FALSE)



#transforming data to plot how physiological and anatomical gsmax values fit 

#new columns with gsmax types
gsmax_type <- c(rep("Physiological", nrow(bdpox_data)), 
                rep("Anatomical (LM)", nrow(bdpox_data)), 
                rep("Anatomical (CM)", nrow(bdpox_data)))
  
phys <- bdpox_data %>%
  select(Individual, Genotype, Physiological_gsmax) %>%
  rename(gsmax = Physiological_gsmax)

LM <- bdpox_data %>%
  select(Individual, Genotype, gsmax_DIC_poles) %>%
  rename(gsmax = gsmax_DIC_poles)

CM <- bdpox_data %>%
  select(Individual, Genotype, gsmax_confocal) %>%
  rename(gsmax = gsmax_confocal)

#binding the data together with types
gsmax_comp <- rbind(phys, LM, CM) %>%
  cbind(., gsmax_type)

#saving the file for comparison
write.csv(gsmax_comp, file = "bdpox_project_gsmax_comparison.csv", quote = FALSE, 
          row.names = FALSE)

####
#plotting the comparison
###

#if needed, load in gsmax_comp data

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/BdPOX_anatomy/output")
gsmax_comp <- read.csv("bdpox_project_gsmax_comparison.csv")

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/BdPOX_anatomy/plots")

gsmax_comp$gsmax_type <- factor(gsmax_comp$gsmax_type, 
                                levels = c("Physiological", "Anatomical (LM)", "Anatomical (CM)")) #gsmax_type as factor

gsmax_comp <- gsmax_comp %>% 
  mutate(Genotype = ifelse(Genotype == "NaN1508 #9", "bdpox", Genotype))

gsmax_comp$Genotype <- factor(gsmax_comp$Genotype,
                              levels = c("WT", "bdpox"))


##statistical test - one-way ANOVA followed by Tukey's test

#filtering out wt_2
gsmax_comp <- gsmax_comp %>%
  filter(Individual != "wt_2")

gsmax_comp_wt <- gsmax_comp %>%
  filter(Genotype == "WT")

#shapiro test for normality
s_wt_p <- gsmax_comp_wt %>%
  filter(gsmax_type == "Physiological") 
shapiro.test(s_wt_p$gsmax)


s_wt_LM <- gsmax_comp_wt %>%
  filter(gsmax_type == "Anatomical (LM)") 
shapiro.test(s_wt_LM$gsmax)

s_wt_CM <- gsmax_comp_wt %>%
  filter(gsmax_type == "Anatomical (CM)") 
shapiro.test(s_wt_CM$gsmax)


gsmax_comp_9 <- gsmax_comp %>%
  filter(Genotype == "bdpox")

s_9_p <- gsmax_comp_9 %>%
  filter(gsmax_type == "Physiological") 
shapiro.test(s_9_p$gsmax)

s_9_LM <- gsmax_comp_9 %>%
  filter(gsmax_type == "Anatomical (LM)") 
shapiro.test(s_9_LM$gsmax)

s_9_CM <- gsmax_comp_9 %>%
  filter(gsmax_type == "Anatomical (CM)") 
shapiro.test(s_9_CM$gsmax)

#t_test
stat_test_LM <- gsmax_comp %>%
  filter(gsmax_type == "Anatomical (LM)") %>%
  t_test(gsmax ~ Genotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat_test_CM <- gsmax_comp %>%
  filter(gsmax_type == "Anatomical (CM)") %>%
  t_test(gsmax ~ Genotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat_test_phys <- gsmax_comp %>%
  filter(gsmax_type == "Physiological") %>%
  t_test(gsmax ~ Genotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

#plot without any significance labels grid by gsmax type
gsmax_plot_nolabels <- ggplot(gsmax_comp, aes(Genotype, gsmax)) +
  facet_grid(.~gsmax_type) +
  geom_boxplot(aes(colour = Genotype, fill = Genotype), alpha = 0.3, 
               lwd = 1,
               outlier.shape = NA) + 
  geom_jitter(aes(colour = Genotype),
              alpha = 1, size  = 3, stroke = 2, shape = 1, height = 0, width = 0.1) + 
  scale_fill_manual(values = c("#8BCF80", "#D95450")) +
  scale_color_manual(values = c("#8BCF80", "#D95450")) +
  theme(strip.text.x = element_text(size = 12, face = "bold"), 
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", colour = "black", size = 1),
        axis.line = element_line(size = 2, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.text.x = element_text(size = 12, face = "bold", colour = "black", angle = 45, margin = unit(c(0.5, 0, 0, 0), "cm")),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold", colour = "black"),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2), expand = c(0,0)) +
  labs(x = element_blank(), 
       y = expression(bold(italic(g)[smax]~(mol~m^-2*s^-1)))) 

gsmax_plot_nolabels

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/anatomical_gsmax_bdpox")

pdf("gsmax_bdpox_plot.pdf", height = 10, width= 6)
print(gsmax_plot_nolabels, device=cairo_pdf)
dev.off()

#comparison of gsmax types

stat_test_wt <- gsmax_comp %>%
  filter(Genotype == "WT") %>%
  aov(gsmax ~ gsmax_type, .) %>%
  tukey_hsd() %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat_test_bdpox <- gsmax_comp %>%
  filter(Genotype == "bdpox") %>%
  aov(gsmax ~ gsmax_type, .) %>%
  tukey_hsd() %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()


gsmax_plot_types <- ggplot(gsmax_comp, aes(gsmax_type, gsmax)) +
  facet_grid(.~Genotype) +
  geom_boxplot(aes(colour = Genotype, fill = Genotype), alpha = 0.3, 
               lwd = 1,
               outlier.shape = NA) + 
  geom_jitter(aes(colour = Genotype),
              alpha = 1, size  = 3, stroke = 2, shape = 1, height = 0, width = 0.1) + 
  scale_fill_manual(values = c("#8BCF80", "#D95450")) +
  scale_color_manual(values = c("#8BCF80", "#D95450")) +
  theme(strip.text.x = element_text(size = 12, face = "bold"), 
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", colour = "black", size = 1),
        axis.line = element_line(size = 2, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.text.x = element_text(size = 12, face = "bold", colour = "black", angle = 45, margin = unit(c(0.75, 0, 0, 0), "cm")),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold", colour = "black"),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2), expand = c(0,0)) +
  labs(x = element_blank(), 
       y = expression(bold(italic(g)[smax]~(mol~m^-2*s^-1)))) 

gsmax_plot_types

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/anatomical_gsmax_bdpox")

pdf("gsmax_plot_types.pdf", height = 10, width= 6)
print(gsmax_plot_types, device=cairo_pdf)
dev.off()

#only wt data

stat_test_wt <- aov(gsmax ~ gsmax_type, gsmax_comp_wt) %>%
  tukey_hsd() %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat_test_wt <- stat_test_wt %>% add_y_position(fun = "max")
  
gsmax_comp_wt_plot <- ggplot(gsmax_comp_wt, aes(gsmax_type, gsmax)) +
  geom_boxplot(aes(colour = gsmax_type, fill = gsmax_type), alpha = 0.2, 
               lwd = 1,
               outlier.shape = NA) + 
  geom_jitter(aes(colour = gsmax_type),
              alpha = 1, size  = 3, stroke = 2, shape = 1, height = 0, width = 0.1) + 
  scale_fill_manual(values = c("#8BCF80", "#D95450", "#70B2A3")) +
  scale_color_manual(values = c("#8BCF80", "#D95450", "#70B2A3")) +
  theme(strip.text.x = element_text(size = 12, face = "bold"), 
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", colour = "black", size = 1),
        axis.line = element_line(size = 2, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.text.x = element_text(size = 22, face = "bold", colour = "black", angle = 45, margin = unit(c(1.75, 0, 0, 0), "cm")),
        axis.text.y = element_text(size = 22, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 22, face = "bold", colour = "black"),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0,0.7), breaks = seq(0, 0.7, 0.1), expand = c(0,0)) +
  labs(x = element_blank(), y = expression(bold(~italic(g)[smax]~(mol~m^-2*s^-1)))) 
  #stat_pvalue_manual(stat_test_wt, label = "p.adj",  y.position = "y.position", hide.ns = FALSE)

gsmax_comp_wt_plot

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/anatomical_formula_optimization")

pdf("gsmax_comp_wt_plot.pdf", height = 10, width= 6)
print(gsmax_comp_wt_plot, device=cairo_pdf)
dev.off()


#summary bdpox

bdpox_data_sum <- bdpox_data %>%
  group_by(Genotype) %>%
  summarize(Mean_GCL_ind = mean(GC_length_ind),
          SD_GCL_ind = sd(GC_length_ind),
          Mean_GCL_DIC = mean(GC_length_DIC),
          SD_GCL_DIC = mean(GC_length_DIC),
          Mean_SD = mean(SD),
          SD_SD = sd(SD),
          Mean_physiological_gsmax = mean(Physiological_gsmax),
          SD_physiological_gsmax = mean(Physiological_gsmax),
          Mean_gsmax_confocal = mean(gsmax_confocal),
          SD_gsmax_confocal = sd(gsmax_confocal),
          Mean_gsmax_DIC = mean(gsmax_DIC_poles),
          SD_gsmax_DIC = mean(gsmax_DIC_poles))

bdpox_data_sum <- bdpox_data_sum %>% 
  mutate(Genotype = ifelse(Genotype == "NaN1508 #9", "bdpox", Genotype))

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/gsmax_bdpox_discussion")

write.csv(bdpox_data_sum, file = "bdpox_project_gsmax_comparison_summary.csv", quote = FALSE, 
          row.names = FALSE)
