#Magdalena Slawinska
#GC length vs pore length for WT, bdpox and BdPOX-sg


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

data <- read.csv("BdPOX_anatomical_data.csv")

data <- data %>% 
  mutate(Genotype = ifelse(Genotype == "NaN1508 #9", "bdpox", Genotype),
         Genotype = ifelse(Genotype == "NaN1508 #3", "BdPOX-sg", Genotype))

data$Genotype <- factor(data$Genotype, levels = c("WT", "bdpox", "BdPOX-sg"))

reg <- lm(Pore_length ~ GC_length, data = data)
reg_eqn <- lm_eqn(reg) #to get coefficients for the equation and R^2; insert them to ggplot

GC_vs_pore_length_plot <- ggplot(data,aes(GC_length, Pore_length)) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              colour="black", size = 2) +
  geom_point(aes(color = Genotype), size = 5, alpha = 0.75, shape = 19) +
  scale_fill_manual(values = c("#8BCF80", "#D95450", "#70B2A3")) +
  scale_color_manual(values = c("#8BCF80", "#D95450", "#70B2A3")) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 22, face = "bold", colour = "black"),
        axis.title = element_text(size = 22, face = "bold"),
        axis.line = element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scale_x_continuous(expand = c(0,0), breaks = seq(24, 34, 2), limits = c(24, 34)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(6, 18, 2), limits = c(6,18)) +
  labs(x = expression(bold(GC~length~(µm))), y = expression(bold(Pore~length~(µm))))+
  annotate(geom = 'text', x = 25,  y = 16, hjust = 0,
           label = "italic(R)^2 == 0.73", 
           parse = TRUE, size = 8)

GC_vs_pore_length_plot

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/anatomical_gsmax_bdpox")

pdf("GC_vs_pore_length_plot.pdf", height = 10, width= 10)
print(GC_vs_pore_length_plot, device=cairo_pdf)
dev.off()