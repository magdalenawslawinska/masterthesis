#Magdalena Slawinska
#anatomical gsmax plots code

#Code for plotting the anatomical data for anatomical gsmax formula optimization

library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(ggpubr)
library(dichromat)
library(RColorBrewer)
library(sjPlot)


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


####
#plots ellipse approximation, rectangle, a*rectangle
####

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/anatomical_paper/figure/output")
data_open_wt <- read.csv("data_open_wt_all.csv")
setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/anatomical_formula_optimization")

#ellipse

ellipse_reg <- lm(Pore_area_polygon_selection ~ Pore_area_ellipse, data = data_open_wt)
ellipse_reg_eqn <- lm_eqn(ellipse_reg) #to get coefficients for the equation and R^2; insert them to ggplot

area_vs_ellipse_plot <- ggplot(data_open_wt,aes(Pore_area_polygon_selection, Pore_area_ellipse)) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              colour="black", size = 2)+
  geom_abline(linetype = "dashed", colour = "#D95450", size = 2) +
  geom_point(colour = "black", size = 5, shape = 19, alpha = 1) +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size = 22, face = "bold", colour = "black"),
        axis.title = element_text(size = 22, face = "bold"),
        axis.line = element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scale_x_continuous(expand = c(0,0), breaks = seq(20, 100, 10), limits = c(20, 100)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(20, 100, 10), limits = c(20, 100)) +
  labs(x = expression(bold(PA[exact]~(µm^2))), y = expression(bold(PA[ellipse]~(µm^2))))+
  annotate(geom = 'text', x = 48,  y = 30, hjust = 0,
           label = "italic(R)^2 == 0.23", 
           parse = TRUE, size = 8)

area_vs_ellipse_plot

pdf("area_vs_ellipse_plot.pdf")
print(area_vs_ellipse_plot, device=cairo_pdf)
dev.off()



#rectangle
rectangle_reg <- lm(Pore_area_polygon_selection ~ Pore_area_rectangle, data = data_open_wt)
rectangle_reg_eqn <- lm_eqn(rectangle_reg)

area_vs_rectangle_plot <- ggplot(data_open_wt, aes(Pore_area_polygon_selection, Pore_area_rectangle)) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              colour="black", size = 2)+
  geom_abline(linetype = "dashed", colour = "#D95450", size = 2) +
  geom_point(colour = "black", size = 5, shape = 19, alpha = 1) +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size = 22, face = "bold", colour = "black"),
        axis.title = element_text(size = 22, face = "bold"),
        axis.line = element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scale_x_continuous(expand = c(0,0), breaks = seq(20, 70, 10), limits = c(20, 70)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(20, 70, 10), limits = c(20, 70)) +
  labs(x = expression(bold(PA[exact]~(µm^2))), y = expression(bold(PA[rectangle]~(µm^2))))+
  annotate(geom = 'text', x = 22,  y = 61, hjust = 0,
           label = "italic(R)^2 == 0.87", 
           parse = TRUE, size = 8) 

area_vs_rectangle_plot

pdf("area_vs_rectangle_plot.pdf")
print(area_vs_rectangle_plot, device=cairo_pdf)
dev.off()

#a*rectangle

approx_reg <- lm(Pore_area_polygon_selection ~ Pore_area_approximation, data = data_open_wt)
approx_reg_eqn <- lm_eqn(approx_reg)

area_vs_approx_plot <- ggplot(data_open_wt,aes(Pore_area_polygon_selection, Pore_area_approximation)) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              colour="black", size = 2)+
  geom_abline(linetype = "dashed", colour = "#D95450", size = 2) +
  geom_point(colour = "black", size = 5, shape = 19, alpha = 1) +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size = 22, face = "bold", colour = "black"),
        axis.title = element_text(size = 22, face = "bold"),
        axis.line = element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scale_x_continuous(expand = c(0,0), breaks = seq(20, 70, 10), limits = c(20, 70)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(20, 70, 10), limits = c(20, 70)) +
  labs(x = expression(bold(PA[exact]~(µm^2))), y = expression(bold(PA[0.9%*%rectangle]~(µm^2))))+
  annotate(geom = 'text', x = 22,  y = 61, hjust = 0,
           label = "italic(R)^2 == 0.87", 
           parse = TRUE, size = 8)

area_vs_approx_plot

pdf("area_vs_approx_plot.pdf")
print(area_vs_approx_plot, device=cairo_pdf)
dev.off()

####
#gsmax plots
####

#ellipse gsmax
ellipse_gsmax_reg <- lm(gsmax_ellipse ~ gsmax_exact, data = data_open_wt)
ellipse_gsmax_reg_eqn <- lm_eqn(ellipse_gsmax_reg)

gsmax_ellipse_plot <- ggplot(data_open_wt,aes(gsmax_exact, gsmax_ellipse)) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              colour="black", size = 2)+
  geom_abline(linetype = "dashed", colour = "#D95450", size = 2) +
  geom_point(colour = "black", size = 5, shape = 19, alpha = 1) +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size = 22, face = "bold", colour = "black"),
        axis.title = element_text(size = 22, face = "bold"),
        axis.line = element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0.2, 0.8, 0.1), limits = c(0.2, 0.80)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.2, 0.8, 0.1), limits = c(0.2, 0.80)) +
  labs(x = expression(bold(italic(g)[smax[exact]]~(mol~m^-2*s^-1))),
       y = expression(bold(italic(g)[smax[ellipse]]~(mol~m^-2*s^-1)))) +
  annotate(geom = 'text', x = 0.35,  y = 0.27, hjust = 0,
           label = "italic(R)^2 == 0.60", 
           parse = TRUE, size = 8) 

gsmax_ellipse_plot

pdf("gsmax_ellipse_plot.pdf")
print(gsmax_ellipse_plot, device=cairo_pdf, height = 10, width = 10)
dev.off()


#rectangle gsmax
rectangle_gsmax_reg <- lm(gsmax_rectangle ~ gsmax_exact, data = data_open_wt)
rectangle_gsmax_reg_eqn <- lm_eqn(rectangle_gsmax_reg)

gsmax_rectangle_plot <- ggplot(data_open_wt,aes(gsmax_exact, gsmax_rectangle)) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              colour="black", size = 2)+
  geom_abline(linetype = "dashed", colour = "#D95450", size = 2) +
  geom_point(colour = "black", size = 5, shape = 19, alpha = 1) +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size = 22, face = "bold", colour = "black"),
        axis.title = element_text(size = 22, face = "bold"),
        axis.line = element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0.2, 0.6, 0.1), limits = c(0.2, 0.60)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.2, 0.6, 0.1), limits = c(0.2, 0.60)) +
  labs(x = expression(bold(italic(g)[smax[exact]]~(mol~m^-2*s^-1))),
       y = expression(bold(italic(g)[smax[rectangle]]~(mol~m^-2*s^-1)))) +
  annotate(geom = 'text', x = 0.21,  y = 0.55, hjust = 0,
           label = "italic(R)^2 == 0.93", 
           parse = TRUE, size = 8) 

gsmax_rectangle_plot

pdf("gsmax_rectangle_plot.pdf")
print(gsmax_rectangle_plot, device=cairo_pdf, height = 10, width = 10)
dev.off()


#a*rectangle gsmax
approx_gsmax_reg <- lm(gsmax_approximation ~ gsmax_exact, data = data_open_wt)
approx_gsmax_reg_eqn <- lm_eqn(approx_gsmax_reg)

gsmax_approx_plot <- ggplot(data_open_wt,aes(gsmax_exact, gsmax_approximation)) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              colour="black", size = 2)+
  geom_abline(linetype = "dashed", colour = "#D95450", size = 2) +
  geom_point(colour = "black", size = 5, shape = 19, alpha = 1) +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size = 22, face = "bold", colour = "black"),
        axis.title = element_text(size = 22, face = "bold"),
        axis.line = element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0.2, 0.6, 0.1), limits = c(0.2, 0.60)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.2, 0.6, 0.1), limits = c(0.2, 0.60)) +
  labs(x = expression(bold(italic(g)[smax[exact]]~(mol~m^-2*s^-1))),
       y = expression(bold(italic(g)[smax[0.9*rectangle]]~(mol~m^-2*s^-1))))+
  annotate(geom = 'text', x = 0.22,  y = 0.55, hjust = 0,
           label = "italic(R)^2 == 0.93", 
           parse = TRUE, size = 8)

gsmax_approx_plot

pdf("gsmax_approx_plot.pdf")
print(gsmax_approx_plot, device=cairo_pdf, height = 10, width = 10)
dev.off()
