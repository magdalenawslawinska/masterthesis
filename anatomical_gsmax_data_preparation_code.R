#Magdalena Slawinska
#anatomical gsmax data preparation

#libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(ggpubr)


#####

#formula optimization: code to generate final tables with calculated values

#####

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/anatomical_paper")

data <- read.csv("anatomical_data.csv")

data <- data %>%
  filter(Genotype == "WT") %>%
  select(-c(Stomatal_density_FM, Stomatal_density_MJ, Stomatal_density_ON, Stomatal_density_ind)) %>%
  rename(Sample_id = ï..Sample_id) %>%
  rename(Stoma_width_poles = GC_width_poles)

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/anatomical_paper/figure/output")


#calculating ratios allowing to estimate pore length and width from GC length (to then use with light microscopy pictures)
#pore length/GC length
#pore width/GC length
#real pore area/pore area rectangle (rectangle: pore length * pore width)

data_open_GCopt_wt <- data_wt %>%
  filter(Treatment == "FUS")  %>%
  mutate(GC_mean_width = (GC_left_width_middle + GC_right_width_middle)/2,
         Pore_area_rectangle_approximation = Pore_length*Pore_width) %>%
  select(Sample_id,  Genotype, Pore_area_polygon_selection, Pore_area_rectangle_approximation,
         Pore_length, GC_length, GC_mean_width, Pore_width) %>%
  mutate(GC_pore_length_ratio = Pore_length/GC_length,
         Polygon_rectangle_ratio = Pore_area_polygon_selection/Pore_area_rectangle_approximation,
         GC_pore_width_ratio = Pore_width/GC_length)

data_open_GCopt_wt_sum <- data_open_GCopt_wt %>%
  summarize(Mean_GC_pore_length_ratio = mean(GC_pore_length_ratio),
            sd_GC_pore_length_ratio= sd(GC_pore_length_ratio),
            Mean_polygon_rectangle_ratio = mean(Polygon_rectangle_ratio),
            sd_polygon_rectangle_ratio= sd(Polygon_rectangle_ratio),
            Mean_GC_pore_width_ratio = mean(GC_pore_width_ratio),
            sd_GC_pore_width_ratio= sd(GC_pore_width_ratio))

#save files with ratios - coefficients
write.csv(data_open_GCopt_wt_sum, file = "data_open_coefficients.csv", quote = FALSE, 
          row.names = FALSE)


#constants for 25 *C
d <- 24.9*10^(-6)
v <- 24.4*10^(-3)
a1 <- 0.9 #coefficient based on the Pore_area_polygon_selection/Pore_area_rectangle_approximation

#final big table with all the morphological data and gsmax calculations
data_wt <-  data %>%
  mutate(GC_mean_width = (GC_left_width_middle + GC_right_width_middle)/2,
         GC_width_poles = Stoma_width_poles/2,
         Pore_area_ellipse = pi*(Pore_length/2)*(Pore_length/4),
         Pore_area_rectangle = Pore_length*Pore_width,
         Pore_area_approximation = a1*Pore_length*Pore_width,
         gsmax_exact = (d*Stomatal_density*Pore_area_polygon_selection)/(v*(Pore_depth+(pi/2)*(sqrt(Pore_area_polygon_selection/pi)))),
         gsmax_ellipse =  (d*Stomatal_density*Pore_area_ellipse)/(v*(Pore_depth+(pi/2)*(sqrt(Pore_area_ellipse/pi)))),
         gsmax_rectangle = (d*Stomatal_density*Pore_area_rectangle)/(v*(Pore_depth+(pi/2)*(sqrt(Pore_area_rectangle/pi)))),
         gsmax_approximation = (d*Stomatal_density*Pore_area_approximation)/(v*(Pore_depth+(pi/2)*(sqrt(Pore_area_approximation/pi)))),
         gsmax_approx_GC_width = (d*Stomatal_density*Pore_area_polygon_selection)/(v*(GC_mean_width+(pi/2)*(sqrt(Pore_area_polygon_selection/pi)))))



data_open_wt <- data_wt %>%
  filter(Treatment == "FUS")

write.csv(data_open_wt, file = "data_open_wt_all.csv", quote = FALSE, 
          row.names = FALSE)

#table with all the calculations and ratios (for the paper)

data_wt_all <-  data %>%
  mutate(GC_mean_width = (GC_left_width_middle + GC_right_width_middle)/2,
         GC_width_poles = Stoma_width_poles/2,
         Pore_area_ellipse = pi*(Pore_length/2)*(Pore_length/4),
         Pore_area_rectangle = Pore_length*Pore_width,
         Pore_area_approximation = a1*Pore_length*Pore_width,
         GC_pore_length_ratio = Pore_length/GC_length,
         Polygon_rectangle_ratio = Pore_area_polygon_selection/Pore_area_rectangle,
         GC_pore_width_ratio = Pore_width/GC_length,
         gsmax_exact = (d*Stomatal_density*Pore_area_polygon_selection)/(v*(Pore_depth+(pi/2)*(sqrt(Pore_area_polygon_selection/pi)))),
         gsmax_ellipse =  (d*Stomatal_density*Pore_area_ellipse)/(v*(Pore_depth+(pi/2)*(sqrt(Pore_area_ellipse/pi)))),
         gsmax_rectangle = (d*Stomatal_density*Pore_area_rectangle)/(v*(Pore_depth+(pi/2)*(sqrt(Pore_area_rectangle/pi)))),
         gsmax_approximation = (d*Stomatal_density*Pore_area_approximation)/(v*(Pore_depth+(pi/2)*(sqrt(Pore_area_approximation/pi)))),
         gsmax_approx_GC_width = (d*Stomatal_density*Pore_area_polygon_selection)/(v*(GC_mean_width+(pi/2)*(sqrt(Pore_area_polygon_selection/pi)))))

data_wt_summary_all <- data_wt_all %>%
  group_by(Treatment) %>%
  summarize(Mean_GC_length = mean(GC_length),
            SD_GC_length = sd(GC_length),
            Mean_pore_area = mean(Pore_area_polygon_selection),
            SD_pore_area = sd(Pore_area_polygon_selection),
            Mean_pore_length = mean(Pore_length),
            SD_pore_length = sd(Pore_length),
            Mean_pore_width = mean(Pore_width),
            SD_pore_width = sd(Pore_width),
            Mean_pore_depth = mean(Pore_depth),
            SD_pore_depth = sd(Pore_depth),
            Mean_GC_width_center = mean(GC_mean_width),
            SD_GC_width_center = sd(GC_mean_width),
            Mean_GC_width_poles = mean(GC_width_poles),
            SD_GC_width_poles = sd(GC_width_poles),
            Mean_GC_pore_length_ratio = mean(GC_pore_length_ratio),
            sd_GC_pore_length_ratio= sd(GC_pore_length_ratio),
            Mean_polygon_rectangle_ratio = mean(Polygon_rectangle_ratio),
            sd_polygon_rectangle_ratio= sd(Polygon_rectangle_ratio),
            Mean_GC_pore_width_ratio = mean(GC_pore_width_ratio),
            sd_GC_pore_width_ratio= sd(GC_pore_width_ratio), 
            Mean_pore_area_ellipse = mean(Pore_area_ellipse),
            SD_pore_area_ellipse = sd(Pore_area_ellipse),
            Mean_pore_area_rectangle = mean(Pore_area_rectangle),
            SD_pore_area_rectangle = sd(Pore_area_rectangle),
            Mean_pore_area_approximation = mean(Pore_area_approximation),
            SD_pore_area_approximation = sd(Pore_area_approximation),
            Mean_gsmax_exact = mean(gsmax_exact),
            SD_gsmax_exact = sd(gsmax_exact),
            Mean_gsmax_ellipse = mean(gsmax_ellipse),
            SD_gsmax_ellipse = sd(gsmax_ellipse),
            Mean_gsmax_rectangle = mean(gsmax_rectangle),
            SD_gsmax_rectangle = sd(gsmax_rectangle),
            Mean_gsmax_approximation = mean(gsmax_approximation),
            SD_gsmax_approximation = sd(gsmax_approximation),
            Mean_gsmax_approx_GC_width = mean(gsmax_approx_GC_width),
            SD_gsmax_approx_GC_width = sd(gsmax_approx_GC_width))

write.csv(data_wt_all, file = "data_wt_all_final_table_paper.csv", quote = FALSE, 
          row.names = FALSE)

write.csv(data_wt_summary_all, file = "data_wt_summary_paper.csv", quote = FALSE, 
          row.names = FALSE)