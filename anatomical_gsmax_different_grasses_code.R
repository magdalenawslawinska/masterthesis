#Magdalena Slawinska
#code for assembling raw Fiji data and calculating gsmax for different grasses

#libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(ggpubr)
library(RColorBrewer)
library(dichromat)
library(FSA)
library(rstatix)

#constants for 25 *C
d <- 24.9*10^(-6)
v <- 24.4*10^(-3)

#coefficients for calculating gsmax from light microscopy data
a <- 0.9 #pore area
b <- 0.44 #pore length
c <- 0.13 #pore width

#reading in SD data and calculating SD  

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/anatomical_gsmax_different_grasses")

SD <- read.csv("SD_grasses.csv")

SD <- SD %>% 
  rename(Species = ï..Species) %>%
  mutate(Dimension1_mm = Dimension1_um/1000, Dimension2_mm = Dimension2_um/1000,
         Area_mm2 = Dimension1_mm*Dimension2_mm, Stomatal_density = Number_of_stomata/Area_mm2)

#reading in Fiji data (stomata length and width)

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/anatomical_gsmax_different_grasses/R_rawdata") #folder with all the Fiji raw results


datafiles <- list.files(path = "C:/Users/Magdalena/Desktop/MAPS/Master_thesis/anatomical_gsmax_different_grasses/R_rawdata") #listing all fines

data <- data.frame() #empty data frame for raw results

for(i in 1:length(datafiles)){
  new_data <- read.csv(datafiles[i])
  data <- rbind(data, new_data)
}


data$Sample_id <- data$Label

#divind Label into Microscopy date, Species, ID, Magnification columns
data <- data %>%
  select(-c(X, Angle)) %>%
  separate(Label, into = c("Microscopy_date", "Species1", "Species2", "ID",
                               "Magnification", "Pic_no"), sep = "\\_") %>%
  mutate(Species = paste(Species1, Species2, sep = " ")) %>%
  select(-c(Species1, Species2, Pic_no))

#dividing data into GC_width_poles and GC_length

data_GCL <- data %>%
  filter(Length >= 20) %>%
  rename(GC_length = Length)
  

data_GCWP <- data %>%
  filter(Length < 20) %>%
  rename(Stoma_width_poles = Length) %>%
  mutate(GC_width_poles = Stoma_width_poles/2)


#calculate mean for individuals

data_GCL <- data_GCL %>%
  group_by(Microscopy_date, Species, ID) %>%
  summarize(GC_length_ind = mean(GC_length))

data_GCWP <- data_GCWP %>%
  group_by(Microscopy_date, Species, ID) %>%
  summarize(GC_width_poles_ind = mean(GC_width_poles))


#merge GCL and GCWP dfs, add SD

data <- merge(data_GCL, data_GCWP)

#add Species_date column to both dfs to use it for merging


data <- merge(data, SD, by = c("Species", "Microscopy_date"))

#get rid of double ID and dimension columns

data <- data %>%
  select(-c(ID.y, Dimension1_mm, Dimension2_um, Dimension1_um, Dimension2_mm, Area_mm2, Number_of_stomata)) %>%
  rename(ID = ID.x) %>%
  rename(SD = Stomatal_density)

#calculate gsmax
data <- data %>%
  mutate(Pore_length = b*GC_length_ind,
         Pore_width = c*GC_length_ind,
         Pore_area = a*Pore_length*Pore_width,
         gsmax_optimized = (d*SD*Pore_area)/(v*((GC_width_poles_ind)+(pi/2)*(sqrt(Pore_area/pi)))),
         gsmax_ellipse = (d*SD*(pi*(Pore_length/2)*(Pore_length/4)))/(v*((GC_width_poles_ind)+(pi/2)*(sqrt((pi*(Pore_length/2)*(Pore_length/4))/pi)))))

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/anatomical_gsmax_different_grasses/R_output")

write.csv(data, file = "gsmax_grasses_output.csv", quote = FALSE, 
          row.names = FALSE)

##read output 

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/anatomical_gsmax_different_grasses/R_output")

data <- read.csv("gsmax_grasses_output.csv")

#read physiological data
setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/anatomical_gsmax_different_grasses")
data_physi <- read.csv("Physi_gsmax_grasses.csv")

data_physi <- data_physi %>%
  rename(Species = ï..Species)

#merge physi and anatomical data
data <- merge(data, data_physi, by = c("Species", "Microscopy_date"))

#filter out rice measured on 20210623
data <- data %>%
  mutate(temp = paste(Species, Microscopy_date, sep =" ")) %>%
  filter(temp != "Oryza sativa 20210623") %>%
  select(-temp)

#gsmax optimized data
data_opt <- data %>%
  select(-c(gsmax_ellipse, Physiological_gsmax)) %>%
  rename(gsmax = gsmax_optimized)

#add column gsmax_type
type_opt <- rep("Anatomical gsmax", nrow(data_opt))
data_opt$gsmax_type <- type_opt

#gsmax Physiological data
data_p <- data %>%
  select(-c(gsmax_ellipse, gsmax_optimized)) %>%
  rename(gsmax = Physiological_gsmax)

#add column gsmax_type
type_physi <- rep("Physiological gsmax", nrow(data_p))
data_p$gsmax_type <- type_physi

#bind anat and physi gsmax
data <- rbind(data_p, data_opt)

#change "Molina caerulea" to "Molinia caerulea"
data <- data %>%
  mutate(Species = ifelse(Species == "Molina caerulea", "Molinia caerulea", Species)) 

data <- data %>%
  mutate(Mdate = Microscopy_date) %>%
  mutate(Mdate = ifelse(Mdate == 20210623, 1, 2)) %>%
  mutate(Species_date = paste(Species, Mdate, sep = "_")) %>%
  select(-c(Mdate))

Species_short <- rep(c("A. donax", "I. cylindrica", "M. sinensis", "M. caerulea", "O. sativa", "S. secundatum"), 2)

data <- cbind(data, Species_short)

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/anatomical_gsmax_different_grasses/R_output")

write.csv(data, file = "gsmax_grasses_plotting_data.csv", quote = FALSE, 
          row.names = FALSE)


#plot according to gsmax type

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/anatomical_gsmax_different_grasses/R_output")

data <- read.csv("gsmax_grasses_plotting_data.csv")

grasses_plot <- ggplot(data,aes(Species_short, gsmax)) +
  geom_point(aes(color = gsmax_type, shape = gsmax_type), size = 7, alpha = 0.7) +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size = 22, face = "bold", colour = "black"),
        axis.text.x = element_text(angle = 45, margin = unit(c(1.75, 0, 0, 0), "cm")),
        axis.title = element_text(size = 22, face = "bold"),
        axis.line = element_line(colour = "black", size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 22, face = "bold"))+
  labs(x = element_blank(),
       y = expression(bold(~italic(g)[smax]~(mol~m^-2*s^-1)))) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1.5, 0.25), limits = c(0, 1.5)) 

  
grasses_plot  

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/anatomical_gsmax_different_grasses")

pdf("grasses_plot.pdf", height = 10, width= 15)
print(grasses_plot, device=cairo_pdf)
dev.off()


  