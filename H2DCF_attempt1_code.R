#code for H2DCF stainings analysis (developmental zone stainings
#leaves taken from plants kept in the greenhouse and transported for the staining to the lab

library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(ggpubr)
library(RColorBrewer)
library(dichromat)
library(FSA)
library(rstatix)

#####
#
#Mean intensity measured on one slice on which the stoma is fully visible
#1508 #9, 1508 #9 - 4 plants, 3 independent growing conditions
#wt - 4 plants, 3 independent growing conditions
#
#####



setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/ROS/H2DCF_april_may_walking/raw_results") #folder with all the Fiji raw results


datafiles <- list.files(path = "C:/Users/Magdalena/Desktop/MAPS/Master_thesis/ROS/H2DCF_april_may_walking/raw_results") #listing all fines

data_H2DCF <- data.frame() #empty data frame for raw results

for(i in 1:length(datafiles)){
  new_data <- read.csv(datafiles[i])
  data_H2DCF <- rbind(data_H2DCF, new_data)
}

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/ROS/H2DCF_april_may_walking")

data_H2DCF$Sample_id <- data_H2DCF$Label

#getting Sample_id and Project out of Label
data_H2DCF <- data_H2DCF %>%
  select(-c(X)) %>%
  separate(Label, into = c("tempA", "tempB", "Sample_id"), sep = "\\ - ") %>%
  separate(Sample_id, into = c("Project", "temp1", "temp2", "temp3", "temp4",
                               "temp5", "temp6", "temp7", "temp8"), sep = "\\_")

#assembling Individual ids from temp data
data_H2DCF_wt <- data_H2DCF %>%
  filter(temp1 == "wt") %>%
  select(-c(tempA, tempB, temp7, temp8)) %>%
  mutate(., Individual = paste(temp1, temp2, temp3, 
                               temp4,temp5, sep = "_"), 
         Genotype = temp1, 
         Stage = temp6) %>%
  select(-c(temp1, temp2, temp3, temp4, temp5, temp6))

data_H2DCF_1508 <- data_H2DCF %>%
  filter(temp1 == "1508") %>%
  select(-c(temp8, tempA, tempB)) %>%
  mutate(., Individual = paste(temp1, temp2, temp3, 
                               temp4,temp5, temp6, sep = "_"), 
         Genotype = paste(temp1, temp2, sep = "_"), 
         Stage = temp7) %>%
  select(-c(temp1, temp2, temp3, temp4, temp5, temp6, temp7))

#combining data for all the genotypes
data_H2DCF_all <- rbind(data_H2DCF_wt, data_H2DCF_1508)

#Genotype as factor, combining substages 6a-e into early, middle and late stage 6 
data_H2DCF_all$Genotype <- factor(data_H2DCF_all$Genotype, 
                                  levels = c("wt", "1508_9", "1508_3"))

levels(data_H2DCF_all$Genotype)<- c("WT", "NaN1508 #9", "NaN1508 #3")

data_H2DCF_all <- data_H2DCF_all %>%
  mutate(Stage = ifelse(Stage %in% c("stage6a", "stage6b", "stage6dab"), "6ab", Stage)) %>%
  mutate(Stage = ifelse(Stage %in% c("stage6d", "stage6e"), "6de", Stage)) %>%
  mutate(Stage = ifelse(Stage %in% c("stage6ab", "6ab"), "early stage 6", Stage)) %>%
  mutate(Stage = ifelse(Stage %in% c("stage6c", "6c"), "middle stage 6", Stage)) %>%
  mutate(Stage = ifelse(Stage %in% c("stage6de", "6de"), "late stage 6", Stage))

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/ROS/H2DCF_april_may_walking/R_output")
write.csv(data_H2DCF_all, file = "data_H2DCF_raw_april_may_walking.csv", quote = FALSE, 
          row.names = FALSE)


#checking normality
#dividing data according to genotypes and stages

#read in the data if you have assembled it from raw results
setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/ROS/H2DCF_april_may_walking/R_output") 
data <- read.csv("data_H2DCF_raw_april_may_walking.csv")

data_wt_a <- data %>%
  filter(Genotype == "WT" & Stage == "early stage 6")

data_wt_b <- data %>%
  filter(Genotype == "WT" & Stage == "middle stage 6")

data_wt_c <- data %>%
  filter(Genotype == "WT" & Stage == "late stage 6")


data_1508_9_a <- data %>%
  filter(Genotype == "NaN1508 #9" & Stage == "early stage 6")

data_1508_9_b <- data %>%
  filter(Genotype == "NaN1508 #9" & Stage == "middle stage 6")

data_1508_9_c <- data %>%
  filter(Genotype == "NaN1508 #9" & Stage == "late stage 6")



data_1508_3_a <- data %>%
  filter(Genotype == "NaN1508 #3" & Stage == "early stage 6")

data_1508_3_b <- data %>%
  filter(Genotype == "NaN1508 #3" & Stage == "middle stage 6")

data_1508_3_c <- data %>%
  filter(Genotype == "NaN1508 #3" & Stage == "late stage 6")

#checking for outliers
#outliers - points outside of the interquantile range: IQR = Q3 - Q1, outlier below [Q1- (1.5)IQR] or above [Q3+(1.5)IQR]

#list of all the data frames (dataset divided according to genotype and stage)
df_list <- list(data_wt_a, data_wt_b, data_wt_c,
                data_1508_9_a, data_1508_9_b, data_1508_9_c,
                data_1508_3_a, data_1508_3_b, data_1508_3_c)

df_list_no <- list() #list to save data frames without the outliers

#loop to find outliers
for(i in 1:length(df_list)){
  df <- df_list[[i]]
  Q <- quantile(df$Mean, probs = c(0.05, 0.95))
  df_no <- subset(df, df$Mean >= Q[1] & 
                    df$Mean <= Q[2])
  df_list_no[[i]] <- df_no
}

data <- bind_rows(df_list_no) 

data <- data %>% 
  mutate(Genotype = ifelse(Genotype == "NaN1508 #9", "bdpox", Genotype),
         Genotype = ifelse(Genotype == "NaN1508 #3", "BdPOX-sg", Genotype))


data$Genotype <- factor(data$Genotype,
                        levels = c("WT", "bdpox", "BdPOX-sg"))

#saving data without outliers
setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/ROS/H2DCF_april_may_walking/R_output")
write.csv(data, file = "data_H2DCF_april_may_walking_no_outliers.csv", quote = FALSE, 
          row.names = FALSE)

#loop for Shaphiro_Wilk normality test
shapiro_p <- list()

for(i in 1:length(df_list_no)){
  df <- df_list_no[[i]]
  normtest <- shapiro.test(df$Mean)
  p <- normtest$p.value
  shapiro_p[[i]] <- p
}

shapiro_p

#data isn't normally distributed --> Kruskal-Wallis
data_a <- bind_rows(df_list_no[[1]], df_list_no[[4]], df_list_no[[7]])
data_b <- bind_rows(df_list_no[[2]], df_list_no[[5]], df_list_no[[8]])
data_c <- bind_rows(df_list_no[[3]], df_list_no[[6]], df_list_no[[9]])

kruskal.test(Mean ~ Genotype, data = data_a)#p-value 0.3.14e-06, difference
kruskal.test(Mean ~ Genotype, data = data_b) #p-value 0.1795 no difference
kruskal.test(Mean ~ Genotype, data = data_c) #p-value 0.4328 no difference

#there is difference between groups --> Dunn's post hoc test

dunnTest(Mean ~ Genotype, data = data_a, method = "bh")
dunnTest(Mean ~ Genotype, data = data_b, method = "bh")
dunnTest(Mean ~ Genotype, data = data_c, method = "bh")



#creating p-values labels from the Dunn's posthoc results (rstatix package)
stat_test <- data %>%
  group_by(Stage) %>%
  dunn_test(Mean ~ Genotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat_test <- stat_test %>% add_y_position(fun = "max", step.increase = 0.3)

#load data if you are plotting later
setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/ROS/H2DCF_april_may_walking/R_output")
data <- read.csv("data_H2DCF_april_may_walking_no_outliers.csv")

#summary data
data_H2DCF_all_sum <-data %>%
  group_by(Genotype, Stage) %>%
  summarize(Mean_mean_intensity = mean(Mean),
            sd_mean_intensity= sd(Mean),
            count = n())

#project as character to use for aes
data$Project <- as.character(data$Project)
data$Project <- factor(data$Project)

#Genotype factor levels
data$Genotype <- factor(data$Genotype,
                        levels = c("WT", "bdpox", "BdPOX-sg"))

#plot
H2DCF_plot <- ggplot(data, aes(Genotype, Mean)) +
  facet_grid(.~factor(Stage, levels = c("early stage 6", "middle stage 6", "late stage 6"))) +
  geom_boxplot(aes(colour = Genotype, fill = Genotype), alpha = 0.3, 
               lwd = 1,
               outlier.shape = NA) + 
  geom_jitter(aes(colour = Genotype), width = 0.1, height = 0.05,
              alpha = 0.8) + 
  scale_fill_manual(values = c("#8BCF80", "#D95450", "#70B2A3")) +
  scale_color_manual(values = c("#8BCF80", "#D95450", "#70B2A3")) +
  theme(strip.text.x = element_text(size = 22, face = "bold"), 
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", colour = "black", size = 1),
        axis.line = element_line(size = 2, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.text.x = element_text(size = 22, face = "bold", colour = "black", 
                                   angle = 45, margin = unit(c(1.2, 0, 0, 0), "cm")),
        axis.text.y = element_text(size = 22, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 22, face = "bold", colour = "black"),
        axis.ticks.length = unit(0.5, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0,320), breaks = seq(0, 300, 50), expand = c(0,0)) +
  labs(x = element_blank(), y = "Mean fluorescence intensity of the guard cells") + 
  stat_compare_means(size = 5, label.y = 300) +
  stat_pvalue_manual(stat_test, label = "p.adj.signif",  y.position = "y.position", hide.ns = TRUE) +
  geom_text(data = data_H2DCF_all_sum, mapping = aes(Genotype, Inf, label = paste("n = ", count)), 
            size = 5, vjust = "inward")

H2DCF_plot

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/H2DCF")

pdf("H2DCF_walking_combined.pdf", height = 10, width= 16)
print(H2DCF_plot, device=cairo_pdf)
dev.off()


######
#all data, but separate boxplots for the day
######

data_gr <- data %>%
  group_by(Project, Stage, Genotype)

data_gr$Stage <- factor(data_gr$Stage, levels = c("early stage 6", "middle stage 6", "late stage 6"))

data_gr_sum <- data_gr %>%
  summarize(Mean_mean_intensity = mean(Mean),
            sd_mean_intensity= sd(Mean),
            count = n())

#creating p-values labels from the Dunn's posthoc results (rstatix package)
stat_test_gr<- data_gr %>%
  group_by(Project, Stage) %>%
  dunn_test(Mean ~ Genotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat_test_gr <- stat_test_gr %>% add_y_position(fun = "max", step.increase = 0.05)

H2DCF_plot3 <- ggplot(data_gr, aes(Genotype, Mean)) +
  facet_grid(Stage~Project) +
  geom_boxplot(aes(colour = Genotype, fill = Genotype), alpha = 0.3, 
               lwd = 0.75,
               outlier.shape = NA) + 
  geom_jitter(aes(colour = Genotype), width = 0.1, height = 0.05,
              alpha = 0.8, size = 3) + 
  scale_fill_manual(values = c("#8BCF80", "#D95450", "#70B2A3")) +
  scale_color_manual(values = c("#8BCF80", "#D95450", "#70B2A3")) +
  theme(strip.text.x = element_text(size = 22, face = "bold"), 
        strip.text.y = element_text(size = 22, face = "bold"), 
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white", colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major.y = element_line(colour = "gray80", size = 0.5),
        panel.grid.major.x = element_blank(),
        axis.line.x =  element_line(size = 1),
        axis.line.y =  element_line(size = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 22, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 22, face = "bold", colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 22, face = "bold"),
        legend.background = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(limits = c(0,375), breaks = seq(0, 350, 50), expand = c(0,0)) +
  labs(x = element_blank(), y = "Mean fluorescence intensity of the guard cells") +
  stat_compare_means(size = 3, label.y = 325) +
  stat_pvalue_manual(stat_test_gr, label = "p.adj.signif",  y.position = "y.position", hide.ns = TRUE) +
  geom_text(data = data_gr_sum, mapping = aes(Genotype, Inf, label = paste("n = ", count)), 
            size = 3, vjust = "inward")

H2DCF_plot3

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/H2DCF")

pdf("H2DCF_may_byday2_walking.pdf", height = 35, width= 15)
print(H2DCF_plot3, device=cairo_pdf)
dev.off()
