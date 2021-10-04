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

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/gsmax_bdpox_discussion")

data_SDSL <- read.csv("bdpox_discussion_data_SLSD.csv")
data_gsmax <- read.csv("bdpox_discussion_data_gsmax.csv")

data_SDSL$season <- factor(data_SDSL$season, levels = c("summer", "autumn", "winter"))

data_s <- data_SDSL %>%
  filter(season == "summer")

data_a <- data_SDSL %>%
  filter(season == "autumn")

data_w <- data_SDSL %>%
  filter(season == "winter")

shapiro.test(data_s$SD) #normal
shapiro.test(data_a$SD) #normal 
shapiro.test(data_w$SD) #normal
shapiro.test(data_s$SL) #normal
shapiro.test(data_a$SL) #normal
shapiro.test(data_w$SL) #normal

stat_test_SD <- data_SDSL %>%
  tukey_hsd(SD ~ season) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

SD_plot <- ggplot(data_SDSL, aes(season, SD)) +
  geom_boxplot(aes(colour = season, fill = season), alpha = 0.3, 
               lwd = 1,
               outlier.shape = NA) + 
  geom_jitter(aes(colour = season),
              alpha = 1, size  = 3, stroke = 2, shape = 1, height = 0, width = 0.1) + 
  scale_fill_manual(values = c("#D95450", "#8BCF80","#70B2A3")) +
  scale_color_manual(values = c("#D95450", "#8BCF80","#70B2A3")) +
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
  scale_y_continuous(limits = c(50, 150), breaks = seq(50, 150, 25), expand = c(0,0)) +
  labs(x = element_blank(), 
       y = expression(bold(Stomatal~density~(stomata/mm^2)))) 

SD_plot

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/gsmax_bdpox_discussion")

pdf("SD_plot.pdf", height = 10, width= 6)
print(SD_plot, device=cairo_pdf)
dev.off()

stat_test_SL <- data_SDSL %>%
  tukey_hsd(SL ~ season) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

SL_plot <- ggplot(data_SDSL, aes(season, SL)) +
  geom_boxplot(aes(colour = season, fill = season), alpha = 0.3, 
               lwd = 1,
               outlier.shape = NA) + 
  geom_jitter(aes(colour = season),
              alpha = 1, size  = 3, stroke = 2, shape = 1, height = 0, width = 0.1) + 
  scale_fill_manual(values = c("#D95450", "#8BCF80","#70B2A3")) +
  scale_color_manual(values = c("#D95450", "#8BCF80","#70B2A3")) +
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
  scale_y_continuous(limits = c(22, 30), breaks = seq(22, 30, 2), expand = c(0,0)) +
  labs(x = element_blank(), 
       y = expression(bold(Stomata~length~(µm)))) 

SL_plot

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/gsmax_bdpox_discussion")

pdf("SL_plot.pdf", height = 10, width= 6)
print(SL_plot, device=cairo_pdf)
dev.off()

#gsmax

data_s <- data_gsmax %>%
  filter(season == "summer")

data_a <- data_gsmax %>%
  filter(season == "autumn")


shapiro.test(data_s$gsmax) #normal
shapiro.test(data_a$gsmax) #normal 

data_gsmax <- data_gsmax %>%
  mutate(season = ifelse(season == "autumn", "autumn/winter", season))

data_gsmax$season <- factor(data_gsmax$season, levels = c("summer", "autumn/winter"))

stat_test_gsmax <- data_gsmax %>%
  t_test(gsmax ~ season) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

gsmax_plot <- ggplot(data_gsmax, aes(season, gsmax)) +
  geom_boxplot(aes(colour = season, fill = season), alpha = 0.3, 
               lwd = 1,
               outlier.shape = NA) + 
  geom_jitter(aes(colour = season),
              alpha = 1, size  = 3, stroke = 2, shape = 1, height = 0, width = 0.1) + 
  scale_fill_manual(values = c("#D95450", "#70B2A3")) +
  scale_color_manual(values = c("#D95450","#70B2A3")) +
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
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2), expand = c(0,0)) +
  labs(x = element_blank(), 
       y = expression(bold(~italic(g)[smax]~(mol~m^-2*s^-1)))) 

gsmax_plot

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/gsmax_bdpox_discussion")

pdf("gsmax_plot.pdf", height = 10, width= 6)
print(gsmax_plot, device=cairo_pdf)
dev.off()

#read in summary from thesis data
setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/gsmax_bdpox_discussion")
bdpox_data_sum <- read.csv("bdpox_project_gsmax_comparison.csv")

Genotype <- rep("WT", nrow(data_SDSL))
data_SDSL <- cbind(data_SDSL, Genotype)


Genotype <- rep("WT", nrow(data_gsmax))
data_gsmax <- cbind(data_gsmax, Genotype)

data_SDSL_sum <- data_SDSL %>%
  group_by(Genotype, season) %>%
  summarize(Mean_SD = mean(SD),
            SD_SD = sd(SD),
            Mean_SL = mean(SL),
            SD_SL = sd(SL))

data_gsmax_sum <- data_gsmax %>%
  group_by(Genotype, season) %>%
  summarize(Mean_gsmax = mean(gsmax),
            SD_gsmax = sd(gsmax))

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/gsmax_bdpox_discussion")

write.csv(data_SDSL_sum, file = "data_SDSL_sum.csv", quote = FALSE, 
          row.names = FALSE)

write.csv(data_gsmax_sum, file = "data_gsmax_sum.csv", quote = FALSE, 
          row.names = FALSE)

#compare 2021 data and Tiago's data

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/BdPOX_anatomy/output")

bdpox_data <- read.csv("bdpox_project_all_data.csv")

bdpox_data <- bdpox_data %>%
  filter(Individual != "wt_2")

bdpox_data_wt <- bdpox_data %>%
  filter(Genotype == "WT")

bdpox_data_wt <- bdpox_data_wt %>%
  select(Individual, Genotype, GC_length_ind, GC_length_DIC, SD,
         Physiological_gsmax, gsmax_confocal, gsmax_DIC_poles)

season <- rep("summer2021", nrow(bdpox_data_wt))
bdpox_data_wt <- cbind(bdpox_data_wt, season)

data_SDSL <- data_SDSL %>%
  rename(Individual = ï..Sample_id)

bdpox_data_wt_SDSL <- bdpox_data_wt %>%
  select(Individual, Genotype, GC_length_DIC, SD, season) %>%
  rename(SL = GC_length_DIC)

SDSL_all <- rbind(data_SDSL, bdpox_data_wt_SDSL)

shapiro.test(bdpox_data_wt$GC_length_DIC)

stat_test_SL21 <- SDSL_all %>%
  tukey_hsd(SL ~ season) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

SL_plot21 <- ggplot(SDSL_all, aes(season, SL)) +
  geom_boxplot(aes(colour = season, fill = season), alpha = 0.3, 
               lwd = 1,
               outlier.shape = NA) + 
  geom_jitter(aes(colour = season),
              alpha = 1, size  = 3, stroke = 2, shape = 1, height = 0, width = 0.1) + 
  scale_fill_manual(values = c("#D95450", "#8BCF80","#70B2A3", "#FFDB6D")) +
  scale_color_manual(values = c("#D95450", "#8BCF80","#70B2A3", "#FFDB6D")) +
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
  scale_y_continuous(limits = c(22, 30), breaks = seq(22, 30, 2), expand = c(0,0)) +
  labs(x = element_blank(), 
       y = expression(bold(Stomata~length~(µm)))) 

SL_plot21

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/gsmax_bdpox_discussion")

pdf("SL_plot21.pdf", height = 10, width= 6)
print(SL_plot21, device=cairo_pdf)
dev.off()

shapiro.test(bdpox_data_wt$SD)

stat_test_SD21 <- SDSL_all %>%
  tukey_hsd(SD ~ season) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

SD_plot21 <- ggplot(SDSL_all, aes(season, SD)) +
  geom_boxplot(aes(colour = season, fill = season), alpha = 0.3, 
               lwd = 1,
               outlier.shape = NA) + 
  geom_jitter(aes(colour = season),
              alpha = 1, size  = 3, stroke = 2, shape = 1, height = 0, width = 0.1) + 
  scale_fill_manual(values = c("#D95450", "#8BCF80","#70B2A3", "#FFDB6D")) +
  scale_color_manual(values = c("#D95450", "#8BCF80","#70B2A3", "#FFDB6D")) +
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
  scale_y_continuous(limits = c(50, 150), breaks = seq(50, 150, 25), expand = c(0,0)) +
  labs(x = element_blank(), 
       y = expression(bold(Stomatal~density~(stomata/mm^2)))) 

SD_plot21

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/gsmax_bdpox_discussion")

pdf("SD_plot21.pdf", height = 10, width= 6)
print(SD_plot21, device=cairo_pdf)
dev.off()

#gsmax

bdpox_data_wt_gsmax <- bdpox_data_wt %>%
  select(Individual, Genotype, Physiological_gsmax, season) %>%
  rename(gsmax = Physiological_gsmax)

data_gsmax <- data_gsmax %>%
  rename(Individual = ï..Sample_id)

gsmax_all <- rbind(bdpox_data_wt_gsmax, data_gsmax)

shapiro.test(bdpox_data_wt_gsmax$gsmax)

stat_test_gsmax21 <- gsmax_all %>%
  tukey_hsd(gsmax ~ season) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

gsmax_all$season <- factor(gsmax_all$season, levels = c("summer", "autumn/winter", "summer2021"))

gsmax_plot21 <- ggplot(gsmax_all, aes(season, gsmax)) +
  geom_boxplot(aes(colour = season, fill = season), alpha = 0.3, 
               lwd = 1,
               outlier.shape = NA) + 
  geom_jitter(aes(colour = season),
              alpha = 1, size  = 3, stroke = 2, shape = 1, height = 0, width = 0.1) + 
  scale_fill_manual(values = c("#D95450", "#70B2A3", "#FFDB6D")) +
  scale_color_manual(values = c("#D95450","#70B2A3", "#FFDB6D")) +
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
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2), expand = c(0,0)) +
  labs(x = element_blank(), 
       y = expression(bold(~italic(g)[smax]~(mol~m^-2*s^-1)))) 

gsmax_plot21

setwd("C:/Users/Magdalena/Desktop/MAPS/Master_thesis/thesis_figures/gsmax_bdpox_discussion")

pdf("gsmax_plot21.pdf", height = 10, width= 6)
print(gsmax_plot21, device=cairo_pdf)
dev.off()