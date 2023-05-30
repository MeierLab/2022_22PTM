library(plyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(Peptides)
library(data.table)
library(iemisc)
library(VennDiagram)
library(ggrepel)
library(ggnewscale)


#Folder with MaxQuant evidence files
WD_Evidence = ""

#Folder for exported figures
WD_Graphs = ""

#Folder for exported files
WD_Files = ""

#color scheme
Colors = c("#941C5A", "#577B80", "#679e1aff", "#00bbd4ff", "#ea8d0bff","#941C5A", "#577B80", "#679e1aff", "#00bbd4ff", "#ea8d0bff","#941C5A")

#### Acquire data, iRT_final_median_cor ####

setwd(WD_Files)

# Read R objects from working directory
iRT_final_median_cor = list.files(pattern = "iRT_final_median_cor") %>%
  map(readRDS)

iRT_final_median_cor = iRT_final_median_cor [[1]]

##### Calculate deviations ######

x = split(iRT_final_median_cor, iRT_final_median_cor$Sequence)

y = lapply(x, function(a) sd(a$CCS_cor)) %>% as.data.frame() %>% t() %>% as.data.frame()
z = lapply(x, function(a) mean(a$CCS_cor)) %>% as.data.frame() %>% t() %>% as.data.frame()

z$SD = y$V1

colnames(z) = c("Mean", "SD")

##### Fig.2A_CCS drifts across measurements + Alignment #####

plot_list = list()

#dotplot CCS vs run number (raw and corrected)
plot1 = ggplot(iRT_final_median_cor, aes(x = CCS, y = ID2, col = Sequence, group = Sequence))  + geom_point(size = 0.5, show.legend = FALSE) +
  xlab('CCS') + ylab('Run Number') + theme_classic() + theme(legend.position = 'top')+
  scale_color_manual(values=Colors) +
  scale_fill_manual (values=Colors) +
  scale_x_continuous(breaks = seq(300, 440, by = 25))+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

plot_list$Uncorrected = plot1

plot2 = ggplot(iRT_final_median_cor, aes(x = CCS_cor, y = ID2, col = Sequence, group = Sequence))  + geom_point(size = 0.5, show.legend = FALSE) +
  xlab('CCS') + ylab('Run Number') + theme_classic() + theme(legend.position = 'top')+
  scale_color_manual(values=Colors) +
  scale_fill_manual (values=Colors)+
  scale_x_continuous(breaks = seq(300, 440, by = 25))+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

plot_list$Corrected = plot2


setwd(WD_Graphs)

pdf(paste("Fig.2A_iCCS_timeline.pdf", sep=''), width = 5.2, height = 2)
for (i in 1:(length(plot_list))) {
  print(plot_list[[i]])
}
dev.off()


#plot RT correction over time
plot_list = list()
#dotplot CCS vs run number (raw and corrected)
plot1 = ggplot(iRT_final_median_cor, aes(x = RT, y = ID2, col = Sequence, group = Sequence))  + geom_point(size = 0.8, show.legend = FALSE) +
  xlab('RT') + ylab('Run Number') + theme_classic() + theme(legend.position = 'top')+
  scale_color_manual(values=Colors) +
  scale_fill_manual (values=Colors) +
  #scale_x_continuous(breaks = seq(300, 440, by = 25))+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

plot_list$Uncorrected = plot1

plot2 = ggplot(iRT_final_median_cor, aes(x = RT_cor, y = ID2, col = Sequence, group = Sequence))  + geom_point(size = 0.8, show.legend = FALSE) +
  xlab('RT') + ylab('Run Number') + theme_classic() + theme(legend.position = 'top')+
  scale_color_manual(values=Colors) +
  scale_fill_manual (values=Colors)+
  #scale_x_continuous(breaks = seq(300, 440, by = 25))+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

plot_list$Corrected = plot2

pdf(paste("Fig.2A_iRT_timeline.pdf", sep=''), width = 5.2, height = 2)
for (i in 1:(length(plot_list))) {
  print(plot_list[[i]])
}
dev.off()



#### Acquire data, CV_replicates_df####
setwd(WD_Files)

Pattern = "CV"
Shorten = c('CV_replicates_df_', ".txt")
Modifications = c('Kmod .', 'Rmod .', 'Ymod .', 'Pmod .', "SpikeMix .")

my.files.path = list.files(pattern = Pattern)

# Remove unwanted characters from file name for naming results
my.files = str_remove_all(my.files.path, paste(Shorten, collapse = "|"))


# read .txt files to data frames
my.data = lapply(my.files.path, 
                 read.csv,
                 header=TRUE, sep="\t")


x = my.data$`C3_ Kmod `
y = split(x, x$ID)

z = my.data$`C3_ Kmod `
n = split(z, z$ID)


names(my.data) = my.files

my.data = rbindlist(my.data)
my.data$ID2 =  gsub("\\_.*","",my.data$ID)

x = split(my.data, my.data$ID)

##### Significance testing ####

library(rstatix)

tTest = t.test(my.data$CV_CCS, my.data$CV_CCS_cor) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance()


tmp = my.data %>%  select(c("CV_CCS", "CV_CCS_cor")) %>% stack()

stat.test = tmp %>%
  t_test(values ~ ind) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

tmp = select(my.data, c("CV_CCS", "CV_CCS_cor")) %>% stack()

med = median(my.data$CV_CCS_cor)
mean = mean(my.data$CV_CCS_cor)

ggplot(tmp, aes(x = ind, y = values)) + geom_boxplot() +
  theme_minimal() + 
  xlab('') + 
  ylab("CV (%)") +
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12)) 



##### Fig.2B_CCS and RT Coefficient of variation ####

plot_list = list()

y = my.data %>% filter(quantile(CV_CCS, 0.99)>CV_CCS)

plot = ggplot(y, aes(x=CV_CCS)) + geom_histogram(binwidth=0.05, color = "#941C5A", fill = "#941C5A") + 
  theme_classic() + 
  ylim(0, 660) + xlim(0, 4.6) + 
  ylab("Count") + xlab("Coefficient of variation (%)")+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

plot_list[["Uncorrected"]] = plot

y = my.data %>% filter(quantile(CV_CCS_cor, 0.99)>CV_CCS_cor)

plot = ggplot(y, aes(x=CV_CCS_cor)) + geom_histogram(binwidth=0.05, color = "#577B80", fill = "#577B80") + 
  theme_classic() + 
  ylim(0, 660) + xlim(0, 4.6) + 
  ylab("Count") + xlab("Coefficient of variation (%)")+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

CV_uncor = median(my.data$CV_CCS)
CV_cor = median(my.data$CV_CCS_cor)

plot_list[["Corrected"]] = plot

setwd(WD_Graphs)

pdf(paste("Fig.2B_CV_CCS_correction.pdf", sep=''), width = 2.6, height = 2)
for (i in 1:(length(plot_list))) {
  print(plot_list[[i]])
}
dev.off()


plot_list = list()
plot = ggplot(my.data, aes(x=CV_RT)) + geom_histogram(binwidth=0.05, color = "#941C5A", fill = "#941C5A") + 
  theme_classic() + 
  ylim(0, 250) + xlim (0, 50) + 
  ylab("Count") + xlab("Coefficient of variation (%)")+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

plot_list[["Uncorrected"]] = plot


plot = ggplot(my.data, aes(x=CV_RT_cor)) + geom_histogram(binwidth=0.05, color = "#577B80", fill = "#577B80") + 
  theme_classic() + 
  ylim(0, 250) + xlim (0, 50) + 
  ylab("Count") + xlab("Coefficient of variation (%)")+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

CV_RT_uncor = median(my.data$CV_RT)
CV_RT_cor = median(my.data$CV_RT_cor)

plot_list[["Corrected"]] = plot

setwd(WD_Graphs)

pdf(paste("Fig.2B_CV_RT_correction.pdf", sep=''), width = 2.6, height = 2)
for (i in 1:(length(plot_list))) {
  print(plot_list[[i]])
}
dev.off()

