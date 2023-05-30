library(tidyverse)
library(magrittr)
library(data.table)
library(iemisc)
library(VennDiagram)
library(ggrepel)
library(ggnewscale)


#Folder for saving graphs
WD_Graphs = 'C:/Users/Andre/OneDrive/Desktop/Paper_plots'

#Folder with previously exported files (Cor_fac_CCS_list, Cor_fac_RT_list)
WD_Files = 'C:/Users/Andre/Data/20220411_21PTM_2/evidence/Exported_files'


Colors = rep(c("#941C5A", "#577B80", "#679e1aff", "#00bbd4ff", "#ea8d0bff"),5)


#### Acquire data Top_int_n3 ####
setwd(WD_Files)

my.files.path = list.files(pattern = "Top_int_n3")

# Cut together data path for naming
Shorten = c('Top_int_n3 ', " .rds")
my.files = str_remove_all(my.files.path, paste(Shorten, collapse = "|"))

# Read R objects from working directory
my.data = list.files(pattern = "Top_int_n3") %>%
  map(readRDS)

names(my.data) = my.files

#Exclude unmodified peptides
data = my.data[1:5]

#Introduce Modification identifier to differentiate identical Sequences with different charges
data_list = list()
for (i in 1:length(data)) {
  tmp = data [[i]]
  name = names(data) [[i]]
  
  tmp_list = list()
  for (j in 1:length(tmp)) {
    
    tmp2 = tmp [[j]]
    name2 = names(tmp) [[j]]
    tmp2$ID2 = name2
    
    name_3 = paste(i,j, sep = "_")
    
    tmp_list [[name_3]] = tmp2
    
  }
  
  
  data_list [[name]] = tmp_list
  
}

##### Supplemental_fig_2; All peptides CCS vs. m.z ####

x = lapply(data_list, function(a) rbindlist(a))
y = x %>% rbindlist()

#split by charge
split = split(y, y$Charge)

#for each charge get unique Sequence for each experiment
z = lapply(split, function(a) distinct(a, Sequence, ID2, .keep_all = TRUE)) %>% rbindlist()

number_of_unique_mod_C_peptides = z


plot = ggplot(z, aes(y = CCS_cor, x = m.z)) + geom_point(size =0.5, alpha=0.2)  + 
  theme_classic()+
  xlab("m/z") + ylab("CCS (A2)") +
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

setwd(WD_Graphs)

pdf("Sup_fig_2_All_unique_peptides.pdf", width = 5, height = 4)
print(plot)
dev.off()



split = split(z, z$ID2)

# Information for ordering
Modifications_full_no_unmod = c("Acetylation", "Biotinylation", "Butyrylation", "Crotonylation", "Dimethylation", "Formylation", "Glutarylation", "Hydroxyisobutyrylation", "Malonylation", "Methylation", "Propionylation", "Succinylation", "Trimethylation", "Ubiquitination", "Hydroxyproline",  "Citrullination", "Dimethylation-asym", "Dimethylation-sym", "Methylation",  "O-GlcNAcylation",  "Nitrotyrosine", "Phosphorylation")
Residue_no_unmod = c(rep("K", 14), rep("P", 1),rep("R", 4),rep("S/T", 1),rep("Y", 2))
Ordering_no_unmod = c(3, 11, 5, 9, 14, 2, 8, 10, 6, 13, 4, 7, 15, 12, 22, 21, 20, 19, 18, 24, 26, 27)

data = rbindlist(my.data, idcol = TRUE)
split = split(data, data$.id)

data = list()
for (i in 1:length(split)) {
  tmp = split[[i]]
  name = names(split) [[i]]
  
  tmp$Ordering = Ordering_no_unmod [[i]]
  tmp$Modifications_full = Modifications_full_no_unmod [[i]]
  tmp$label = paste("(", tmp$Residue, ") ", tmp$Modifications_full, sep = "" )
  data [[name]] = tmp
}


##### Supplemental_fig_1, Histograms (Andromeda_score, Modifications probabilities ####
data = rbindlist(data)
data$Probabilities_Extract = data$Probabilities_Extract %>% as.numeric()

data = transform(data, label = reorder(label, Ordering) )


plot1 = ggplot(data, aes(Score)) + 
  #geom_vline(xintercept = 40, color = "red") +
  geom_histogram(binwidth = 10, center = 45, color = "#577B80", fill = "#577B80") + theme_classic() +
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 11)) + 
  #xlim (0, 400) +
  xlab("Andromeda score") + ylab("Count")

plot2 = ggplot(data, aes(Probabilities_Extract)) + 
  #geom_vline(xintercept = 40, color = "red") +
  geom_histogram(binwidth = 0.01, color = "#577B80", fill = "#577B80") + theme_classic() +
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 11)) + 
  xlab("Modification Probabilities") + ylab("Count")

plot_list = list(plot1, plot2)

setwd('C:/Users/Andre/OneDrive/Desktop/Paper_plots')

pdf("Sup_fig_1_Histogram_top_int_n3_C2C3.pdf", width = 3.3, height = 2)
for (i in 1:(length(plot_list))) {
  print(plot_list[[i]])
}
dev.off()

##### Supplemental_fig_1, Box plots (Andromeda_score, Modifications probabilities ####

split = split(data, data$Modifications)

x_med = ddply(data, .(label), summarise, med = median(Score))

y = split(data, data$label)
x_med$x_n = sapply(y, NROW)


plot1 = ggplot(data, aes(x=label, y= Score)) + 
  geom_boxplot(aes(fill = Residue), size = 0.3, outlier.size = 0.3, show.legend = FALSE) +
  xlab('') + ylab('Andromeda Score') +
  scale_x_discrete(guide = guide_axis(angle = 40))+ 
  theme_classic()+
  scale_color_manual(values=Colors) +
  scale_fill_manual (values=Colors)+
  geom_text(data = x_med, aes(x = label, y= 410, label = round(x_n, digits = 1)), position = position_dodge(width = 1), size = 2) +
  geom_text(data = x_med, aes(x = label, y= 460, label = round(med, digits = 1)), position = position_dodge(width = 1), size = 2) +
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 11))



plot2 = ggplot(data, aes(x=label, y= Probabilities_Extract)) + 
  geom_boxplot(aes(fill = Residue), size = 0.3, outlier.size = 0.3, show.legend = FALSE) +
  xlab('') + ylab('Modification Prob.') +
  scale_x_discrete(guide = guide_axis(angle = 40))+ 
  theme_classic()+
  scale_color_manual(values=Colors) +
  scale_fill_manual (values=Colors)+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 11), panel.grid.major.y = element_line(colour = "grey", size = (0.2)))



score_median = median(data$Score)
score_min = min(data$Score)

prob_median = median(data$Probabilities_Extract)
prob_min = min(data$Probabilities_Extract)


plot_list = list()

plot_list [[1]] = plot1
plot_list [[2]] = plot2


pdf(paste("Sup_fig_1_box_Andromeda_Probability.pdf", sep=''), width = 7, height = 3)
for (i in 1:(length(plot_list))) {
  print(plot_list[[i]])
}
dev.off()

