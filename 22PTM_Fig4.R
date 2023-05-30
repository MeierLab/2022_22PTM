library(tidyverse)
library(magrittr)
library(data.table)
library(iemisc)
library(VennDiagram)
library(ggrepel)
library(ggnewscale)
library(rstatix)

#Folder with MaxQuant evidence files
WD_Evidence = 'C:/Users/Andre/Data/20220411_21PTM_2/evidence'

#Folder for saving graphs
WD_Graphs = 'C:/Users/Andre/OneDrive/Desktop/Paper_plots'

#Folder with previously exported files (Cor_fac_CCS_list, Cor_fac_RT_list)
WD_Files = 'C:/Users/Andre/Data/20220411_21PTM_2/evidence/Exported_files'


Modifications = c('Kmod .', 'Rmod .', 'Ymod .', 'Pmod .', "Spike .")
Colors_2 = rep(c("#941C5A", "#577B80", "#679e1aff", "#00bbd4ff", "#ea8d0bff"),5)
white_black = c("white", "grey", "black")

#Parameters for labeling and ordering
Modifications_full = c("Acetyl", "Biotinyl", "Butyryl", "Crotonyl", "Dimethyl", "Formyl", "Glutaryl", "Hydroxyisobutyryl", "Malonyl", "Methyl", "Propionyl", "Succinyl", "Trimethyl", "GlyGly", "Unmodified", "Hydroxyproline", "Unmodified",  "Citrullin", "Dimethyl-asym", "Dimethyl-sym", "Methyl", "Unmodified", "O-GlcNAcyl", "Unmodified", "Nitrotyrosine", "Phosphoryl", "Unmodified")
Residue = c(rep("K", 15), rep("P", 2),rep("R", 5),rep("S/T", 2),rep("Y", 3))
Ordering = c(3, 11, 5, 9, 14, 2, 8, 10, 6, 13, 4, 7, 15, 12, 1, 22, 21, 20, 19, 18, 17, 16, 24, 23, 26, 27, 25)


#### Acquire data Matched_list_C2_C3 ####
setwd(WD_Files)

my.files.path = list.files(pattern = "Matched_list_C2_C3")

Shorten = c('Matched_list_C2_C3 ', " .rds")
my.files = str_remove_all(my.files.path, paste(Shorten, collapse = "|"))

# Read R objects from working directory
my.data = list.files(pattern = "Matched_list_C2_C3") %>%
  map(readRDS)

names(my.data) = my.files

my.data = flatten(my.data)

# Export Matched_peptides.tsv
x = rbindlist(my.data)
setwd(WD_Graphs)
write.table(x, file='Matched_peptides.tsv', quote=FALSE, row.names = FALSE, col.names = TRUE, sep='\t')

data_list = list()
for (i in 1:length(my.data)) {
  
  data = my.data [[i]]
  
  name = names(my.data) [[i]]
  
  x = Residue [[i]]
  
  data$Modifications_full = paste( "(", x, ") ", Modifications_full [[i]] ,sep = "")
  data$mod = Residue [[i]]
  data$order = Ordering [[i]]
  
  data_list [[name]] = data
}

my.data = data_list[!grepl("unmod", names(data_list))] 


##### Fig.4b/c CCS shift vs mz and RT dotplot #####

delta_list = list()
delta_median = list()
for (i in 1:length(my.data)) {
  data = my.data [[i]]
  name = names(my.data) [[i]]
  
  data$Charge = data$Charge %>% as.numeric()
  
  data$dCCS = ((data[,"CCS_cor_mod"]-data[,"CCS_cor_unmod"])/data[,"CCS_cor_unmod"])*100
  data$dCCS_absolute = data[,"CCS_cor_mod"]-data[,"CCS_cor_unmod"]
  data$dRT = (data[,"RT_cor_mod"]-data[,"RT_cor_unmod"])
  data$dM = (data[,"m.z_mod"]-data[,"m.z_unmod"])/data[,"m.z_unmod"]*100
  data$dM_absolute = (data[,"m.z_mod"]-data[,"m.z_unmod"])*data[,"Charge"]
  
  dCCS = median(data$dCCS)
  dCCS_absolute = median(data$dCCS_absolute)
  dRT = median(data$dRT)
  dM = median(data$dM)
  dM_absolute = median(data$dM_absolute)
  mod = data$mod [[1]]
  Modifications_full = data$Modifications_full [[1]]
  
  delta_list[[name]] = data
  delta_median [[name]] = c(dCCS, dCCS_absolute, dRT, dM, dM_absolute, mod, Modifications_full)
}

delta_list = rbindlist(delta_list)

x = as.data.frame(delta_median) %>% t() %>% as.data.frame()
colnames(x) = c("dCCS", "dCCS_absolute", "dRT", "dM", "dM_absolute", "mod", "Modifications_full")

x$dCCS = x$dCCS %>% as.numeric()
x$dCCS_absolute = x$dCCS_absolute %>% as.numeric()
x$dRT = x$dRT %>% as.numeric()
x$dM = x$dM %>% as.numeric()
x$dM_absolute = x$dM_absolute %>% as.numeric()

Shorten = c('\\(K\\)', "\\(R\\) ", "\\(S/T\\) ", "\\(P\\) ", "\\(Y\\) ", "ation")
#Shorten = c('(K) ', "(R) ", "(S/T) ", "(P) ", "(Y) ")
x$Modifications = str_remove_all(x$Modifications_full, paste(Shorten, collapse = "|"))

plot_list = list()

plot1 = ggplot(x, aes(y = dCCS, x = dM))  + 
  geom_vline(xintercept=0, color = "grey") +
  geom_hline(yintercept=0, color = "grey") +
  #geom_smooth(method='lm', formula= y~x, color = "black", show.legend = FALSE, size = 0.5) +
  geom_point(aes(color = mod, fill = mod), alpha=1, shape=21, color="black", stroke = 0.2, size = 2.5, show.legend = FALSE)+
  geom_text_repel(aes(label = Modifications, color = mod), size = 3, show.legend = FALSE) +
  xlab(expression(Delta * "M (%)")) +
  ylab(expression(Delta * "CCS (%)")) +
  theme_classic() +
  scale_color_manual(values=Colors_2) +
  scale_fill_manual (values=Colors_2)+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

plot2 = ggplot(x, aes(y = dCCS_absolute, x = dM_absolute))  + 
  geom_vline(xintercept=0, color = "grey") +
  geom_hline(yintercept=0, color = "grey") +
  #geom_smooth(method='lm', formula= y~x, color = "black", show.legend = FALSE, size = 0.5) +
  geom_point(aes(color = mod, fill = mod), alpha=1, shape=21, color="black", stroke = 0.2, size = 2.5, show.legend = FALSE)+
  geom_text_repel(aes(label = Modifications, color = mod), size = 3, show.legend = FALSE) +
  xlab(expression(Delta * "M (Da)")) +
  ylab(expression(Delta * "CCS (" ~ ring(A) ^ 2 * ")")) +
  theme_classic() +
  scale_color_manual(values=Colors_2) +
  scale_fill_manual (values=Colors_2)+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

plot3 = ggplot(x, aes(y = dCCS, x = dM_absolute))  + 
  geom_vline(xintercept=0, color = "grey") +
  geom_hline(yintercept=0, color = "grey") +
  #geom_smooth(method='lm', formula= y~x, color = "black", show.legend = FALSE, size = 0.5) +
  geom_point(aes(color = mod, fill = mod), alpha=1, shape=21, color="black", stroke = 0.2, size = 2.5, show.legend = FALSE)+
  geom_text_repel(aes(label = Modifications, color = mod), size = 3, show.legend = FALSE) +
  xlab(expression(Delta * "M (Da)")) +
  ylab(expression(Delta * "CCS (%)")) +
  theme_classic() +
  scale_color_manual(values=Colors_2) +
  scale_fill_manual (values=Colors_2)+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))


tmp = lm(x$dM~x$dCCS) 
z = summary(tmp) %>% print()

z = cor(x$dM, x$dCCS, method="pearson")

plot4 = ggplot(x, aes(y = dCCS, x = dRT, color = mod, fill = mod))  + 
  geom_vline(xintercept=0, color = "grey") +
  geom_hline(yintercept=0, color = "grey") +
  geom_point(alpha=1, shape=21, color="black", stroke = 0.2, size = 2.5, show.legend = FALSE)+
  geom_text_repel(aes(label = Modifications, color = mod), size = 3, show.legend = FALSE) +
  xlab(expression(Delta * "RT (min)")) +
  ylab(expression(Delta * "CCS (%)")) +
  theme_classic() +
  scale_color_manual(values=Colors_2) +
  scale_fill_manual (values=Colors_2)+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

plot_list = list(plot1, plot2, plot3, plot4)

setwd(WD_Graphs)

pdf("Fig.4bc_CCS_vs_mz_RT_dotplot.pdf", width = 4, height = 3.7)
for (i in 1:(length(plot_list))) {
  print(plot_list[[i]])
}
dev.off()


##### Correlation in dotplot ####
tmp = lm(x$dM~x$dCCS) 
z = summary(tmp)
R2_dotplot = z$adj.r.squared

y = x[-c(2, 20),]
tmp = lm(y$dM~y$dCCS) 
z = summary(tmp)
R2_dotplot_reduced = z$adj.r.squared


##### One sample testing #####

tmp = delta_list

one_sample_wilcow = tmp %>%
  group_by(Modifications_full) %>%
  wilcox_test(dCCS ~ 1, mu = 0) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

one_sample_t_test = tmp %>%
  group_by(Modifications_full) %>%
  t_test( dCCS ~ 1, mu = 0) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")



##### Fig.4a CCS shift boxplot ####

x = delta_list
x$Charge = x$Charge %>% as.character()

split = split(x, x$Modifications_full)

data_Median_list = list()
for (i in 1:length(split)) {
  data = split [[i]]
  name = names(split) [[i]]
  
  tmp = median(data$dCCS)
  data_Median_list [[name]] = tmp
}

Ordering_no_unmod = c(3, 11, 5, 9, 14, 2, 8, 10, 6, 13, 4, 7, 15, 12, 22, 21, 20, 19, 18, 24, 26, 27)

data_Median  = cbind(data_Median = as.numeric(data_Median_list), Ordering_no_unmod = as.numeric(Ordering_no_unmod)) %>% as.data.frame()
data_Median$Modifications_full = names(data_Median_list)

y = transform(x, Modifications_full = reorder(Modifications_full, order))

plot = ggplot(y, aes(y = dCCS, x = Modifications_full)) + 
  geom_hline(yintercept = 0, size = 0.3, color = "grey") + 
  geom_hline(yintercept = 0.69, linetype = "dashed", size = 0.6, color = "black") + 
  geom_hline(yintercept = -0.69, linetype = "dashed", size = 0.6, color = "black") +
  geom_boxplot(aes(fill = Charge),size = 0.3, outlier.size = 0.3, show.legend = FALSE) + 
  theme_classic() + 
  xlab('') + 
  ylab(expression(Delta * "CCS (%)")) +
  scale_x_discrete(guide = guide_axis(angle = 40))+ 
  scale_color_manual(values=white_black) +
  scale_fill_manual (values=white_black)+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12)) + 
  geom_text(data = data_Median, aes(x = Modifications_full, y= 16, label = round(data_Median, digits = 1)), position = position_dodge(width = 1), size = 3) +
  geom_text(data = one_sample_wilcow, aes(x = Modifications_full, y= 19, label = p.adj.signif), position = position_dodge(width = 1), size = 3)


pdf("Fig.4a_dCCS_boxplot.pdf", width = 8, height = 3.6)
print(plot)
dev.off()


plot = ggplot(y) + geom_bar(mapping = aes(x=Modifications_full, fill = as.character(Charge)), width = 0.7, colour="black", show.legend = FALSE, position = "dodge") + 
  theme_classic() + 
  xlab('') + ylab('N (pairs)') + ggtitle("") +
  scale_x_discrete(guide = guide_axis(angle = 40))+ 
  scale_fill_manual (values=white_black)+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12)) +
  theme(axis.text.x=element_blank(), panel.grid.major.y = element_line(colour = "grey", size = (0.2)))

pdf("Fig.4a_n_peptide_pairs.pdf", width = 8, height = 1.3)
print(plot)
dev.off()

x = split(delta_list, delta_list$Modifications_full)
x = my.data$Ymod_Phospho
x = split(x, x$Charge)
