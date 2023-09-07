library(tidyverse)
library(magrittr)
library(data.table)
library(iemisc)
library(VennDiagram)
library(ggrepel)
library(ggnewscale)
library(heatmaply)
library(orca)


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


### Acquire data Matched_list_C2_C3 ####
setwd(WD_Files)

my.files.path = list.files(pattern = "Matched_list_C2_C3")

Shorten = c('Matched_list_C2_C3 ', " .rds")
my.files = str_remove_all(my.files.path, paste(Shorten, collapse = "|"))

# Read R objects from working directory
my.data = list.files(pattern = "Matched_list_C2_C3") %>%
  map(readRDS)

names(my.data) = my.files

my.data = flatten(my.data)
number_of_matched_peptides = rbindlist(my.data) %>% filter(Modifications != "Unmodified")


##### Calculate delta values ####

delta_list = list()
delta_median = list()
for (i in 1:length(my.data)) {
  data = my.data [[i]]
  name = names(my.data) [[i]]
  
  x = Residue[[i]]
  
  data$Charge = data$Charge %>% as.numeric()
  
  data$dCCS_absolute = data[,"CCS_cor_mod"]-data[,"CCS_cor_unmod"]
  data$dM = (data[,"m.z_mod"]-data[,"m.z_unmod"])*data[,"Charge"]
  
  data$Modifications_full = paste( "(", x, ") ", Modifications_full [[i]] ,sep = "")
  data$mod = Residue [[i]]
  data$order = Ordering [[i]]
  
  delta_list[[name]] = data
}

delta = delta_list[!grepl("unmod", names(delta_list))] 

delta_list = rbindlist(delta_list)


##### Export Matched_peptides.tsv ####
setwd(WD_Graphs)

x = delta_list %>% select(-c(mod, order, ID, Modifications_full)) %>% filter(Modifications != "Unmodified")
write.table(x, file='Matched_peptides.tsv', quote=FALSE, row.names = FALSE, col.names = TRUE, sep='\t')


##### Calculate residuals of linear regression (mod/unmod, CCS/mz) ####

residuals_list_CCS = list()
coefficients_list = list()

for (i in 1:length(delta)) {
  data = delta [[i]]
  name = names(delta) [[i]]
  
  split = split(data, data$Charge)
  
  residuals_list = list()
  coefficients = list()
  for (j in 1:2) {
    
    Charge = c("C2", "C3") [[j]]
    matched = split[[j]]
    
    # Create linear regression CCS/mz, extract Residuals (Distance of every point to the regression line)
    tmp = lm(matched$CCS_cor_unmod~matched$m.z_unmod) 
    residual = matched
    residual$residuals = tmp$residuals
    
    # Create linear regression between dCCS and the residuals of CCS/mz
    tmp2 = lm(residual$dCCS_absolute~residual$residuals)
    x = tmp2$coefficients
    
    # Create linear regression CCS/CCS, extract Residuals (Distance of every point to the regression line)
    tmp = lm(matched$CCS_cor_unmod~matched$CCS_cor_mod) 
    residual$residuals_CCS_CCS = tmp$residuals
    residuals_list [[j]] = residual
    
    
    coefficients [[Charge]] = x
    
  }
  
  residuals_list = rbindlist(residuals_list)
  
  residuals_list_CCS [[name]] = residuals_list
  
  x = as.data.frame(coefficients)
  x$ID = row.names(x)
  
  coefficients_list [[name]] = x
  
}

residuals_CCS = rbindlist(residuals_list_CCS)

residuals_CCS$Charge = as.character(residuals_CCS$Charge)


x = rbindlist(coefficients_list, idcol = TRUE)
y = split(x, x$ID)

z = y [[1]]
n = stack(z[,2:3])
n$ID = rep(z$.id)
names(n)[names(n) == 'values'] <- 'Intercept'

z = y [[2]]
m = stack(z[,2:3])
m$ID = rep(z$.id)
names(m)[names(m) == 'values'] <- 'Abline'

z = cbind(n, m)
z = z[,-c(5, 6)]

Modifications_full = c("Acetylation", "Biotinylation", "Butyrylation", "Crotonylation", "Dimethylation", "Formylation", "Glutarylation", "Hydroxyisobutyrylation", "Malonylation", "Methylation", "Propionylation", "Succinylation", "Trimethylation", "GlyGly", "Hydroxyproline",  "Citrullination", "Dimethylation-asym", "Dimethylation-sym", "Methylation",  "O-GlcNAcylation",  "Nitrotyrosine", "Phosphorylation")
Modifications_with_blank = c("Acetylation", "Biotinylation", "Butyrylation", "Crotonylation", "Dimethylation", "Formylation", "Glutarylation", "Hydroxyisobutyrylation", "Malonylation", "Methylation", "Propionylation", "Succinylation", "Trimethylation", "GlyGly", "Hydroxyproline",  "Citrullination", "Dimethylation-asym", "Dimethylation-sym", "Methylation",  "O-GlcNAcylation",  "Nitrotyrosine", "Phosphorylation", rep("", 22))

Residue = c(rep("K", 14), rep("P", 1),rep("R", 4),rep("S/T", 1),rep("Y", 2))

z$Modifications = Modifications_with_blank #rep(Modifications_full, 2)
z$Residue = rep(Residue, 2)
z$Charge = z$ind


##### Fig.5d_slope_overview ####

Ordering2 = c(3, 11, 5, 9, 14, 2, 8, 10, 6, 13, 4, 7, 15, 12,  22, 20, 19, 18, 17,  21, 23, 26)

z_lineplot = z
tmp = rep(z_lineplot$Modifications[1:22], 2)
z_lineplot$Modifications_full = paste("(", z_lineplot$Residue, ") ", tmp, sep = "")
z_lineplot$order = rep(Ordering2,2)

z_lineplot = transform(z_lineplot, Modifications_full = reorder(Modifications_full, order) )


plot = z_lineplot %>%
  ggplot( aes(x=Modifications_full, y=Abline)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_line(aes(group = ind, color = ind), show.legend = FALSE) +
  geom_point(aes(fill = ind), shape=21, color="black", size=4, show.legend = FALSE) +
  scale_color_manual(values  = c("black", "grey")) +
  scale_fill_manual(values  = c("black", "grey")) +
  theme_classic() +
  #scale_y_continuous(trans = "reverse")+
  scale_x_discrete(guide = guide_axis(angle = 40)) +
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

setwd(WD_Graphs)

pdf(paste("Fig.5d_slope_overview.pdf", sep=''), width = 8, height = 3)
print(plot)
dev.off()

z_lineplot$Modifications_short = sapply(z_lineplot$Modifications_full, function(a) gsub("ation$", "", a))
z_lineplot = transform(z_lineplot, Modifications_short = reorder(Modifications_short, order) )

plot = z_lineplot %>%
  ggplot( aes(x=Modifications_short, y=Abline)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_line(aes(group = ind, color = ind), show.legend = FALSE) +
  geom_point(aes(fill = ind), shape=21, color="black", size=4, show.legend = FALSE) +
  scale_color_manual(values  = c("black", "grey")) +
  scale_fill_manual(values  = c("black", "grey")) +
  theme_classic() +
  #scale_y_continuous(trans = "reverse")+
  #scale_x_discrete(guide = guide_axis(angle = 40)) +
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12)) +
  coord_flip()+
  scale_x_discrete(limits=rev)

setwd(WD_Graphs)

pdf(paste("Fig.5d_slope_overview_poster.pdf", sep=''), width = 3.7, height = 5.2)
print(plot)
dev.off()


##### Fig.5c_CCS_vs_mz Gradient ####

setwd(WD_Graphs)

plot_list = list()
for (i in 1:length(delta)) {
  name = names(delta) [[i]]
  data = delta [[i]]
  data$Charge = data$Charge %>% as.character()
  
  plot = ggplot(data=data, aes(x=m.z_unmod, y=CCS_cor_unmod, color = dCCS_absolute, fill=dCCS_absolute))+
    geom_smooth(data, mapping = aes(y = CCS_cor_unmod, x = m.z_unmod, group = Charge, color = Charge),method='lm', formula= y~x, color = "black", size = 0.3, show.legend = FALSE)+
    #geom_point(shape=21, color="black", stroke = 0.2, size=1.5) +
    geom_point(shape=21, color="black", stroke = 0.2, size=2, show.legend = TRUE) +
    scale_fill_gradientn(colours = c(low="#006AF5", mid = "white", high="#E00724"),
                         values = scales::rescale(c(min(data$dCCS_absolute), median(data$dCCS_absolute), max(data$dCCS_absolute)))) +
    guides(fill = guide_colourbar( ticks.colour = "black", frame.colour = "black")) + 
    #scale_fill_gradient2("dCCS", low="#006AF5", mid = "white", high="#E00724") +
    scale_color_gradient2("dCCS", low="#006AF5", mid = "white", high="#E00724") +
    xlab('m/z') + 
    ylab(expression("CCS (" ~ ring(A) ^ 2 * ")")) +
    theme_classic() +
    ggtitle (data$Modifications_full[[1]]) +
    theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12),plot.title = element_text(size=14)) +
    theme(legend.title=element_blank())
  
  
  plot_list [[name]] = plot
  
}

pdf(paste("Fig.5c_CCS_mz_blot_gradient.pdf", sep=''), width = 4, height = 2.5)
for (i in 1:(length(plot_list))) {
  print(plot_list[[i]])
}
dev.off()


##### Fig.5.Supplementary_dCCS_vs_residuals ####

plot_list2 = list()
plot_list3 = list()
for (i in 1:length(residuals_list_CCS)) {
  
  data = residuals_list_CCS [[i]] 
  name = names(residuals_list_CCS) [[i]]
  
  data$Charge = data$Charge %>% as.character()
  selected = data %>% filter(Charge == 3)
  
  x_limit = c(min(data$residuals)-2, max(data$residuals)+2)
  y_limit = c(min(data$dCCS_absolute)-2, max(data$dCCS_absolute)+2)
  
  plot = ggplot(selected, aes(y = dCCS_absolute, x = residuals))  + 
    geom_vline(xintercept=0, linetype="dashed", color = "grey") +
    geom_hline(yintercept=0, linetype="dashed", color = "grey") +
    geom_smooth(method='lm', formula= y~x, color = "black", show.legend = FALSE, size = 0.5) +
    geom_point(aes(fill = dCCS_absolute), shape=21, color="black", stroke = 0.2, size=2, show.legend = TRUE) +
    #ylab('dCCS (A2)') + 
    ylab(expression(Delta * "CCS (" ~ ring(A) ^ 2 * ")")) +
    xlab('Residuals') + 
    theme_classic() +
    scale_y_continuous(limits = y_limit)+
    scale_x_continuous(limits = x_limit)+
    scale_fill_gradientn(colours = c(low="#006AF5", mid = "white", high="#E00724"),
                         values = scales::rescale(c(min(data$dCCS_absolute), median(data$dCCS_absolute), max(data$dCCS_absolute)))) +
    guides(fill = guide_colourbar( ticks.colour = "black", frame.colour = "black")) + 
    ggtitle (data$Modifications_full[[1]]) +
    theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12),plot.title = element_text(size=14)) +
    theme(legend.title=element_blank())
  
  
  plot_list2 [[name]] = plot
  
  # Plot C2 data with axis limits of C3
  selected = data %>% filter(Charge == 2)
  
  plot = ggplot(selected, aes(y = dCCS_absolute, x = residuals))  + 
    geom_vline(xintercept=0, linetype="dashed", color = "grey") +
    geom_hline(yintercept=0, linetype="dashed", color = "grey") +
    geom_smooth(method='lm', formula= y~x, color = "black", show.legend = FALSE, size = 0.5) +
    geom_point(aes(fill = dCCS_absolute), shape=21, color="black", stroke = 0.2, size=2, show.legend = TRUE) +
    ylab(expression(Delta * "CCS (" ~ ring(A) ^ 2 * ")")) +
    xlab('Residuals') + 
    theme_classic() +
    scale_y_continuous(limits = y_limit)+
    scale_x_continuous(limits = x_limit)+
    scale_fill_gradientn(colours = c(low="#006AF5", mid = "white", high="#E00724"),
                         values = scales::rescale(c(min(data$dCCS_absolute), median(data$dCCS_absolute), max(data$dCCS_absolute)))) +
    guides(fill = guide_colourbar( ticks.colour = "black", frame.colour = "black")) + 
    ggtitle (data$Modifications_full[[1]]) +
    theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12),plot.title = element_text(size=14)) +
    theme(legend.title=element_blank())
  
  plot_list3 [[name]] = plot
  
  
}

pdf(paste("Fig.5.Sup_dCCS_vs_residuals_C3.pdf", sep=''), width = 4, height = 2.5)
for (i in 1:(length(plot_list2))) {
  print(plot_list2[[i]])
}
dev.off()

pdf(paste("Fig.5.Sup_dCCS_vs_residuals_C2.pdf", sep=''), width = 4, height = 2.5)
for (i in 1:(length(plot_list2))) {
  print(plot_list3[[i]])
}
dev.off()


##### Fig.5.Supplementary_dCCS_vs_residuals_no_color ####

plot_list2 = list()
plot_list3 = list()
for (i in 1:length(residuals_list_CCS)) {
  
  data = residuals_list_CCS [[i]] 
  name = names(residuals_list_CCS) [[i]]
  
  data$Charge = data$Charge %>% as.character()
  selected = data %>% filter(Charge == 3)
  
  x_limit = c(min(data$residuals)-2, max(data$residuals)+2)
  y_limit = c(min(data$dCCS_absolute)-2, max(data$dCCS_absolute)+2)
  
  plot = ggplot(selected, aes(y = dCCS_absolute, x = residuals))  + 
    geom_vline(xintercept=0, linetype="dashed", color = "grey") +
    geom_hline(yintercept=0, linetype="dashed", color = "grey") +
    geom_smooth(method='lm', formula= y~x, color = "black", show.legend = FALSE, size = 0.5) +
    geom_point(shape=21, color="black", fill = "black", stroke = 0.2, size=1, show.legend = FALSE, alpha = 0.5) +
    ylab(expression(Delta * "CCS (" ~ ring(A) ^ 2 * ")")) +
    xlab('Residuals') + 
    theme_classic() +
    scale_y_continuous(limits = y_limit)+
    scale_x_continuous(limits = x_limit)+
    guides(fill = guide_colourbar( ticks.colour = "black", frame.colour = "black")) + 
    ggtitle (paste(data$Modifications_full[[1]], ", Charge 3+", sep = "")) +
    theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12),plot.title = element_text(size=14)) +
    theme(legend.title=element_blank())
  
  
  plot_list2 [[name]] = plot
  
  # Plot C2 data with axis limits of C3
  selected = data %>% filter(Charge == 2)
  
  plot = ggplot(selected, aes(y = dCCS_absolute, x = residuals))  + 
    geom_vline(xintercept=0, linetype="dashed", color = "grey") +
    geom_hline(yintercept=0, linetype="dashed", color = "grey") +
    geom_smooth(method='lm', formula= y~x, color = "black", show.legend = FALSE, size = 0.5) +
    geom_point(shape=21, color="black", fill = "black", stroke = 0.2, size=1, show.legend = FALSE, alpha = 0.5) +
    ylab(expression(Delta * "CCS (" ~ ring(A) ^ 2 * ")")) +
    xlab('Residuals') + 
    theme_classic() +
    scale_y_continuous(limits = y_limit)+
    scale_x_continuous(limits = x_limit)+
    guides(fill = guide_colourbar( ticks.colour = "black", frame.colour = "black")) + 
    ggtitle (paste(data$Modifications_full[[1]], ", Charge 2+", sep = "")) +
    theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12),plot.title = element_text(size=14)) +
    theme(legend.title=element_blank())
  
  plot_list3 [[name]] = plot
  
  
}



pdf(paste("Fig.5.Sup_dCCS_vs_residuals_C3_black.pdf", sep=''), width = 3, height = 2.5)
for (i in 1:(length(plot_list2))) {
  print(plot_list2[[i]])
}
dev.off()

pdf(paste("Fig.5.Sup_dCCS_vs_residuals_C2_black.pdf", sep=''), width = 3, height = 2.5)
for (i in 1:(length(plot_list2))) {
  print(plot_list3[[i]])
}
dev.off()



#### Fig.5A Line Plot C2 Acylation ####

data_position = delta_list
data_position = data_position %>% filter(Charge ==2)

data_position$mod_length = sub("\\(.*", "", data_position$Modified.sequence)
data_position$mod_position = nchar(data_position$mod_length)-1
data_position$mod_position = data_position$mod_position %>% as.character()

z = data_position
z = transform(z, Modifications_full = reorder(Modifications_full, order) )
n = select(z, c("Sequence", "Modifications_full", "mod_position", "dCCS_absolute"))
m = split(n, n$Modifications_full)
o = m[2:15]
merged = o %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("mod_position", "Sequence")), .) 
merged = merged[complete.cases(merged), ]
merged2 = select(merged, c(1,4,6,8,10,12,14,16,18,20,22,24,26,28,30))
names(merged2) = c("Sequence", "Formyl", "Acetyl", "Propionyl", "Butyryl", "Malonyl", "Succinyl", "Glutaryl", "Crotonyl", "Hydroxyisobutyryl", "Biotinyl", "Glygly", "Methyl", "Dimethyl", "Trimethyl")
merged2 = transform(merged2, Sequence = reorder(Sequence, Formyl) )

library(GGally)

plot = ggparcoord(merged2, scale = "globalminmax",
                  columns = c(2:15), groupColumn = "Sequence",
                  showPoints = TRUE, 
                  alphaLines = 0.5)+
  theme_classic() + 
  scale_x_discrete(guide = guide_axis(angle = 40))+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12)) +
  ylab(expression(Delta * "CCS (" ~ ring(A) ^ 2 * ")")) +
  theme(legend.position="none") 

pdf("Fig5_LinePlot_dCCS.pdf", width = 8, height = 3.5)
print(plot)
dev.off()

pdf("Fig5_LinePlot_dCCS_poster.pdf", width = 6, height = 3.5)
print(plot)
dev.off()


#### Fig.5A Line Plot C3 Acylation, Supplementary ####

data_position = delta_list
data_position = data_position %>% filter(Charge ==3)

data_position$mod_length = sub("\\(.*", "", data_position$Modified.sequence)
data_position$mod_position = nchar(data_position$mod_length)-1
data_position$mod_position = data_position$mod_position %>% as.character()

z = data_position
z = transform(z, Modifications_full = reorder(Modifications_full, order) )
n = select(z, c("Sequence", "Modifications_full", "mod_position", "dCCS_absolute"))
m = split(n, n$Modifications_full)
o = m[2:15]
# Remove Propionyl and Biotinyl, as intersection is low
o = o[names(o) %in% c("(K) Propionyl", "(K) Biotinyl") == FALSE]

merged = o %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("mod_position", "Sequence")), .) 
merged = merged[complete.cases(merged), ]
merged2 = select(merged, c(1,4,6,8,10,12,14,16,18,20,22,24,26))
names(merged2) = c("Sequence", "Formyl", "Acetyl", "Butyryl", "Malonyl", "Succinyl", "Glutaryl", "Crotonyl", "Hydroxyisobutyryl",  "Glygly", "Methyl", "Dimethyl", "Trimethyl")
merged2 = transform(merged2, Sequence = reorder(Sequence, Formyl) )

library(GGally)

plot = ggparcoord(merged2, scale = "globalminmax",
                  columns = c(2:13), groupColumn = "Sequence",
                  showPoints = TRUE, 
                  alphaLines = 0.5)+
  theme_classic() + 
  scale_x_discrete(guide = guide_axis(angle = 40))+
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12)) +
  ylab(expression(Delta * "CCS (" ~ ring(A) ^ 2 * ")")) +
  theme(legend.position="none") 

pdf("Fig5_LinePlot_dCCS_C3.pdf", width = 8, height = 3.5)
print(plot)
dev.off()

#### Fig.5B Heatmap Succinyl  ####

data_position = delta_list
split = split(data_position, data_position$Charge)

C2 = split$`2` %>% filter(Modifications == "Succinyl")
C2$Modifications_full = paste(C2$Modifications_full, " C", C2$Charge, sep = "")
C3 = split$`3` %>% filter(Modifications == "Succinyl")
C3$Modifications_full = paste(C3$Modifications_full, " C", C3$Charge, sep = "")

data_position = rbind(C2, C3)

split = split(data_position, data_position$Modifications_full)
name = c("Sequence", names(split))

# Extract modification site to make sure that charge 2 and 3 data have the same site
data_position$mod_length = sub("\\(.*", "", data_position$Modified.sequence)
data_position$mod_position = nchar(data_position$mod_length)-1
data_position$mod_position = data_position$mod_position %>% as.character()

z = data_position
z = transform(z, Modifications_full = reorder(Modifications_full, order) )
n = select(z, c("Sequence", "Modifications_full", "mod_position", "dCCS_absolute"))
m = split(n, n$Modifications_full)
o = m

# Merge the data frames by Sequence and modification position
merged = o %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("mod_position", "Sequence")), .) 
merged = merged[complete.cases(merged), ]
merged2 = select(merged, c(1,4,6))
names(merged2) = name
merged2 = transform(merged2, Sequence = reorder(Sequence, `(K) Butyryl C2`) )

# Transform data for plotting 
tmp = merged2$Sequence %>% as.character()
merged2 = merged2[,-1] %>% as.matrix()
row.names(merged2) = tmp


p = heatmaply(
  merged2,
  file = "Fig.5B_heatmap_Succinyl_C2_C3.pdf",
  scale_fill_gradient_fun = scale_fill_gradient2(
    low = "#006AF5", 
    high = "#E00724"),
  fontsize_row = 7,
  fontsize_column = 12,
  width = 500, 
  height = 1000
)

#### Fig.5B Heatmap Trimethyl  ####

data_position = delta_list
split = split(data_position, data_position$Charge)

C2 = split$`2` %>% filter(Modifications == "Trimethyl (K)")
C2$Modifications_full = paste(C2$Modifications_full, " C", C2$Charge, sep = "")
C3 = split$`3` %>% filter(Modifications == "Trimethyl (K)")
C3$Modifications_full = paste(C3$Modifications_full, " C", C3$Charge, sep = "")

data_position = rbind(C2, C3)

split = split(data_position, data_position$Modifications_full)
name = c("Sequence", names(split))


data_position$mod_length = sub("\\(.*", "", data_position$Modified.sequence)
data_position$mod_position = nchar(data_position$mod_length)-1
data_position$mod_position = data_position$mod_position %>% as.character()

z = data_position
z = transform(z, Modifications_full = reorder(Modifications_full, order) )
n = select(z, c("Sequence", "Modifications_full", "mod_position", "dCCS_absolute"))
m = split(n, n$Modifications_full)
o = m
merged = o %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("mod_position", "Sequence")), .) 
merged = merged[complete.cases(merged), ]
merged2 = select(merged, c(1,4,6))
names(merged2) = name
merged2 = transform(merged2, Sequence = reorder(Sequence, `(K) Trimethyl C2`) )

tmp = merged2$Sequence %>% as.character()
merged2 = merged2[,-1] %>% as.matrix()
row.names(merged2) = tmp


p = heatmaply(
  merged2,
  file = "Fig.5B_heatmap_Trimethyl_C2_C3.pdf",
  scale_fill_gradient_fun = scale_fill_gradient2(
    low = "#006AF5", 
    high = "#E00724"),
  fontsize_row = 7,
  fontsize_column = 12,
  width = 500, 
  height = 1000
)

##### Calculate R2 of modified and unmodified peptide clouds ####

x = my.data %>% rbindlist() %>% filter(Modifications != "Unmodified")

y = filter(x, x$Charge == 2)

unmod = lm(y$CCS_unmod_cor~y$m.z_unmod)
tmp = summary(unmod)
R2_C2_unmod = tmp$adj.r.squared

mod = lm(y$CCS_mod~y$m.z_mod)
tmp = summary(mod)
R2_C2_mod = tmp$adj.r.squared


y = filter(x, x$Charge == 3)

unmod = lm(y$CCS_unmod~y$m.z_unmod)
tmp = summary(unmod)
R2_C3_unmod = tmp$adj.r.squared

mod = lm(y$CCS_mod~y$m.z_mod)
tmp = summary(mod)
R2_C3_mod = tmp$adj.r.squared

#### Do modifications "normalize" the data? #####

#Charge 3

C2 = delta_list %>% filter(Charge == 3)
split = split(C2, C2$Modifications_full)

Function_unmod = function(x){ Coefficients_unmod[2] *log(x)+ Coefficients_unmod[1] }
Function_mod = function(x){ Coefficients_mod[2] *log(x)+ Coefficients_mod[1] }

list = list()
R2_list_C3 = list()
R2_change_list_C3 = list()

for (i in 1:length(split)) {
  
  name = names(split) [[i]]
  a = split [[i]]
  
  Regr_log_unmod = lm(a$CCS_cor_unmod~log(a$m.z_unmod))
  #Regr_log_unmod = lm(a$CCS_cor_unmod~a$m.z_unmod)
  Coefficients_unmod = Regr_log_unmod$coefficients
  R2_unmod = summary(Regr_log_unmod)$adj.r.squared
  
  Regr_log_mod = lm(a$CCS_cor_mod~log(a$m.z_mod))
  #Regr_log_mod = lm(a$CCS_cor_mod~a$m.z_mod)
  Coefficients_mod = Regr_log_mod$coefficients
  R2_mod = summary(Regr_log_mod)$adj.r.squared
  
  a$Res_unmod = a$CCS_cor_unmod - Function_unmod(a$m.z_unmod)
  a$Res_mod = a$CCS_cor_mod - Function_mod(a$m.z_mod)
  
  list[[name]] = a
  R2_list_C3 [[name]] = c(R2_unmod,R2_mod, a$order[[2]])
  
  R2_change_list_C3 [[name]] = c(R2_mod - R2_unmod, a$order[[2]], "C3")
  
}

# plot violin plot of residues
tmp2 = rbindlist(list)

tmp = rbindlist(list) %>% select(c(Sequence, Modifications_full, Res_mod, Res_unmod, order)) %>%
  pivot_longer(., cols = c(Res_mod,Res_unmod), names_to = "Var", values_to = "Val") %>% 
  arrange(desc(order)) %>%
  mutate(Modifications_full=fct_reorder(Modifications_full,order)) %>%
  filter(!grepl("Unmodified", Modifications_full))

plot = ggplot(tmp, aes(x = Modifications_full, y = Val, fill = fct_rev(Var))) +
  theme_classic() +
  geom_hline(yintercept = 0, size = 0.3, color = "grey") + 
  geom_boxplot(size = 0.3, outlier.size = 0.3, show.legend = FALSE) + 
  scale_fill_manual(values  = c("white", "grey")) +
  xlab("") + 
  ylab(expression("Residual (" ~ ring(A) ^ 2 * ")")) +
  scale_x_discrete(guide = guide_axis(angle = 40)) +
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

setwd(WD_Graphs)

pdf(paste("Residuals_mod_vs_unmod_boxplot_C3.pdf", sep=''), width = 8, height = 3.5)
print(plot)
dev.off()



#Charge 2
C2 = delta_list %>% filter(Charge == 2)
split = split(C2, C2$Modifications_full)

list = list()
R2_list_C2 = list()
R2_change_list_C2 = list()

for (i in 1:length(split)) {
  
  name = names(split) [[i]]
  a = split [[i]]
  
  Regr_log_unmod = lm(a$CCS_cor_unmod~log(a$m.z_unmod))
  #Regr_log_unmod = lm(a$CCS_cor_unmod~a$m.z_unmod)
  Coefficients_unmod = Regr_log_unmod$coefficients
  R2_unmod = summary(Regr_log_unmod)$adj.r.squared
  
  Regr_log_mod = lm(a$CCS_cor_mod~log(a$m.z_mod))
  #Regr_log_mod = lm(a$CCS_cor_mod~a$m.z_mod)
  Coefficients_mod = Regr_log_mod$coefficients
  R2_mod = summary(Regr_log_mod)$adj.r.squared
  
  a$Res_unmod = a$CCS_cor_unmod - Function_unmod(a$m.z_unmod)
  a$Res_mod = a$CCS_cor_mod - Function_mod(a$m.z_mod)
  
  list[[name]] = a
  R2_list_C2 [[name]] = c(R2_unmod,R2_mod, a$order[[2]])
  
  R2_change_list_C2 [[name]] = c(R2_mod - R2_unmod, a$order[[2]], "C2")
  
}

x = R2_change_list_C3 %>% as_tibble() %>% t() %>% as.data.frame()
x$Modifications_full = row.names(y)

y = R2_change_list_C2 %>% as_tibble() %>% t() %>% as.data.frame()
y$Modifications_full = row.names(y)


tmp2 = rbindlist(list)

tmp = rbindlist(list) %>% select(c(Sequence, Modifications_full, Res_mod, Res_unmod, order)) %>%
  pivot_longer(., cols = c(Res_mod,Res_unmod), names_to = "Var", values_to = "Val") %>% 
  arrange(desc(order)) %>%
  mutate(Modifications_full=fct_reorder(Modifications_full,order)) %>%
  filter(!grepl("Unmodified", Modifications_full))

plot = ggplot(tmp, aes(x = Modifications_full, y = Val, fill = fct_rev(Var))) +
  theme_classic() +
  geom_hline(yintercept = 0, size = 0.3, color = "grey") + 
  geom_boxplot(size = 0.3, outlier.size = 0.3, show.legend = FALSE) + 
  scale_fill_manual(values  = c("white", "grey")) +
  xlab("") + 
  ylab(expression("Residual (" ~ ring(A) ^ 2 * ")")) +
  scale_x_discrete(guide = guide_axis(angle = 40)) +
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12))

setwd(WD_Graphs)

pdf(paste("Residuals_mod_vs_unmod_boxplot_C2.pdf", sep=''), width = 8, height = 3.5)
print(plot)
dev.off()



