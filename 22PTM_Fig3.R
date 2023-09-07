library(tidyverse)
library(magrittr)
library(Peptides)
library(data.table)
library(iemisc)
library(VennDiagram)
library(ggrepel)

#Folder of MaxQuant evidence files
WD_Evidence = 'C:/Users/Andre/Data/20220411_21PTM_2/evidences'

WD_Graphs = 'C:/Users/Andre/OneDrive/Desktop/Paper_plots'

#Folder of previously exported files (Cor_fac_CCS_list, Cor_fac_RT_list)
WD_Files = 'C:/Users/Andre/Data/20220411_21PTM_2/evidence/Exported_files'

Colors = rep(c("#941C5A", "#577B80", "#679e1aff", "#00bbd4ff", "#ea8d0bff","#941C5A"), 2)
Colors_black_white = c("white", "grey", "black")

#Parameters for labeling and ordering
Modifications_short = c("Acetyl", "Biotinyl", "Butyryl", "Crotonyl", "Dimethyl", "Formyl", "Glutaryl", "Hydroxyisobutyryl", "Malonyl", "Methyl", "Propionyl", "Succinyl", "Trimethyl", "GlyGly", "Unmodified", "Hydroxyproline", "Unmodified",  "Citrullinyl", "Dimethyl-asym", "Dimethyl-sym", "Methyl", "Unmodified", "O-GlcNAcyl", "Unmodified", "Nitrotyrosine", "Phosphoryl", "Unmodified")
Residue = c(rep("K", 15), rep("P", 2),rep("R", 5),rep("S/T", 2),rep("Y", 3))
Ordering = c(3, 11, 5, 9, 14, 2, 8, 10, 6, 13, 4, 7, 15, 12, 1, 22, 21, 20, 19, 18, 17, 16, 24, 23, 26, 27, 25)

#### Acquire Correction Factors ####

setwd(WD_Files)

# Read R objects from working directory
Cor_fac_K0_list = list.files(pattern = "Cor_fac_K0_list") %>%
  map(readRDS) %>% flatten()
Cor_fac_RT_list = list.files(pattern = "Cor_fac_RT_list") %>%
  map(readRDS) %>% flatten()

#### Acquire Evidence files ####

setwd(WD_Evidence)

my.files.path = list.files(pattern = "evidence")

# Remove unwanted characters from file name for naming results
my.files.identified = str_remove_all(my.files.path, '.txt')
my.files.identified = str_remove_all(my.files.identified, 'evidence_')

cut_names = c('Kmod_', 'Rmod_', 'Ymod_', 'Pmod_', "SpikeMix_")
my.files = str_remove_all(my.files.identified, paste(cut_names, collapse = "|"))


# read .txt files to data frames
my.data = lapply(my.files.path, 
                 read.csv,
                 header=TRUE, sep="\t")

names(my.data) = my.files
my.data.identified = my.data
names(my.data.identified) = my.files.identified

# get number of evidence entries
x = lapply(my.data, function(a) select(a, c(Sequence, Modifications, Intensity, Score)))
x = rbindlist(x)

test = my.data.identified$Kmod_unmod
test2 = split(test, test$Missed.cleavages)

#### Clean data #####

my.data.clean = list()

for (i in 1:(length(my.data.identified))){
  
  data = my.data.identified [[i]]
  name = names(my.data.identified) [[i]]
  
  #if all reverse rows were "na", an error would occur -> change "na" in Reverse to "0"
  data$Reverse[is.na(data$Reverse)] = '0'  
  
  data = filter(data, Potential.contaminant != "+") %>%
    filter(Reverse != "+")  %>%
    filter(Charge != "1")  %>%
    filter(Proteins != 'sp|PTM_Trainingskit_QC|PTM_Trainingskit_QC') %>%
    filter(!grepl("iRT", Proteins)) %>%
    filter(!grepl("2", Modifications)) %>%
    filter(!grepl("3", Modifications)) %>%
    filter(!grepl("Oxidation", Modifications))
  
  #keep only modified peptides for modified samples
  if (grepl( "unmod", name) != TRUE) {
    data = filter(data, data$Modifications != "Unmodified")
    # Extract position of Modification
    data$mod_length = sub("\\(.*", "", data$Modified.sequence)
    data$mod_position = nchar(data$mod_length)-1
    data$rel_position = data$mod_position/data$Length
    
    # filter out peptides where Modification has been localized at C-terminus
    data = filter(data, rel_position != 1) 
    
  }
  #keep only unmodified peptides for unmodified samples
  else {
    data = filter(data, data$Modifications == "Unmodified")
  }
  
  # unmodified Kmod and Rmod peptide should have at least one missed cleavage
  if (name == "Kmod_unmod"|| name == "Rmod_unmod") {
    data = filter(data, data$Missed.cleavages != 0)
  }
  
  my.data.clean [[name]] = data
}

my.data.clean = my.data.clean[ - which(names(my.data.clean) == "iRT")]


##### B_Number_top_int_peptides ####

top_int_list = list()
for (i in 1:length(my.data.clean)) {
  
  data = my.data.clean [[i]]
  name = names(my.data.clean) [[i]]
  
  data = data %>% group_by(Sequence) %>% top_n(1, Intensity)
  
  data = select(data, c(Sequence, Charge, Intensity, CCS, m.z, Raw.file, Retention.time, X1.K0))
  
  data$Modifications_short = paste("(", Residue [[i]], ") ", Modifications_short [[i]], sep = "")
  data$Residue = Residue [[i]]
  data$Ordering = Ordering [[i]]
  
  top_int_list[[name]] = data
  
}

x = rbindlist(top_int_list)

y = transform(x, Modifications_short = reorder(Modifications_short, Ordering) ) 
y$Charge = y$Charge %>% as.character()

plot = ggplot(y, aes(x = Modifications_short, fill = Charge)) + 
  geom_bar(show.legend = FALSE, color = "black",position = position_dodge(preserve = "single"))  +
  scale_fill_manual(values = Colors_black_white, drop = FALSE) + theme_classic() + 
  scale_x_discrete(guide = guide_axis(angle = 40), labels=labels_ordered$Modifications_short)+
  xlab('') + ylab('n') + #ggtitle("C2_C3") +
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 11), panel.grid.major.y = element_line(colour = "grey", size = (0.2)))

setwd(WD_Graphs)


pdf(paste("Fig3_Number_top_int_peptides.pdf", sep=''), width = 7.5, height = 2)
print(plot)
dev.off()

x = top_int_list %>% rbindlist()
y = split(x, x$Charge)


##### B_Charge distribution in CCS/mz clouds ####

CCS.correction.function = function(a){
  
  output = list()
  
  for (i in 1:(length(a))){
    
    data = a [[i]]
    name = names(a) [[i]]
    
    #Get Identifier from Raw.file name
    data$ID = substr(data$Raw.file, nchar(data$Raw.file) - 3, nchar(data$Raw.file))
    
    data_split = split(data, data$ID)
    
    #Correct CCS values according to specific run ID
    data_corrected = list()
    for (j in 1:length(data_split)) {
      x = data_split[[j]]
      name2 = names(data_split)[[j]]
      x = x %>% group_by(Sequence) %>% top_n(1, Intensity)
      
      #Get correction factor for specific ID
      
      Cor_fac_K0 = Cor_fac_K0_list [[name2]]
      Cor_fac_RT = Cor_fac_RT_list [[name2]]
      
      x$K0_cor = x$X1.K0 - Cor_fac_K0
      
      Charge = x$Charge
      Mass = x$m.z * Charge
      Red_mass = sqrt(1/305*(28+(Mass + Charge * 1.00727647))/(28*(Mass + Charge * 1.00727647)))
      
      x$CCS_cor = 18500 * Charge *x$K0_cor * Red_mass
      
      x$RT_cor = x$Retention.time - Cor_fac_RT
      
      data_corrected [[name2]] = x
      
    }
    tmp = rbindlist(data_corrected)
    
    output [[name]] = tmp
    
  }
  return(output)
}

top_int_corrected = CCS.correction.function(top_int_list)


data_split = rbindlist(top_int_corrected)

data_split = split(data_split, data_split$Residue)

plot_list = list()
for (i in 1:length(data_split)) {
  
  data = data_split [[i]]
  name = names(data_split) [[i]]
  
  split = split(data, data$Modifications_short)
  
  plot_list_2 = list()
  plot_list_3 = list()
  for (j in 1:(length(split)-1)) {
    
    data_2 = split [[j]]
    name_2 = names(split) [[j]]
    
    x = split[[length(split)]]
    
    data_2 = rbind(x, data_2)
    
    data_2$Charge = data_2$Charge %>% as.character()
    
    data_2_mod = filter(data_2, !grepl("Unmodified",  Modifications_short))
    
    #plot modified data
    plot = ggplot(data_2_mod, aes(x = m.z, y = CCS_cor)) + 
      geom_point(alpha=1, shape=21, color="black", fill ="#E00724", stroke = 0.1, size = 1, show.legend = FALSE)+
      xlab('m/z') +
      ylab(expression("CCS (" ~ ring(A) ^ 2 * ")")) +
      theme_classic() +
      ggtitle (name_2) +
      scale_color_manual(values="#E00724") +
      scale_fill_manual (values="#E00724")+
      theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12), plot.title = element_text(size=14)) + 
      scale_x_continuous( limits = c(350,1110)) +
      scale_y_continuous( limits = c(320,690))
    
    plot_list_2 [[name_2]] = plot
    
    data_2_unmod = filter(data_2, grepl("Unmodified",  Modifications_short))
    
    #plot unmodified data (same plot repeated multiple times for convenience)
    plot2 = ggplot(data_2_unmod, aes(x = m.z, y = CCS_cor))  + 
      geom_point(alpha=1, shape=21, color="black", fill ="#4373B2", stroke = 0.1, size = 1, show.legend = FALSE)+
      xlab('m/z') + 
      ylab(expression("CCS (" ~ ring(A) ^ 2 * ")")) +
      theme_classic() +
      ggtitle (name_2) +
      theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 12), plot.title = element_text(size=14))+ 
      scale_x_continuous( limits = c(350,1110)) +
      scale_y_continuous( limits = c(320,690))
    
    plot_list_3 [[name_2]] = plot2
    
  }
  plot_list = c(plot_list, plot_list_2, plot_list_3)
}

setwd(WD_Graphs)

pdf(paste("Fig3_Charge_distribution_Modification.pdf", sep=''), width = 2.75, height = 2.5)
for (i in 1:(length(plot_list))) {
  print(plot_list[[i]])
}
dev.off()

##### B_Relative Charge intensity ####

Int_list = list()
for (i in 1:(length(my.data.clean))){
  
  data = my.data.clean [[i]]
  name = names(my.data.clean) [[i]]
  
  tmp = split(data, data$Charge)
  
  # Calculate summed intensity of each charge
  Intensity_sums = sapply(tmp, function(a) sum(a$Intensity))
  
  # Calculate relative intensity of each charge
  tmp2 = Intensity_sums/sum(Intensity_sums)*100
  
  # For some modifications no C4 exists 
  if (length(tmp2) ==2) {
    tmp2 [[3]] = 0
  }
  
  Int_list [[name]] = tmp2
}

z = Int_list %>% as.data.frame()
y = stack(z)
y$ID = rep(c(2:4)) %>% as.character()


# Create information for sorting of data
OrderingX3 = vector()
Modifications_fullX3 = vector()
ModificationsX3 = vector()

for (i in 1:length(Ordering)) {
  tmp = rep(Ordering [[i]],3)
  OrderingX3 = c(OrderingX3, tmp)
  
  tmp = paste("(", Residue [[i]], ") ", Modifications_short [[i]], sep = "")
  tmp = rep(tmp,3)
  Modifications_fullX3 = c(Modifications_fullX3, tmp)
  
  tmp2 = rep(Modifications_short[[i]], 3)
  ModificationsX3 = c(ModificationsX3, tmp2)
}

y$OrderingX3 = OrderingX3
y$Modifications_fullX3 = Modifications_fullX3
y$ModificationsX3 = ModificationsX3
y$Charge = rep(c("2","3","4"))

# Sort data
n = transform(y, Modifications_fullX3 = reorder(Modifications_fullX3, OrderingX3) )

# Create labels for plot
df = cbind(Ordering, Modifications_short) %>% as.data.frame()
df$Ordering = as.numeric(df$Ordering)
labels_ordered = df[order(df$Ordering, decreasing = FALSE), ]

plot = ggplot(n, aes(x = Modifications_fullX3, y = values, fill = Charge)) + 
  geom_col(position = "dodge", show.legend = FALSE, color = "black")  +
  scale_fill_manual(values = Colors_black_white) + theme_classic() + 
  xlab('') + ylab('Rel. Int. (%)') +
  scale_x_discrete(guide = guide_axis(angle = 40), labels=labels_ordered$Modifications_short)+ 
  theme(text = element_text(size=14, color = "black"), axis.text = element_text(colour = "black", size = 11), panel.grid.major.y = element_line(colour = "grey", size = (0.2)))

setwd(WD_Graphs)

pdf(paste("Fig3_Charge_distribution_BW.pdf", sep=''), width = 7.5, height = 2)
print(plot)
dev.off()

