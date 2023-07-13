library(tidyverse)
library(magrittr)
library(data.table)
library(iemisc)
library(VennDiagram)
library(ggrepel)
library(xlsx)

#Folder with MaxQuant evidence files
WD_Evidence = 'C:/Users/Andre/Data/20220411_21PTM_2/evidences'

#Folder for exported figures
WD_Graphs = 'C:/Users/Andre/OneDrive/Desktop/Paper_plots'

#Folder for exported files
WD_Files = 'C:/Users/Andre/Data/20220411_21PTM_2/evidence/Exported_files'

#color scheme
Colors = c("#941C5A", "#577B80", "#679e1aff", "#00bbd4ff", "#ea8d0bff","#941C5A", "#577B80", "#679e1aff", "#00bbd4ff", "#ea8d0bff","#941C5A")

#### Acquire data ####
setwd(WD_Evidence)

Pattern = "evidence"
cut_names = c('Kmod_', 'Rmod_', 'Ymod_', 'Pmod_', "SpikeMix_")

my.files.path = list.files(pattern = Pattern)

# Remove unwanted characters from file name for naming results
my.files.identified = str_remove_all(my.files.path, '.txt')
my.files.identified = str_remove_all(my.files.identified, 'evidence_')
my.files = str_remove_all(my.files.identified, paste(cut_names, collapse = "|"))


# read .txt files to data frames
my.data = lapply(my.files.path, 
                 read.csv,
                 header=TRUE, sep="\t")
names(my.data) = my.files
my.data.identified = my.data
names(my.data.identified) = my.files.identified

#### Get number of PSMs and unique peptide sequences ####

# Reduce columns to bind modified and unmodified data frames
x = lapply(my.data, function(a) select(a, c(Sequence, Modifications, Intensity, Score, Charge)))
x = rbindlist(x)
# Total number of evidence entries
A_Peptide_spectrum_matches = nrow(x)

# Remove peptides with multiple modifications and Oxidation
split = split(x, x$Modifications)
`%notlike%` <- Negate(`%like%`)
z = split[names(split) %notlike% "2|3|4|Oxidation"]

# Get highest intensity peptides
y = lapply(z, function(x) group_by(x, Sequence) %>% top_n(1, Intensity) %>% ungroup() %>% distinct())

# Count unique peptides (unmodified or singly modified, without Oxidation)
A_Unique_peptide_sequences = y %>% lapply(function(x) nrow(x)) %>% as.data.frame() %>% sum()

x = rbindlist(y)
x = split(x, x$Charge)

#### Prepare data for iRT analysis ####

#loop for cleaning data.
iRT_list = list()

for (i in 1:(length(my.data))){
  
  data = my.data [[i]]
  name = paste(my.files.identified [i])
  
  #if all reverse rows were "na", an error would occur -> change "na" in Reverse to "0"
  data$Reverse[is.na(data$Reverse)] = '0'  
  
  iRT = filter(data, Potential.contaminant != "+") %>%
    filter(Reverse != "+") %>%
    filter(grepl("iRT", Proteins)) %>%
    filter(grepl("Unmod", Modifications))
  iRT_list [[name]] = iRT
  
}

##### Allocate IDs to MS runs, Reduce data frames ####

iRT_all_list = list()

for (i in 1:(length(iRT_list))){
  
  data = iRT_list [[i]]
  name = names(iRT_list) [[i]]
  
  #extract run-number from raw.file name to later sort time dependent
  data$ID = substr(data$Raw.file, nchar(data$Raw.file) - 3, nchar(data$Raw.file))
  
  # Reduce data frame, as modified and unmodified have different number of columns
  data_df = data.frame(Sequence = data$Sequence,
                       CCS = data$CCS,
                       RT = data$Retention.time,
                       Intensity = data$Intensity,
                       Raw.file = data$Raw.file,
                       ID = data$ID,
                       charge = data$Charge,
                       m.z = data$m.z,
                       Modified.sequence = data$Modified.sequence,
                       name = name,
                       X1.K0 = data$X1.K0,
                       Score = data$Score)

  
  #split into technical replicates, get individual top intensity peptides
  iRT_split = split(data_df, data$Raw.file)
  
  
  ##### Take top intensity peptides for each run; Choose correct RT peak for DGLDAASYYAPVR #####
  iRT_split_list = list()
  
  for (j in 1:(length(iRT_split))){
    
    df = iRT_split [[j]]
    ID = df$ID [[2]]
    
    n = iRT_split [[j]] %>% group_by(Sequence) %>% top_n(1, Intensity) %>% 
      ungroup %>% distinct(Sequence, .keep_all = TRUE)

    # DGLDAASYYAPVR has two RT-peaks. For some measurements (ID < 2200), the wrong peak is selected
    # From the two highest intensity entries take the one with lower RT
    
    if(ID <= 2200){
    
    n_DGLDAA = df %>% filter(Sequence == "DGLDAASYYAPVR") %>% top_n(2, Intensity) %>% slice(which.min(RT))
    n = n %>% filter(Sequence != "DGLDAASYYAPVR")
    n = rbind(n_DGLDAA, n)
    }
    
    iRT_split_list [[j]] = n
  }
  
  iRT_all_list [[name]] = rbindlist(iRT_split_list)
  
}

iRT_all = rbindlist(iRT_all_list)  



##### Split distinct iRT peptides of all runs by sequence ####

iRT_seq_split = split(iRT_all, iRT_all$Sequence)

iRT_seq = rbindlist(iRT_seq_split) 

numbering = as.data.frame(iRT_seq$ID)
numbering = distinct(numbering, numbering[1])
colnames(numbering) = "ID"
numbering = as.data.frame(numbering[order(numbering$ID),])
colnames(numbering) = "ID"
numbering$ID2 = c(1:nrow(numbering))

iRT_final = merge(iRT_seq, numbering, by.x = "ID", by.y = "ID")


##### Create correction factor for iRT peptides and apply it to iRT peptides ####

iRT = iRT_all_list$iRT 

split_Sequence = split(iRT, iRT$Sequence)


#Get mean CCS and RT values for all iRT peptides
iCCS_mean_list = vector()
iCCS_CV_list = vector()
iCCS_sd_list = vector()
iRT_mean_list = vector()
iRT_CV_list = vector()
iRT_sd_list = vector()
im.z_mean_list = vector()


for (i in 1:length(split_Sequence)) {
  
  tmp = split_Sequence[[i]]
  name = names(split_Sequence) [[i]]
  
  #CCS_median = median(tmp$CCS)
  CCS_mean = mean(tmp$CCS)
  CCS_sd = sd(tmp$CCS)
  CCS_CV = cv(tmp$CCS)
  RT_mean = mean(tmp$RT)
  RT_sd = sd(tmp$RT)
  RT_CV = cv(tmp$RT)
  m.z_mean = mean(tmp$m.z)
  
  iCCS_mean_list [[name]] = CCS_mean
  iCCS_sd_list [[name]] = CCS_sd
  iCCS_CV_list [[name]] = CCS_CV  
  iRT_mean_list [[name]] = RT_mean
  iRT_sd_list [[name]] = RT_sd
  iRT_CV_list [[name]] = RT_CV
  im.z_mean_list [[name]] = m.z_mean
}

iCCS_mean_df = data.frame(Sequence = names(iCCS_mean_list),
                            m.z = im.z_mean_list,
                            CCS = iCCS_mean_list,
                            CCS_sd = iCCS_sd_list,
                            RT = iRT_mean_list,
                            RT_sd = iRT_sd_list)

all_by_ID = split(iRT_final, iRT_final$ID)

iCCS_final_median_cor_list = list()
iRT_final_median_cor_list = list()
Cor_fac_CCS_list = list()
Cor_fac_RT_list = list()

for (i in 1:length(all_by_ID)) {
  data = all_by_ID [[i]] 
  name = names(all_by_ID) [[i]]
  
  #select only peptides which were identified in a specific run
  iRT = iCCS_mean_df[which(iCCS_mean_df$Sequence %in% data$Sequence),]
  data = data[which(data$Sequence %in% iRT$Sequence),]
  
  #Calculate Correction for CCS
  tmp = data
  tmp$CF = tmp$CCS - iRT$CCS
  tmp$Sequence_iRT = iRT$Sequence
  tmp$CCS_iRT = iRT$CCS
  #substract median Correction factor from CCS
  tmp$CCS_cor = tmp$CCS - median(tmp$CF)

  
  iCCS_final_median_cor_list [[name]] = tmp
  Cor_fac_CCS_list [[name]] = median(tmp$CF)
  
  #Calculate Correction for RT
  tmp = data
  tmp$CF_RT = tmp$RT - iRT$RT
  tmp$Sequence_iRT = iRT$Sequence
  tmp$RT_iRT = iRT$RT
  #substract median Correction factor from RT
  tmp$RT_cor = tmp$RT - median(tmp$CF_RT)
  
  iRT_final_median_cor_list [[name]] = tmp
  Cor_fac_RT_list [[name]] = median(tmp$CF)
  
}

iRT_final_median_cor = rbindlist(iCCS_final_median_cor_list)
iRT_RT_final_median_cor = rbindlist(iRT_final_median_cor_list)
tmp = select(iRT_RT_final_median_cor, c(CF_RT, RT_cor, RT_iRT))

iRT_final_median_cor = cbind(iRT_final_median_cor, tmp)

iRT_final_median_cor_export = select(iRT_final_median_cor, c(Sequence, ID, CCS, CCS_cor, RT, RT_cor, m.z, charge, Raw.file))



##### Export files and Data frames regarding iRT peptides ####
setwd(WD_Files)

# Supl.Table.2
write.xlsx(iCCS_mean_df, "iCCS_mean_df.xlsx", col.names = TRUE, row.names = TRUE, append = FALSE)

# For Fig.3
saveRDS(Cor_fac_CCS_list, paste("Cor_fac_CCS_list.rds"))
saveRDS(Cor_fac_RT_list, paste("Cor_fac_RT_list.rds"))

# For Fig.2A
saveRDS(iRT_final_median_cor, paste("iRT_final_median_cor.rds"))

# Export all used iRT peptides
setwd(WD_Graphs)
write.xlsx(iRT_final_median_cor_export, "iRT_final_median_cor.xlsx", col.names = TRUE, row.names = TRUE, append = FALSE)

#### Median CV of iRT CCS values ####

x = split(iRT_final_median_cor, iRT_final_median_cor$Sequence)

CV_iRT = lapply(x, function(a) cv(a$CCS)) %>% as.data.frame() %>% t() %>% median()
CV_iRT_cor = lapply(x, function(a) cv(a$CCS_cor)) %>% as.data.frame() %>% t() %>% median()

#### Export global iRT CCS deviation ####

x = iCCS_final_median_cor_list %>% rbindlist()
x = split(x, x$Sequence)

y = lapply(x, function(a) as.list(c(mean(a$CCS_cor), sd(a$CCS_cor), cv(a$CCS_cor)))) %>% rbindlist() %>% as.data.frame()
colnames(y) = c("mean", "sd", "cv")
rownames(y) = names(x)

setwd(WD_Graphs)
write.xlsx(y, "iRT_global_deviation.xlsx", col.names = TRUE, row.names = TRUE, append = FALSE)


#### Prepare PTM data ####

#SpikeMix, Rmod, Kmod, Pmod, Ymod, iRT (Partition of peptide controls)

Modification_list = c("SpikeMix", "Rmod", "Kmod", "Pmod", "Ymod")
Residue_list = c("S/T", "R", "K", "P", "Y")

for (t in 1:length(Modification_list)) {
  
  Modification = Modification_list [[t]]
  

PTM_list = list()
names(my.data.identified) = my.files.identified

my.data_selected = my.data.identified[grep(Modification, names(my.data.identified))]


#loop for cleaning data
my.data.mod = list()
my.data.unmod = list()

my.data.mod.C3 = list()
my.data.unmod.C3 = list()

my.data.mod.C4 = list()
my.data.unmod.C4 = list()

for (i in 1:(length(my.data_selected))){
  
  data = my.data_selected [[i]]
  name = names(my.data_selected) [[i]]
  
  data$Residue = Residue_list [[t]]
  
  
  #if all reverse rows were "na", an error would occur -> change "na" in Reverse to "0"
  data$Reverse[is.na(data$Reverse)] = '0'  
  
  data = filter(data, Potential.contaminant != "+") %>%
    filter(Reverse != "+")  %>%
    filter(Proteins != 'sp|PTM_Trainingskit_QC|PTM_Trainingskit_QC') %>%
    filter(!grepl("iRT", Proteins)) %>%
    filter(!grepl("Oxidation", Modifications)) 
  
  # This column interferes with extracting modification probabilities for the specific Modifications
  data = select(data, -c("Oxidation..M..Probabilities"))
  
  unmod = filter(data, Modifications == "Unmodified")
  
  tmp = filter(unmod, Charge == "2")
  my.data.unmod [[name]] = tmp
  
  tmp = filter(unmod, Charge == "3")
  my.data.unmod.C3 [[name]] = tmp
  
  if ("4" %in% unmod$Charge) {
  tmp = filter(unmod, Charge == "4")
  my.data.unmod.C4 [[name]] = tmp
  }
  
  mod = filter(data, Modifications != "Unmodified") %>%
    filter(!grepl("Oxidation", Modifications)) %>%
    filter(!grepl("2", Modifications)) %>%
    filter(!grepl("3", Modifications)) #exclude doubly modified peptides
  
  #Extract highest number out of Probability Sequence
  #Column name of Probabilities contains modification information -> Col 5
  if (grepl( "unmod", name) != TRUE) {
  mod$Probabilities = mod[,5]
  x = sapply(mod$Probabilities, function(x) str_match_all(x, pattern = "(?<=\\().+?(?=\\))"))
  y = lapply(x, function(x) max(as.numeric(x)))
  mod$Probabilities_Extract = y
  }
  
    #mod = filter(mod, Probabilities_Extract >= 0.99)
 
  # Extract position of Modification
  mod$mod_length = sub("\\(.*", "", mod$Modified.sequence)
  mod$mod_position = nchar(mod$mod_length)-1
  mod$rel_position = mod$mod_position/mod$Length
  
  # filter out peptides where Modification has been localized at C-terminus
  mod = filter(mod, rel_position != 1) 
  
  
  tmp = filter(mod, Charge == "2")
  my.data.mod [[name]] = tmp
  
  tmp = filter(mod, Charge == "3")
  my.data.mod.C3 [[name]] = tmp
  
  if ("4" %in% mod$Charge) {
    tmp = filter(mod, Charge == "4")
    my.data.mod.C4 [[name]] = tmp
  }


  
}


my.data.mod = my.data.mod[grep("unmod", names(my.data.mod), invert = TRUE)] 
my.data.unmod = my.data.unmod[grep("unmod", names(my.data.unmod))]

my.data.mod.C3 = my.data.mod.C3[grep("unmod", names(my.data.mod.C3), invert = TRUE)] 
my.data.unmod.C3 = my.data.unmod.C3[grep("unmod", names(my.data.unmod.C3))]

my.data.mod.C4 = my.data.mod.C4[grep("unmod", names(my.data.mod.C4), invert = TRUE)]
my.data.unmod.C4 = my.data.unmod.C4[grep("unmod", names(my.data.unmod.C4))]

###### Correct CCS values with iRT data by ID (run number) ####


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
      x = x %>% group_by(Sequence) %>% top_n(1, Intensity) %>% ungroup() %>% distinct(Sequence, .keep_all = TRUE)
      
      #Get correction factor for specific ID
      
      Cor_fac_CCS = Cor_fac_CCS_list [[name2]]
      Cor_fac_RT = Cor_fac_RT_list [[name2]]
      
      x$CCS_cor = x$CCS - Cor_fac_CCS
      x$RT_cor = x$Retention.time - Cor_fac_RT
      
      data_corrected [[name2]] = x
      
    }
    tmp = rbindlist(data_corrected)
    
    output [[name]] = tmp
    
  }
  return(output)
}


  my.data.mod.cor = CCS.correction.function(my.data.mod)
  my.data.unmod.cor = CCS.correction.function(my.data.unmod)
  
  my.data.mod.cor.C3 = CCS.correction.function(my.data.mod.C3)
  my.data.unmod.cor.C3 = CCS.correction.function(my.data.unmod.C3)
  
  # For some samples there is no Charge 4 modified peptides -> create empty list
  if (length(my.data.mod.C4) != 0) {
    my.data.mod.cor.C4 = CCS.correction.function(my.data.mod.C4)
  } else{my.data.mod.cor.C4 = list()}
  
    my.data.unmod.cor.C4 = CCS.correction.function(my.data.unmod.C4)

  
#### Export top_int_n3 (For suppl. Fig.1 and 2) ####
  
  #unmod and mod have different number of columns -> preselection needed
  select.funtion = function(a){(lapply(a, function(x) select(x, c(Sequence, Modified.sequence, Charge, Score, Intensity, ID, Residue, CCS, m.z, Mass, Retention.time, CCS_cor, RT_cor, Modifications, Probabilities_Extract))))}
  
  x = select.funtion(my.data.mod.cor)
  y = select.funtion(my.data.mod.cor.C3)
  z = select.funtion(my.data.mod.cor.C4)
  
  Export_top_int_n3 = c(x, y, z)
  Export_top_int_unmod_n3 = c(my.data.unmod.cor, my.data.unmod.cor.C3, my.data.unmod.cor.C4)
  
  setwd(WD_Files)

  saveRDS(Export_top_int_n3, paste("Top_int_n3", Modification, ".rds"))
  saveRDS(Export_top_int_unmod_n3, paste("Top_int_n3_unmod", Modification, ".rds"))
  

###### Match 3 technical replicates and calculate CV (For Fig.2B) ####

unmod_name = paste(Modification, "_unmod", sep = "")
  
mod.and.unmod = my.data.mod.cor
mod.and.unmod[[unmod_name]] = my.data.unmod.cor [[1]]

mod.and.unmod.C3 = my.data.mod.cor.C3
mod.and.unmod.C3[[unmod_name]] = my.data.unmod.cor.C3 [[1]]

  
peptide.matching.3rep.function = function(a){
  
  matched_list = list()
  
  for (i in 1:(length(a))){
    
    n = a [[i]]
    name = names(a) [[i]]
    
    split.n = split(n, n$ID)
    
    x = split.n [[1]]
    y = split.n [[2]]
    z = split.n [[3]]
    
    matched_x = x[which(x$Sequence %in% y$Sequence),]
    matched_x = matched_x[which(matched_x$Sequence %in% z$Sequence),]
    
    matched_y = y[which(y$Sequence %in% matched_x$Sequence),]
    matched_z = z[which(z$Sequence %in% matched_x$Sequence),]
    
    # create data frames from matched data
    matched = data.frame(Sequence = matched_x$Sequence,
                         CCS_x = matched_x$CCS,
                         CCS_y = matched_y$CCS, 
                         CCS_z = matched_z$CCS, 
                         CCS_cor_x = matched_x$CCS_cor,
                         CCS_cor_y = matched_y$CCS_cor, 
                         CCS_cor_z = matched_z$CCS_cor,
                         RT_x = matched_x$Retention.time,
                         RT_y = matched_y$Retention.time,
                         RT_z = matched_z$Retention.time,
                         RT_cor_x = matched_x$RT_cor,
                         RT_cor_y = matched_y$RT_cor,
                         RT_cor_z = matched_z$RT_cor,
                         Charge = matched_x$Charge,
                         Modified.sequence = matched_x$Modified.sequence,
                         Modifications = matched_x$Modifications,
                         ID = rep(name))
    
    
    matched$CV_CCS = apply(matched[2:4], MARGIN = 1, FUN = cv)
    matched$CV_CCS_cor = apply(matched[5:7], MARGIN = 1, FUN = cv)
    
    matched$CV_RT = apply(matched[8:10], MARGIN = 1, FUN = cv)
    matched$CV_RT_cor = apply(matched[11:13], MARGIN = 1, FUN = cv)
    
    matched_list [[name]] = matched
    
    }
  return(matched_list) 
  }

CV_replicates_df = peptide.matching.3rep.function(mod.and.unmod)
CV_replicates_df = rbindlist(CV_replicates_df)

CV_replicates_df_C3 = peptide.matching.3rep.function(mod.and.unmod.C3)
CV_replicates_df_C3 = rbindlist(CV_replicates_df_C3)

setwd(WD_Files)

# For Fig.2B
#Includes only peptides matched between all 3 technical replicates  
write.table(CV_replicates_df,paste("CV_replicates_df_", Modification,".txt"),sep="\t",row.names=TRUE)
write.table(CV_replicates_df_C3,paste("CV_replicates_df_C3_", Modification,".txt"),sep="\t",row.names=TRUE)


###### Get Overlap between technical replicates -> Calculate mean ####
       #mean of n3, n2 and also single values are included for matching

Calculate.mean.function = function(a, Column) {
  mean_list = list()
  sd_list = list()
  
  for (j in 1:length(a)) {
    
    tmp = a[[j]]
    tmp =  tmp[, ..Column] [[Column]]
    mean = mean(tmp)
    sd = sd(tmp)
    
    name2 = names(a) [[j]]
    
    mean_list[[name2]] = mean
    sd_list[[name2]] = sd
    
  }
  return(list(mean_list, sd_list))
}
  
  
Calculate.mean.function.total = function(a) {
  
  n_all_list = list()
  n_in2_list = list() 
  n_in1_list = list()
  
  #Calculate mean of peptides overlapping between technical replicates
  mean_list = list()
  for (i in 1:(length(a))){
    n = a [[i]]
    name = names(a) [[i]]
    
    #split into 3 technical replicates
    split.n = split(n, n$ID)
    
    x = split.n [[1]]$Sequence
    y = split.n [[2]]$Sequence
    z = split.n [[3]]$Sequence
    
    #Calculate overlap between technical replicates (Package VennDiagram)
    overlap = calculate.overlap(x = list(
      "x" = x,
      "y" = y,
      "z" = z
    ))
    
    #Clarify/ annotate the overlaps
    names(overlap) <- c("a123", "a12", "a13", "a23", "a1", "a2", "a3")
    
    n_all = length(overlap$a123)
    n_in2 = length(c(overlap$a12, overlap$a13, overlap$a23))
    n_in1 = length(c(overlap$a1, overlap$a2, overlap$a3))
    
    n_all_list [[name]] = n_all
    n_in2_list [[name]] = n_in2
    n_in1_list [[name]] = n_in1
    
    # Ions present in 3 technical replicates
    n3 = n[which(n$Sequence %in% overlap$a123),]
    n3 = split(n3, n3$Sequence)
    
    n2 = n[which(n$Sequence %in% c(overlap$a12, overlap$a13, overlap$a23)),]
    n2 = split(n2, n2$Sequence)
    
    n1 = n[which(n$Sequence %in% c(overlap$a1, overlap$a2, overlap$a3)),]

    
    #loop for calculating mean of different parameters
    Columns = c("CCS", "CCS_cor", "m.z", "Retention.time", "RT_cor")
    
    mean.pre.list = list()
    sd.pre.list = list()
    for (j in 1:length(Columns)) {
      Col = Columns [[j]]
    
    #Ions with n3 overlap
    mean = Calculate.mean.function(n3, Col)  
    mean_3 = mean[[1]]
    sd_3 = mean[[2]]
    
    #Ions with n2 overlap
    mean = Calculate.mean.function(n2, Col)
    mean_2 = mean[[1]]
    sd_2 = mean[[2]]
    
    #Ions detected in only one measurement
    N1 = n1[,..Col] [[1]] %>% as.list()
    names(N1) = n1$Sequence
  
    #put n1, n2 and n3 into one table
    mean = c(mean_2, mean_3, N1) %>% as.data.frame() %>% t() %>% as.data.frame()
    mean$sd = as.numeric(c(sd_2, sd_3, rep(0, length(N1))))
    mean$Sequence = rownames(mean) 
    
    #Get number of overlapping technical replicates, but only once during the loop (in the cas of CCS)
    if (Col == "CCS") {
    mean$N = c(rep(2, length(mean_2)), rep(3, length(mean_3)), rep(1, length(N1)))
    }
    
    mean.pre.list[[Col]] = mean
    
    }
    
    df = mean.pre.list %>% reduce(full_join, by='Sequence', .init)
    names(df) = c("CCS", "sd_CCS", "Sequence", "N", "CCS_cor", "sd_CCS_cor",  "m.z", "sd_m.z", "Retention.time", "sd_Retention.time",  "RT_cor", "sd_RT_cor")

    #not in correct order, yet!
    mean_list [[name]] = df 
    
  }
  
  return.overlap = cbind(n_all_list, n_in2_list, n_in1_list) %>% as.data.frame()
  return(list(return.overlap, mean_list))
  
}

#Returns:
  #Overlap between 3 technical replicates (retur.overlap)
  #mean between technical replicates (mean_list)

x = Calculate.mean.function.total(mod.and.unmod)
mean_list = x[[2]]
return.overlap = x[[1]]

x = Calculate.mean.function.total(mod.and.unmod.C3)
mean_list_C3 = x[[2]]
return.overlap_C3 = x[[1]]


# Bring Charge, Modifications, Modified Sequence back in

add_mod_Seq_Function = function(mean_list, mod.and.unmod) {

tmp_list = list()

for (i in 1:length(mean_list)) {
  
  name = names(mean_list) [[i]]
  x = mean_list [[i]]
  x = x[order(x$Sequence),]
  
  y = mod.and.unmod [[i]]
  y = y[order(y$Sequence),]

  z = y %>% group_by(Sequence) %>% top_n(1, Intensity) %>% ungroup() %>% distinct(Sequence, .keep_all = TRUE)
  
  x$Modified.sequence = z$Modified.sequence
  x$Modifications = z$Modifications
  x$Charge = z$Charge
  x$Score = z$Score
  
  tmp_list [[name]] = x
}
return(tmp_list)
}

mean_list = add_mod_Seq_Function(mean_list, mod.and.unmod)
mean_list_C3 = add_mod_Seq_Function(mean_list_C3, mod.and.unmod.C3)

#### Export mean_list and Overlap ####
x = rbindlist(mean_list, idcol = TRUE)

setwd(WD_Files)

saveRDS(x, paste("Mean_list", Modification, ".rds"))

saveRDS(return.overlap, paste("Overlap", Modification, ".rds"))
saveRDS(return.overlap_C3, paste("Overlap", Modification, "_C3.rds"))


#### Match peptides ####

  peptide.matching.function.2 = function(a, Charge){
    
    matched_list = list()
    
    for (i in 1:(length(a))){
    
      x = a [[i]]
      y = a [[length(a)]]
      
      name = names(a) [[i]]
      
      #match mod/ unmod by sequence, order by Sequence
      matched_x = x[which(x$Sequence %in% y$Sequence),]
      matched_x = matched_x[order(matched_x$Sequence),]
      
      matched_y = y[which(y$Sequence %in% x$Sequence),]
      matched_y = matched_y[order(matched_y$Sequence),]
      
      # create data frames from matched data
      matched = data.frame(Sequence = matched_x$Sequence,
                           N_rep_mod = matched_x$N,
                           N_rep_unmod = matched_y$N,
                           CCS_mod = matched_x$CCS,
                           CCS_mod_sd = matched_x$sd_CCS,
                           CCS_unmod = matched_y$CCS,
                           CCS_unmod_sd = matched_y$sd_CCS,
                           CCS_cor_mod = matched_x$CCS_cor,
                           CCS_cor_mod_sd = matched_x$sd_CCS_cor, 
                           CCS_cor_unmod = matched_y$CCS_cor,
                           CCS_cor_unmod_sd = matched_y$sd_CCS_cor,
                           m.z_mod = matched_x$m.z, 
                           m.z_unmod = matched_y$m.z,
                           RT_mod = matched_x$Retention.time, 
                           RT_unmod = matched_y$Retention.time,
                           RT_cor_mod = matched_x$RT_cor,
                           RT_cor_mod_sd = matched_x$sd_RT_cor, 
                           RT_cor_unmod = matched_y$RT_cor,
                           RT_cor_unmod_sd = matched_y$sd_RT_cor,
                           Charge = Charge,
                           Modified.sequence = matched_x$Modified.sequence,
                           Modifications = matched_x$Modifications,
                           ID = rep(name))
      
      
      matched_list [[name]] = matched
      
    }
    return(matched_list)
  }
  

  matched_list = peptide.matching.function.2(mean_list, Charge = 2)
  matched_list_C3 = peptide.matching.function.2(mean_list_C3, Charge = 3)
  
  matched_list_C2_C3 = list()
  for (i in 1:length(matched_list)) {
    x = matched_list [[i]]
    y = matched_list_C3 [[i]]
    name = names(matched_list) [[i]]
    
    tmp = rbind(x,y)
    
    matched_list_C2_C3 [[name]] = tmp
    
  }
  

#### Export Data Frames regarding matched data (Fig.2B, Fig.4) ####

# For Fig.4
saveRDS(matched_list_C2_C3, paste("Matched_list_C2_C3", Modification, ".rds"))

}
