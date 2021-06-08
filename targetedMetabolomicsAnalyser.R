#load libraries
library(dplyr)
library(reshape2) #pivoting
library(readxl)
library(ggplot2)
library(stringr)
library(data.table) #newly added
#library('ggpmisc')
#install.packages(sigr)

##################################################################################################################################################
#FUNCTIONS

plotQC_histo <- function(QCdf, binwidth = 0.1) {
  ggplot(QCdf, aes(RSD)) +
    geom_histogram(binwidth) +
    ggtitle(deparse(substitute(QCdf)))
}

meanRTs <- function(retentionTimeDF) {
  retentionTimeDF %>%
    group_by(Precursor.Ion.Name) %>%
    summarize(Min.Start.Time = mean(Min.Start.Time),
              Max.End.Time = mean(Max.End.Time))
}

#returns pmol wide table from nM skyline output
transformTable <- function(df, nMtoPmolFactor = 0.2) {
  #calc conc to numeric
  df$Calculated.Concentration <- as.numeric(df$Calculated.Concentration) 
  #filter for C12 data
  df <- df %>%
    select(c( "Precursor.Ion.Name","File.Name", "Calculated.Concentration", "Isotope.Label.Type")) %>%
    filter(Isotope.Label.Type == "light") %>%
    select(c( "Precursor.Ion.Name","File.Name", "Calculated.Concentration"))
  #calculate absolute pmol from nM
  df$Calculated.Concentration <- nMtoPmolFactor * df$Calculated.Concentration
  #remove .raw from Filename
  df$File.Name <- str_replace_all(df$File.Name, ".raw", "")
  #wide format 
  df <- reshape2::dcast(df, File.Name ~ Precursor.Ion.Name, fun.aggregate = mean)
  #remove blanks
  metabolites_wide_pmol <- df %>%
    filter(!str_detect(df$File.Name, "Blank"))
  return(metabolites_wide_pmol)
}

# clean this up
transformTableArea <- function(df, nMtoPmolFactor = 0.2) {
  #calc conc to numeric
  df$Area <- as.numeric(df$Area) 
  #filter for C12 data
  df <- df %>%
    select(c( "Precursor.Ion.Name","File.Name", "Area", "Isotope.Label.Type")) %>%
    filter(Isotope.Label.Type == "light") %>%
    select(c( "Precursor.Ion.Name","File.Name", "Area"))
  #calculate absolute pmol from nM
  df$Calculated.Concentration <- nMtoPmolFactor * df$Area
  #remove .raw from Filename
  df$File.Name <- str_replace_all(df$File.Name, ".raw", "")
  #wide format 
  df <- reshape2::dcast(df, File.Name ~ Precursor.Ion.Name, fun.aggregate = mean)
  #remove blanks
  metabolites_wide_pmol <- df %>%
    filter(!str_detect(df$File.Name, "Blank"))
  return(metabolites_wide_pmol)
}



#append polarity to the name of the metabolites
add_polarity <- function(df, polarity) {
  names(df) <- paste(names(df), polarity, sep = "")
  names(df)[1] <- "File.Name" #but not to first column name
  return(df)
}

#analyse RSD of poolQCs
analyseQCs <- function(df, QCnameSTR) {
  #select QCs only 
  QCs <- df %>%
    filter(str_detect(df$File.Name, QCnameSTR))
  
  #melt df to aggregate column-wise for RSD
  molten_QC <- reshape2::melt(QCs, id = c("File.Name"), variable.name = "metabolite", value.name = "pmol")
  
  QCsRSD <- molten_QC %>%
    group_by(metabolite) %>% 
    summarize(pmol_mean = mean(pmol, na.rm = TRUE), RSD=sd(pmol, na.rm = TRUE)/abs(mean(pmol, na.rm = TRUE)), N=n()) 
  
  
  
  return(QCsRSD)
}

#filter for metabolites wich are not reproducible in QCs
removeHighRSD_metabolites <- function(df, levelRSD = 0.6, QCnameSTR) {
  #select QCs only 
  QCs <- df %>%
    filter(str_detect(df$File.Name, QCnameSTR))
  
  #melt df to aggregate column-wise for RSD
  molten_QC <- reshape2::melt(QCs, id = c("File.Name"), variable.name = "metabolite", value.name = "pmol")
  
  QCsRSD <- molten_QC %>%
    group_by(metabolite) %>% 
    summarize(pmol_mean = mean(pmol, na.rm = TRUE), RSD=sd(pmol, na.rm = TRUE)/mean(pmol, na.rm = TRUE)) %>%
    filter(RSD < levelRSD & RSD > 0)
  
  cols_to_keep <- names(df) %in% QCsRSD$metabolite
  cols_to_keep <- names(df)[cols_to_keep]
  
  df_cleaned <- subset(df, select = c("File.Name", cols_to_keep))
  return(df_cleaned)
}

#compares mean pmol QC with LOD values
removeBelowLODs <- function(df, QCs_df, LODs_df) {
  #convert to data tables
  QCs_df <- as.data.table(QCs_df)
  LOD_df <- as.data.table(LODs_df)
  
  #join on metabolites
  QC_LOD <- QCs_df[LODs_df, on = .(metabolite)]
  QC_LOD <- na.omit(QC_LOD)
  
  #keep only above LOD metabolites or no LOD determined
  QC_LOD <- QC_LOD[pmol_mean > LOD | is.na(LOD)]
  cols_to_keep <- names(df) %in% QC_LOD$metabolite
  cols_to_keep <- names(df)[cols_to_keep]
  
  df_cleaned <- subset(df, select = c("File.Name", cols_to_keep))
  return(df_cleaned)
}

#order columns alphabetically
orderColumnsByAlphabet <- function(df) {
  df_ordered <- df[, order(colnames(df))]
  return(df_ordered)
}

makeMetaboAnalystCompatible <- function(df) {
  #first column Sample/File.Name & second Group
  columns <- names(df)
  remove <- c("File.Name", "Group")
  columns <- columns[!columns %in% remove]
  #add the rest of table without File.Name and Group
  columns <- c("File.Name", "Group", columns)
  compatibleDf <- df[,columns]
  
  #filter out NA groups
  compatibleDf <- compatibleDf %>%
    filter(Group != "NA" || Group != "x_x")
  
  return(compatibleDf)
}

#replace negative values with NA
replaceNegativeWithNA <- function(df) {
  df[,3:ncol(df)][df[,3:ncol(df)] <= 0] <- NA
  return(df)
}

extractFirstLettersAsLabel <- function(df) {
  #extract group labels from FileName to a vector
  labels <- substring(df$File.Name, 1, 1)
  
  newLabels <- c()
  for (label in labels) {
    if (label == "C") {
      label <- "Cntrl"
    } else if (label == "H") {
      label <- "H2O2"
    } else if (label == "I") {
      label <- "IL1B"
    } else if (label == "q") {
      label <- "poolQC"
    } 
    else {
      label <- "NA"
    }
    newLabels <- c(newLabels, label)
  }
  return(newLabels)
}

removeMetabolites <- function(metabolitesDF, metabolitesToRemove) {
  reducedMetabolitesDF <- metabolitesDF[, !(names(metabolitesDF) %in% metabolitesToRemove)]
  return(reducedMetabolitesDF)
}

extractSpecialLabel <- function(df, parts=3) {
  #extract group labels from FileName to a vector
  labels <- substring(df$File.Name, 1, 1)
  
  #extract group labels from FileName and add to dataframe
  labels <- strsplit(df$File.Name, "_")
  
  newLabels <- c()
  for (label in labels) {
    if (parts == 3) {
      label <- paste(label[1], "_", label[2], "_", label[3], sep = "")
    } else if (parts == 2) {
      label <- paste(label[1], "_", label[2], sep = "")
    }
    newLabels <- c(newLabels, label)
  }
  
  return(newLabels)
}

getNumUniqueMetabolites <- function(df) {
  df <- df[,3:ncol(df)]
  nameMetabolites <- c()
  metabolites <- strsplit(names(df), "_")
  for (metabolite in metabolites) {
    nameMetabolites <- c(nameMetabolites, metabolite)
  }
  uniqueMetabolites <- unique(nameMetabolites)
  print(uniqueMetabolites)
  return(length(uniqueMetabolites))
}

extractWithoutThirdLabel <- function(df) {
  #extract group labels from FileName to a vector
  labels <- substring(df$File.Name, 1, 1)
  
  #extract group labels from FileName and add to dataframe
  labels <- strsplit(df$File.Name, "_")
  
  newLabels <- c()
  for (label in labels) {
    label <- paste(label[1], "_", label[2], "_", label[4], sep = "")
    newLabels <- c(newLabels, label)
  }
  
  return(newLabels)
}


addGroupColumn <- function(pmolDF) {
  newLabels <- extractSpecialLabel(pmolDF)
  pmolDF$Group <- newLabels
  return(pmolDF)
}

normalizeByProtein <- function(pmolDF, proteinDF) {
  pmol_with_protein <- inner_join(pmolDF, proteinDF, by="File.Name")
  pmol_with_protein[,2:ncol(pmol_with_protein)] <-
    pmol_with_protein[,2:ncol(pmol_with_protein)]/pmol_with_protein[,ncol(pmol_with_protein)]
  pmol_PER_protein <- pmol_with_protein
  pmol_PER_protein$protein <- NULL #remove protein column as no longer needed
  return(pmol_PER_protein)
}
