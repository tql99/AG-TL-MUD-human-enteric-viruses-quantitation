# This script is written to perform the equivalent of SAS's PROG GLM (General Linear Models, not to be confused with Generalized Linear Models).
# It can also give estimate pairwise difference.
# Also, it can perform a PCA and create corresponding biplots.
# In addition, it is capable of creating coefficient matrices and melting them down to a 3-column format that can be processed by other applications.
# The datasets used in this script are from Audrey's qPCR runs of samples from her project.
# Area(s) of improvement: Less code repetition (especially GLM & PDIFF).


install.packages("sasLM") # Contains the GLM function which performs the core function of SAS PROC GLM. Also contains PDIFF for pairwise difference.
install.packages("readxl") # Imports Excel files (xls/xlsx) into R w/o external dependencies.
install.packages("dplyr") # Data cleaning.
install.packages("rstudioapi") # Needed to set working directory using "getActiveDocumentContext()" function.
install_github("vqv/ggbiplot") # Visualize PCA. If can't install package, load library(devtools) first.
install.packages("Hmisc") # Needed for (Spearman's) correlation analyses.
install.packages("reshape2") # Melts correlation matrix to visualize it.
install.packages("report") # Contains funtion cite_packages() that lists every used package's citation.

library(sasLM)
library(readxl)
library(dplyr)
library(rstudioapi)
library(devtools) 
library(ggbiplot)
library(Hmisc)
library(reshape2)
library(report)


setwd(dirname(getActiveDocumentContext()$path)) # IMPORTANT - Set working directory to whatever folder that's containing this script.


# Some aesthetics presets:
boxplot_theme <- theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 3.5), #hjust = 0.5 to center horizontally.
                       plot.margin = margin(1, 1, 1, 1, "cm"),
                       panel.background = element_rect(fill = 'white'),
                       panel.grid.major = element_line(color = 'grey95'),
                       panel.grid.major.x = element_blank(),
                       axis.title.x = element_text(size = 14, vjust = -1),
                       axis.title.y = element_text(size = 14, vjust = 2.5),
                       axis.text.x = element_text(size = 10),
                       axis.text.y = element_text(size = 10),
                       axis.line = element_line(color = 'black'),
                       legend.title = element_text(size = 14, face = "bold"),
                       legend.text = element_text(size = 10.5),
                       legend.key = element_blank())

# Custom function(s) to reduce code repetition:
input.cleaning <- function(x){
  # A function to simplify the initial data cleaning process.
  x <- subset(x, select = c('Sample', 'Time', 'GC per mL or g sample', 'log(GC per mL or g sample)', 
                       'GC per ng DNA', 'log(GC per ng DNA)',
                       'GC per ng RNA', 'log(GC per ng RNA)'))
  rename(x, c('sample' = 'Sample',
              'time' = 'Time', 
              'GCpermLorg' = 'GC per mL or g sample', 
              'logGCpermLorg' = 'log(GC per mL or g sample)',
              'GCperngDNA' = 'GC per ng DNA',
              'logGCperngDNA' = 'log(GC per ng DNA)',
              'GCperngRNA' = 'GC per ng RNA',
              'logGCperngRNA' = 'log(GC per ng RNA)')) 
}

# Information for references:
citation() # Info on how to cite R.
RStudio.Version() # Contain "copy-and-paste" reference entry for RStudio.
citation("sasLM") # How to cite package in quotation marks. Apply to other packages as needed.
cite_packages() # Get references for all packages used.
# NOTE for reshape2: Info in "cite_package()" is for reshape. Do manual searches for reshape2 info.
sessionInfo() # Get version number for all packages if cite_packages() didn't help. Just reference the packages in cite_packages(), though.
################################################################################

# Adeno:
Adeno <- read_excel('Excel_analyses/Adeno.xls') %>%
  input.cleaning()
dir.create('R_analyses_output/Adeno') # Create folder to contain output files.

## Effluents:
Adeno_eff <- Adeno %>%
  subset(sample == 'Effluents')

Adeno_eff_GCpermLorg_GLM <- GLM(logGCpermLorg~time,Adeno_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_eff_GCpermLorg_GLM, file = 'R_analyses_output/Adeno/Adeno_eff_GCpermLorg_GLM.txt') # Export results.
Adeno_eff_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,Adeno_eff) %>% # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_eff_GCpermLorg_PDIFF, file = 'R_analyses_output/Adeno/Adeno_eff_GCpermLorg_PDIFF.txt') # Export results.

Adeno_eff_GCperngDNA_GLM <- GLM(logGCperngDNA~time,Adeno_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_eff_GCperngDNA_GLM, file = 'R_analyses_output/Adeno/Adeno_eff_GCperngDNA_GLM.txt') # Export results.
Adeno_eff_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,Adeno_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_eff_GCperngDNA_PDIFF, file = 'R_analyses_output/Adeno/Adeno_eff_GCperngDNA_PDIFF.txt') # Export results.

Adeno_eff_GCperngRNA_GLM <- GLM(logGCperngRNA~time,Adeno_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_eff_GCperngRNA_GLM, file = 'R_analyses_output/Adeno/Adeno_eff_GCperngRNA_GLM.txt') # Export results.
Adeno_eff_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,Adeno_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_eff_GCperngRNA_PDIFF, file = 'R_analyses_output/Adeno/Adeno_eff_GCperngRNA_PDIFF.txt') # Export results.



## Activated Sludge:
Adeno_as <- Adeno %>%
  subset(sample == 'Activated Sludge')

Adeno_as_GCpermLorg_GLM <- GLM(logGCpermLorg~time,Adeno_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_as_GCpermLorg_GLM, file = 'R_analyses_output/Adeno/Adeno_as_GCpermLorg_GLM.txt') # Export results.
Adeno_as_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,Adeno_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_as_GCpermLorg_PDIFF, file = 'R_analyses_output/Adeno/Adeno_as_GCpermLorg_PDIFF.txt') # Export results.

Adeno_as_GCperngDNA_GLM <- GLM(logGCperngDNA~time,Adeno_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_as_GCperngDNA_GLM, file = 'R_analyses_output/Adeno/Adeno_as_GCperngDNA_GLM.txt') # Export results.
Adeno_as_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,Adeno_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_as_GCperngDNA_PDIFF, file = 'R_analyses_output/Adeno/Adeno_as_GCperngDNA_PDIFF.txt') # Export results.

Adeno_as_GCperngRNA_GLM <- GLM(logGCperngRNA~time,Adeno_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_as_GCperngRNA_GLM, file = 'R_analyses_output/Adeno/Adeno_as_GCperngRNA_GLM.txt') # Export results.
Adeno_as_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,Adeno_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_as_GCperngRNA_PDIFF, file = 'R_analyses_output/Adeno/Adeno_as_GCperngRNA_PDIFF.txt') # Export results.


## Raw Sewage:
Adeno_rs <- Adeno %>%
  subset(sample == 'Raw Sewage')

Adeno_rs_GCpermLorg_GLM <- GLM(logGCpermLorg~time,Adeno_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_rs_GCpermLorg_GLM, file = 'R_analyses_output/Adeno/Adeno_rs_GCpermLorg_GLM.txt') # Export results.
Adeno_rs_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,Adeno_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_rs_GCpermLorg_PDIFF, file = 'R_analyses_output/Adeno/Adeno_rs_GCpermLorg_PDIFF.txt') # Export results.

Adeno_rs_GCperngDNA_GLM <- GLM(logGCperngDNA~time,Adeno_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_rs_GCperngDNA_GLM, file = 'R_analyses_output/Adeno/Adeno_rs_GCperngDNA_GLM.txt') # Export results.
Adeno_rs_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,Adeno_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_rs_GCperngDNA_PDIFF, file = 'R_analyses_output/Adeno/Adeno_rs_GCperngDNA_PDIFF.txt') # Export results.

Adeno_rs_GCperngRNA_GLM <- GLM(logGCperngRNA~time,Adeno_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_rs_GCperngRNA_GLM, file = 'R_analyses_output/Adeno/Adeno_rs_GCperngRNA_GLM.txt') # Export results.
Adeno_rs_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,Adeno_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_rs_GCperngRNA_PDIFF, file = 'R_analyses_output/Adeno/Adeno_rs_GCperngRNA_PDIFF.txt') # Export results.


## Sludge Cake:
Adeno_sc <- Adeno %>%
  subset(sample == 'Sludge Cake')

Adeno_sc_GCpermLorg_GLM <- GLM(logGCpermLorg~time,Adeno_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_sc_GCpermLorg_GLM, file = 'R_analyses_output/Adeno/Adeno_sc_GCpermLorg_GLM.txt') # Export results.
Adeno_sc_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,Adeno_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_sc_GCpermLorg_PDIFF, file = 'R_analyses_output/Adeno/Adeno_sc_GCpermLorg_PDIFF.txt') # Export results.

Adeno_sc_GCperngDNA_GLM <- GLM(logGCperngDNA~time,Adeno_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_sc_GCperngDNA_GLM, file = 'R_analyses_output/Adeno/Adeno_sc_GCperngDNA_GLM.txt') # Export results.
Adeno_sc_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,Adeno_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_sc_GCperngDNA_PDIFF, file = 'R_analyses_output/Adeno/Adeno_sc_GCperngDNA_PDIFF.txt') # Export results.

Adeno_sc_GCperngRNA_GLM <- GLM(logGCperngRNA~time,Adeno_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_sc_GCperngRNA_GLM, file = 'R_analyses_output/Adeno/Adeno_sc_GCperngRNA_GLM.txt') # Export results.
Adeno_sc_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,Adeno_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_sc_GCperngRNA_PDIFF, file = 'R_analyses_output/Adeno/Adeno_sc_GCperngRNA_PDIFF.txt') # Export results.


## GLM & PDIFF of treatment averages:
Adeno_GCpermLorg_treatment_GLM <- GLM(logGCpermLorg~sample,Adeno,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_GCpermLorg_treatment_GLM, file = 'R_analyses_output/Adeno/Adeno_GCpermLorg_treatment_GLM.txt') # Export results.
Adeno_GCpermLorg_treatment_PDIFF <- PDIFF(logGCpermLorg~sample,Adeno) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_GCpermLorg_treatment_PDIFF, file = 'R_analyses_output/Adeno/Adeno_GCpermLorg_treatment_PDIFF.txt') # Export results.

Adeno_GCperngDNA_treatment_GLM <- GLM(logGCperngDNA~sample,Adeno,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_GCperngDNA_treatment_GLM, file = 'R_analyses_output/Adeno/Adeno_GCperngDNA_treatment_GLM.txt') # Export results.
Adeno_GCperngDNA_treatment_PDIFF <- PDIFF(logGCperngDNA~sample,Adeno) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_GCperngDNA_treatment_PDIFF, file = 'R_analyses_output/Adeno/Adeno_GCperngDNA_treatment_PDIFF.txt') # Export results.

Adeno_GCperngRNA_treatment_GLM <- GLM(logGCperngRNA~sample,Adeno,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Adeno_GCperngRNA_treatment_GLM, file = 'R_analyses_output/Adeno/Adeno_GCperngRNA_treatment_GLM.txt') # Export results.
Adeno_GCperngRNA_treatment_PDIFF <- PDIFF(logGCperngRNA~sample,Adeno) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Adeno_GCperngRNA_treatment_PDIFF, file = 'R_analyses_output/Adeno/Adeno_GCperngRNA_treatment_PDIFF.txt') # Export results.


# ## Potential - Boxplots:
# qPCR_Adeno_GCpermL$site <- factor(qPCR_Adeno_GCpermL$site,levels = c(1:11,"Negative control")) # Reorder sites from 1-10-11-2-3-4...-Neg to 1-2-3...-10-11-Neg on x-axis.
# pdf("R_analyses_output/qPCR_Adeno/qPCR_Adeno_GCpermL_boxplot.pdf", width = 13, height = 8) # Create the pdf output file.
# qPCR_Adeno_GCpermL_boxplot <- ggplot(qPCR_Adeno_GCpermL, aes(x = site, y = GCpermL, fill = site)) + 
#   geom_boxplot() +
#   geom_jitter(color = 'black', size = 1.5, width = 0.2, height = 0) +
#   scale_fill_discrete(limits = c(1:11,"Negative control")) + # Reorder sites in legend.
#   labs(title = "GC per mL sample by site, Adeno", x = "Site", y = "GC per mL sample", fill = "Site") +
#   boxplot_theme
# print(qPCR_Adeno_GCpermL_boxplot)
# dev.off() # Tells R you're done adding content to the pdf file.
# 
# qPCR_Adeno_GCperngDNA$site <- factor(qPCR_Adeno_GCperngDNA$site,levels = c(1:11,"Negative control")) # Reorder sites from 1-10-11-2-3-4...-Neg to 1-2-3...-10-11-Neg on x-axis.
# pdf("R_analyses_output/qPCR_Adeno/qPCR_Adeno_GCperngDNA_boxplot.pdf", width = 13, height = 8) # Create the pdf output file.
# qPCR_Adeno_GCperngDNA_boxplot <- ggplot(qPCR_Adeno_GCperngDNA, aes(x = site, y = logGCperngDNA, fill = site)) + 
#   geom_boxplot() +
#   geom_jitter(color = 'black', size = 1.5, width = 0.2, height = 0) +
#   scale_fill_discrete(limits = c(1:11,"Negative control")) + # Reorder sites in legend.
#   labs(title = "Normalized GC per ng by site, Adeno", x = "Site", y = "log(GC per ng)", fill = "Site") +
#   boxplot_theme
# print(qPCR_Adeno_GCperngDNA_boxplot)
# dev.off() # Tells R you're done adding content to the pdf file.



# CrAss:
CrAss <- read_excel('Excel_analyses/CrAss.xls') %>%
  input.cleaning()
dir.create('R_analyses_output/CrAss') # Create folder to contain output files.

## Effluents:
CrAss_eff <- CrAss %>%
  subset(sample == 'Effluents')

CrAss_eff_GCpermLorg_GLM <- GLM(logGCpermLorg~time,CrAss_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_eff_GCpermLorg_GLM, file = 'R_analyses_output/CrAss/CrAss_eff_GCpermLorg_GLM.txt') # Export results.
CrAss_eff_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,CrAss_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_eff_GCpermLorg_PDIFF, file = 'R_analyses_output/CrAss/CrAss_eff_GCpermLorg_PDIFF.txt') # Export results.

CrAss_eff_GCperngDNA_GLM <- GLM(logGCperngDNA~time,CrAss_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_eff_GCperngDNA_GLM, file = 'R_analyses_output/CrAss/CrAss_eff_GCperngDNA_GLM.txt') # Export results.
CrAss_eff_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,CrAss_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_eff_GCperngDNA_PDIFF, file = 'R_analyses_output/CrAss/CrAss_eff_GCperngDNA_PDIFF.txt') # Export results.

CrAss_eff_GCperngRNA_GLM <- GLM(logGCperngRNA~time,CrAss_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_eff_GCperngRNA_GLM, file = 'R_analyses_output/CrAss/CrAss_eff_GCperngRNA_GLM.txt') # Export results.
CrAss_eff_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,CrAss_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_eff_GCperngRNA_PDIFF, file = 'R_analyses_output/CrAss/CrAss_eff_GCperngRNA_PDIFF.txt') # Export results.


## Activated Sludge:
CrAss_as <- CrAss %>%
  subset(sample == 'Activated Sludge')

CrAss_as_GCpermLorg_GLM <- GLM(logGCpermLorg~time,CrAss_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_as_GCpermLorg_GLM, file = 'R_analyses_output/CrAss/CrAss_as_GCpermLorg_GLM.txt') # Export results.
CrAss_as_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,CrAss_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_as_GCpermLorg_PDIFF, file = 'R_analyses_output/CrAss/CrAss_as_GCpermLorg_PDIFF.txt') # Export results.

CrAss_as_GCperngDNA_GLM <- GLM(logGCperngDNA~time,CrAss_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_as_GCperngDNA_GLM, file = 'R_analyses_output/CrAss/CrAss_as_GCperngDNA_GLM.txt') # Export results.
CrAss_as_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,CrAss_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_as_GCperngDNA_PDIFF, file = 'R_analyses_output/CrAss/CrAss_as_GCperngDNA_PDIFF.txt') # Export results.

CrAss_as_GCperngRNA_GLM <- GLM(logGCperngRNA~time,CrAss_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_as_GCperngRNA_GLM, file = 'R_analyses_output/CrAss/CrAss_as_GCperngRNA_GLM.txt') # Export results.
CrAss_as_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,CrAss_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_as_GCperngRNA_PDIFF, file = 'R_analyses_output/CrAss/CrAss_as_GCperngRNA_PDIFF.txt') # Export results.


## Raw Sewage:
CrAss_rs <- CrAss %>%
  subset(sample == 'Raw Sewage')

CrAss_rs_GCpermLorg_GLM <- GLM(logGCpermLorg~time,CrAss_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_rs_GCpermLorg_GLM, file = 'R_analyses_output/CrAss/CrAss_rs_GCpermLorg_GLM.txt') # Export results.
CrAss_rs_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,CrAss_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_rs_GCpermLorg_PDIFF, file = 'R_analyses_output/CrAss/CrAss_rs_GCpermLorg_PDIFF.txt') # Export results.

CrAss_rs_GCperngDNA_GLM <- GLM(logGCperngDNA~time,CrAss_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_rs_GCperngDNA_GLM, file = 'R_analyses_output/CrAss/CrAss_rs_GCperngDNA_GLM.txt') # Export results.
CrAss_rs_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,CrAss_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_rs_GCperngDNA_PDIFF, file = 'R_analyses_output/CrAss/CrAss_rs_GCperngDNA_PDIFF.txt') # Export results.

CrAss_rs_GCperngRNA_GLM <- GLM(logGCperngRNA~time,CrAss_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_rs_GCperngRNA_GLM, file = 'R_analyses_output/CrAss/CrAss_rs_GCperngRNA_GLM.txt') # Export results.
CrAss_rs_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,CrAss_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_rs_GCperngRNA_PDIFF, file = 'R_analyses_output/CrAss/CrAss_rs_GCperngRNA_PDIFF.txt') # Export results.


## Sludge Cake:
CrAss_sc <- CrAss %>%
  subset(sample == 'Sludge Cake')

CrAss_sc_GCpermLorg_GLM <- GLM(logGCpermLorg~time,CrAss_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_sc_GCpermLorg_GLM, file = 'R_analyses_output/CrAss/CrAss_sc_GCpermLorg_GLM.txt') # Export results.
CrAss_sc_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,CrAss_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_sc_GCpermLorg_PDIFF, file = 'R_analyses_output/CrAss/CrAss_sc_GCpermLorg_PDIFF.txt') # Export results.

CrAss_sc_GCperngDNA_GLM <- GLM(logGCperngDNA~time,CrAss_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_sc_GCperngDNA_GLM, file = 'R_analyses_output/CrAss/CrAss_sc_GCperngDNA_GLM.txt') # Export results.
CrAss_sc_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,CrAss_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_sc_GCperngDNA_PDIFF, file = 'R_analyses_output/CrAss/CrAss_sc_GCperngDNA_PDIFF.txt') # Export results.

CrAss_sc_GCperngRNA_GLM <- GLM(logGCperngRNA~time,CrAss_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_sc_GCperngRNA_GLM, file = 'R_analyses_output/CrAss/CrAss_sc_GCperngRNA_GLM.txt') # Export results.
CrAss_sc_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,CrAss_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_sc_GCperngRNA_PDIFF, file = 'R_analyses_output/CrAss/CrAss_sc_GCperngRNA_PDIFF.txt') # Export results.


## GLM & PDIFF of treatment averages:
CrAss_GCpermLorg_treatment_GLM <- GLM(logGCpermLorg~sample,CrAss,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_GCpermLorg_treatment_GLM, file = 'R_analyses_output/CrAss/CrAss_GCpermLorg_treatment_GLM.txt') # Export results.
CrAss_GCpermLorg_treatment_PDIFF <- PDIFF(logGCpermLorg~sample,CrAss) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_GCpermLorg_treatment_PDIFF, file = 'R_analyses_output/CrAss/CrAss_GCpermLorg_treatment_PDIFF.txt') # Export results.

CrAss_GCperngDNA_treatment_GLM <- GLM(logGCperngDNA~sample,CrAss,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_GCperngDNA_treatment_GLM, file = 'R_analyses_output/CrAss/CrAss_GCperngDNA_treatment_GLM.txt') # Export results.
CrAss_GCperngDNA_treatment_PDIFF <- PDIFF(logGCperngDNA~sample,CrAss) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_GCperngDNA_treatment_PDIFF, file = 'R_analyses_output/CrAss/CrAss_GCperngDNA_treatment_PDIFF.txt') # Export results.

CrAss_GCperngRNA_treatment_GLM <- GLM(logGCperngRNA~sample,CrAss,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(CrAss_GCperngRNA_treatment_GLM, file = 'R_analyses_output/CrAss/CrAss_GCperngRNA_treatment_GLM.txt') # Export results.
CrAss_GCperngRNA_treatment_PDIFF <- PDIFF(logGCperngRNA~sample,CrAss) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(CrAss_GCperngRNA_treatment_PDIFF, file = 'R_analyses_output/CrAss/CrAss_GCperngRNA_treatment_PDIFF.txt') # Export results.



# GI:
GI <- read_excel('Excel_analyses/GI.xls') %>%
  input.cleaning()
dir.create('R_analyses_output/GI') # Create folder to contain output files.

## Effluents:
GI_eff <- GI %>%
  subset(sample == 'Effluents')

GI_eff_GCpermLorg_GLM <- GLM(logGCpermLorg~time,GI_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_eff_GCpermLorg_GLM, file = 'R_analyses_output/GI/GI_eff_GCpermLorg_GLM.txt') # Export results.
GI_eff_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,GI_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_eff_GCpermLorg_PDIFF, file = 'R_analyses_output/GI/GI_eff_GCpermLorg_PDIFF.txt') # Export results.

GI_eff_GCperngDNA_GLM <- GLM(logGCperngDNA~time,GI_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_eff_GCperngDNA_GLM, file = 'R_analyses_output/GI/GI_eff_GCperngDNA_GLM.txt') # Export results.
GI_eff_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,GI_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_eff_GCperngDNA_PDIFF, file = 'R_analyses_output/GI/GI_eff_GCperngDNA_PDIFF.txt') # Export results.

GI_eff_GCperngRNA_GLM <- GLM(logGCperngRNA~time,GI_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_eff_GCperngRNA_GLM, file = 'R_analyses_output/GI/GI_eff_GCperngRNA_GLM.txt') # Export results.
GI_eff_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,GI_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_eff_GCperngRNA_PDIFF, file = 'R_analyses_output/GI/GI_eff_GCperngRNA_PDIFF.txt') # Export results.


## Activated Sludge:
GI_as <- GI %>%
  subset(sample == 'Activated Sludge')

GI_as_GCpermLorg_GLM <- GLM(logGCpermLorg~time,GI_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_as_GCpermLorg_GLM, file = 'R_analyses_output/GI/GI_as_GCpermLorg_GLM.txt') # Export results.
GI_as_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,GI_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_as_GCpermLorg_PDIFF, file = 'R_analyses_output/GI/GI_as_GCpermLorg_PDIFF.txt') # Export results.

GI_as_GCperngDNA_GLM <- GLM(logGCperngDNA~time,GI_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_as_GCperngDNA_GLM, file = 'R_analyses_output/GI/GI_as_GCperngDNA_GLM.txt') # Export results.
GI_as_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,GI_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_as_GCperngDNA_PDIFF, file = 'R_analyses_output/GI/GI_as_GCperngDNA_PDIFF.txt') # Export results.

GI_as_GCperngRNA_GLM <- GLM(logGCperngRNA~time,GI_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_as_GCperngRNA_GLM, file = 'R_analyses_output/GI/GI_as_GCperngRNA_GLM.txt') # Export results.
GI_as_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,GI_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_as_GCperngRNA_PDIFF, file = 'R_analyses_output/GI/GI_as_GCperngRNA_PDIFF.txt') # Export results.


## Raw Sewage:
GI_rs <- GI %>%
  subset(sample == 'Raw Sewage')

GI_rs_GCpermLorg_GLM <- GLM(logGCpermLorg~time,GI_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_rs_GCpermLorg_GLM, file = 'R_analyses_output/GI/GI_rs_GCpermLorg_GLM.txt') # Export results.
GI_rs_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,GI_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_rs_GCpermLorg_PDIFF, file = 'R_analyses_output/GI/GI_rs_GCpermLorg_PDIFF.txt') # Export results.

GI_rs_GCperngDNA_GLM <- GLM(logGCperngDNA~time,GI_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_rs_GCperngDNA_GLM, file = 'R_analyses_output/GI/GI_rs_GCperngDNA_GLM.txt') # Export results.
GI_rs_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,GI_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_rs_GCperngDNA_PDIFF, file = 'R_analyses_output/GI/GI_rs_GCperngDNA_PDIFF.txt') # Export results.

GI_rs_GCperngRNA_GLM <- GLM(logGCperngRNA~time,GI_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_rs_GCperngRNA_GLM, file = 'R_analyses_output/GI/GI_rs_GCperngRNA_GLM.txt') # Export results.
GI_rs_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,GI_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_rs_GCperngRNA_PDIFF, file = 'R_analyses_output/GI/GI_rs_GCperngRNA_PDIFF.txt') # Export results.


## Sludge Cake:
GI_sc <- GI %>%
  subset(sample == 'Sludge Cake')

GI_sc_GCpermLorg_GLM <- GLM(logGCpermLorg~time,GI_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_sc_GCpermLorg_GLM, file = 'R_analyses_output/GI/GI_sc_GCpermLorg_GLM.txt') # Export results.
GI_sc_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,GI_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_sc_GCpermLorg_PDIFF, file = 'R_analyses_output/GI/GI_sc_GCpermLorg_PDIFF.txt') # Export results.

GI_sc_GCperngDNA_GLM <- GLM(logGCperngDNA~time,GI_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_sc_GCperngDNA_GLM, file = 'R_analyses_output/GI/GI_sc_GCperngDNA_GLM.txt') # Export results.
GI_sc_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,GI_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_sc_GCperngDNA_PDIFF, file = 'R_analyses_output/GI/GI_sc_GCperngDNA_PDIFF.txt') # Export results.

GI_sc_GCperngRNA_GLM <- GLM(logGCperngRNA~time,GI_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_sc_GCperngRNA_GLM, file = 'R_analyses_output/GI/GI_sc_GCperngRNA_GLM.txt') # Export results.
GI_sc_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,GI_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_sc_GCperngRNA_PDIFF, file = 'R_analyses_output/GI/GI_sc_GCperngRNA_PDIFF.txt') # Export results.


## GLM & PDIFF of treatment averages:
GI_GCpermLorg_treatment_GLM <- GLM(logGCpermLorg~sample,GI,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_GCpermLorg_treatment_GLM, file = 'R_analyses_output/GI/GI_GCpermLorg_treatment_GLM.txt') # Export results.
GI_GCpermLorg_treatment_PDIFF <- PDIFF(logGCpermLorg~sample,GI) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_GCpermLorg_treatment_PDIFF, file = 'R_analyses_output/GI/GI_GCpermLorg_treatment_PDIFF.txt') # Export results.

GI_GCperngDNA_treatment_GLM <- GLM(logGCperngDNA~sample,GI,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_GCperngDNA_treatment_GLM, file = 'R_analyses_output/GI/GI_GCperngDNA_treatment_GLM.txt') # Export results.
GI_GCperngDNA_treatment_PDIFF <- PDIFF(logGCperngDNA~sample,GI) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_GCperngDNA_treatment_PDIFF, file = 'R_analyses_output/GI/GI_GCperngDNA_treatment_PDIFF.txt') # Export results.

GI_GCperngRNA_treatment_GLM <- GLM(logGCperngRNA~sample,GI,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GI_GCperngRNA_treatment_GLM, file = 'R_analyses_output/GI/GI_GCperngRNA_treatment_GLM.txt') # Export results.
GI_GCperngRNA_treatment_PDIFF <- PDIFF(logGCperngRNA~sample,GI) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GI_GCperngRNA_treatment_PDIFF, file = 'R_analyses_output/GI/GI_GCperngRNA_treatment_PDIFF.txt') # Export results.



# GII:
GII <- read_excel('Excel_analyses/GII.xls') %>%
  input.cleaning()
dir.create('R_analyses_output/GII') # Create folder to contain output files.

## Effluents:
GII_eff <- GII %>%
  subset(sample == 'Effluents') 

GII_eff_GCpermLorg_GLM <- GLM(logGCpermLorg~time,GII_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_eff_GCpermLorg_GLM, file = 'R_analyses_output/GII/GII_eff_GCpermLorg_GLM.txt') # Export results.
GII_eff_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,GII_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_eff_GCpermLorg_PDIFF, file = 'R_analyses_output/GII/GII_eff_GCpermLorg_PDIFF.txt') # Export results.

GII_eff_GCperngDNA_GLM <- GLM(logGCperngDNA~time,GII_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_eff_GCperngDNA_GLM, file = 'R_analyses_output/GII/GII_eff_GCperngDNA_GLM.txt') # Export results.
GII_eff_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,GII_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_eff_GCperngDNA_PDIFF, file = 'R_analyses_output/GII/GII_eff_GCperngDNA_PDIFF.txt') # Export results.

GII_eff_GCperngRNA_GLM <- GLM(logGCperngRNA~time,GII_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_eff_GCperngRNA_GLM, file = 'R_analyses_output/GII/GII_eff_GCperngRNA_GLM.txt') # Export results.
GII_eff_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,GII_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_eff_GCperngRNA_PDIFF, file = 'R_analyses_output/GII/GII_eff_GCperngRNA_PDIFF.txt') # Export results.


## Activated Sludge:
GII_as <- GII %>%
  subset(sample == 'Activated Sludge')

GII_as_GCpermLorg_GLM <- GLM(logGCpermLorg~time,GII_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_as_GCpermLorg_GLM, file = 'R_analyses_output/GII/GII_as_GCpermLorg_GLM.txt') # Export results.
GII_as_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,GII_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_as_GCpermLorg_PDIFF, file = 'R_analyses_output/GII/GII_as_GCpermLorg_PDIFF.txt') # Export results.

GII_as_GCperngDNA_GLM <- GLM(logGCperngDNA~time,GII_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_as_GCperngDNA_GLM, file = 'R_analyses_output/GII/GII_as_GCperngDNA_GLM.txt') # Export results.
GII_as_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,GII_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_as_GCperngDNA_PDIFF, file = 'R_analyses_output/GII/GII_as_GCperngDNA_PDIFF.txt') # Export results.

GII_as_GCperngRNA_GLM <- GLM(logGCperngRNA~time,GII_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_as_GCperngRNA_GLM, file = 'R_analyses_output/GII/GII_as_GCperngRNA_GLM.txt') # Export results.
GII_as_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,GII_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_as_GCperngRNA_PDIFF, file = 'R_analyses_output/GII/GII_as_GCperngRNA_PDIFF.txt') # Export results.


## Raw Sewage:
GII_rs <- GII %>%
  subset(sample == 'Raw Sewage') 

GII_rs_GCpermLorg_GLM <- GLM(logGCpermLorg~time,GII_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_rs_GCpermLorg_GLM, file = 'R_analyses_output/GII/GII_rs_GCpermLorg_GLM.txt') # Export results.
GII_rs_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,GII_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_rs_GCpermLorg_PDIFF, file = 'R_analyses_output/GII/GII_rs_GCpermLorg_PDIFF.txt') # Export results.

GII_rs_GCperngDNA_GLM <- GLM(logGCperngDNA~time,GII_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_rs_GCperngDNA_GLM, file = 'R_analyses_output/GII/GII_rs_GCperngDNA_GLM.txt') # Export results.
GII_rs_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,GII_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_rs_GCperngDNA_PDIFF, file = 'R_analyses_output/GII/GII_rs_GCperngDNA_PDIFF.txt') # Export results.

GII_rs_GCperngRNA_GLM <- GLM(logGCperngRNA~time,GII_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_rs_GCperngRNA_GLM, file = 'R_analyses_output/GII/GII_rs_GCperngRNA_GLM.txt') # Export results.
GII_rs_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,GII_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_rs_GCperngRNA_PDIFF, file = 'R_analyses_output/GII/GII_rs_GCperngRNA_PDIFF.txt') # Export results.


## Sludge Cake:
GII_sc <- GII %>%
  subset(sample == 'Sludge Cake') 

GII_sc_GCpermLorg_GLM <- GLM(logGCpermLorg~time,GII_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_sc_GCpermLorg_GLM, file = 'R_analyses_output/GII/GII_sc_GCpermLorg_GLM.txt') # Export results.
GII_sc_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,GII_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_sc_GCpermLorg_PDIFF, file = 'R_analyses_output/GII/GII_sc_GCpermLorg_PDIFF.txt') # Export results.

GII_sc_GCperngDNA_GLM <- GLM(logGCperngDNA~time,GII_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_sc_GCperngDNA_GLM, file = 'R_analyses_output/GII/GII_sc_GCperngDNA_GLM.txt') # Export results.
GII_sc_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,GII_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_sc_GCperngDNA_PDIFF, file = 'R_analyses_output/GII/GII_sc_GCperngDNA_PDIFF.txt') # Export results.

GII_sc_GCperngRNA_GLM <- GLM(logGCperngRNA~time,GII_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_sc_GCperngRNA_GLM, file = 'R_analyses_output/GII/GII_sc_GCperngRNA_GLM.txt') # Export results.
GII_sc_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,GII_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_sc_GCperngRNA_PDIFF, file = 'R_analyses_output/GII/GII_sc_GCperngRNA_PDIFF.txt') # Export results.


## GLM & PDIFF of treatment averages:
GII_GCpermLorg_treatment_GLM <- GLM(logGCpermLorg~sample,GII,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_GCpermLorg_treatment_GLM, file = 'R_analyses_output/GII/GII_GCpermLorg_treatment_GLM.txt') # Export results.
GII_GCpermLorg_treatment_PDIFF <- PDIFF(logGCpermLorg~sample,GII) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_GCpermLorg_treatment_PDIFF, file = 'R_analyses_output/GII/GII_GCpermLorg_treatment_PDIFF.txt') # Export results.

GII_GCperngDNA_treatment_GLM <- GLM(logGCperngDNA~sample,GII,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_GCperngDNA_treatment_GLM, file = 'R_analyses_output/GII/GII_GCperngDNA_treatment_GLM.txt') # Export results.
GII_GCperngDNA_treatment_PDIFF <- PDIFF(logGCperngDNA~sample,GII) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_GCperngDNA_treatment_PDIFF, file = 'R_analyses_output/GII/GII_GCperngDNA_treatment_PDIFF.txt') # Export results.

GII_GCperngRNA_treatment_GLM <- GLM(logGCperngRNA~sample,GII,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(GII_GCperngRNA_treatment_GLM, file = 'R_analyses_output/GII/GII_GCperngRNA_treatment_GLM.txt') # Export results.
GII_GCperngRNA_treatment_PDIFF <- PDIFF(logGCperngRNA~sample,GII) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(GII_GCperngRNA_treatment_PDIFF, file = 'R_analyses_output/GII/GII_GCperngRNA_treatment_PDIFF.txt') # Export results.



# PMMV:
PMMV <- read_excel('Excel_analyses/PMMV.xls') %>%
  input.cleaning()
dir.create('R_analyses_output/PMMV') # Create folder to contain output files.

## Effluents:
PMMV_eff <- PMMV %>%
  subset(sample == 'Effluents')

PMMV_eff_GCpermLorg_GLM <- GLM(logGCpermLorg~time,PMMV_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_eff_GCpermLorg_GLM, file = 'R_analyses_output/PMMV/PMMV_eff_GCpermLorg_GLM.txt') # Export results.
PMMV_eff_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,PMMV_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_eff_GCpermLorg_PDIFF, file = 'R_analyses_output/PMMV/PMMV_eff_GCpermLorg_PDIFF.txt') # Export results.

PMMV_eff_GCperngDNA_GLM <- GLM(logGCperngDNA~time,PMMV_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_eff_GCperngDNA_GLM, file = 'R_analyses_output/PMMV/PMMV_eff_GCperngDNA_GLM.txt') # Export results.
PMMV_eff_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,PMMV_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_eff_GCperngDNA_PDIFF, file = 'R_analyses_output/PMMV/PMMV_eff_GCperngDNA_PDIFF.txt') # Export results.

PMMV_eff_GCperngRNA_GLM <- GLM(logGCperngRNA~time,PMMV_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_eff_GCperngRNA_GLM, file = 'R_analyses_output/PMMV/PMMV_eff_GCperngRNA_GLM.txt') # Export results.
PMMV_eff_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,PMMV_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_eff_GCperngRNA_PDIFF, file = 'R_analyses_output/PMMV/PMMV_eff_GCperngRNA_PDIFF.txt') # Export results.


## Activated Sludge:
PMMV_as <- PMMV %>%
  subset(sample == 'Activated Sludge') 

PMMV_as_GCpermLorg_GLM <- GLM(logGCpermLorg~time,PMMV_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_as_GCpermLorg_GLM, file = 'R_analyses_output/PMMV/PMMV_as_GCpermLorg_GLM.txt') # Export results.
PMMV_as_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,PMMV_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_as_GCpermLorg_PDIFF, file = 'R_analyses_output/PMMV/PMMV_as_GCpermLorg_PDIFF.txt') # Export results.

PMMV_as_GCperngDNA_GLM <- GLM(logGCperngDNA~time,PMMV_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_as_GCperngDNA_GLM, file = 'R_analyses_output/PMMV/PMMV_as_GCperngDNA_GLM.txt') # Export results.
PMMV_as_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,PMMV_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_as_GCperngDNA_PDIFF, file = 'R_analyses_output/PMMV/PMMV_as_GCperngDNA_PDIFF.txt') # Export results.

PMMV_as_GCperngRNA_GLM <- GLM(logGCperngRNA~time,PMMV_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_as_GCperngRNA_GLM, file = 'R_analyses_output/PMMV/PMMV_as_GCperngRNA_GLM.txt') # Export results.
PMMV_as_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,PMMV_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_as_GCperngRNA_PDIFF, file = 'R_analyses_output/PMMV/PMMV_as_GCperngRNA_PDIFF.txt') # Export results.


## Raw Sewage:
PMMV_rs <- PMMV %>%
  subset(sample == 'Raw Sewage')

PMMV_rs_GCpermLorg_GLM <- GLM(logGCpermLorg~time,PMMV_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_rs_GCpermLorg_GLM, file = 'R_analyses_output/PMMV/PMMV_rs_GCpermLorg_GLM.txt') # Export results.
PMMV_rs_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,PMMV_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_rs_GCpermLorg_PDIFF, file = 'R_analyses_output/PMMV/PMMV_rs_GCpermLorg_PDIFF.txt') # Export results.

PMMV_rs_GCperngDNA_GLM <- GLM(logGCperngDNA~time,PMMV_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_rs_GCperngDNA_GLM, file = 'R_analyses_output/PMMV/PMMV_rs_GCperngDNA_GLM.txt') # Export results.
PMMV_rs_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,PMMV_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_rs_GCperngDNA_PDIFF, file = 'R_analyses_output/PMMV/PMMV_rs_GCperngDNA_PDIFF.txt') # Export results.

PMMV_rs_GCperngRNA_GLM <- GLM(logGCperngRNA~time,PMMV_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_rs_GCperngRNA_GLM, file = 'R_analyses_output/PMMV/PMMV_rs_GCperngRNA_GLM.txt') # Export results.
PMMV_rs_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,PMMV_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_rs_GCperngRNA_PDIFF, file = 'R_analyses_output/PMMV/PMMV_rs_GCperngRNA_PDIFF.txt') # Export results.


## Sludge Cake:
PMMV_sc <- PMMV %>%
  subset(sample == 'Sludge Cake')

PMMV_sc_GCpermLorg_GLM <- GLM(logGCpermLorg~time,PMMV_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_sc_GCpermLorg_GLM, file = 'R_analyses_output/PMMV/PMMV_sc_GCpermLorg_GLM.txt') # Export results.
PMMV_sc_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,PMMV_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_sc_GCpermLorg_PDIFF, file = 'R_analyses_output/PMMV/PMMV_sc_GCpermLorg_PDIFF.txt') # Export results.

PMMV_sc_GCperngDNA_GLM <- GLM(logGCperngDNA~time,PMMV_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_sc_GCperngDNA_GLM, file = 'R_analyses_output/PMMV/PMMV_sc_GCperngDNA_GLM.txt') # Export results.
PMMV_sc_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,PMMV_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_sc_GCperngDNA_PDIFF, file = 'R_analyses_output/PMMV/PMMV_sc_GCperngDNA_PDIFF.txt') # Export results.

PMMV_sc_GCperngRNA_GLM <- GLM(logGCperngRNA~time,PMMV_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_sc_GCperngRNA_GLM, file = 'R_analyses_output/PMMV/PMMV_sc_GCperngRNA_GLM.txt') # Export results.
PMMV_sc_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,PMMV_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_sc_GCperngRNA_PDIFF, file = 'R_analyses_output/PMMV/PMMV_sc_GCperngRNA_PDIFF.txt') # Export results.


## GLM & PDIFF of treatment averages:
PMMV_GCpermLorg_treatment_GLM <- GLM(logGCpermLorg~sample,PMMV,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_GCpermLorg_treatment_GLM, file = 'R_analyses_output/PMMV/PMMV_GCpermLorg_treatment_GLM.txt') # Export results.
PMMV_GCpermLorg_treatment_PDIFF <- PDIFF(logGCpermLorg~sample,PMMV) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_GCpermLorg_treatment_PDIFF, file = 'R_analyses_output/PMMV/PMMV_GCpermLorg_treatment_PDIFF.txt') # Export results.

PMMV_GCperngDNA_treatment_GLM <- GLM(logGCperngDNA~sample,PMMV,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_GCperngDNA_treatment_GLM, file = 'R_analyses_output/PMMV/PMMV_GCperngDNA_treatment_GLM.txt') # Export results.
PMMV_GCperngDNA_treatment_PDIFF <- PDIFF(logGCperngDNA~sample,PMMV) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_GCperngDNA_treatment_PDIFF, file = 'R_analyses_output/PMMV/PMMV_GCperngDNA_treatment_PDIFF.txt') # Export results.

PMMV_GCperngRNA_treatment_GLM <- GLM(logGCperngRNA~sample,PMMV,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(PMMV_GCperngRNA_treatment_GLM, file = 'R_analyses_output/PMMV/PMMV_GCperngRNA_treatment_GLM.txt') # Export results.
PMMV_GCperngRNA_treatment_PDIFF <- PDIFF(logGCperngRNA~sample,PMMV) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(PMMV_GCperngRNA_treatment_PDIFF, file = 'R_analyses_output/PMMV/PMMV_GCperngRNA_treatment_PDIFF.txt') # Export results.



# Rota:
Rota <- read_excel('Excel_analyses/Rota.xls') %>%
  input.cleaning()
dir.create('R_analyses_output/Rota') # Create folder to contain output files.

## Effluents:
Rota_eff <- Rota %>%
  subset(sample == 'Effluents')

Rota_eff_GCpermLorg_GLM <- GLM(logGCpermLorg~time,Rota_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_eff_GCpermLorg_GLM, file = 'R_analyses_output/Rota/Rota_eff_GCpermLorg_GLM.txt') # Export results.
Rota_eff_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,Rota_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_eff_GCpermLorg_PDIFF, file = 'R_analyses_output/Rota/Rota_eff_GCpermLorg_PDIFF.txt') # Export results.

Rota_eff_GCperngDNA_GLM <- GLM(logGCperngDNA~time,Rota_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_eff_GCperngDNA_GLM, file = 'R_analyses_output/Rota/Rota_eff_GCperngDNA_GLM.txt') # Export results.
Rota_eff_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,Rota_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_eff_GCperngDNA_PDIFF, file = 'R_analyses_output/Rota/Rota_eff_GCperngDNA_PDIFF.txt') # Export results.

Rota_eff_GCperngRNA_GLM <- GLM(logGCperngRNA~time,Rota_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_eff_GCperngRNA_GLM, file = 'R_analyses_output/Rota/Rota_eff_GCperngRNA_GLM.txt') # Export results.
Rota_eff_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,Rota_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_eff_GCperngRNA_PDIFF, file = 'R_analyses_output/Rota/Rota_eff_GCperngRNA_PDIFF.txt') # Export results.


## Activated Sludge:
Rota_as <- Rota %>%
  subset(sample == 'Activated Sludge')

Rota_as_GCpermLorg_GLM <- GLM(logGCpermLorg~time,Rota_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_as_GCpermLorg_GLM, file = 'R_analyses_output/Rota/Rota_as_GCpermLorg_GLM.txt') # Export results.
Rota_as_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,Rota_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_as_GCpermLorg_PDIFF, file = 'R_analyses_output/Rota/Rota_as_GCpermLorg_PDIFF.txt') # Export results.

Rota_as_GCperngDNA_GLM <- GLM(logGCperngDNA~time,Rota_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_as_GCperngDNA_GLM, file = 'R_analyses_output/Rota/Rota_as_GCperngDNA_GLM.txt') # Export results.
Rota_as_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,Rota_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_as_GCperngDNA_PDIFF, file = 'R_analyses_output/Rota/Rota_as_GCperngDNA_PDIFF.txt') # Export results.

Rota_as_GCperngRNA_GLM <- GLM(logGCperngRNA~time,Rota_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_as_GCperngRNA_GLM, file = 'R_analyses_output/Rota/Rota_as_GCperngRNA_GLM.txt') # Export results.
Rota_as_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,Rota_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_as_GCperngRNA_PDIFF, file = 'R_analyses_output/Rota/Rota_as_GCperngRNA_PDIFF.txt') # Export results.


## Raw Sewage:
Rota_rs <- Rota %>%
  subset(sample == 'Raw Sewage')

Rota_rs_GCpermLorg_GLM <- GLM(logGCpermLorg~time,Rota_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_rs_GCpermLorg_GLM, file = 'R_analyses_output/Rota/Rota_rs_GCpermLorg_GLM.txt') # Export results.
Rota_rs_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,Rota_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_rs_GCpermLorg_PDIFF, file = 'R_analyses_output/Rota/Rota_rs_GCpermLorg_PDIFF.txt') # Export results.

Rota_rs_GCperngDNA_GLM <- GLM(logGCperngDNA~time,Rota_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_rs_GCperngDNA_GLM, file = 'R_analyses_output/Rota/Rota_rs_GCperngDNA_GLM.txt') # Export results.
Rota_rs_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,Rota_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_rs_GCperngDNA_PDIFF, file = 'R_analyses_output/Rota/Rota_rs_GCperngDNA_PDIFF.txt') # Export results.

Rota_rs_GCperngRNA_GLM <- GLM(logGCperngRNA~time,Rota_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_rs_GCperngRNA_GLM, file = 'R_analyses_output/Rota/Rota_rs_GCperngRNA_GLM.txt') # Export results.
Rota_rs_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,Rota_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_rs_GCperngRNA_PDIFF, file = 'R_analyses_output/Rota/Rota_rs_GCperngRNA_PDIFF.txt') # Export results.


## Sludge Cake:
Rota_sc <- Rota %>%
  subset(sample == 'Sludge Cake')

Rota_sc_GCpermLorg_GLM <- GLM(logGCpermLorg~time,Rota_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_sc_GCpermLorg_GLM, file = 'R_analyses_output/Rota/Rota_sc_GCpermLorg_GLM.txt') # Export results.
Rota_sc_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,Rota_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_sc_GCpermLorg_PDIFF, file = 'R_analyses_output/Rota/Rota_sc_GCpermLorg_PDIFF.txt') # Export results.

Rota_sc_GCperngDNA_GLM <- GLM(logGCperngDNA~time,Rota_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_sc_GCperngDNA_GLM, file = 'R_analyses_output/Rota/Rota_sc_GCperngDNA_GLM.txt') # Export results.
Rota_sc_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,Rota_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_sc_GCperngDNA_PDIFF, file = 'R_analyses_output/Rota/Rota_sc_GCperngDNA_PDIFF.txt') # Export results.

Rota_sc_GCperngRNA_GLM <- GLM(logGCperngRNA~time,Rota_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_sc_GCperngRNA_GLM, file = 'R_analyses_output/Rota/Rota_sc_GCperngRNA_GLM.txt') # Export results.
Rota_sc_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,Rota_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_sc_GCperngRNA_PDIFF, file = 'R_analyses_output/Rota/Rota_sc_GCperngRNA_PDIFF.txt') # Export results.


## GLM & PDIFF of treatment averages:
Rota_GCpermLorg_treatment_GLM <- GLM(logGCpermLorg~sample,Rota,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_GCpermLorg_treatment_GLM, file = 'R_analyses_output/Rota/Rota_GCpermLorg_treatment_GLM.txt') # Export results.
Rota_GCpermLorg_treatment_PDIFF <- PDIFF(logGCpermLorg~sample,Rota) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_GCpermLorg_treatment_PDIFF, file = 'R_analyses_output/Rota/Rota_GCpermLorg_treatment_PDIFF.txt') # Export results.

Rota_GCperngDNA_treatment_GLM <- GLM(logGCperngDNA~sample,Rota,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_GCperngDNA_treatment_GLM, file = 'R_analyses_output/Rota/Rota_GCperngDNA_treatment_GLM.txt') # Export results.
Rota_GCperngDNA_treatment_PDIFF <- PDIFF(logGCperngDNA~sample,Rota) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_GCperngDNA_treatment_PDIFF, file = 'R_analyses_output/Rota/Rota_GCperngDNA_treatment_PDIFF.txt') # Export results.

Rota_GCperngRNA_treatment_GLM <- GLM(logGCperngRNA~sample,Rota,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(Rota_GCperngRNA_treatment_GLM, file = 'R_analyses_output/Rota/Rota_GCperngRNA_treatment_GLM.txt') # Export results.
Rota_GCperngRNA_treatment_PDIFF <- PDIFF(logGCperngRNA~sample,Rota) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(Rota_GCperngRNA_treatment_PDIFF, file = 'R_analyses_output/Rota/Rota_GCperngRNA_treatment_PDIFF.txt') # Export results.



# uidA:
uidA <- read_excel('Excel_analyses/uidA.xls') %>%
  input.cleaning()
dir.create('R_analyses_output/uidA') # Create folder to contain output files.

## Effluents:
uidA_eff <- uidA %>%
  subset(sample == 'Effluents') 

uidA_eff_GCpermLorg_GLM <- GLM(logGCpermLorg~time,uidA_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_eff_GCpermLorg_GLM, file = 'R_analyses_output/uidA/uidA_eff_GCpermLorg_GLM.txt') # Export results.
uidA_eff_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,uidA_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_eff_GCpermLorg_PDIFF, file = 'R_analyses_output/uidA/uidA_eff_GCpermLorg_PDIFF.txt') # Export results.

uidA_eff_GCperngDNA_GLM <- GLM(logGCperngDNA~time,uidA_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_eff_GCperngDNA_GLM, file = 'R_analyses_output/uidA/uidA_eff_GCperngDNA_GLM.txt') # Export results.
uidA_eff_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,uidA_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_eff_GCperngDNA_PDIFF, file = 'R_analyses_output/uidA/uidA_eff_GCperngDNA_PDIFF.txt') # Export results.

uidA_eff_GCperngRNA_GLM <- GLM(logGCperngRNA~time,uidA_eff,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_eff_GCperngRNA_GLM, file = 'R_analyses_output/uidA/uidA_eff_GCperngRNA_GLM.txt') # Export results.
uidA_eff_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,uidA_eff) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_eff_GCperngRNA_PDIFF, file = 'R_analyses_output/uidA/uidA_eff_GCperngRNA_PDIFF.txt') # Export results.


## Activated Sludge:
uidA_as <- uidA %>%
  subset(sample == 'Activated Sludge') 

uidA_as_GCpermLorg_GLM <- GLM(logGCpermLorg~time,uidA_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_as_GCpermLorg_GLM, file = 'R_analyses_output/uidA/uidA_as_GCpermLorg_GLM.txt') # Export results.
uidA_as_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,uidA_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_as_GCpermLorg_PDIFF, file = 'R_analyses_output/uidA/uidA_as_GCpermLorg_PDIFF.txt') # Export results.

uidA_as_GCperngDNA_GLM <- GLM(logGCperngDNA~time,uidA_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_as_GCperngDNA_GLM, file = 'R_analyses_output/uidA/uidA_as_GCperngDNA_GLM.txt') # Export results.
uidA_as_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,uidA_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_as_GCperngDNA_PDIFF, file = 'R_analyses_output/uidA/uidA_as_GCperngDNA_PDIFF.txt') # Export results.

uidA_as_GCperngRNA_GLM <- GLM(logGCperngRNA~time,uidA_as,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_as_GCperngRNA_GLM, file = 'R_analyses_output/uidA/uidA_as_GCperngRNA_GLM.txt') # Export results.
uidA_as_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,uidA_as) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_as_GCperngRNA_PDIFF, file = 'R_analyses_output/uidA/uidA_as_GCperngRNA_PDIFF.txt') # Export results.


## Raw Sewage:
uidA_rs <- uidA %>%
  subset(sample == 'Raw Sewage') 

uidA_rs_GCpermLorg_GLM <- GLM(logGCpermLorg~time,uidA_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_rs_GCpermLorg_GLM, file = 'R_analyses_output/uidA/uidA_rs_GCpermLorg_GLM.txt') # Export results.
uidA_rs_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,uidA_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_rs_GCpermLorg_PDIFF, file = 'R_analyses_output/uidA/uidA_rs_GCpermLorg_PDIFF.txt') # Export results.

uidA_rs_GCperngDNA_GLM <- GLM(logGCperngDNA~time,uidA_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_rs_GCperngDNA_GLM, file = 'R_analyses_output/uidA/uidA_rs_GCperngDNA_GLM.txt') # Export results.
uidA_rs_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,uidA_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_rs_GCperngDNA_PDIFF, file = 'R_analyses_output/uidA/uidA_rs_GCperngDNA_PDIFF.txt') # Export results.

uidA_rs_GCperngRNA_GLM <- GLM(logGCperngRNA~time,uidA_rs,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_rs_GCperngRNA_GLM, file = 'R_analyses_output/uidA/uidA_rs_GCperngRNA_GLM.txt') # Export results.
uidA_rs_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,uidA_rs) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_rs_GCperngRNA_PDIFF, file = 'R_analyses_output/uidA/uidA_rs_GCperngRNA_PDIFF.txt') # Export results.


## Sludge Cake:
uidA_sc <- uidA %>%
  subset(sample == 'Sludge Cake')

uidA_sc_GCpermLorg_GLM <- GLM(logGCpermLorg~time,uidA_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_sc_GCpermLorg_GLM, file = 'R_analyses_output/uidA/uidA_sc_GCpermLorg_GLM.txt') # Export results.
uidA_sc_GCpermLorg_PDIFF <- PDIFF(logGCpermLorg~time,uidA_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_sc_GCpermLorg_PDIFF, file = 'R_analyses_output/uidA/uidA_sc_GCpermLorg_PDIFF.txt') # Export results.

uidA_sc_GCperngDNA_GLM <- GLM(logGCperngDNA~time,uidA_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_sc_GCperngDNA_GLM, file = 'R_analyses_output/uidA/uidA_sc_GCperngDNA_GLM.txt') # Export results.
uidA_sc_GCperngDNA_PDIFF <- PDIFF(logGCperngDNA~time,uidA_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_sc_GCperngDNA_PDIFF, file = 'R_analyses_output/uidA/uidA_sc_GCperngDNA_PDIFF.txt') # Export results.

uidA_sc_GCperngRNA_GLM <- GLM(logGCperngRNA~time,uidA_sc,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_sc_GCperngRNA_GLM, file = 'R_analyses_output/uidA/uidA_sc_GCperngRNA_GLM.txt') # Export results.
uidA_sc_GCperngRNA_PDIFF <- PDIFF(logGCperngRNA~time,uidA_sc) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_sc_GCperngRNA_PDIFF, file = 'R_analyses_output/uidA/uidA_sc_GCperngRNA_PDIFF.txt') # Export results.


## GLM & PDIFF of treatment averages:
uidA_GCpermLorg_treatment_GLM <- GLM(logGCpermLorg~sample,uidA,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_GCpermLorg_treatment_GLM, file = 'R_analyses_output/uidA/uidA_GCpermLorg_treatment_GLM.txt') # Export results.
uidA_GCpermLorg_treatment_PDIFF <- PDIFF(logGCpermLorg~sample,uidA) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_GCpermLorg_treatment_PDIFF, file = 'R_analyses_output/uidA/uidA_GCpermLorg_treatment_PDIFF.txt') # Export results.

uidA_GCperngDNA_treatment_GLM <- GLM(logGCperngDNA~sample,uidA,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_GCperngDNA_treatment_GLM, file = 'R_analyses_output/uidA/uidA_GCperngDNA_treatment_GLM.txt') # Export results.
uidA_GCperngDNA_treatment_PDIFF <- PDIFF(logGCperngDNA~sample,uidA) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_GCperngDNA_treatment_PDIFF, file = 'R_analyses_output/uidA/uidA_GCperngDNA_treatment_PDIFF.txt') # Export results.

uidA_GCperngRNA_treatment_GLM <- GLM(logGCperngRNA~sample,uidA,lsm=T) # GC per mL. Set lsm=T to get Least Square Means table.
capture.output(uidA_GCperngRNA_treatment_GLM, file = 'R_analyses_output/uidA/uidA_GCperngRNA_treatment_GLM.txt') # Export results.
uidA_GCperngRNA_treatment_PDIFF <- PDIFF(logGCperngRNA~sample,uidA) # Estimate pairwise difference. Set argument PLOT=T for diffogram.
capture.output(uidA_GCperngRNA_treatment_PDIFF, file = 'R_analyses_output/uidA/uidA_GCperngRNA_treatment_PDIFF.txt') # Export results.



# PCA:
dir.create('R_analyses_output/PCA') # Create folder to contain output files.
eff_PCA_in <- read_excel('Excel_analyses/eff_PCA_in.xlsx')
eff_PCA_out <- prcomp(eff_PCA_in[,c(2:ncol(eff_PCA_in))], center = T, scale = T)
eff_PCA_out_summary <- summary(eff_PCA_out) # Get importance of components.
capture.output(eff_PCA_out_summary, file = 'R_analyses_output/PCA/eff_PCA_out_summary.txt')


## Biplots:
### PC1 vs PC2:
pdf('R_analyses_output/PCA/eff_PCA_out_biplot_12.pdf', width = 13, height = 8) # Create the pdf output file.
eff_PCA_out_biplot_12 <- ggbiplot(eff_PCA_out, ellipse = T, groups = as.factor(unique(eff_PCA_in$Time))) +
  geom_text(label = eff_PCA_in$Time, color = 'black', size = 4, position = position_jitter(width = .1, height = .1)) +
  geom_point(shape = 21, colour = 'black', fill = 'black')
print(eff_PCA_out_biplot_12)
dev.off() # Tell R you're done adding content to the pdf file.

pdf('R_analyses_output/PCA/eff_PCA_out_biplot_noAxes_12.pdf', width = 13, height = 8) # Create the pdf output file.
eff_PCA_out_biplot_noAxes_12 <- ggbiplot(eff_PCA_out, ellipse = T, groups = as.factor(unique(eff_PCA_in$Time)), var.axes = F) +
  geom_text(label = eff_PCA_in$Time, color = 'black', size = 4, position = position_jitter(width = .1, height = .1)) +
  geom_point(shape = 21, colour = 'black', fill = 'black') 
print(eff_PCA_out_biplot_noAxes_12)
dev.off() # Tell R you're done adding content to the pdf file.

### PC2 vs PC3:
pdf('R_analyses_output/PCA/eff_PCA_out_biplot_23.pdf', width = 13, height = 8) # Create the pdf output file.
eff_PCA_out_biplot_23 <- ggbiplot(eff_PCA_out, choices = c(2,3), ellipse = T, groups = as.factor(unique(eff_PCA_in$Time))) +
  geom_text(label = eff_PCA_in$Time, color = 'black', size = 4, position = position_jitter(width = .1, height = .1)) +
  geom_point(shape = 21, colour = 'black', fill = 'black')
print(eff_PCA_out_biplot_23)
dev.off() # Tell R you're done adding content to the pdf file.

pdf('R_analyses_output/PCA/eff_PCA_out_biplot_noAxes_23.pdf', width = 13, height = 8) # Create the pdf output file.
eff_PCA_out_biplot_noAxes_23 <- ggbiplot(eff_PCA_out, choices = c(2,3), ellipse = T, groups = as.factor(unique(eff_PCA_in$Time)), var.axes = F) +
  geom_text(label = eff_PCA_in$Time, color = 'black', size = 4, position = position_jitter(width = .1, height = .1)) +
  geom_point(shape = 21, colour = 'black', fill = 'black') 
print(eff_PCA_out_biplot_noAxes_23)
dev.off() # Tell R you're done adding content to the pdf file.

### PC1 vs PC3:
pdf('R_analyses_output/PCA/eff_PCA_out_biplot_13.pdf', width = 13, height = 8) # Create the pdf output file.
eff_PCA_out_biplot_13 <- ggbiplot(eff_PCA_out, choices = c(1,3), ellipse = T, groups = as.factor(unique(eff_PCA_in$Time))) +
  geom_text(label = eff_PCA_in$Time, color = 'black', size = 4, position = position_jitter(width = .1, height = .1)) +
  geom_point(shape = 21, colour = 'black', fill = 'black')
print(eff_PCA_out_biplot_13)
dev.off() # Tell R you're done adding content to the pdf file.

pdf('R_analyses_output/PCA/eff_PCA_out_biplot_noAxes_13.pdf', width = 13, height = 8) # Create the pdf output file.
eff_PCA_out_biplot_noAxes_13 <- ggbiplot(eff_PCA_out, choices = c(1,3), ellipse = T, groups = as.factor(unique(eff_PCA_in$Time)), var.axes = F) +
  geom_text(label = eff_PCA_in$Time, color = 'black', size = 4, position = position_jitter(width = .1, height = .1)) +
  geom_point(shape = 21, colour = 'black', fill = 'black') 
print(eff_PCA_out_biplot_noAxes_13)
dev.off() # Tell R you're done adding content to the pdf file.



# Spearman's correlation:
dir.create('R_analyses_output/Spearman') # Create folder to contain output files.
eff_Spearman_in <- eff_PCA_in_noGCperng[,2:ncol(eff_PCA_in_noGCperng)]
eff_Spearman_out <- rcorr(as.matrix(eff_Spearman_in), type = "spearman") # Specify "type" argument for correlation type ("pearson", "spearman").

eff_Spearman_out_coeff <- eff_Spearman_out$r # Extract the correlation coefficients.
eff_Spearman_out_coeff_melted <- melt(eff_Spearman_out_coeff) # Creates an output w/ 3 cols (var1, var2, val) that can be used as input for other data-handling applications, e.g., Tableau.
write.csv(eff_Spearman_out_coeff_melted, "R_analyses_output/Spearman/eff_Spearman_out_coeff_melted.csv", row.names= F)

eff_Spearman_out_P <- eff_Spearman_out$P # Extract p-values.
eff_Spearman_out_P_melted <- melt(eff_Spearman_out_P) # Creates an output w/ 3 cols (var1, var2, val) that can be used as input for other data-handling applications, e.g., Tableau.
write.csv(eff_Spearman_out_P_melted, "R_analyses_output/Spearman/eff_Spearman_out_P_melted.csv", row.names= F)


