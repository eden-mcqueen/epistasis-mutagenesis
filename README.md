# epistasis-mutagenesis
Data and code for 2025 Epistasis Mutagenesis paper

General Information

Repository contains files associated with epistasis yeast mutagenesis project in Wittkopp laboratory at the University of Michigan.
Project origination and majority of code was written by former Wittkopp laboratory graduate student Dr. Bing Yang. Additions, modifications, and cuts were made by Dr. E. W. McQueen for the final analysis and figures for publication.


Included files


Python scripts:

Create_template.py
Creates template file (“TEMPLATE.csv”) for importing and analyzing flow data.

R scripts:

Clean.Data.R
Uses “TEMPLATE.csv”
Reading in fcs files and transforming data

Cleaning.Functions.2.R
Functions used in “Clean.Data.R”

ANALYSIS.R
Base scripts for analyzing and scoring mutants from transformed data from “Clean.Data.R”. Further cleaning of the data by removing small count samples, etc.
Calculations for the medians, expression estimates, etc.
This code was run at separate times for different portions of the data, so the output file names include ascending numbers. 
Outputs are the "STRAIN.ESTIMATES.txt" files, below.

Analysis.1.9.R
Script to run final analysis for manuscript figures and tables using processed data files.

Processed data files:
STRAIN.ESTIMATES.1.txt
STRAIN.ESTIMATES.2.txt
STRAIN.ESTIMATES.3.txt
STRAIN.ESTIMATES.4.txt
STRAIN.ESTIMATES.5.txt
STRAIN.ESTIMATES.6NEW.txt
SUMMARY.TRANS.txt
raw.dist.mv.sorted.csv  
raw.dist.mc.sorted.csv
Note: SUMMARY.TRANS.txt uses previously collected data (Metzger et al., 2016).
raw.dist.mv.sorted.csv and raw.dist.mc.sorted.csv are the specific set of bootstrap permutations used in the manuscript.


External data

Flow files for the raw data used to generate STRAIN.ESTIMATES.1.txt and STRAIN.ESTIMATES.2.txt will be made publicly available. Unfortunately, due to a hard drive failure, the raw flow data for STRAIN.ESTIMATES.3.txt, STRAIN.ESTIMATES.4.txt, STRAIN.ESTIMATES.5.txt and STRAIN.ESTIMATES.6NEW.txt were lost.  

Questions can be directed to the principal investigator Patricia Wittkopp:

Patricia J. Wittkopp, Ph.D. (she/her)
LSA Associate Dean for Natural Sciences 
Deborah E. Goldberg Distinguished University Professor of Ecology & Evolutionary Biology and Molecular, Cellular & Developmental Biology
Arthur F. Thurnau Professor
University of Michigan
Ann Arbor, MI 48109-1085
tel: 734.763.1548 (office); 734.647.5483 (lab)
https://sites.lsa.umich.edu/wittkopp-lab/

References:
Metzger, Brian PH, Fabien Duveau, David C. Yuan, Stephen Tryban, Bing Yang, and Patricia J. Wittkopp. "Contrasting frequencies and effects of cis-and trans-regulatory mutations affecting gene expression." Molecular biology and evolution 33, no. 5 (2016): 1131-1146.
