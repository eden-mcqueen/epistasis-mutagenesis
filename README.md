# epistasis-mutagenesis
Data and code for 2026 Epistasis Mutagenesis paper

General Information

Repository contains files associated with the epistasis yeast mutagenesis project in the Wittkopp laboratory at the University of Michigan.
Project origination and majority of code was written by former Wittkopp laboratory graduate student Dr. Bing Yang. Additions, modifications, and cuts were made by Dr. E. W. McQueen for the final analysis and figures for publication.


Included files


Python scripts:

1. Create_template.py
Creates template file (“TEMPLATE.csv”) for importing and analyzing flow data.

R scripts:

1. Clean.Data.R
Uses “TEMPLATE.csv”
Reading in fcs files and transforming data

2. Cleaning.Functions.2.R
Functions used in “Clean.Data.R”

3. ANALYSIS.R
Base scripts for analyzing and scoring mutants from transformed data from “Clean.Data.R”. Further cleaning of the data by removing small count samples, etc.
Calculations for the medians, expression estimates, etc.
This code was run at separate times for different portions of the data, so the output file names include ascending numbers. 
Outputs are the "STRAIN.ESTIMATES.txt" files, below.

4. Compile.Data.R
Takes STRAIN.ESTIMATES.txt files and compiles them into one file for analysis called “SUMMARY.PER.ISOLATE.csv”.

5. Analysis.1.13.R
Script to run final analysis for manuscript figures and tables using processed data files.

Processed data files:
1. STRAIN.ESTIMATES.1.txt
2. STRAIN.ESTIMATES.2.txt
3. STRAIN.ESTIMATES.3.txt
4. STRAIN.ESTIMATES.4.txt
5. STRAIN.ESTIMATES.5.txt
6. STRAIN.ESTIMATES.6NEW.txt
7. SUMMARY.TRANS.txt
8. SUMMARY.PER.ISOLATE.csv
9. raw.dist.mv.sorted.csv  
10. raw.dist.mc.sorted.csv
Note: Item 7 uses previously collected data (Metzger et al., 2016).
Items 9 and 10 contain the specific set of bootstrap permutations used in the manuscript.

Columns in the STRAIN.ESTIMATES.X.txt files (items 1-6 above) are defined below:

STRAIN
Identifier for the genotype or experimental group (e.g., mutant, control, wild-type).

PLATE
Identifier for the experimental plate (batch in which the sample was processed/measured).

POSITION
The specific well/location on the plate where the sample was processed.

REPLICATE
Number of valid, non-outlier measurements for the STRAIN+PLATE+POSITION group after outlier removal.

YFP.MEDIAN.R2WT.MEAN
Mean (across replicates) of plate-corrected, background-subtracted median YFP fluorescence, normalized to wild-type control for each group. This value is labeled median.mean in SUMMARY.PER.ISOLATE.csv.

YFP.MEDIAN.R2WT.SD
Standard deviation among replicates of the above median value, showing reproducibility relative to wild-type. This value is labeled median.sd in SUMMARY.PER.ISOLATE.csv.

YFP.SD.R2WT.MEAN
Mean (across replicates) of plate-corrected, background-subtracted YFP SD, normalized to wild-type control. This value is labeled sd.mean in SUMMARY.PER.ISOLATE.csv.

YFP.SD.R2WT.SD
Standard deviation among replicates of the above YFP SD value. This value is labled sd.sd in SUMMARY.PER.ISOLATE.csv.


External data

Flow cytometry files containing the raw data used to generate STRAIN.ESTIMATES.1.txt and STRAIN.ESTIMATES.2.txt are available on University of Michigan’s Deep Blue Data Repository. Unfortunately, flow cytometry files containing the raw data used to generate STRAIN.ESTIMATES.3.txt, STRAIN.ESTIMATES.4.txt, STRAIN.ESTIMATES.5.txt and STRAIN.ESTIMATES.6NEW.txt were lost to a hard drive failure.  

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
