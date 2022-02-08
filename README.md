# Field-goose
Data and code for a project testing the effects of simulated grazing disturbance on eelgrass genotypic diversity

Metadata for project: Distrubance reduces colonization
Authors: Kollars and Stachowicz 

# File: "Kollars_feedback_usats.csv" 
--------------------------------------
This file has the raw microsat genotypes for each sample 

Variables included in the dataset:

population: designation for each plot at each specific time point   

time: 1 = initial sampling, 5 = final sampling

block: experimental block designation. Note that in the MS, we report nine blocks but block numbering in the dataset starts at "2". We originally had 10 blocks and removed the first block from the study due to human sampling error. 

plot: there are four plots in a block labeled as plot A, B, C, D.  A and B are on the shoreward side and C and D are on the channel side of the transect

trt: four levels - 0X: not clipped; 1X: clipped once per season; 2X: clipped twice per season; 4X: clipped 4 times in the first year but not clipped in the second year of the experiment

sample: Up to 36 samples were collected per plot

name: unique name for each sample

X/Y: each sample was collected from a grid of 6 X 6 points. The X/Y coordinates designate the location the sample is from.   

The rest of the columns are headers for each locus with the "_1" and "_2" labels indicated "Allele 1" and "Allele 2", respectively.
Values listed for each locus are fragment lengths with "NA" indicating missing data due to PCR failure and "0" indicating that the sample was not collected in the field. Note that loci 8 and 11 were removed from analysis due to high frequency of PCR failure (a potential indication of null alleles).   

Zm01 - CL559Contig1
Zm02 - CL32Contig2
Zm03 - ZMC12075
Zm04 - ZMC13053
Zm05 - CL172Contig1
Zm06 - ZMC19062
Zm07 - ZMC19017 
Zm08 - ZosmarCT-3
Zm09 - ZosmarCT-12
Zm10 - ZosmarCT-19
Zm11 - ZosmarGA-3

# File: "Kollars_feedback_mlgs"  
--------------------------------
This file contains all the same information as "Kollars_feedback_usats" but with an additional column, "mlg" that has the assigned clonal id for each sample. 

# File: "Kollars_feedback_clones"
---------------------------------
This file is a modified output file of running "Kollars_feedback_usats" through the R package "Rclone" and generating clonal summary stats. The modifications included adding plot identifying information. 

Variables included in the dataset: 

pop: designation for each plot at each specific time point   

time: T1 = initial sampling, T5 = final sampling 

block: experimental block designation. Note that in the MS, we report nine blocks but block numbering in the dataset starts at "2". We originally had 10 blocks and removed the first block from the study due to human sampling error. 

plot: there are four plots in a block labeled as plot A, B, C, D.  A and B are on the shoreward side and C and D are on the channel side of the transect

trt: four levels - 0X: not clipped; 1X: clipped once per season; 2X: clipped twice per season; 4X: clipped 4 times in the first year but not clipped in the second year of the experiment

N:total number of samples (maximum number of samples collected was 36)

G: total number of genotypes

R: Genotypic richness calculated as (G-1)/(N-1)

H..: Shannon-Weiner Index 

J.: Pielou's J (evenness)		
		
D:Simpson's unibased D

V:Simpson's unbiased V

H:Reciprocal of Simpson's index	

# File: "Kollars_feedback_shoots"
--------------------------------
data for vegetative and flowering shoot counts throughout the course of the experiment

Variables included in the dataset: 

block: experimental block designation. Note that in the MS, we report nine blocks but block numbering in the dataset starts at "2". We originally had 10 blocks and removed the first block from the study due to human sampling error. 

plot: there are four plots in a block labeled as plot A, B, C, D.  A and B are on the shoreward side and C and D are on the channel side of the transect

trt: four levels - 0X: not clipped; 1X: clipped once per season; 2X: clipped twice per season; 4X: clipped 4 times in the first year but not clipped in the second year of the experiment

columns listed as dates: inidivudal sampling dates for vegetative shoot counts

"no.flowers.###" where ### = "march2018, april2018, may2018, or june2018": data for number of flowering shoots in each plot. Note that data from march2018 was not included in the analysis. 

# File: "Kollars_feedback_genet_assignment.code"
--------------------------------
R code used to work with the microsatellite data and assign genets

# File: "Kollars_feedback_stats_figure.code"
--------------------------------
R code used to conduct statistical analysis and generate figures 
