# DuganFinal-

This code analyzes mass spectrometry data generated using Skyline,
a program that extracts and organizes raw data taken directly from
#the mass spectrometer. The csv file provided contains mass spec data
for proteins isolated from the liver of a Sirt5 knock-out mouse.
Sirt5 is an enzyme that acts as a deacetylase, desuccinylase, and demalonylase
to make post-translational modifications on lysine residues of proteins.
Sirt5 has been shown to be important for mitochondrial metabolism
and is relevant for understanding metabolic diseases.
The csv file contains the fragmentation data for WT and Sirt5 KO proteins
and shows the differences in post-translational modifications in the
modifed sequence column. In this code, an interface will prompt the user to enter 
the number of custom experimental conditions* (2), WT and KO, 
and the code creates a tabular statistical report that includes the average
peak area, standard deviation, and coefficient of variation.

*The two experimental conditions will always be KO and WT  

DuganFinal.py is formatted for Jupyter notebook, while AveragePeakArea.py is formatted with 
a sys argument to be run in command line. The arguments are as follows: 
AveragePeakArea.py, Dugan_MSdata.csv. KO, WT 



