# The root microbiome code

This repository hosts the code for "Characterizing both bacteria and fungi improves understanding of the Arabidopsis root
microbiome" by Bergelson, Mittelstrass, and Horton

Author: Matt Horton (horton.matthew.w@gmail.com)

This package contains R-scripts for analyzing the root microbiome data described in Bergelson et al. (2018). The code should also work for other types of data, provided those data are in the 'OTU table' file format of Qiime (version 1.3). The code will eventually be extended to process data in the ASV file format.

The code requires the environmental variable 'ROOT_MICROBIOME_DIR' to be set to the data and code stored in the file 'root_microbes.tgz'. For example, in bash, one could use the command:
'''
***export*** ROOT_MICROBIOME_DIR='/Volumes/projects/root_microbes/'
'''

Thanks,
Matt
