## Analyzing the Arabidopsis root microbiome

Authors: Matthew W. Horton (horton.matthew.w@gmail.com), 

This package contains code (R-scripts) for analyzing plant-microbial communities; although the code was originally developed for [*Arabidopsis thaliana*](https://www.nature.com/articles/s41598-018-37208-z), it can be applied to other organisms provided that the data are in the OTU-table file format of Qiime (version 1.3). The code will eventually be extended to process data in the ASV file format used by newer versions of QIIME. Suggestions and code contributions are welcomed.

The current version is 1.0

### The main dependencies are:
* R
* The R-package vegan

Unless otherwise noted, the R-scripts were written using R version 3.4.3.

### Environmental variables need to be set

To avoid using hardcoded paths, we use R's ability to load environmental variables (using the command syntax Sys.getenv('variable_name')).
The variables that need to be loaded, by the user, can be done so using the following bash code:

```console
*export ROOT_MICROBIOME_DIR='<your_directory>/root_microbes/'*
```

Thanks,
Matt
