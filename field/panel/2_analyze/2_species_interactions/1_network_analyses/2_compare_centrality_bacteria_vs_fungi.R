## This R-script uses the estimate of centrality calculated in:
##  1_estimate_centrality.R
## 
## to determine whether centrality metrics differ among kingdoms (16S vs ITS)
## these differences are plotted in boxplots.
## 
## Author: matt.horton
###############################################################################

rm(list=ls());

require(igraph); require(RColorBrewer);
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.base_methods.R"));

###############################################################################
## hard coded variables
###############################################################################
organ <- "root";
otuCutoff <- 97;

listColors <- list();
listColors[["Bacteria"]] <- "tomato1";
listColors[["Fungi"]] <- "cadetblue2"; 

minimumNumberOfReads <- 400;
sizeOfCommunity <- c(100); 
minimumRhoThreshold <- 0;
normalizeCentralityMeasures <- TRUE;

simpleCap <- function(x) {
	s <- strsplit(x, " ")[[1]]
	paste(toupper(substring(s, 1,1)), substring(s, 2),
			sep="", collapse=" ")
}

###############################################################################
## iterate over the samples of size n_ ... and read the centrality report
###############################################################################
setwd(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "structure/", organ, "/otus", otuCutoff, "/"));

{
	## for the choosefile function:
	cat("Choose the file of pvalues.\n");
	for( j in 1:length(sizeOfCommunity)){

		sizeOfCommunity_j <- sizeOfCommunity[j];
		cat("Working with the community at size:", sizeOfCommunity_j, "\n");
		setwd(paste0("n_", sizeOfCommunity_j));

		cat("Using rho threshold:", minimumRhoThreshold, "\n");
		suffixMatch <- paste0("_minRho", minimumRhoThreshold, ifelse( normalizeCentralityMeasures, "_norml", ""), "_centrality");

		fileName <- chooseFile(c(paste0(organ, ".", minimumNumberOfReads), suffixMatch, "\\.txt$")); ## maybe choose this?
		fileName <- fileName$filename;
		centrality <- read.table(fileName, header=T, sep="\t", row.names=1); ## alt: sort the rhoThresholds, but this forces this to be clean & clear

		###############################################################################
		## we know the centrality of every node... so now
		## we can report these in boxplots or ???
		###############################################################################
		measures <- c( "degree", "between" );
		boxplotFileName <- paste0("boxplot_", gsub(".txt$", ".eps", fileName));

		setEPS(paper="special", horizontal=TRUE);
		postscript(boxplotFileName);

		par(mfrow=c(1, 4));
		for( m in 1:length(measures)){
			measure_m <- measures[m];
			plot(centrality[,measure_m] ~ as.factor(centrality[,"kingdom"]), ylab=simpleCap(measure_m), xlab="", col=unlist(listColors[levels(factor(centrality[,"kingdom"]))]), cex.axis=0.9);
		}

		dev.off();

		## the results from the following code indicate that - among the most heavily sequenced taxa - the fungi have, on average, more keystone/effectual taxa than the bacteria do
		for( m in 1:length(measures)){
			measure_m <- measures[m];
			cat("--------------------------------------\n");
			cat("Analyzing metric:", measure_m, "\n");
			lm0 <- lm(get(measure_m) ~ 1, data=na.exclude(centrality)); 
			lm1 <- lm(get(measure_m) ~ as.factor(kingdom), data=na.exclude(centrality));
			print(anova(lm0, lm1, test="F"));
		}

		## by taxa???
		setwd("../");
	}
}

