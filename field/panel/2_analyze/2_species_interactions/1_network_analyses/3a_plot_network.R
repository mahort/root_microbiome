## This R-script plots the network, using igraph, in a davidson-harel layout.
## 
## Author: matt.horton
###############################################################################

rm(list=ls());

require(igraph); require(RColorBrewer);
source(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "/code/hdr.base_methods.R"));

###############################################################################
## the user can indicate variables here... 
###############################################################################
organ <- "root";
otuCutoff <- 97;

plottingGranularity <- 6;

minimumNumberOfReads <- 400;
sizeOfCommunity <- c(100);
normalizeCentralityMeasure <- TRUE;

minimumRho <- 0; ##

{
	###############################################################################
	## the following function can be used to relocate the label; in the plot
	## function, one needs to specify these new locations using the argument:
	## ... , vertex.label.degree=newLabelLocations, ...
	###############################################################################
	radian.rescale <- function(x, start=0, direction=2) {
		c.rotate <- function(x) (x + start) %% (2 * pi) * direction;
		c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)));
	}

	getColors <- function(theValues, colorRange){

		minimumValue <- min(theValues);
		maximumValue <- max(theValues);
		
		breaks <- seq(minimumValue - 0.03, maximumValue + 0.03, len=(1 + length(colorRange)));
		
		correspondingSlices <- cut(theValues, breaks);
		print(levels(correspondingSlices));
		
		assignedColors <- colorRange[match(correspondingSlices, levels(correspondingSlices))];
		
		return(list(colors=assignedColors, slices=levels(correspondingSlices), keycolors=colorRange));
	}

	setwd(paste0(Sys.getenv('ROOT_MICROBIOME_DIR'), "structure/", organ, "/otus", otuCutoff, "/"));

	## for the choosefile function:
	cat("Choose the graph to plot.\n");
	for( j in 1:length(sizeOfCommunity)){

		sizeOfCommunity_j <- sizeOfCommunity[j];
		cat("Working with the community at size:", sizeOfCommunity_j, "\n");
		setwd(paste0("n_", sizeOfCommunity_j));

		###############################################################################
		## the graph knows the (1.) degree and (2.) betweenness centrality of every node
		###############################################################################
		fileName <- chooseFile(c(paste0(organ, ".", minimumNumberOfReads), "gml$", ifelse(normalizeCentralityMeasure, "norml", ""))); 
		fileName <- fileName$filename;		
		graph <- read.graph( fileName, "gml" );

		E(graph)$color <- getColors(E(graph)$cor, brewer.pal(n=plottingGranularity, "RdBu"))$colors; #ifelse( E(graph)$cor[k] > 0, "seagreen4", "firebrick4");
		E(graph)$width <- 0.75;
		E(graph)$curved <- TRUE;
		
		## 
		dhLayout <- layout.davidson.harel(graph);

		measures <- c("degree"); #, "between"); 
		for( m in 1:length(measures)){
			measure_m <- measures[m];

			epsOutputFileName <- gsub(".gml", paste("_g", plottingGranularity, "_rho", minimumRho, "_dh.", measure_m, "_centrality.sc", ".eps", sep=""), fileName);
			postscript(epsOutputFileName, width=7, height=7);
			plot(graph, layout=dhLayout, vertex.size=(7*get.vertex.attribute(graph, measure_m))/max(get.vertex.attribute(graph, measure_m)),
					vertex.label=NA, vertex.label.cex=0.5, vertex.label.family="Helvetica"); ##
			dev.off();
		}

		setwd("../");
	}
}

