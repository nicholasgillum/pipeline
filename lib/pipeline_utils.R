library(limma)
library(maDB)
library(RColorBrewer)
library('yaml')
library('outliers')

getConditions <- function(Slides) {
	candidates = Slides$targets$Name
	dye_mask = grepl("_DF", candidates,ignore.case=TRUE)
	return(matrix(candidates[!dye_mask],ncol=1))
}

getDyeSwaps <- function(Slides) {
	candidates = Slides$targets$Name
	dye_mask = grepl("_DF", candidates,ignore.case=TRUE)
	return(matrix(candidates[dye_mask],ncol=1))
}

divide <-function(name, Slides,config, conditions, dyeswap) {
	#Find the spots with the given gene name
	gene_mask = Slides$genes$GeneName == name;
	tmp = Slides[gene_mask,]
	
	#Get intensities
	output <- apply(conditions, 1, collector, 
				dye_swaps, tmp, config$setup$exp_color, name)
	
	return(output)
}

collector<-function(condition, dye_swaps, fSlides, exp_dye,gname) {
	
	expr <- fSlides$M[,condition]
	df_mask = grepl(paste("^",condition,sep="false"), dye_swaps, ignore.case=TRUE)
	
	if(sum(df_mask) == 1) {
		name <- dye_swaps[df_mask]
		expr_df <- fSlides$M[,name]
	} else if(sum(df_mask) == 0) {
		expr_df <- c()
	} else {
		stop(paste("Error! ambiguous condition name: ", condition))
	}
	
	if(exp_dye == "R"){
		output = unlist(c(expr, -1*expr_df))
	} else if(exp_dye == "G") {
		output = unlist(c(-1*expr, expr_df))
	} else {
		stop(paste("Unrecognized exp dye color (not R or G): ", exp_dye))
	}
	
	return(collapse(output, paste(condition,"-",gname)));
}

collapse <- function(dat, name) {
	cleaned <- rm_outliers(dat, name)
	return(mean(cleaned))
}

rm_outliers  <-function(dat, name) {
	punk <- outlier(dat)
	outlier_mask = dat==punk
	test=dixon.test(dat)
	if(test$p.value < config$dixon$cutoff) {
		cat("Removed outlier ", punk, " from dataset ", 
				name, "(p=", sprintf('%0.5e',test$p.value), ')\n')
		return(dat[!outlier_mask])
	} else {
		return(dat)
	}
}
