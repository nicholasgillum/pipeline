library(limma)
library(maDB)
library(RColorBrewer)
library('yaml')
library('outliers')

getConditions <- function(Slides) {
	candidates = Slides$targets$Name
	dye_mask = grepl("_DF", candidates,ignore.case=TRUE)
#	fwd_mask = grepl("forward", Slides$targets$SlideType, ignore.case=TRUE)
	
	# NDF/fwd and DF/rev
	good_mask <- !dye_mask#((!dye_mask) & fwd_mask) | (dye_mask & (!fwd_mask))
	
	return(matrix(candidates[good_mask],ncol=1))
}
getForwardMask<-function(Slides,candidates) {
	return(grepl("forward", Slides$targets$SlideType))
}

getType<-function(config) {
	if(config$algorithms$median) {
		return("agilent.median");
	} else {
		return("agilent")
	}
}

getUniqueGeneNames<-function(config) {
	if(config$algorithms$median) {
		return(matrix(unique(Slides.norm$genes$ProbeName),ncol=1))
	} else {
		return(matrix(unique(Slides.norm$genes$GeneName),ncol=1))
	}
}

getDyeSwaps <- function(Slides) {
	candidates = Slides$targets$Name
	dye_mask = grepl("_DF", candidates,ignore.case=TRUE)
	
	fwd_mask = grepl("forward", Slides$targets$SlideType, ignore.case=TRUE)
	
	# NDF/fwd and DF/rev
	good_mask <- dye_mask#((!dye_mask) & fwd_mask) | (dye_mask & (!fwd_mask))
	
	return(matrix(candidates[!good_mask],ncol=1))
}

divide <-function(name, Slides,config, conditions, dyeswap, fwdmask) {
	#Find the spots with the given gene name
	if(config$algorithms$median) {
		gene_mask <- Slides$genes$ProbeName == name;
	} else {
		gene_mask <- Slides$genes$GeneName == name;
	}
	
	tmp<-Slides[gene_mask,]
	
	
	#Get intensities
	output <- apply(conditions, 1, collector, conditions,
				dye_swaps, tmp, name,fwdmask)
	
	return(output)
}

collector<-function(condition, conditions, dye_swaps, fSlides,gname,fwdmask) {
	fwd <- sum((condition == conditions) & fwdmask) >0
	
	
	expr <- fSlides$M[, condition]
	df_mask = grepl(paste("^",condition,"_DF",sep="false"),dye_swaps, ignore.case=TRUE)
	if(sum(df_mask) == 1) {
		name <- dye_swaps[df_mask]
		expr_df <- fSlides$M[,name]
	} else if(sum(df_mask) == 0) {
		expr_df <- c()
	} else {
		stop(paste("Error! ambiguous condition name: ", condition))
	}
	
	if(fwd) {
		output <- unlist(c(expr, -1*expr_df))
	} else {
		output <- unlist(c(-1*expr, expr_df))
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
