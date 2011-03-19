# TODO: Add comment
# 
# Author: nag9
###############################################################################
library(limma)
library(maDB)
library(RColorBrewer)
library('yaml')

config<-yaml.load_file("config/config.yaml")
source(paste(config$files$library,"/pipeline_utils.R",sep=""))
source(paste(config$files$library,"/latex_utils.R",sep=""))

targets <-readTargets(config$files$metadata)
targets$FileName <- gsub(" ", "", targets$FileName)
Slides.raw<-read.maimages(files=targets$FileName,source=getType(config), 
		path=config$files$rawDir, names=targets$Name, verbose=TRUE)
Slides.raw$targets <-targets
config$count<-length(targets$FileName)
config$targets<-targets
config$vsn<-(config$algorithms$betweenArray=="vsn")


Slides.carma <- read.maimages(files = c(
	"US22502698_253010810006_S01_GE2-v5_95_Feb07_1_1.txt",
	"US22502698_253010810006_S01_GE2-v5_95_Feb07_2_1.txt",
	"US22502698_253010810004_S01_GE2-v5_95_Feb07_1_1.txt",
	"US22502698_253010810004_S01_GE2-v5_95_Feb07_2_1.txt"), 
	names = c("-40 ura ", "-40 ura DF", "60 ura", "60 ura DF"),
	source = "agilent", sep = "\t", 
	columns = list(Rf = "rMeanSignal", Gf = "gMeanSignal",
					Rb = "rBGMedianSignal", Gb = "gBGMedianSignal"),
			path=config$files$rawDir)
	
Slides.raw <-backgroundCorrect(Slides.raw, method=config$algorithms$bgCorrection)
Slides.carma <- backgroundCorrect(Slides.carma, method="normexp")

Slides.norm <- normalizeWithinArrays(Slides.raw, layout=Slides.raw$printer,
		method=config$algorithms$withinArray)
rm(Slides.raw)
gc()

Method <- "printtiploess"
if (is.null(Slides.carma$printer)) {
			Method <- "loess"
}
Slides.carma <- normalizeWithinArrays(Slides.carma, layout = Slides.carma$printer,
	method = Method)
gc()

Slides.norm<-normalizeBetweenArrays(Slides.norm, method=config$algorithms$betweenArray)
Slides.carma <-normalizeBetweenArrays(Slides.carma, method = "quantile")

Eset <- newMadbSet(Slides.carma)
Eset <- average(Eset, average.which = c(1, 2, 2, 1, 3, 4, 4, 3), 
	method = "mean", array.names = c("40 ura A Red", "40 ura A Green",
		"-40 ura A Red", "-40 ura A Green"), average.genes = TRUE,
		exclude.flagged = TRUE)


control_mask<-(Slides.norm$genes$ControlType != 0)
Slides.norm<- Slides.norm[!control_mask,]
cat("Removed ", sum(control_mask), " control probes.\n")
gene_names <- getUniqueGeneNames(config)
cat(length(gene_names), " genes were found. ")
