## Quantifying copy number aberration from shallow whole-genome sequencing 


* This document contains R code necessary to reproduce analyses in the manuscript: Mouliere, Chandrananda, Piskorz and Moore et al. **Boosting non-invasive cancer genomics with fragment size analysis of plasma cell-free DNA** (Manuscript submitted). 
* The sections below explain how copy number aberration (CNA) was quantified from shallow whole-genome sequencing (sWGS, <0.5x) data generated from the Illumina HiSeq platform.


![Analysis workflow](/images/tMAD_figure.png)

## Install CNAclinic software

*  CNAclinic is a bioinformatics pipeline in R that can facilitate end-to-end analysis of CNAs from sWGS (Manuscript under review).
* For more information see https://github.com/sdchandra/CNAclinic

```R

install.packages("devtools") 
library(devtools)

install_github("sdchandra/CNAclinic", build_vignettes = TRUE))
library(CNAclinic)

```

## Downsample sWGS data to have e.g. 10 Million reads

* This can be done using samtools within the R environment
* samtools must be available in your PATH

```R

###########################################
# 3 Arguments to change
###########################################

total_reads_needed <- 10000000

bam_path <- "/path/to/bamfiles"

downsample_dir <- "/path/to/write/new/bamfiles"

###########################################
###########################################


command <- paste0("ls ", bam_path, "/*.bam")
bamfiles <- as.character(system(command, intern=TRUE, wait=TRUE))

bamfiles <- unlist(lapply(bamfiles, function(x){
    y <- strsplit(x, "/")[[1]]
    y[length(y)]
    
}))

random_seed <- 8
total_reads_label <- paste0(round(total_reads_needed/1000000, 2), "Million")


for(t in 1:length(total_reads_needed)){
    
    total_reads <- total_reads_needed[t]

    for(i in 1:length(bamfiles)){
    
        bamfile <- paste0(bam_path, "/", bamfiles[i])
    
        if(file.exists(bamfile)){
        
            # Get the read count total from BAM
            command <- paste0("samtools view -c -q 20 ", bamfile)
            
            BAM_total <- as.numeric(system(command, intern=TRUE, wait=TRUE))
        
            if(total_reads >= BAM_total){
                
             
                stop(paste0("BAM has lower read count than ", 
                	total_reads_label[t], 
                	bamfiles[i]))         
                     
            }else{
                
                sampling_prop <- 8 + round(total_reads/BAM_total, 3)
                
                command <- 
                    paste0("samtools view -s ", 
                           sampling_prop, 
                           " -q 20 -b ", bamfile, " > ", 
                           downsample_dir, 
                           "/samp_", total_reads_label[t], "_", bamfiles[i])
                
                system(command, intern=FALSE, wait=TRUE)
            }        
        
        }else{
        
           stop(paste0("File does not exists : ", bamfile))
        
        }
    }
}


 
```


## Create an empirical blacklist from control data


```R


library(CNAclinic)

###########################################
# Arguments to change
###########################################

binSizes <- c(30) # Reads binned into 30 Kbp windows

total_reads_needed <- 10000000

downsample_dir <- "/path/to/sampled/bamfiles"

###########################################
###########################################

total_reads_label <- paste0(round(total_reads_needed/1000000, 2), "Million")

for(b in 1:length(binSizes)){
    
    binSize <- binSizes[b]
    
    for(t in 1:length(total_reads_label)){
        
        bamfiles <- Sys.glob(paste0(downsample_dir,
                "/samp_", total_reads_label[t], "_*.bam"))
        
        
        userMadeBins <- QDNAseq::getBinAnnotations(binSize=binSize,
                                                   genome="hg19")
        
        
        readCounts <- QDNAseq::binReadCounts(bins=userMadeBins,
                                             bamfiles=bamfiles,
                                             cache=TRUE,
                                             pairedEnds=TRUE)
        
        ctrl <- readCounts
        
        readCounts <- QDNAseq::applyFilters(readCounts, residual=FALSE, 
                                            blacklist=FALSE,
                                            mappability=FALSE, 
                                            bases=FALSE,
                                            chromosomes = c("X", "Y", "MT", "M"))
        
        userMadeBins$residual <- QDNAseq::iterateResiduals(readCounts)
        
        chromosomes = c("X", "Y", "MT", "M")
        
        # Create a residual filter from cfDNA controls
        condition <- rep(TRUE, times=nrow(readCounts))
        condition <- !(Biobase::fData(readCounts)$chromosome %in% chromosomes)
        condition <- condition & !is.na(Biobase::fData(readCounts)$gc)
        
        residuals <- userMadeBins$residual
        cutoff <- TRUE * matrixStats::madDiff(residuals, na.rm=TRUE)
        residualsMissing <- aggregate(residuals,
                            by=list(chromosome=Biobase::fData(readCounts)$chromosome),
                            function(x) all(is.na(x)))
        chromosomesWithResidualsMissing <-
                residualsMissing$chromosome[residualsMissing$x]
        chromosomesToInclude <-
                setdiff(chromosomesWithResidualsMissing, chromosomes)
        if (length(chromosomesToInclude) > 0) {
            message("Note: Residual filter missing for chromosomes: ",
            paste(chromosomesToInclude, collapse=", "))
                residuals[Biobase::fData(readCounts)$chromosome %in% chromosomesToInclude] <- 0
        }

        # If FALSE, filter genomic bin from analysis, if TRUE keep bin
        condition <- condition & !is.na(residuals)

        saveRDS(condition, 
                file=paste0("control_blacklist_", total_reads_label[t], 
                            "_", binSize, "Kbp"))   
        
    }
}

```

## Calculate t-MAD from test samples

* creates a .csv file with t-MAD values per sample 

```
###########################################
# Arguments to change
###########################################

downsample_dir <- "/path/to/sampled/bamfiles"

total_reads_needed <- 10000000 

selected_binSize <- 30  

downsampled_control <- "path/to/sampled/control/samp_10Million_XXX.bam"

control_name <- "CONTROL_XXX"

path_to_blacklist <- "path/to/saved/blacklist"


###########################################
# Check the code below for necessary changes
###########################################

library(CNAclinic)

total_reads_label <- paste0(round(total_reads_needed/1000000, 2), "Million")

controlSample <- c(control_name, "none")

segType <- c("CH")


for(t in 1:length(total_reads_label)){
    
    if(length(total_reads_label) != length(selected_binSize))
        stop("Check variables")
    
    binSize <- selected_binSize[t]
    
    cfDNA_blacklist <- readRDS(paste0(path_to_blacklist, 
                                      "control_blacklist_", 
                                      total_reads_label[t], 
                                      "_", binSize, "Kbp"))
    
    bamfiles <- Sys.glob(paste0(downsample_dir,
                                "/samp_", total_reads_label[t], "_*.bam"))
    
    bamnames <- unlist(lapply(strsplit(bamfiles, "/"), function(x){ x[[length(x)]]}))
    
    bamfiles <- c(downsampled_control, bamfiles)
    
    bamnames <-  c(control_name, bamnames)
    
    for(j in 1:length(controlSample)){
        
        # Doing only median normalization for these samples        
        control <- controlSample[j]
        
        processedData <- NULL
        
        if(control == "none"){
            
            outfile_suffix <- paste0(total_reads_label[t], "_",
                                     binSize, "_control_",
                                     control, "_none")
            
            processedData <- processForSegmentation(
                                bamfiles=bamfiles,
                                binSize=binSize,
                                chromosomesFilter=c("X", "Y", "MT", "M"),
                                cache=FALSE,
                                isPaired=TRUE,
                                saveCountData=FALSE)


            # Run segmentation
            if(segType == "CH"){
                CNAData <- runSegmentation(processedData, genome="hg19",
                                           segmentType=c("CBS", "HMM"),
                                           summaryMethod="mean")
                
                saveRDS(CNAData, file=paste0("CNAData_", outfile_suffix))
            }
        }else if(control == control_name){
            
            if(length(control_name) != 1)
                stop("Check control specified")
            
            which_control <- which(bamnames == control_name)
            refSamples <- rep(control_name, length(bamnames))
            
            # drop the control sample after normalizing test samples 
            # as we will have log2R=0 if it is normalized by its own bin counts
            refSamples[which_control] <- 'drop'
            
            outfile_suffix <- paste0(study, "_", libPrep[i], "_",
                                     total_reads_label[t], "_",
                                     binSize, "_control_",
                                     control, "_noMedianNorm")
            
            processedData <- processForSegmentation(cache=FALSE,
                                bamfiles=bamfiles,
                                bamnames=bamnames,
                                refSamples=refSamples,
                                binSize=binSize,
                                isPaired=TRUE,
                                skipMedianNormalization=TRUE,
                                saveCountData=FALSE)
            
            

            # Run segmentation
            if(segType == "CH"){
                CNAData <- runSegmentation(processedData, genome="hg19",
                                           segmentType=c("CBS", "HMM"),
                                           summaryMethod="mean")
                
                saveRDS(CNAData, file=paste0("CNAData_", outfile_suffix))
            }
        }
        
        sampleNames <- unlist(lapply(sampleNames(CNAData), 
                        function(x){strsplit(x, ".vs.")[[1]][1]}))
        
        sampleNames(CNAData) <- sampleNames
        
        stopifnot(length(cfDNA_blacklist) == length(usebin(CNAData)))
        
        segData <-  segSummary(CNAData)
        
        segData <- segData[cfDNA_blacklist & usebin(CNAData), ]
        
        segData[abs(segData) >=5 ] <- NA
        
        tMAD <- apply(segData, 2,  function(x){ 
            abs(mad(x = x, center=0, na.rm = TRUE))
        })
        
        outData <- data.frame(sampleNames=sampleNames(CNAData),
                              tMAD,
                              stringsAsFactors=FALSE)

        write.csv(outData, file=paste0("tMAD_", outfile_suffix, ".csv"), 
                  quote=F, row.names=F)
        
    }
}

```




 

