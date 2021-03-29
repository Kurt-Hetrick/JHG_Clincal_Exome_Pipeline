

# Eric Minikel
# script to run ExomeDepth
# how to run:
# runExomeDepth.r -b bamlist.txt -o /output/path/ -v > runExomeDepthOutput.txt

start_time = Sys.time()
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings
# sessionInfo()
.libPaths("/mnt/linuxtools/clinical/R/3.5.3/lib64/R/library")
library(optparse) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
library(ExomeDepth)
##################################################
### this part as replacement of QCpipeline library
#library(QCpipeline)
readConfig <- function(file, ...) {
  config.table <- read.table(file, as.is=TRUE, ...)
  if (any(duplicated(config.table[, 1]))) stop("duplicated parameters in config file are not allowed!")
  config <- config.table[,2]
  names(config) <- config.table[,1]
  # recode tabs
  config[config %in% "\\t"] <- "\t"
  
  return(config)
}

setConfigDefaults <- function(config, required, optional, default) {
  stopifnot(length(optional) == length(default))
  
  config.params <- names(config)
  found.params <- intersect(config.params, c(required, optional))
  if (length(found.params) > 0) {
    message("found parameters: ", paste(found.params, collapse=", "))
  }
  
  # if required params not in config, stop
  missing.params <- setdiff(required, config.params)
  if (length(missing.params) > 0) {
    stop("missing required parameters: ", paste(missing.params, collapse=", "))
  }
  
  # if not in config, set default value
  set.params <- setdiff(optional, config.params)
  if (length(set.params) > 0) {
    config[set.params] <- default[match(set.params, optional)]
    message("using default values: ", paste(set.params, collapse=", "))
  }
  
  # note unsed params in config
  extra.params <- setdiff(config.params, c(required, optional))
  if (length(extra.params) > 0) {
    message("unused parameters: ", paste(extra.params, collapse=", "))
  }
  
  # return config with default values set
  return(config)
}

##########end of read.config function#######################################

## read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check for type
if (length(args) < 4) stop("missing bam")
name <- args[2]
sex <- args[3]
bam <- args[4]

# check config and set defaults
if (sex %in% c("female", "F", "Female", "FEMALE")) {
  required <- c("ref_female")
  refcounts <- config["ref_female"]
  cat(paste("female referecence panel selected for sample ",name,"\n",sep=''),file=stdout())
} else if (sex %in% c("male", "M", "Male", "MALE")) {
  required <- c("ref_male")
  refcounts <- config["ref_male"]
  cat(paste("male referecence panel selected for sample ",name,"\n",sep=''),file=stdout())
} else {
  required <- c("ref_all")
  refcounts <- config["ref_all"]
  cat(paste("all referecence panel selected for sample ",name,"\n",sep=''),file=stdout())
}

## check config and set defaults
optional <- c("bed_file", "out_dir")
default <- c(NA, "./")
config <- setConfigDefaults(config, required, optional, default)
print(config)

if (file.exists(config["bed_file"])) {
  # read bam list directly into a vector (note use of $V1)
  bed = read.table(config["bed_file"],header=TRUE, stringsAsFactors = FALSE)
  colnames(bed) = c("chromosome", "start", "end","name")
} else {
  cat("You need to specify a valid bed file using -f.\n",file=stderr())
  stop()
}

if (file.exists(refcounts)) {
  # refcounts
  load(refcounts)
} else {
  cat("You need to specify the name of reference set using -r.\n",file=stderr())
  stop()
}

# read output directory
# note right now if not specified, I stop execution. 
# an alternative is to create the dir.
# see http://stackoverflow.com/questions/4216753/check-existence-of-directory-and-create-if-doesnt-exist
if (file.exists(config["out_dir"])) {
  wk_dir <- paste(config["out_dir"], "/", name, sep="")
  system(paste("mkdir -p ", wk_dir, sep=""))
  setwd(wk_dir)
} else {
  cat("You need to specify a valid output directory using -o.\n",file=stderr())
#    cat(paste("The directory you specified was '",opt$outdir,"'.",sep=''),file=stderr())
  stop()
}

# if (opt$verbose) {
     cat(paste("Read BAM list from ",name,"\n",sep=''),file=stdout())
# }
# get the read count if it exist in .rda
counts_file <- grep(".+\\.rda$", list.files(), value = TRUE)
# file.rename(counts, "counts.rda")
# if there is a unique rda file, using the existing rda
if (length(counts_file) == 1) {
  load(counts_file)
  cat(paste("Use the existing counts for this sample ", counts_file, "\n", sep=''),file=stderr())
} else {
  # count the study sample
  counts = getBamCounts(bed.frame = bed, bam.files = bam)
#  save(counts,file=paste("counts", ".rda", sep=""))
  cat(paste("Calculated counts\n",sep=''),file=stdout())
}


#####
# # re-load the counts, 'ref_counts'
ref_countdf = as.data.frame(refcounts)
ref_countmat = as.matrix(ref_countdf[,6:dim(ref_countdf)[2], drop = FALSE])
# counts is an S4 object.
# you need to cast it to a data frame (for bin length)
# AND to a matrix (for reference.count)
# and for some reason you can't cast S4 directly to matrix, only via df
countdf = as.data.frame(counts)
countmat = as.matrix(countdf[,6:dim(countdf)[2], drop = FALSE]) # remove cols 1-5 metadata

# for plotting 
min_pos <- aggregate(start ~ chromosome, bed, function(x) min(x))
max_pos <- aggregate(end ~ chromosome, bed, function(x) max(x))

# keep in mind the min_pos and max_pos order is not numeric
chr_order <- c("1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
# reorder chr in chr_order
min_pos <- min_pos[match(chr_order, min_pos$chromosome),]
max_pos <- max_pos[match(chr_order, max_pos$chromosome),]

# for (i in 1:dim(countmat)[2]) {
i <- 1
    sample_name = colnames(countmat)[i]
    reference_list = select.reference.set(test.counts = countmat[,i], 
        reference.count = ref_countmat,
        bin.length=(countdf$end-countdf$start)/1000,
        n.bins.reduced = 10000)
    reference_set = apply(
        X = as.matrix(ref_countmat[, reference_list$reference.choice]), 
        MAR=1, FUN=sum)
    all_exons = new('ExomeDepth', test=countmat[,i], 
        reference=reference_set,
        formula = 'cbind(test,reference) ~ 1')
    all_exons = CallCNVs(x = all_exons, transition.probability=10^-4,
        chromosome=countdf$space, start=countdf$start,
        end=countdf$end, name=countdf$names)
   
     pdf(paste(sample_name, "_all_chr_CNV_plot_cnt30.pdf", sep=""))
    for (k in 1:dim(min_pos)[1]-2){
      title = paste("chr_", k, "_",sample_name, sep="")
      plot (all_exons,
            sequence = as.character(k),
            xlim = c(min_pos$start[k], max_pos$end[k]),
            count.threshold = 30,
            main = title,
            cex.lab = 0.8,
            with.gene = FALSE)
    }

    # plot chr X
     title = paste("chr_", "X", "_",sample_name, sep="")
     plot (all_exons,
        sequence = "X",
        xlim = c(min_pos$start[23], max_pos$end[23]),
        count.threshold = 30,
        main = title,
        cex.lab = 0.8,
        with.gene = FALSE)
		
	 # plot chr Y
     title = paste("chr_", "Y", "_",sample_name, sep="")
     plot (all_exons,
        sequence = "Y",
        xlim = c(min_pos$start[23], max_pos$end[23]),
        count.threshold = 30,
        main = title,
        cex.lab = 0.8,
        with.gene = FALSE)	
    dev.off()

    # add total number of calls and cnv size columns HERE
    calls <- all_exons@CNV.calls
    calls$sample <- sample_name
    calls$num.calls <- dim(calls)[1]
    # reorder the columns for annotSV annotation, index of "type" will affect command of annotSV
    calls <- calls[c("chromosome","start","end","type","sample", "num.calls", "id","BF","reads.expected","reads.observed","reads.ratio","start.p", "end.p", "nexons")]
    # table(calls$chromosome)
    exomeDepth_output <- paste(sample_name,".bed",sep='')
    write.table(calls, file=exomeDepth_output, 
        sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
#    if (opt$verbose) {
        cat(paste("Total # CNVs called for ",sample_name," is: ", dim(calls)[1], "\n",sep=''),file=stdout())
        cat(paste("Test_sample: ",sample_name,"\n",sep=''),file=stdout())
        cat(paste(reference_list[[1]], sep=' '), file=stdout())
        cat(paste("\n",sep=''),file=stdout())
#    }

   # compute the percentage
   cnv_call_pct_file <- paste(sample_name,".call.pct.by.chr.txt",sep='')
   code_dir <- config["script_dir"]
   bed_pct_file <- config["bed_pct"]
   cmd <- paste("bash ",code_dir, "/scripts/compute.chr.cover.pct.by.sample.sh ",  exomeDepth_output," ", bed_pct_file, " >", cnv_call_pct_file, sep="")
   system(cmd)
   # # set the threshold of cnv percentage to 5% for each chr
   # cnv_pct_threshold <- 3
    pct_chr <- read.table(cnv_call_pct_file, header=FALSE)
    colnames(pct_chr) <- c("chr", "pct")
    pct_chr_out_file <- paste(sample_name,".call.pct.by.chr.txt",sep='')
    write.table(pct_chr, file=pct_chr_out_file, row.names=FALSE, quote=FALSE, sep="\t")
   # pct_chr_out <- pct_chr[pct_chr$pct>cnv_pct_threshold,]
   # annote with annotSV
   # set the path, it's currently in the runall.sh script, can also be added to user's .bash_profile
   # system(paste("export ANNOTSV=",config["annotSV_dir"], sep=''))
   # system("export ANNOTSV=/mnt/research/statgen/AnnotSV/AnnotSV_2.0")
   annotSV_temp <- paste(sample_name,".annotSV.temp.tsv",sep='')
   cmd <- paste("$ANNOTSV/bin/AnnotSV --SVinputFile ", exomeDepth_output," -outputFile ", annotSV_temp, " -typeOfAnnotation split -outputDir ./ -SVinputInfo 1 -svtBEDcol 4", sep="")
   system(cmd)

   annotSV <- read.table(annotSV_temp, header=TRUE, stringsAsFactors = FALSE, sep="\t")
# NEED TO CHANGE TO TSV INSTEAD OF CSV FORMAT
   colname <- names(annotSV)
   colname[7:16] <- c("sample", "num.calls", "id","BF","reads.expected","reads.observed","reads.ratio","start.p","end.p","nexons")
   #Now set the column name in data frame
   colnames(annotSV) <- colname
   # extract and rearrange columns
   annotSV_out <- annotSV[, c("sample", "SV.chrom","SV.start","SV.end","SV.length","Gene.name", "SV.type", "num.calls" ,
                            "id","BF","reads.expected","reads.observed", "reads.ratio", "nexons","NM","location",
                            "DGV_GAIN_IDs","DGV_GAIN_n_samples_with_SV", "DGV_GAIN_n_samples_tested",  "DGV_GAIN_Frequency",
                            "DGV_LOSS_IDs","DGV_LOSS_n_samples_with_SV","DGV_LOSS_n_samples_tested",  "DGV_LOSS_Frequency",
                            "DDD_SV","DDD_DUP_n_samples_with_SV",  "DDD_DUP_Frequency", "DDD_DEL_n_samples_with_SV","DDD_DEL_Frequency", 
                            "X1000g_event", "X1000g_AF","X1000g_max_AF","promoters","dbVar_event", "dbVar_variant",  "dbVar_status", 
                            "ACMG",  "HI_CGscore", "TriS_CGscore", "DDD_status", "DDD_mode","DDD_consequence", "DDD_disease", 
                            "DDD_pmids","HI_DDDpercent","synZ_ExAC", "misZ_ExAC", "pLI_ExAC", "delZ_ExAC", "dupZ_ExAC","cnvZ_ExAC",
                            "Mim.Number","Phenotypes","Inheritance","AnnotSV.ranking")]

   annotSV_out_file <- paste(sample_name,".annotSV.result.tsv",sep='')
   write.table(annotSV_out, file=annotSV_out_file, row.names=FALSE, quote=FALSE, sep="\t")
   # extract genes of interest
  if (file.exists(config["gene_list_file"])) {
    geneList <- read.csv(config["gene_list_file"],as.is=TRUE,header=TRUE)
    zoom <- annotSV_out[annotSV$Gene.name %in% geneList$Gene,]
    zoom_file <- paste(sample_name,".annotSV.result.zoom.", config["gene_list_name"],".tsv",sep='')
    write.table(zoom, file=zoom_file, row.names=FALSE, quote=FALSE, sep="\t")
  } else {
    cat("No subset result output as no zoom gene list provided\n",file=stdout())
  }
   
   # library(plyr)
   # if (!empty(pct_chr_out)) {
   #   pct_chr_out_file <- paste(sample_name,".call.pct.cutoff.", cnv_pct_threshold, ".by.chr.txt",sep='')
   #   write.table(pct_chr_out, file=pct_chr_out_file, row.names=FALSE, quote=FALSE, sep="\t")
   #   # # plot the intensities for big chr anomaly
   #   #     pdf(paste(sample_name, "_all_chr_CNV_plot_cnt30.pdf", sep=""))
   #   #     for (k in 1:dim(min_pos)[1]){
   #   #       title = paste("chr_", k, "_",sample_name, sep="")
   #   #       plot (all_exons,
   #   #             sequence = k,
   #   #             xlim = c(min_pos$start[k], max_pos$end[k]),
   #   #             count.threshold = 30,
   #   #             main = title,
   #   #             cex.lab = 0.8,
   #   #             with.gene = FALSE)
   #   #     }
   #   #     # plot chr X
   #   #      title = paste("chr_", "X", "_",sample_name, sep="")
   #   #      plot (all_exons,
   #   #         sequence = "X",
   #   #         xlim = c(min_pos$start[23], max_pos$end[23]),
   #   #         count.threshold = 30,
   #   #         main = title,
   #   #         cex.lab = 0.8,
   #   #         with.gene = FALSE)
   #   #     dev.off()
   # } else {
   #   # only remove exomeDepth results when there is no big chr anomaly found
   #   file.remove(exomeDepth_output)
   # }
   # remove the intermediate files
   # file.remove(cnv_call_pct_file)
   file.remove(annotSV_temp)
   save(counts,file=paste(sample_name, ".rda", sep=""))
#   file.rename(counts, paste(sample_name, ".rda", sep=''))
# }

duration = format(Sys.time() - start_time)

# if(opt$verbose) {
    cat(paste("Completed execution in ",duration,"\n",sep=''),file=stdout())
# }

#     all columns
#     [1] "AnnotSV.ID"                 "SV.chrom"                   "SV.start"                   "SV.end"                     "SV.length"                  "SV.type"                   
#     [7] "sample"                     "num.calls"                  "id"                         "BF"                         "reads.expected"             "reads.observed"            
#     [13] "reads.ratio"                "start.p"                    "end.p"                      "nexons"                     "AnnotSV.type"               "Gene.name"                 
#     [19] "NM"                         "CDS.length"                 "tx.length"                  "location"                   "intersectStart"             "intersectEnd"              
#     [25] "DGV_GAIN_IDs"               "DGV_GAIN_n_samples_with_SV" "DGV_GAIN_n_samples_tested"  "DGV_GAIN_Frequency"         "DGV_LOSS_IDs"               "DGV_LOSS_n_samples_with_SV"
#     [31] "DGV_LOSS_n_samples_tested"  "DGV_LOSS_Frequency"         "DDD_SV"                     "DDD_DUP_n_samples_with_SV"  "DDD_DUP_Frequency"          "DDD_DEL_n_samples_with_SV" 
#     [37] "DDD_DEL_Frequency"          "X1000g_event"               "X1000g_AF"                  "X1000g_max_AF"              "promoters"                  "dbVar_event"               
#     [43] "dbVar_variant"              "dbVar_status"               "TADcoordinates"             "ENCODEexperiments"          "GCcontent_left"             "GCcontent_right"           
#     [49] "Repeats_coord_left"         "Repeats_type_left"          "Repeats_coord_right"        "Repeats_type_right"         "ACMG"                       "HI_CGscore"                
#     [55] "TriS_CGscore"               "DDD_status"                 "DDD_mode"                   "DDD_consequence"            "DDD_disease"                "DDD_pmids"                 
#     [61] "HI_DDDpercent"              "synZ_ExAC"                  "misZ_ExAC"                  "pLI_ExAC"                   "delZ_ExAC"                  "dupZ_ExAC"                 
#     [67] "cnvZ_ExAC"                  "morbidGenes"                "morbidGenesCandidates"      "Mim.Number"                 "Phenotypes"                 "Inheritance"               
#     [73] "AnnotSV.ranking" 

