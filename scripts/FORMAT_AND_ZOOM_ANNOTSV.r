# Eric Minikel
# script to run ExomeDepth
# how to run:
# runExomeDepth.r -b bamlist.txt -o /output/path/ -v > runExomeDepthOutput.txt

start_time = Sys.time()

require(optparse) # http://cran.r-project.org/web/packages/optparse/optparse.pdf

options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("-a", "--annotSVtemp"), action="store", default='', 
    type='character', help="interim annotSV output"),
  make_option(c("-f", "--zoomfile"), action="store", default='', 
    type='character', help="zoom gene list if given"),
  make_option(c("-z", "--zoomname"), action="store", default='', 
    type='character', help="path to reference panel rda file"),
  make_option(c("-n", "--smTag"), action="store", default='', 
    type='character', help="sample name"),
  make_option(c("-o", "--outdir"), action="store", default='./',
    type='character', help="Output directory [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
    help="Print verbose output [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
    help="Do not print verbose output (this is the default)")
)
opt = parse_args(OptionParser(option_list=option_list))

print(opt,file=stderr())

##############################################################################################
# set working directory to output directory ($CORE_PATH/$PROJECT/$FAMILY/$SAMPLE/CNV_OUTPUT) #
##############################################################################################

  setwd(opt$outdir)

####

   annotSV <- read.table(opt$annotSVtemp, header=TRUE, stringsAsFactors = FALSE, sep="\t")
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

   annotSV_out_file <- paste(opt$smTag, ".annotSV.result.tsv", sep='')
   write.table(annotSV_out, file=annotSV_out_file, row.names=FALSE, quote=FALSE, sep="\t")

  # extract genes of interest
  
  if (file.exists(opt$zoomfile)) {
    geneList <- read.csv(opt$zoomfile,as.is=TRUE,header=TRUE)
    zoom <- annotSV_out[annotSV$Gene.name %in% geneList$Gene,]
    zoom_file <- paste(opt$smTag,".annotSV.result.zoom.", opt$zoomname , ".tsv" ,sep='')
    write.table(zoom, file=zoom_file, row.names=FALSE, quote=FALSE, sep="\t")
  } else {
    cat("No subset result output as no zoom gene list provided\n",file=stdout())
  }
