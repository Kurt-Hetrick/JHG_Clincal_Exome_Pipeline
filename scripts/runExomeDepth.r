# Eric Minikel
# script to run ExomeDepth
# how to run:
# runExomeDepth.r -b bamlist.txt -o /output/path/ -v > runExomeDepthOutput.txt

start_time = Sys.time()

require(optparse) # http://cran.r-project.org/web/packages/optparse/optparse.pdf
require(ExomeDepth)

options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("-b", "--bamfile"), action="store", default='', 
              type='character', help="Path to input bam file"),
  make_option(c("-r", "--refcounts"), action="store", default='', 
              type='character', help="path to reference panel rda file"),
  # $exomeDEPTH_BED
   make_option(c("-f", "--bedfile"), action="store", default='', 
               type='character', help="bed file"),
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

############################################################################################
# BED FILE is the full bed file, $exomeDEPTH_BED (that peng modified that exomeDepth uses) #
############################################################################################

  bed = read.table(opt$bedfile,header=TRUE, stringsAsFactors = FALSE)
  colnames(bed) = c("chromosome", "start", "end", "name")

##################################################################################################
# load in the reference panel rda file. patient gender determines which one is used ##############
# $REF_PANEL_FEMALE_READ_COUNT_RDA is used if patient is female ##################################
# $REF_PANEL_MALE_READ_COUNT_RDA is used if patient is male ######################################
# $REF_PANEL_ALL_READ_COUNT_RDA is used if patient gender is unknown #############################
# the three above are generalized/recoded as a $REF_PANEL_READ_COUNT_RDA depending on the gender #
##################################################################################################

  load(opt$refcounts)

##############################################################################################
# set working directory to output directory ($CORE_PATH/$PROJECT/$FAMILY/$SAMPLE/CNV_OUTPUT) #
##############################################################################################

  setwd(opt$outdir)

####################################################################################
# keep, but probaby won't use. this is for re-running exomeDepth ###################
# might be helpful while setting developing or when changes are added to save time #
####################################################################################

  # get the read count if it exist in .rda
  counts_file <- grep(".+\\.rda$", list.files(), value = TRUE)
  # file.rename(counts, "counts.rda")

#################################
# keep, same reasoning as above #
#################################

  # if there is a unique rda file, using the existing rda
  if (length(counts_file) == 1) {
    load(counts_file)
    cat(paste("Use the existing counts for this sample ", counts_file, "\n", sep=''),file=stderr())
  } else {

    # count the study sample
    counts = getBamCounts(bed.frame = bed, bam.files = opt$bamfile)
  #  save(counts,file=paste("counts", ".rda", sep=""))
    cat(paste("Calculated counts\n",sep=''),file=stdout())
  }

########################################################

# re-load the counts, 'ref_counts'

  ref_countdf = as.data.frame(refcounts)
  ref_countmat = as.matrix(ref_countdf[,6:dim(ref_countdf)[2], drop = FALSE])

# counts is an S4 object.
# you need to cast it to a data frame (for bin length)
# AND to a matrix (for reference.count)
# and for some reason you can't cast S4 directly to matrix, only via df
# reformats patient read count data into r data frame and then into a matrix

  countdf = as.data.frame(counts)
  countmat = as.matrix(countdf[,5:dim(countdf)[2], drop = FALSE]) # remove cols 1-5 metadata

#################################################################
#### this is creating the rda for the patient sample. #####
############################################################

   save(counts,file=paste(opt$smTag, ".rda", sep=""))

# for plotting

  min_pos <- aggregate(start ~ chromosome, bed, function(x) min(x))
  max_pos <- aggregate(end ~ chromosome, bed, function(x) max(x))

# keep in mind the min_pos and max_pos order is not numeric

  chr_order <- c("1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

# reorder chr in chr_order

  min_pos <- min_pos[match(chr_order, min_pos$chromosome),]
  max_pos <- max_pos[match(chr_order, max_pos$chromosome),]

# for (i in 1:dim(countmat)[2]) {
#  loop through each chromosome

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
      chromosome=countdf$chromosome, start=countdf$start,
      end=countdf$end, name=countdf$exon)
   
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

duration = format(Sys.time() - start_time)

# if(opt$verbose) {
    cat(paste("Completed execution in ",duration,"\n",sep=''),file=stdout())
# }
