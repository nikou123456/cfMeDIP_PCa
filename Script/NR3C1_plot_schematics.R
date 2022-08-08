library(Gviz);
library(GenomicRanges);
library(BSgenome.Hsapiens.UCSC.hg19);
setwd('/Users/sujunchen/Document/labs/2020/5mc');
mysnp <- data.frame(chr = 'chr5', start = 142782300, end = 142782600);
mysnp <- makeGRangesFromDataFrame(mysnp);
atrack <- AnnotationTrack(mysnp, name = "SNP");
mygene <- read.table('gene_NR3C1.gtf', sep = '\t');
mygene <- mygene[grep('exon', mygene$V9), ];
mygene <- mygene[grepl('exon|utr', mygene$V3), ];
mygene$exon <- gsub('.*exon_id |exon_version.*|; |.*ccds_id |havana_transcript', '', mygene$V9);
mygene$transcript <- gsub('.*transcript_name |transcript_source.*|; ', '', mygene$V9);
mygene <- data.frame(chromosome = 'chr5', start = mygene$V4, end = mygene$V5, width = mygene$V5 - mygene$V4, 
	strand = '-', feature = mygene$V3, gene = 'NR3C1', exon = mygene$exon, transcript = mygene$transcript);
mygene$feature <- gsub('exon', 'protein_coding', '')
chr <- as.character(unique(seqnames(mysnp)));
gen <- 'hg19'
gtrack <- GenomeAxisTrack();
itrack <- IdeogramTrack(genome = gen, chromosome = chr);
grtrack <- GeneRegionTrack(mygene, genome = gen, chromosome = chr, 
                           name = "Gene Model",
                           transcriptAnnotation = "transcript",
                           background.title = "brown")
###
pdf(paste0(Sys.Date(), '_gene_NR3C1_transcript.pdf'), width = 4, height = 3);
plotTracks(list(itrack, gtrack, atrack, grtrack), col = NULL);
dev.off();

### plot regions
mysnp <- r1;
atrack <- AnnotationTrack(mysnp, name = "V011");
chr <- 'chr22';
gen <- 'hg19';
gtrack <- GenomeAxisTrack();
itrack <- IdeogramTrack(genome = gen, chromosome = chr);
###
pdf(paste0(Sys.Date(), '_region_V011.pdf'), width = 10, height = 6);
plotTracks(list(itrack, gtrack, atrack), col = NULL);
dev.off();

mysnp <- r1;
atrack <- AnnotationTrack(mysnp, name = "V011");
chr <- 'chr8';
gen <- 'hg19';
gtrack <- GenomeAxisTrack();
itrack <- IdeogramTrack(genome = gen, chromosome = chr);
###
pdf(paste0(Sys.Date(), '_region_V048.pdf'), width = 10, height = 6);
plotTracks(list(itrack, gtrack, atrack), col = NULL);
dev.off();

