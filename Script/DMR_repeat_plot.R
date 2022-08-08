library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(GenomicRanges);
library(pheatmap);
library(dendextend);
setwd('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/clean');
### enrichment of downregulated regions in repeat
dmr <- readRDS("/cluster/projects/hansengroup/wye/5mc/Paired_End/BWA/Sort_Bamfiles_redo/DMR_DESeq2/DMR_VPC67_vs_CPC30_with_age_control_result_deseq.Rdata");
region <- readRDS('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/2021-02-05_hg19_region_337420.rds');
rep <- read.table('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/region_337420_ol_repeat.bed');
dmr$rep <- ifelse(rownames(dmr)%in%rep$V6, 'y', 'n');
mydn <- as.matrix(table(dmr$rep, ifelse(dmr$status=='Down regulated', 'y', 'n')));
fisher.test(mydn)$p.value;
fisher.test(mydn)$estimate;
(mydn[2,2]/mydn[1,2])/(mydn[2,1]/mydn[1,1])
myup <- as.matrix(table(dmr$rep, ifelse(dmr$status=='Up regulated', 'y', 'n')));
fisher.test(myup)$p.value;
fisher.test(myup)$estimate;
pdf(generate.filename('contingency_table', 'dn_reg', 'pdf'));
pheatmap(mydn, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE, 
	color = colorRampPalette(c('white', 'blue'))(100), 
	display_numbers = TRUE, numbers_format = 's', number_color = 'black', fontsize = 50);
dev.off();
#### enrichment of downregulated peaks in repeat
samp <- readRDS('/cluster/projects/hansengroup/sujunc/methylation/5mC//sample/2021-02-03_clinical_vpc_all_clean.rds')
rpkm.vc <- readRDS('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/2021-01-26_peak_vpc_cpc_filtered.rds');
gr.peak <- readRDS('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/2021-01-26_grange_peaks.rds');

peri <- read.table('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/hg19_centromere_1mb.bed');
colnames(peri) <- c('chr', 'start', 'end');
peri <- makeGRangesFromDataFrame(peri);
ov <- findOverlaps(gr.peak, peri);
ov <- as.data.frame(ov)
colnames(ov) <- c("from","to")
ov <- cbind(ov,gr.peak[ov$from, ],peri[ov$to, ]);
ov.peri.peak <- ov;
####
rep <- read.table('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/clean/peak_838649_ol_repeat.bed');
dmrD <- readRDS('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/2021-01-26_diff_peaks_filtered.rds');
dmr <- dmrD;
dmr$rep <- ifelse(rownames(dmr)%in%rep$V6, 'y', 'n');
dmr$status <- 'non';
dmr[dmr$fdr<0.05&dmr$lfc>0, ]$status <- 'up';
dmr[dmr$fdr<0.05&dmr$lfc<0, ]$status <- 'dn';
dmr$peri <- ifelse(rownames(dmr)%in%ov.peri.peak$id, 'y', 'n');
dmrD <- dmr;
mydn <- as.matrix(table(dmr$rep, ifelse(dmr$status=='dn', 'y', 'n')));
fisher.test(mydn)$p.value;
fisher.test(mydn)$estimate;
(mydn[2,2]/mydn[1,2])/(mydn[2,1]/mydn[1,1])
myup <- as.matrix(table(dmr$rep, ifelse(dmr$status=='up', 'y', 'n')));
fisher.test(myup)$p.value;
fisher.test(myup)$estimate;
### 
plot_repeatType <- function(dmr, rep, myname){
	mytype <- as.data.frame(table(droplevels(rep[rep$V12%in%rownames(dmr[dmr$status=='dn', ]), ]$V4)));
	mytype <- mytype[order(-mytype$Freq), ];
	rownames(mytype) <- mytype$Var1;
	mytype$pos <- seq(nrow(mytype));
	to.plot <- mytype[1:20, ];
	to.plot$pos <- seq(nrow(to.plot), 1)
	create.barplot(
		file = generate.filename('RepType_peak_dn', myname, 'pdf'),
		formula = pos~Freq,
		data = to.plot,
		stack = FALSE,
		plot.horizontal = TRUE,
		yaxis.lab = rev(rownames(to.plot)),
		ylab.label = 'Type',
		xlab.label = 'Frequncy',
		xlimits = c(0, max(to.plot)),
		#xat = seq(0, 100, 20),
		#yaxis.rot = 90,
		style = 'Nature',
		width = 6,
		height = 8,
		border.col = 'white',
		#box.ratio = 8
		);
};
rep <- read.table('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/clean/repeat_ol_peak_838649.bed');
plot_repeatType(dmrD, rep, 'discovery');
####
dmr <- dmrD;
mdat <- rpkm.vc[rownames(rpkm.vc)%in%rownames(dmr[dmr$peri=='y'&dmr$fdr<0.05, ]), ];
mdat <- data.frame(t(scale(t(mdat))));
ann_col <- data.frame(type = ifelse(grepl('V', colnames(mdat)), 'VPC', 'CPC'));
rownames(ann_col) <- colnames(mdat);
nbreaks <- unique(c(seq(-10, -1, 1), seq(-1, 1, 0.05), seq(1, 10, 1)));
ann_col <- ann_col[order.hclust(hclust(dist(t(mdat)))), ,drop = FALSE];
ann_col$cl <- seq(nrow(ann_col));
ann_col <- ann_col[order(ann_col$type, ann_col$cl), ];
ann_col <- ann_col[, 'type', drop = FALSE];
cl.row <- order.hclust(hclust(dist((mdat))));
mdat <- mdat[cl.row, rownames(ann_col)];
pdf(paste0(Sys.Date(), '_pericentr_diff_peak_', 'cpc_vpc', '.pdf'), width = 10);
pheatmap(mdat, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames = FALSE,
        color = colorRampPalette(c('blue', 'white', 'red'))(length(nbreaks)), breaks = nbreaks,
        annotation_col = ann_col)
dev.off();

tiff(paste0(Sys.Date(), '_pericentr_diff_peak_', 'cpc_vpc', '.tiff'), width = 600);
pheatmap(mdat, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames = FALSE,
        color = colorRampPalette(c('blue', 'white', 'red'))(length(nbreaks)), breaks = nbreaks,
        annotation_col = ann_col, legend = FALSE, annotation_legend = FALSE)
dev.off();

pdf(paste0(Sys.Date(), '_pericentr_diff_peak_', 'cpc_vpc_cl', '.pdf'), width = 10);
pheatmap(mdat, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames = FALSE,
        color = colorRampPalette(c('blue', 'white', 'red'))(length(nbreaks)), breaks = nbreaks,
        annotation_col = ann_col)
dev.off();
####
to.plot <- data.frame(score = colSums(mdat>0));
to.plot$mean.scale <- colMeans(mdat)
to.plot$group <- ifelse(grepl('V', rownames(to.plot)), 'VPC', 'CPC')
pval1 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='CPC', ]$mean.scale, to.plot[to.plot$group=='VPC', ]$mean.scale)$p.value);
create.boxplot(
	formula = mean.scale~group,
	data = to.plot,
	add.stripplot = TRUE,
	xlab.label = "",
	ylab.label = 'Mean',
	add.text = TRUE,
	main = '',
	main.cex = 1.5,
	text.x = c(1.5),
	text.y = max(to.plot$mean.scale),
	text.label = c(pval1),
	filename = generate.filename('peir_diff_cf', 'mean', 'pdf'),
	style = 'Nature',
	width = 5,
	height = 4
    );
write.csv(to.plot[, -1], 'data_F2J.csv');
###
pval1 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='CPC', ]$score, to.plot[to.plot$group=='VPC', ]$score)$p.value);
create.boxplot(
	formula = score~group,
	data = to.plot,
	add.stripplot = TRUE,
	xlab.label = "",
	ylab.label = 'Score',
	add.text = TRUE,
	main = '',
	main.cex = 1.5,
	text.x = c(1.5),
	text.y = max(to.plot$score),
	text.label = c(pval1),
	filename = generate.filename('peir_diff_cf', 'score', 'pdf'),
	style = 'Nature',
	width = 5,
	height = 4
    );

pval1 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='Benign', ]$mean.scale, to.plot[to.plot$group=='Tumor', ]$mean.scale)$p.value);
create.boxplot(
	formula = mean.scale~group,
	data = to.plot,
	add.stripplot = TRUE,
	xlab.label = "",
	ylab.label = 'Mean',
	add.text = TRUE,
	main = '',
	main.cex = 1.5,
	text.x = c(1.5),
	text.y = max(to.plot$mean.scale),
	text.label = c(pval1),
	filename = generate.filename('peir_diff_cf', 'Mean', 'pdf'),
	style = 'Nature',
	width = 5,
	height = 4
    );
samp <- readRDS('/cluster/projects/hansengroup/sujunc/methylation/5mC/sample/2021-02-03_clinical_vpc_all_clean.rds');
samp1 <- readRDS('/cluster/projects/hansengroup/wye/5mc/OICR/Sort_Bamfiles_redo/Wiggle/txt/Data_QC_Report_Full_addMEDIPS_20210414.Rdata');
samp$age <- samp1[match(samp$sampleID, samp1$simpleName), ]$age;
to.plot$ctDNA <- samp[rownames(to.plot), ]$ctDNA;
to.plot <- na.omit(to.plot);
to.plot$age <- samp[rownames(to.plot), ]$age;
create.scatterplot(
	formula = mean.scale~ctDNA,
	data = to.plot,
	ylab.label = 'Mean',
	xlab.label = 'ctDNA %',
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	style = 'Nature',
	xaxis.cex = 2,
	yaxis.cex = 2,
	filename = generate.filename('peri_mean', 'ctDNA_cf', 'pdf'),
	legend = list(
 	inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot$mean.scale,
                y = to.plot$ctDNA,
                label.items = c('pearson','pearson.p'),
                alpha.background = 0,
                key.cex = 1
                )
            ),
        x = 0.7,
        y = 0.98,
        corner = c(0,1)
        )
    )
	);

####
peak.ch <- readRDS('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/clean/Changhai_hg19_peaks/2021-06-02_changhai_hg19_peaks.rds');
peak.uc <- readRDS('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/clean/Felix_hg19_peaks/2021-06-02_felix_hg19_peaks.rds');
mdat <- peak.ch;
mdat <- mdat[rownames(mdat)%in%rownames(dmr[dmr$peri=='y'&dmr$fdr<0.05, ]), ];
mdat <- na.omit(data.frame(t(scale(t(mdat)))));
ann_col <- data.frame(type = ifelse(grepl('N', colnames(mdat)), 'Benign', 'Tumor'));
rownames(ann_col) <- colnames(mdat);
nbreaks <- unique(c(seq(-20, -1, 1), seq(-1, 1, 0.05), seq(1, 20, 1)));
ann_col <- ann_col[order.hclust(hclust(dist(t(mdat)))), ,drop = FALSE];
ann_col$cl <- seq(nrow(ann_col));
ann_col <- ann_col[order(ann_col$type, ann_col$cl), ];
ann_col <- ann_col[, 'type', drop = FALSE];
cl.row <- order.hclust(hclust(dist((mdat))));
mdat <- mdat[cl.row, rownames(ann_col)];
pdf(paste0(Sys.Date(), '_pericentr_diff_peak_', 'ch', '.pdf'), width = 10);
pheatmap(mdat, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames = FALSE,
        color = colorRampPalette(c('blue', 'white', 'red'))(length(nbreaks)), breaks = nbreaks,
        annotation_col = ann_col)
dev.off();

colnames(ann_col) <- ' '
tiff(paste0(Sys.Date(), '_pericentr_diff_peak_', 'ch', '.tiff'), width = 600);
pheatmap(mdat, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames = FALSE,
        color = colorRampPalette(c('blue', 'white', 'red'))(length(nbreaks)), breaks = nbreaks,
        annotation_col = ann_col, legend = FALSE, annotation_legend = FALSE)
dev.off();

to.plot <- data.frame(score = colSums(mdat>0));
to.plot$mean.scale <- colMeans(mdat)
to.plot$group <- ifelse(grepl('N', rownames(to.plot)), 'Benign', 'Tumor');
pval1 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='Benign', ]$score, to.plot[to.plot$group=='Tumor', ]$score)$p.value);
create.boxplot(
	formula = score~group,
	data = to.plot,
	add.stripplot = TRUE,
	xlab.label = "",
	ylab.label = 'Score',
	add.text = TRUE,
	main = '',
	main.cex = 1.5,
	text.x = c(1.5),
	text.y = max(to.plot$score),
	text.label = c(pval1),
	filename = generate.filename('peir_diff_ch', 'score', 'pdf'),
	style = 'Nature',
	width = 5,
	height = 4
    );

####
mdat <- cbind(peak.ch, peak.uc);
mdat <- mdat[rownames(mdat)%in%rownames(dmr[dmr$peri=='y'&dmr$fdr<0.05, ]), ];
#mdat <- mdat[rownames(mdat)%in%rownames(dmr[dmr$peri=='y'&dmr$fdr<0.05&dmr$status=='dn', ]), ];
mdat <- na.omit(data.frame(t(scale(t(mdat)))));
ann_col <- data.frame(type = ifelse(grepl('N', colnames(mdat)), 'Benign', 'Primary'));
rownames(ann_col) <- colnames(mdat);
ann_col$type <- as.vector(ann_col$type);
ann_col[grepl('DTB', rownames(ann_col)), ] <- 'Mets';
ann_col$type <- factor(ann_col$type, levels = c('Benign', 'Primary', 'Mets'));
nbreaks <- unique(c(seq(-23, -1, 1), seq(-1, 1, 0.05), seq(1, 23, 1)));
ann_col <- ann_col[order.hclust(hclust(dist(t(mdat)))), ,drop = FALSE];
ann_col$cl <- seq(nrow(ann_col));
ann_col <- ann_col[order(ann_col$type, ann_col$cl), ];
ann_col <- ann_col[, 'type', drop = FALSE];
cl.row <- order.hclust(hclust(dist((mdat))));
mdat <- mdat[cl.row, rownames(ann_col)];
pdf(paste0(Sys.Date(), '_pericentr_diff_peak_', 'ch_uc', '.pdf'), width = 10);
pheatmap(mdat, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames = FALSE,
        color = colorRampPalette(c('blue', 'white', 'red'))(length(nbreaks)), breaks = nbreaks,
        annotation_col = ann_col)
dev.off();

colnames(ann_col) <- ' '
tiff(paste0(Sys.Date(), '_pericentr_diff_peak_', 'ch', '.tiff'), width = 600);
pheatmap(mdat, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames = FALSE,
        color = colorRampPalette(c('blue', 'white', 'red'))(length(nbreaks)), breaks = nbreaks,
        annotation_col = ann_col, legend = FALSE, annotation_legend = FALSE)
dev.off();

to.plot <- data.frame(score = colSums(mdat>0));
to.plot$mean.scale <- colMeans(mdat)
to.plot$group <- ann_col[rownames(to.plot), ];
pval1 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='Benign', ]$score, to.plot[to.plot$group=='Primary', ]$score)$p.value);
pval2 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='Primary', ]$score, to.plot[to.plot$group=='Mets', ]$score)$p.value);
create.boxplot(
	formula = score~group,
	data = to.plot,
	add.stripplot = TRUE,
	xlab.label = "",
	ylab.label = 'Score',
	add.text = TRUE,
	main = '',
	main.cex = 1.5,
	text.x = c(1.5, 2.5),
	text.y = max(to.plot$score),
	text.label = c(pval1, pval2),
	filename = generate.filename('peir_diff_ch_uc', 'score', 'pdf'),
	style = 'Nature',
	width = 5,
	height = 4
    );

pval1 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='Benign', ]$mean.scale, to.plot[to.plot$group=='Primary', ]$mean.scale)$p.value);
pval2 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='Primary', ]$mean.scale, to.plot[to.plot$group=='Mets', ]$mean.scale)$p.value);
create.boxplot(
	formula = mean.scale~group,
	data = to.plot,
	add.stripplot = TRUE,
	xlab.label = "",
	ylab.label = 'Mean',
	add.text = TRUE,
	main = '',
	main.cex = 1.5,
	text.x = c(1.5, 2.5),
	text.y = max(to.plot$mean.scale),
	text.label = c(pval1, pval2),
	filename = generate.filename('peir_diff_ch_uc', 'Mean', 'pdf'),
	style = 'Nature',
	width = 5,
	height = 4
    );
write.csv(to.plot[, -1], 'data_F2K.csv')

mdat <- cbind(peak.ch, peak.uc);
mdat <- mdat[rownames(mdat)%in%rownames(dmr[dmr$peri=='y'&dmr$fdr<0.05, ]), ];
group <- ann_col[, 1];
names(group) <- rownames(ann_col);
mydiff <- data.frame(tissue = apply(mdat+0.1, 1, function(x) log2(mean(x[group=='Mets'])/mean(x[group=='Primary']))));
mydiff$cf <- dmr[rownames(mydiff), ]$lfc;
create.scatterplot(
	formula = cf~tissue,
	data = mydiff,
	#ylab.label = '',
	#xlab.label = 'ctDNA %',
	xaxis.fontface = 'plain', 
	yaxis.fontface = 'plain',
	style = 'Nature',
	xaxis.cex = 2,
	yaxis.cex = 2,
	filename = generate.filename('peri_lfc', 'cf_tissue', 'pdf'),
	legend = list(
    	inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = mydiff$cf,
                y = mydiff$tissue,
                label.items = c('pearson','pearson.p'),
                alpha.background = 0,
                key.cex = 1
                )
            ),
        x = 0.7,
        y = 0.98,
        corner = c(0,1)
        )
    )
	);

