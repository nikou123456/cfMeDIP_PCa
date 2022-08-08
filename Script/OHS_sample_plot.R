##R/3.6.1
library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(readxl);
setwd('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis');
select <- data.frame(read_excel('../sample/20211125_Prediciton_probability_matrix_based_on_methylation_feature.xlsx'));
samp <- read.table('/cluster/projects/hansengroup/sujunc/methylation/5mC/data/philip/predx.prad.sampleinfo.txt', sep = '\t', header = TRUE);
samp <- samp[samp$GRP_Id%in%gsub('-', '_', gsub('.*_', '', select$...1)), ];
samp <- samp[order(samp$Time_to_Diagnosis), ];
samp$pos <- seq(nrow(samp));
samp$col <- 'grey';
samp[samp$Diagnosis.Stage%in%c('I', 'II', 'IIA', 'IIB'), ]$col <- 'chartreuse4';
samp[samp$Diagnosis.Stage%in%c('III', 'IV'), ]$col <- 'darkorchid4';
samp$min <- 0;

create.segplot(
	formula = pos ~ min + Time_to_Diagnosis,
	data = samp,
	plot.horizontal = TRUE,
	centers = samp$Time_to_Diagnosis,
	segments.col = samp$col,
	ylab.label = 'Time to diagnosis (days)',
	xlab.label = 'Patient',
	xlab.cex = 1.5, ylab.cex = 1.5,
	xaxis.cex = 1.5, yaxis.cex = 1.5,
	xaxis.fontface = 'plain', yaxis.fontface = 'plain',
	xaxis.rot = 90,
	#xlimits = c(0, 1.0), xat = seq(0, 1.0, 0.2),
	#ylimits = c(0.5, nrow(plot.data)+0.5), yat = seq(1, nrow(plot.data), 1),
	#yaxis.lab = rownames(plot.data),
	#resolution = 600,
	style = 'Nature',
	filename = generate.filename('Time_to_Diagnosis', 'OHS30', 'pdf'),
	width = 4,
	height = 6,
	key = list(
	  text = list(
	  lab = c('Unkown', 'I/II', 'III/IV'),
	  cex = 1
	  ),
	points = list(
	  pch = 20,
	  col = c('grey', 'chartreuse4', 'darkorchid4')
	  ),
	x = 0.65, 
	y = 0.25
	)
	);
