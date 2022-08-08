library(BoutrosLab.plotting.general);
library(tidyverse)
library(GenomicRanges)
library(cowplot)
setwd('/cluster/projects/hansengroup/sujunc/methylation/5mC/analysis/clean/fragmentation');
df.fr3 <- readRDS('bins_5mbcompartments_all.rds');
samp <- readRDS('/cluster/projects/hansengroup/wye/5mc/OICR/Sort_Bamfiles_redo/Wiggle/txt/Data_QC_Report_Full_addMEDIPS_20210414.Rdata');
to.plot <- data.frame(samp = df.fr3$id, group = df.fr3$group, bin = df.fr3$bin, ratio = df.fr3$ratio.corrected2);
to.plot <- to.plot[!to.plot$group%in%c('Benign', 'HC'), ];
to.plot$group <- droplevels(to.plot$group);
p1 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='Barrier', ]$ratio, to.plot[to.plot$group=='CPC', ]$ratio)$p.value)
p2 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='OCIR', ]$ratio, to.plot[to.plot$group=='CPC', ]$ratio)$p.value)
p3 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='VPC', ]$ratio, to.plot[to.plot$group=='CPC', ]$ratio)$p.value)
p4 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='WCDT', ]$ratio, to.plot[to.plot$group=='CPC', ]$ratio)$p.value)
create.boxplot(
  data = to.plot, 
  formula = ratio~group,
  add.stripplot = TRUE,
  xlab.label = 'Group',
  ylab.label = 'Ratio',
  add.text = TRUE,
  text.x = c(1,3,4,5),
  text.y = 0.8,
  text.label = c(p1, p2, p3, p4),
  filename = generate.filename('ratio', 'all', 'pdf'),
  style = 'Nature'
  );
to.plot$ctDNA <- samp[match(to.plot$samp, samp$simpleName), ]$ctDNA;

df.fr3 <- df.fr3 %>% group_by(id) %>% mutate(ratio.mean = mean(ratio.corrected2, na.rm = TRUE));
to.plot <- data.frame(samp = df.fr3$id, group = df.fr3$group, mean = df.fr3$ratio.mean);
to.plot <- to.plot[!to.plot$group%in%c('Benign', 'HC', 'OCIR'), ];
to.plot$group <- droplevels(to.plot$group);
to.plot$group <- factor(to.plot$group, levels = c('CPC', 'Barrier', 'VPC', 'WCDT'));
to.plot <- unique(to.plot);
p1 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='Barrier', ]$mean, to.plot[to.plot$group=='CPC', ]$mean)$p.value);
#p2 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='OCIR', ]$mean, to.plot[to.plot$group=='CPC', ]$mean)$p.value)
p3 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='VPC', ]$mean, to.plot[to.plot$group=='CPC', ]$mean)$p.value);
p4 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='WCDT', ]$mean, to.plot[to.plot$group=='CPC', ]$mean)$p.value);
mycol <- c('#1b9e77', '#d95f02', '#7570b3', '#a6761d');
names(mycol) <- levels(to.plot$group);
to.plot$col <- mycol[to.plot$group];
create.boxplot(
  data = to.plot, 
  formula = mean~group,
  add.stripplot = TRUE,
  ylimits = c(0, max(to.plot$mean)+0.05),
  xlab.label = '',
  ylab.label = 'Mean ratio',
  add.text = TRUE,
  border.col = mycol,
  points.col = to.plot$col,
  text.x = c(2,3,4),
  text.y = 0.5,
  text.label = c(p1, p3, p4),
  filename = generate.filename('ratio', 'mean', 'pdf'),
  style = 'Nature'
  );
####
df.fr3 <- df.fr3 %>% group_by(id) %>% mutate(sd = sd(ratio.corrected2, na.rm = TRUE));
to.plot <- data.frame(samp = df.fr3$id, group = df.fr3$group, mean = df.fr3$sd);
to.plot <- to.plot[!to.plot$group%in%c('Benign', 'HC', 'OCIR'), ];
to.plot$group <- droplevels(to.plot$group);
to.plot$group <- factor(to.plot$group, levels = c('CPC', 'Barrier', 'VPC', 'WCDT'));
to.plot <- unique(to.plot);
p1 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='Barrier', ]$mean, to.plot[to.plot$group=='CPC', ]$mean)$p.value);
#p2 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='OCIR', ]$mean, to.plot[to.plot$group=='CPC', ]$mean)$p.value)
p3 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='VPC', ]$mean, to.plot[to.plot$group=='CPC', ]$mean)$p.value);
p4 <- scientific.notation(wilcox.test(to.plot[to.plot$group=='WCDT', ]$mean, to.plot[to.plot$group=='CPC', ]$mean)$p.value);
mycol <- c('#1b9e77', '#d95f02', '#7570b3', '#a6761d');
names(mycol) <- levels(to.plot$group);
to.plot$col <- mycol[to.plot$group];
create.boxplot(
  data = to.plot, 
  formula = mean~group,
  add.stripplot = TRUE,
  ylimits = c(0, max(to.plot$mean)+0.005),
  xlab.label = '',
  ylab.label = 'Ratio standard deviation',
  add.text = TRUE,
  border.col = mycol,
  points.col = to.plot$col,
  text.x = c(2,3,4),
  text.y = 0.035,
  text.label = c(p1, p3, p4),
  filename = generate.filename('ratio', 'sd', 'pdf'),
  style = 'Nature'
  );
####
to.plot$ctDNA <- samp[match(to.plot$samp, samp$simpleName), ]$ctDNA;
to.plot <- na.omit(to.plot)
create.scatterplot(
  formula = mean~ctDNA,
  data = to.plot,
  ylab.label = 'SD',
  xlab.label = 'ctDNA %',
  xaxis.fontface = 'plain', 
  yaxis.fontface = 'plain',
  style = 'Nature',
  xaxis.cex = 2,
  yaxis.cex = 2,
  filename = generate.filename('SD', 'ctDNA', 'pdf'),
  legend = list(
    inside = list(
        fun = draw.key,
        args = list(
            key = get.corr.key(
                x = to.plot$mean,
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

