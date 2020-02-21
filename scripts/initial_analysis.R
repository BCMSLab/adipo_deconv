library(CAMTHC)
library(reshape2)
library(tidyverse)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)
library("RColorBrewer")

# load data
load('data/adipo_counts.rda')

#n_samples <- table(adipo_counts$time) > 5
#ind <- adipo_counts$time %in% names(n_samples)[n_samples]
ind <- TRUE

# subset the data
mat <- assay(adipo_counts)[, ind]
time <- adipo_counts$time[ind]
stage <- adipo_counts$stage

# log transformation
mat2 <- log2(mat + 1)

# calculate pca
pc <- prcomp(t(mat2))

# plot pca
png(filename = 'manuscript/pca.png',
    height = 13, width = 13, units = 'cm', res = 300)

col <- c('gray', 'blue', 'black', 'red')[as.factor(stage)]

plot(pc$x[, 1], pc$x[, 2],
     pch = 19,
     col = col,
     cex = .7,
     xlab = 'PC One', ylab = 'PC Two')

legend('top',
       legend = c(0, 1, 2, 3),
       col = c('gray', 'blue', 'black', 'red'),
       pch = 16,
       horiz = TRUE,
       xpd = TRUE,
       box.lwd = 0,
       inset = c(0, -.2),
       title = 'Differentiation Stage')

dev.off()

# set a seed
set.seed(111)

# run deconvolution
rCAM <- CAM(mat,
            K = 2:5,
            thres.low = 0.30,
            thres.high = 0.95,
            cores = 4)

# plot mdl
plot(MDL(rCAM), data.term = TRUE)

# extract estimate matrices and markers
n <- 3
aset <- Amat(rCAM, n)
sset <- Smat(rCAM, n)
markers <- MGsforA(rCAM, n)

aset %>%
    melt() %>%
    as_tibble() %>%
    mutate(time = rep(time, n)) %>%
    ggplot(aes(x = (time), y = value, group = as.factor(Var2), color = as.factor(Var2))) +
    geom_jitter() +
    geom_smooth()

# plot fraction of cells by time
png(filename = 'manuscript/fraction_time.png', height = 13, width = 13, units = 'cm', res = 300)

splines <- apply(aset, 2, function(y) smooth.spline(time, y, spar = .8))

plot(c(min(time), max(time)), 0:1,
     type = 'n',
     xlab = 'Time Point (hour)',
     ylab = 'Population Fraction')

points(time, aset[, 1], col = 'darkgreen', pch = 19, cex = .5)
lines(splines[[1]], lwd = 2, col = 'darkgreen')

points(time, aset[, 2], col = 'darkred', pch = 19, cex = .5)
lines(splines[[2]], lwd = 2, col = 'darkred')

points(time, aset[, 3], col = 'black', pch = 19, cex = .5)
lines(splines[[3]], lwd = 2, col = 'black')

legend('top',
       legend = c(1, 2, 3),
       col = c('darkgreen', 'darkred', 'black'),
       pch = 16,
       horiz = TRUE,
       xpd = TRUE,
       box.lwd = 0,
       inset = c(0, -.2),
       title = 'Sub-population Group')

dev.off()

# simplexplot
png(filename = 'manuscript/simplexplot.png', height = 13, width = 13, units = 'cm', res = 300)

debCAM::simplexplot(mat,
                    aset,
                    markers,
                    pch = 19,
                    mg.pch = 19,
                    mg.col = c('darkgreen', 'darkred', 'black'),
                    cex = .6,
                    mg.cex = .7,
                    xlab = 'Dimension One',
                    ylab = 'Dimension Two')

legend('top',
       legend = c(1, 2, 3),
       col = c('darkgreen', 'darkred', 'black'),
       pch = 16,
       horiz = TRUE,
       xpd = TRUE,
       box.lwd = 0,
       inset = c(0, -.2),
       title = 'Sub-population Markers')
dev.off()

# perform enrichment analysis on the markers
# extract annotations
symbols <- keys(org.Mm.eg.db, 'SYMBOL')
term2gene <- AnnotationDbi::select(org.Mm.eg.db,
                    symbols,
                    'GO',
                    'SYMBOL') %>%
    dplyr::select(term = GO, gene = SYMBOL) %>%
    unique()

terms <- as.data.frame(GOTERM)
terms <- tibble(ID = terms$go_id,
                term = terms$Term,
                definition = terms$Definition,
                ontology = terms$Ontology)

markers_symbols <- map(markers, intersect, y = unique(term2gene$gene))

# run enrichment
comp <- map_df(markers_symbols,
            function(x) {
                enricher(x,
                         TERM2GENE = term2gene,
                         pAdjustMethod = 'fdr')@result
            }, .id = 'population')

# merge annotation and enrichemnt
inner_join(terms, comp) %>%
    dplyr::filter(p.adjust < .2, Count > 2) %>%
    unique() %>%
    group_by(population, ontology) %>%
    dplyr::slice(1:10) %>%
    write_csv('manuscript/marker_enrichment.csv')

# plot the adipogenic markers
adipo_markers <- list('Adipogenic' = c('Pparg', 'Cebpb', 'Cebpa'),
                      'Lipogenic' = c('Lpl', 'Acly', 'Fasn'))

adipo_mat <- sset[rownames(sset) %in% unlist(adipo_markers),]

png(filename = 'manuscript/adipogenic_markers.png', height = 13, width = 13, units = 'cm', res = 300)

heatmap(adipo_mat)

dev.off()

aset2 <- aset + .001
sset2 <- sset + .001

markers_stats <- MGstatistic(mat,
                             aset,
                             boot.alpha = 0.05,
                             nboot = 1000,
                             cores = 4)

(markers_stats %>%
    rownames_to_column('gene') %>%
    dplyr::filter(gene %in% unlist(markers), OVE.FC != Inf, !grepl('Rik', gene)) %>%
    group_by(idx) %>%
    slice(1:5) %>%
    ggplot(aes(x = gene, y = log2(OVE.FC))) +
    geom_col() +
    scale_y_continuous(limits = c(0,6.2), expand = c(0, .1)) +
    facet_wrap(~idx, scales = 'free_x') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0, 'cm')) +
    labs(x = '', y = 'Fold-change (log_2)')) %>%
    ggsave(plot = .,
           filename = 'manuscript/markers_fc.png',
           height = 7, width = 10, units = 'cm')

# primary adipocytes
primary_se <- read_rds('primary_adipocyte.rds')

primary_mat <- 2 ^ na.omit(exprs(primary_se))
group <- primary_se$group
patient <- primary_se$patient
unique(patient)

new_markers <- map(markers, function(x) intersect(toupper(x), rownames(primary_mat)))

rre <- redoASest(primary_mat, new_markers, maxIter = 10)

df <- as.data.frame(rre$Aest) %>%
    mutate(group = factor(group, levels = c('NDf', 'Df'))) %>%
    gather(pop, frac, -group)
ave <- df %>% group_by(pop, group) %>% summarise(frac = mean(frac))

(df %>% 
    ggplot() +
    geom_jitter(aes(x = group, y = frac, color = pop, group = pop),
                width = .1) +
    geom_point(data = ave, 
               aes(x = group, y = frac, color = pop, group = pop)) +
    geom_line(data = ave, aes(x = group, y = frac, color = pop, group = pop)) +
    theme_bw() +
    labs(x = '', y = 'Population Fraction', color = 'Sub-population Group') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'top') +
        scale_color_manual(values = c('darkgreen', 'darkred', 'black'))) %>%
    ggsave(plot = .,
           filename = 'manuscript/primary_adipocytes.png',
           height = 8, width = 8, units = 'cm')
    
