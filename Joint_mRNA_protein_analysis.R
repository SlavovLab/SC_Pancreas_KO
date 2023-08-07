#Source functions for integration
source('/Users/andrewleduc/Desktop/Current_projects/Senescense/Trachea/Joint_mRNA_prot_functions.R')

# Read in protein data

convert_mouse <- Proc_fasta('/Users/andrewleduc/Desktop/Current_projects/Senescense/Mouse.fasta')

protein_norm_imp <- read.prot_gene('/Users/andrewleduc/Desktop/Current_projects/Juan/prot_imp.csv',convert_mouse)
protein_norm_noimp <- read.prot_gene('/Users/andrewleduc/Desktop/Current_projects/Juan/prot_noimp.csv',convert_mouse)

prot_meta <- read.csv('/Users/andrewleduc/Desktop/Current_projects/Senescense/Trachea/Prot_3_7/prot_meta.csv',row.names = 1)

#  Read in mRNA data

Panc <- readRDS("/Users/andrewleduc/Desktop/Current_projects/Senescense/Trachea/Trachea_3_7.rds")
Panc <-juan


#### Filtering mRNA and protein data correctly + correlation vector analysis ####

counts <- as.matrix(Panc@assays$RNA@counts)
counts <- counts[intersect(rownames(protein_norm_imp),rownames(counts)),]

counts_filt <- counts
counts_filt[counts_filt==0] <- NA

counts_filt <- filt.mat.cr(counts_filt,.98,.98)#Filter_missingness(counts_filt,clust)
protein_norm_noimp <- filt.mat.cr(protein_norm_noimp,.95,.95)#Filter_missingness(protein_norm_noimp,clust)

intersect_genes <- intersect(rownames(protein_norm_noimp),rownames(counts_filt))

# filter for intersected
protein_norm_noimp <- protein_norm_noimp[intersect_genes,]
protein_norm_imp <- protein_norm_imp[intersect_genes,]
counts <- counts[intersect_genes,]

# make new scaledata
Panc <- ScaleData(Panc, do.scale = TRUE ,features = intersect_genes)
scale_data <- Panc@assays[["RNA"]]@scale.data
scale_data <- scale_data[intersect_genes,]


imputed_cor_vect <- GetCorVect(protein_norm_imp,scale_data)
imputed_cor_vect_null <- GetCorVect(protein_norm_imp[sample(rownames(protein_norm_imp)),],scale_data)
unimputed_cor_vect <- GetCorVect(protein_norm_noimp,scale_data)


ggplot(imputed_cor_vect, aes(x = cors)) + geom_histogram()  +
  ylab('# of Genes') + xlab('Correlations') + ggtitle('Correlation vector comparison')
ggplot(unimputed_cor_vect, aes(x = cors)) + geom_histogram()
ggplot(imputed_cor_vect_null, aes(x = cors)) + geom_histogram()

unimputed_cor_vect_good <- unimputed_cor_vect %>% filter(cors > .25)
imputed_cor_vect_good <- imputed_cor_vect %>% filter(cors > .2)

keep_corV <- intersect(unimputed_cor_vect_good$genes,imputed_cor_vect_good$genes)

protein_norm_imp_filt <- protein_norm_imp[keep_corV,]
scale_dat_filt <- scale_data[ keep_corV,]
counts_filt <- counts[keep_corV,]


#### Integration plots showing mixing, cell type amounts per clust, and age per clust ####
counts_filt <- counts_filt[,exclude_immune]
integrated <- Integrate_liger(counts_filt,protein_norm_imp_filt,10)
integrated[[1]]
umap_liger <-integrated[[2]]
ggplot(umap_liger, aes(x = tsne1,y = tsne2, color = Dataset)) + geom_point()
umap_liger_prot <- umap_liger %>% filter(Dataset == 'prot')
umap_liger_rna <- umap_liger %>% filter(Dataset == 'mRNA')


umap_liger_prot <- umap_liger_prot[rownames(prot_meta),]
prot_meta$assign <- umap_liger_prot$Cluster


dot_plot <- theme_classic()+theme(plot.title = element_text(hjust = .5,size = 24),
                                  axis.title.x = element_text(size = 20),
                                  axis.title.y = element_text(size = 20),
                                  axis.text.x = element_text(size = 14),
                                  axis.text.y = element_text(size = 14),
                                  legend.title = element_text(size = 14),
                                  legend.position = c(.7, .2),
                                  legend.text = element_text(size = 14),
                                  legend.title.align=0.9)


prot_plot <- prot_no_imp['P55095',]
prot_meta$prot <- prot_plot

# Display 
ggscatter(prot_meta,color = 'prot', x ="UMAP_1", y = "UMAP_2" ,  size = 5, alpha=0.5) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red",name = '')

ggplot(prot_meta, aes(x = UMAP_1,y = UMAP_2, color = Cell_type)) + geom_point() + dot_plot

Ppy <- ggplot(prot_meta, aes(x = UMAP_1,y = UMAP_2, color = prot)) + geom_point() + dot_plot+
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                      high = "red",name = 'log2(Protein Abs.)')+ggtitle('Ppy')

Ins1 <- ggplot(prot_meta, aes(x = UMAP_1,y = UMAP_2, color = prot)) + geom_point() + dot_plot+
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red",name = 'log2(Protein Abs.)')+ggtitle('Ins1')

Gcg <- ggplot(prot_meta, aes(x = UMAP_1,y = UMAP_2, color = prot)) + geom_point() + dot_plot+
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red",name = 'log2(Protein Abs.)')+ggtitle('Gcg')


(Ppy + Gcg) / (Ins1+Ins2)

prot_meta$Cell_type[prot_meta$UMAP_1 > 5] <- 'Beta 2'
prot_meta$Cell_type[prot_meta$UMAP_1 < 5 & prot_meta$UMAP_2 > -5] <- 'Beta 1'

prot_meta$Cell_type[prot_meta$Cell_type == 'Endothelial'] <- 'Alpha'

#init rna umap
um_plot_mrna <- as.data.frame(Panc@reductions[["umap"]]@cell.embeddings)
um_plot_mrna <- um_plot_mrna %>% filter(rownames(um_plot_mrna) %in% exclude_immune)
um_plot_mrna$cond <- Panc@meta.data[["orig.ident"]]
um_plot_mrna$cell <- colnames(counts_filt) #== rownames(um_plot_mrna) 
umap_liger_rna <- umap_liger_rna[um_plot_mrna$cell,]
um_plot_mrna$assign <- umap_liger_rna$Cluster

ggplot(um_plot_mrna, aes(x = UMAP_1,y = UMAP_2, color = assign)) + geom_point()

# Map annotations to clusters

clusters <- c(0,1,4,5,6)
celltypes <- c('Basal','Fibroblast1','Fibroblast2','Goblet','Immune','Chondrocytes','Ciliated','Endothelial')
prot_meta$Cell_type <- NA
prot_meta$Cell_type[prot_meta$assign == 0] <- 'Basal 1'
prot_meta$Cell_type[prot_meta$assign == 1] <- 'Basal 2'
prot_meta$Cell_type[prot_meta$assign == 4] <- 'Delta'
prot_meta$Cell_type[prot_meta$assign == 5] <- 'Alpha'
prot_meta$Cell_type[prot_meta$assign == 6] <- 'Endothelial'







