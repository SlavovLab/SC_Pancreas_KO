# Functions in Github folder
source("/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/My Drive/MS/Collaborators/Juan_ap_al/Knock_out_Study/sc_pipeline-master 2/functions_parameters.R")


# Reference channel number (1-11, or 1-16)
ref_channel<-2

# Add your cell type labels, must match those used in experimental design
your_labels<-c("wt","ko","neg")
your_control_label<-"neg"

# Import ------------------------------------------------------------------

# Load raw data evidence file can be found on MASSIVE
ev<-read.delim('/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Collaborators/Juan_ap_al/Knock_out_Study/Panc_MQ_search/evidence.txt')


# Meta data are in github folder
design<-read.csv("/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/My Drive/MS/Collaborators/Juan_ap_al/Knock_out_Study/CR_KO_annotations.csv")

# Also in github folder, needed for code but no important info
batch<-read.csv("/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/My Drive/MS/Collaborators/Juan_ap_al/Knock_out_Study/sc_pipeline-master 2/batch.csv")












# Attach batch data to protein data
ev[,colnames(batch)[-1]]<-NA
for(X in batch$set){

  ev$lcbatch[ev$Raw.file==X] <- as.character(batch$lcbatch[batch$set%in%X])
  ev$sortday[ev$Raw.file==X] <- as.character(batch$sortday[batch$set%in%X])
  ev$digest[ev$Raw.file==X] <- as.character(batch$digest[batch$set%in%X])

}


  


  
parse_row<-grep("|",ev$Leading.razor.protein, fixed=T)

if(length(parse_row)>0){
  split_prot<-str_split(ev$Leading.razor.protein[parse_row], pattern = fixed("|"))
  split_prot2<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  ev$Leading.razor.protein[parse_row]<-split_prot2
}
# Create unique peptide+charge column:
ev$modseq<-paste0(ev$Modified.sequence,ev$Charge)

# Add X in front of experiment names because R doesn't like column names starting with numbers
ev$Raw.file<-paste0("X",ev$Raw.file)
design$Set<-paste0("X",design$Set)


# Which columns hold the TMT Reporter ion (RI) data
ri.index<-which(colnames(ev)%in%paste0("Reporter.intensity.",1:18))

# Make sure all runs are described in design, if not, print and remove them:
not.described<-unique(ev$Raw.file)[ !unique(ev$Raw.file) %in% paste0(design$Set) ]
ev<-ev[!ev$Raw.file%in%not.described,]

# Filter out reverse hits, contaminants, and contaminated spectra...
ev<-ev[-which(ev$Reverse=="+"),]
if(length(grep("REV", ev$Leading.razor.protein))>0){ ev<-ev[-grep("REV", ev$Leading.razor.protein),] }
if(length(grep("CON", ev$Leading.razor.protein))>0){ ev<-ev[-grep("CON", ev$Leading.razor.protein),] }
if(length(which(ev$Potential.contaminant=="+"))>0){ ev<-ev[-which(ev$Potential.contaminant=="+"),] }
ev<-ev[!is.na(ev$PIF),]
ev<-ev[ev$PIF>0.75,]

# Remove peptides that are more the 10% the intensity of the carrier in the single cell runs (only)
ev<-as.data.frame(ev)
ev$mrri<-0
ev$mrri <- rowMeans(ev[, ri.index[4:length(ri.index)]] / ev[, ri.index[1]], na.rm = T)
ev<-ev[ev$mrri < 0.1, ]


# Filter by PEP or FDR: CHOOSE ONE

ev<-ev[ev$PEP<0.07, ]



# Normalize single cell runs to normalization channel
ev<-as.data.frame(ev)
ev[, ri.index] <- ev[, ri.index] / ev[, ri.index[ref_channel]]


# Organize data into a more convenient data structure:
# Create empty data frame
ev.melt<-melt(ev[0, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(ev)[ri.index]) ],
              id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest"))

colnames(ev.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")


# Record mapping of cell type to Channel:
ct.v<-c()
qt.v<-c()

# Create a unique ID string
unique.id.numeric<-1:length(ri.index)
unique.id<-paste0("i",unique.id.numeric)

RI_keep<-ri.index

# Give each sample a unique identifier
for(X in unique(ev$Raw.file)){

  # Subset data by X'th experiment
  ev.t<-ev[ev$Raw.file%in%X, ]

  # Name the RI columns by what sample type they are: carrier, single cell, unused, etc...
  colnames(ev.t)[ri.index]<-paste0(as.character(unlist(design[design$Set==X,-1])),"-", unique.id)

    # Melt it! and combine with other experimental sets
    ev.t.melt<-melt(ev.t[, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(ev.t)[RI_keep]) ],
                    id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest"));

    # Record mapping of cell type to Channel:
    ct.v<-c(ct.v, unique.id[which(ri.index%in%RI_keep)] )
    qt.v<-c(qt.v, colnames(ev)[RI_keep] )

    colnames(ev.t.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")

    ev.melt<-rbind(ev.melt, ev.t.melt)

  # Update unique ID string
  unique.id.numeric<-unique.id.numeric + length(ri.index)
  unique.id<-paste0("i", unique.id.numeric)

}

c2q<-data.frame(ct.v, qt.v); colnames(c2q)<-c("id","channel")

# Grab the unique number associate to each and every cell, carrier channel, and empty channel
ev.melt$id<-unlist(strsplit(as.character(ev.melt$celltype),"-"))[seq(2,length(unlist(strsplit(as.character(ev.melt$celltype),"-"))),2)]
ev.melt$celltype<-unlist(strsplit(as.character(ev.melt$celltype),"-"))[seq(1,length(unlist(strsplit(as.character(ev.melt$celltype),"-"))),2)]
ev.melt$id<-as.factor(ev.melt$id)

# Remove duplicate observations of peptides from a single experiment
ev.melt<-remove.duplicates(ev.melt,c("sequence","id") )
ev.melt<-ev.melt[!is.na(ev.melt$protein), ]

# Create additional meta data matrices
ev.melt.uniqueID<-remove.duplicates(ev.melt,"id")

c2q$id <- as.character(c2q$id)
ev.melt.uniqueID$id <- as.character(ev.melt.uniqueID$id)
ev.melt.uniqueID <- left_join(ev.melt.uniqueID,c2q, by = 'id')


ev.melt.pep<-remove.duplicates(ev.melt, c("sequence","protein") )

# Create data frame of peptides x cells, populated by quantitation
ev.unmelt<-dcast(ev.melt, sequence ~ id, value.var = "quantitation", fill=NA)

# Also create matrix of same shape
ev.matrix<-as.matrix(ev.unmelt[,-1]); row.names(ev.matrix)<-ev.unmelt$sequence

# Replace all 0s with NA
ev.matrix[ev.matrix==0]<-NA
ev.matrix[ev.matrix==Inf]<-NA
ev.matrix[ev.matrix==-Inf]<-NA



# Divide matrix into single cells (including intentional blanks) and carriers
sc_cols<-unique(ev.melt$id[(ev.melt$celltype%in%c(your_labels))])
ev.matrix.sc<-ev.matrix[, colnames(ev.matrix)%in%sc_cols]




# Filter single cells ----------------------------------------------------------------------


sc.melt<-ev.melt

xd<-as_tibble( sc.melt )

xd <- xd %>% group_by(id) %>% mutate(med_per_c = median(quantitation, na.rm=T)); length(unique(xd$id))

length(unique(xd$id))

xd$quantitation[(xd$quantitation)==Inf]<-NA
xd$quantitation[(xd$quantitation)==0]<-NA

xd <- xd %>% mutate_if(is.factor, as.character)

xd1 <- xd %>%
  group_by(id) %>%
  mutate(norm_q1 = quantitation / median(quantitation, na.rm=T))

xd2 <- xd1 %>%
  group_by(sequence, Raw.file) %>%
  mutate(norm_q = quantitation / mean(norm_q1, na.rm=T))

xd3<- xd2 %>%
  filter(celltype%in%c(your_labels))

xd4<- xd3 %>%
  group_by(protein, id) %>%
  mutate(cvq = cv(norm_q))

xd5<- xd4 %>%
  group_by(protein, id) %>%
  mutate(cvn = cvna(norm_q))

xd6<- xd5 %>%
  filter(cvn > 2)


xd7<-xd6 %>% group_by(id) %>% mutate(cvm=median(cvq, na.rm=T))

xdf<-xd7


# Filter out variable wells and controls
cvPar <- 0.43
sc_kept<-unique( xdf$id[xdf$celltype!=your_control_label & xdf$cvm < cvPar])
sc0_kept<-unique( xdf$id[xdf$celltype==your_control_label & xdf$cvm > cvPar])

# Which wells to keep
keep_these<-unique( xdf$id)

sc_total<-unique( xdf$id[xdf$celltype!=your_control_label])
sc0_total<-unique( xdf$id[xdf$celltype==your_control_label])
scrate<-round(length(sc_kept) / length(sc_total),2)*100

ev.matrix.sc.f<-ev.matrix.sc[,colnames(ev.matrix.sc)%in%sc_kept]; dim(ev.matrix.sc.f)
ev.matrix.sc.f[ev.matrix.sc.f==Inf]<-NA
ev.matrix.sc.f[ev.matrix.sc.f==-Inf]<-NA
ev.matrix.sc.f[ev.matrix.sc.f==0]<-NA

xdf$control<-"sc"
xdf$control[xdf$celltype==your_control_label]<-"ctl"

my_col3<-c( "black", "purple2")

# Plot!
ggplot(data=xdf, aes(x=cvm)) + geom_density(aes(fill=control, alpha=0.5), adjust=4) + theme_pubr() +
  scale_fill_manual(values=my_col3[c(1,2)]) +
  xlab("Quantification variability") + ylab("Density") + rremove("y.ticks") + rremove("y.text") +
  font("xylab", size=35) +
  font("x.text", size=30) +
  coord_cartesian(xlim=c(0,1))+
  #xlim(c(-0.15, 0.35)) +
  # annotate("text", x=0.27, y= 14, label=paste0(scrate,"% single cells passed"), size=8, color=my_col3[c(2)])+
  # annotate("text", x=0.27, y= 12.5, label=paste0(sc0rate,"% control wells passed"), size=8, color=my_col3[c(1)])+
  annotate("text", x=0.172, y= 14, label=paste0(length(sc_kept)," single cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.6, y= 12, label=paste0(length(sc0_kept)," control wells"), size=10, color=my_col3[c(1)])+
  annotate("text", x=0.6, y= 14, label=paste0(length(sc_total) -length(sc_kept)," single cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.165, y= 12, label=paste0(length(sc0_total) - length(sc0_kept)," control wells"), size=10, color=my_col3[c(1)])+
  #annotate("text", x=0.25, y= 3, label="Macrophage-like", size=6) +
  rremove("legend") + geom_vline(xintercept=0.43, lty=2, size=2, color="gray50")




# Data transformations ----------------------------------------------------



# Perform normalizations / transformations in multiple steps with visual sanity checks:
b.t<-"FD"
xlim.t<-c(-2,2)
par(mfrow=c(3,3))

# Original data, normalized to reference channel, filtered for failed wells:
t0<-ev.matrix.sc.f

hist(c(t0), breaks=b.t, xlim=xlim.t)

# Column then row normalize by median or mean (see source functions):
t1<-cr_norm(t0)
hist(c(t1), breaks=b.t, xlim=xlim.t)


# Filter for missing data:
t2<-filt.mat.rc(t1, .95, .95)
hist(c(t2), breaks=b.t, xlim=xlim.t)


# Log2 transform:
t3<-log2(t2)
t3[t3==Inf]<-NA
t3[t3==-Inf]<-NA
t3[t3==0]<-NA
hist(c(t3), breaks=b.t, xlim=xlim.t)


# # Collapse to protein level by median:
t3m<-data.frame(t3)
t3m$pep<-rownames(t3)
t3m$prot <- ev.melt.pep$protein[match(t3m$pep, ev.melt.pep$sequence)]
t3m<-melt(t3m, variable.names = c("pep", "prot"))
colnames(t3m) <-c("pep","prot","id","quantitation")
t3m2<- t3m %>% group_by(prot, id) %>% dplyr::summarise(qp = median(quantitation, na.rm=T))
t4m<-dcast(t3m2, prot ~ id, value.var = "qp", fill=NA)
t4<-as.matrix(t4m[,-1]); row.names(t4)<-t4m[,1]
hist(c(t4), breaks=b.t, xlim=xlim.t)

# Re-column and row normalize:
t4b<-cr_norm_log(t4)
hist(c(t4b), breaks=b.t, xlim=xlim.t)

# Assign to a final variable name:
ev.matrix.sc.f.n<-t4b

hist(colSums(is.na(t4b)==F))

## Impute single celldata
imp.input<-ev.matrix.sc.f.n
sc.imp <- hknn(imp.input, k.t)
t5<-sc.imp
sum(is.na(sc.imp))
dim(sc.imp)

sc.imp[(is.na(sc.imp))]<-0



# Define the batches and model:
batch.covs <- ev.melt.uniqueID$Raw.file[match(colnames(sc.imp), ev.melt.uniqueID$id)]
#RI Bias Batches
batch.covsRI <- ev.melt.uniqueID$channel[match(colnames(sc.imp), ev.melt.uniqueID$id)]

#Cell Type Info Labels
mod<-data.frame(ev.melt.uniqueID$celltype[match(colnames(sc.imp), ev.melt.uniqueID$id)]); colnames(mod)<-"celltype"
mod<-model.matrix(~as.factor(celltype), data=mod)
library(limma)
matrix.sc.batch <- removeBatchEffect(sc.imp, batch=factor(batch.covs), batch2=factor(batch.covsRI),design=mod)


t6<-cr_norm_log(matrix.sc.batch)

t7 <- t6
t7[is.na(t4b)==T] <- NA

# visual sanity checks post-imputation:
hist(c(t5), breaks=b.t, xlim=xlim.t)
hist(c(t6), breaks=b.t, xlim=xlim.t)

par(mfrow=c(1,1))





# PCA ------------------------------------------------------------------------

mat.sc.imp <- t6

X.m <-  mat.sc.imp

write.csv(t6,'/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Collaborators/Juan_ap_al/Knock_out_Study/Protein_imputed.csv')
write.csv(t7,'/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Collaborators/Juan_ap_al/Knock_out_Study/Protein_uniputed.csv')
write.csv(t3,'/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Collaborators/Juan_ap_al/Knock_out_Study/peptide_uniputed.csv')
write.csv(ev.melt.uniqueID,'/Users/andrewleduc/Desktop/meta.csv')


pca.imp.cor <- cor(X.m, use = 'pairwise.complete.obs',method = c('pearson'))

# PCA
sc.pca<-eigen(pca.imp.cor)
scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)

# Percent of variance explained by each principle component
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")

# Map meta data
pca.melt <- melt(scx); colnames(pca.melt)<-c("id","pc","value")
#ev.melt.uniqueID_lim<-ev.melt.uniqueID[,c(4:6,8,10:11)]

pca.melt <- left_join(pca.melt,ev.melt.uniqueID,by= 'id')
pca.display <- dcast(pca.melt, id + celltype + channel ~ pc, value.var = "value", fill=NA)

PCx<-"PC1"
PCy<-"PC2"


# Display 
ggscatter(pca.display,color = 'celltype', x =PCx, y = PCy ,  size = 5, alpha=0.5) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  annotate("text", x=0.05-0.02, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  annotate("text", x=0.062-0.03, y=-0.11, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8)




beta_cell <- pca.display %>% filter(PC2 < .05)
beta_cell <- pca.display %>% filter(PC1 < .01)



mat.sc.imp_ko <- as.data.frame(t7_filt) %>% dplyr::select(ko$id)
mat.sc.imp_wt <- as.data.frame(t7_filt) %>% dplyr::select(wt$id)


prots <- (rownames(mat.sc.imp_ko))
df <- as.data.frame(prots)
diff <- (rowMedians(as.matrix(mat.sc.imp_wt),na.rm = T)-rowMedians(as.matrix(mat.sc.imp_ko),na.rm = T))
df <- as.data.frame(diff)
df$time <- as.numeric(rownames(mat.sc.imp_wt))
df2 <- df
df2$time <- df2$time + 2*pi
df <- rbind(df,df2)
df$FC <- 'FC'
ggplot(df,aes(y = diff,x = FC)) +ggbeeswarm::geom_quasirandom()

ggplot(df,aes(x = time,y = diff)) + geom_point(size = 3) + ylim(c(-.5,.5)) + ylab('(average WT abundance) - (average KO abundance)') + 
  xlab('Cycle Time') #+ geom_vline(xintercept=2*pi, lty=2, size=1, color="gray50")
  #geom_vline(xintercept=7.23, lty=2, size=2, color="gray50")+
  ggtitle('WT cells show higher abundance of late phase proteins compared to KO')

diff_all <- (rowMeans(as.matrix(mat.sc.imp_wt_all),na.rm = T)-rowMeans(as.matrix(mat.sc.imp_ko_all),na.rm = T))
hist(diff_all,50)






library(Seurat)

# Need to go to paper to download seurat object in rna seq data
load(file ='/Users/andrewleduc/Downloads/intobj.rda')
intobj.lst <- intobj.lst[[1]]
intobj.lst = UpdateSeuratObject(object = intobj.lst)

DimPlot(intobj.lst,label = TRUE)

juan <- CreateSeuratObject(counts = intobj.lst@assays[["RNA"]]@counts, project = "juan")
rm(intobj.lst)

juan <- NormalizeData(juan, normalization.method = "LogNormalize", scale.factor = 10000)
juan <- FindVariableFeatures(juan, selection.method = "vst", nfeatures = 4000)
juan <- ScaleData(juan, features = juan@assays[["RNA"]]@var.features)
juan <- RunPCA(juan, features = VariableFeatures(object = juan))
DimPlot(juan, reduction = "pca")

juan <- FindNeighbors(juan, dims = 1:10)
juan <- FindClusters(juan, resolution = .25)

juan <- RunUMAP(juan, dims = 1:10)

DimPlot(juan, reduction = "umap")

FeaturePlot(juan,features = "Pyy")

length(intersect(convert_mouse$split_prot, rownames(t6)))
nrow(t6)



prot_umap <- CreateSeuratObject(counts = t6, project = "prot_mat")

prot_umap <- NormalizeData(prot_umap, normalization.method = "LogNormalize", scale.factor = 10000)
prot_umap@assays$RNA@scale.data <- t6

all.genes <- rownames(t6)
prot_umap <- ScaleData(prot_umap, features = all.genes)
prot_umap

prot_umap@assays$RNA@scale.data <- t6

prot_umap <- RunPCA(prot_umap, features = all.genes)

DimHeatmap(prot_umap, dims = 1:6, cells = 200, balanced = TRUE)

#DimPlot(prot_umap, reduction = "pca")

prot_umap <- FindNeighbors(prot_umap, dims = 1:6)
prot_umap <- FindClusters(prot_umap, resolution = 1.5)

prot_umap <- RunUMAP(prot_umap, dims = 1:4,n.neighbors = 15)


DimPlot(prot_umap, reduction = "umap") + ggtitle('UMAP colored by Louvain cluster')

prot_meta <- (as.data.frame(prot_umap@reductions[["umap"]]@cell.embeddings))
prot_meta$clust <- prot_umap@meta.data$seurat_clusters
prot_meta$Condition <- pca.display$celltype

t7_plot <- t7['P55095',]
prot_meta$prot <- t7_plot

# Display 
ggscatter(prot_meta,color = 'prot', x ="UMAP_1", y = "UMAP_2" ,  size = 5, alpha=0.5) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red",name = '')

DimPlot(juan, reduction = "umap") # 8, 9, 12, 13
FeaturePlot(juan,features = "Arf1") # Ppy and Sst for delta

                                      # Acot13 Hadh for beta 2

exclude_immune <- colnames(counts_filt)[!juan@meta.data$seurat_clusters %in% c(8, 9, 12, 13)]


prot_meta_c1 <- prot_meta %>% filter(clust %in% c(4))
prot_meta_c2 <- prot_meta %>% filter(!clust %in% c(4))
p1 <- t7[,rownames(prot_meta_c1)]
p2 <- t7[,rownames(prot_meta_c2)]

prot <- c()
dif <- c()
pval <- c()

for(i in 1:nrow(t7)){
  if(sum(is.na(p1[i,])==F) > 4){
    if(sum(is.na(p2[i,])==F) > 4){
      pval <- c(pval,t.test(p1[i,],p2[i,])$p.value)
      dif <- c(dif,median(p1[i,],na.rm=T)-median(p2[i,],na.rm=T))
      prot <- c(prot,rownames(p2)[i])
    }
  }
}

df <- as.data.frame(prot)
df$dif <- dif
df$pval <- pval



