library(Seurat)
library(dplyr)
library(UpSetR)
library(rlist)

set.seed(12345)

# Set date and create out folder
getwd()
Sys.Date()
main_dir <- "~/Dropbox/Lab data/PCLS/PCLS_061921/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

getwd()

# Read in Seurat object
pcls <- readRDS("~/pcls_merged_annotated.rds")

# Running DE analysis using negative binominal test
pcls_list <- SplitObject(pcls, split.by = 'celltype')



pcls_vs_fresh <- lapply(pcls_list, function(xx){
  print(unique(xx@meta.data$celltype))
  if(length(unique(xx@meta.data$diagnosis)) > 1) {
    FindMarkers(xx, group.by = "diagnosis", ident.1 = "PCLS", ident.2 = "Fresh Control", test.use = "negbinom")
  } 
  else{
    return(NULL)
  } 
})

for(i in 1:length(pcls_list)){
  write.csv(pcls_vs_fresh[[i]], paste(gsub("/", "", unique(pcls_list[[i]]@meta.data$celltype)), "_pcls_vs_fresh", ".csv"), sep =",", quote = F)
}


# Myeloid cells
myeloid <- list.subset (pcls_vs_fresh, c("Alveolar macrophage", 'cDC', 'Monocyte', 'Monocyte-derived macrophage', 'moDC', 'Proliferating macrophage'))

onion <- lapply(myeloid, function(xx){ row.names(xx[xx$p_val_adj <= .05,])})
onion <- unique(unlist(onion))
onion2 <- lapply(myeloid, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .05,])})
upset_d_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_d_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]], onion2[[5]], onion2[[6]])
upset_d_vs_c <- as.data.frame(upset_d_vs_c)

row.names(upset_d_vs_c) <- onion
colnames(upset_d_vs_c) <- names(myeloid)

upset_d_vs_c[upset_d_vs_c == T] <- 1
upset_d_vs_c[upset_d_vs_c == F] <- 0

pdf("20210726_myeloid_upset.pdf", height = 5, width = 15)
upset(upset_d_vs_c, nsets = 25, text.scale = 1, show.numbers = 'yes', nintersects = NA, point.size = 1, line.size = 0.5, mb.ratio = c(0.7, 0.3))
dev.off()

write.csv(upset_d_vs_c, file = '~/Dropbox/Lab data/PCLS/PCLS_061921/20210719/myeloid_pcls_vs_fresh_upset.csv')


# Lymphoid cells
lymphoid <- list.subset (pcls_vs_fresh, c("NK", 'NKT', 'CD8', 'CD4', 'Proliferating T-cells'))

onion <- lapply(lymphoid, function(xx){ row.names(xx[xx$p_val_adj <= .05,])})
onion <- unique(unlist(onion))
onion2 <- lapply(lymphoid, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .05,])})
upset_d_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_d_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]], onion2[[5]])
upset_d_vs_c <- as.data.frame(upset_d_vs_c)

row.names(upset_d_vs_c) <- onion
colnames(upset_d_vs_c) <- names(lymphoid)

upset_d_vs_c[upset_d_vs_c == T] <- 1
upset_d_vs_c[upset_d_vs_c == F] <- 0

pdf("20210728_lymphoid_upset.pdf", height = 5, width = 15)
upset(upset_d_vs_c, nsets = 25, text.scale = 1, show.numbers = 'yes', nintersects = NA, point.size = 1, line.size = 0.5, mb.ratio = c(0.7, 0.3))
dev.off()


write.csv(upset_d_vs_c, file = '~/Dropbox/Lab data/PCLS/PCLS_061921/20210719/lymphoid_pcls_vs_fresh_upset.csv')


