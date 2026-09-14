#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1]) else normalizePath(getwd())
out <- file.path(root, "04-reviewer-reanalysis", "inputs")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
source <- "/data32TB/shared/ppmi-foundin/ibm-aspera-gene-data/processed/SCRN/iPSCsDopaALL_integratedAfterBroadCellType.RDS"
say <- function(...) {cat(format(Sys.time()), paste0(...), "\n"); flush.console()}
write_gz <- function(x, path) {con <- gzfile(path, open="wt"); on.exit(close(con)); write.table(x, con, sep="\t", quote=FALSE, row.names=FALSE)}
say("Reading current source object ", source)
obj <- readRDS(source)
meta <- obj@meta.data
cnt <- obj@assays$RNA@counts
dat <- obj@assays$RNA@data
stopifnot(identical(colnames(cnt), rownames(meta)), identical(colnames(dat), rownames(meta)))
say("Loaded ", ncol(obj), " cells; ", nrow(cnt), " RNA genes")
cache <- read.delim(file.path(root,"02-reanalysis","expr_SNCA_autophagy.tsv.gz"), check.names=FALSE)
stopifnot(identical(as.character(cache$cell), colnames(dat)))
cached_genes <- setdiff(colnames(cache),"cell")
actual <- t(as.matrix(dat[cached_genes,,drop=FALSE]))
maxdiff <- max(abs(actual - as.matrix(cache[,cached_genes,drop=FALSE])))
stopifnot(maxdiff < 1e-10)
say("All cached RNA values matched; max absolute difference ", maxdiff)
rm(cache,actual); gc()
extra <- c("PRKN","PARK2","MAP2","NEFM","NEFL","DCX","SOX2","MKI67","RBFOX3")
genes <- unique(c(cached_genes,intersect(extra,rownames(dat))))
write_gz(data.frame(cell=colnames(dat),t(as.matrix(dat[genes,,drop=FALSE])),check.names=FALSE), file.path(out,"rna_panel_logexpr.tsv.gz"))
write_gz(data.frame(cell=colnames(cnt),t(as.matrix(cnt[genes,,drop=FALSE])),check.names=FALSE), file.path(out,"rna_panel_counts.tsv.gz"))
meta$cell <- rownames(meta)
meta$subtype <- meta$DESCRP_CAT
isna <- is.na(meta$DESCRP_CAT) | meta$DESCRP_CAT=="na"
for(g in c("HC","PD","SWEDD")) meta$subtype[isna & meta$ENROLL_CAT==g] <- paste0("na_",if(g=="PD")"iPD"else g)
write_gz(meta, file.path(out,"cell_metadata_full.tsv.gz"))
if("SCT" %in% names(obj@assays)) {
  sct <- obj@assays$SCT@data
  sg <- intersect(genes,rownames(sct))
  stopifnot(identical(colnames(sct),colnames(dat)))
  write_gz(data.frame(cell=colnames(sct),t(as.matrix(sct[sg,,drop=FALSE])),check.names=FALSE),file.path(out,"sct_panel_logexpr.tsv.gz"))
  say("SCT panel genes ", length(sg))
}
key <- paste(meta$SampleID,meta$CellType,sep="|")
fac <- factor(key)
ind <- sparse.model.matrix(~0+fac)
colnames(ind) <- levels(fac)
pb <- cnt %*% ind
nc <- as.integer(Matrix::colSums(ind))
first <- match(levels(fac),key)
sm <- meta[first,c("SampleID","PPMI_ID","CellType","subtype","genetic_sex","BATCH"),drop=FALSE]
sm$group <- levels(fac); sm$n_cells <- nc
stopifnot(identical(paste(sm$SampleID,sm$CellType,sep="|"),colnames(pb)))
write_gz(data.frame(gene=rownames(pb),as.matrix(pb),check.names=FALSE),file.path(out,"sample_pseudobulk_counts.tsv.gz"))
write.table(sm,file.path(out,"sample_pseudobulk_metadata.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(gene=rownames(cnt)),file.path(out,"rna_gene_names.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
sink(file.path(out,"source_validation.txt"))
cat("Source",source,"\nSource_bytes",file.info(source)$size,"\nSource_mtime",as.character(file.info(source)$mtime),"\n")
cat("Cells",ncol(obj),"\nGenes",nrow(cnt),"\nCached_genes_compared",length(cached_genes),"\nCached_values_max_abs_diff",maxdiff,"\n")
cat("Panel_genes",paste(genes,collapse=","),"\nSample_celltype_groups",ncol(pb),"\nTotal_counts",sum(pb),"\n")
print(sessionInfo()); sink()
say("Extraction complete; total counts ",sum(pb))
