#!/usr/bin/env Rscript
# Independent-donor pseudobulk DE, continuous-statistic gene sets, and corrected legacy ORA.
suppressPackageStartupMessages({library(limma);library(data.table)})
root <- normalizePath(getwd())
here <- file.path(root,"04-reviewer-reanalysis")
out <- file.path(here,"results")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(20260913)
say <- function(...) {cat(format(Sys.time()),paste0(...),"\n");flush.console()}
read_counts <- function(path) {
  d <- as.data.frame(fread(cmd=paste("gzip -dc",shQuote(path))),check.names=FALSE)
  rn <- d[[1]];d[[1]] <- NULL
  x <- as.matrix(d);rownames(x)<-rn;x
}
independent <- function(x) {
  keep<-integer()
  for(j in seq_len(ncol(x))) if(qr(x[,c(keep,j),drop=FALSE])$rank>length(keep)) keep<-c(keep,j)
  x[,keep,drop=FALSE]
}
make_design <- function(md,adjust_batch=TRUE) {
  cov<-cbind(Intercept=1,female=as.integer(md$genetic_sex==2))
  if(adjust_batch) cov<-cbind(cov,as.matrix(md[,paste0("batch_fraction_",2:5),drop=FALSE]))
  cov<-independent(cov)
  x<-cbind(cov,iPD=as.integer(md$group=="iPD"))
  stopifnot(qr(x)$rank==ncol(x),nrow(x)-ncol(x)>=4)
  x
}
median_ratio_lib <- function(y) {
  positive<-rowSums(y>0)==ncol(y)
  stopifnot(sum(positive)>100)
  geom<-exp(rowMeans(log(y[positive,,drop=FALSE])))
  sf<-apply(sweep(y[positive,,drop=FALSE],1,geom,"/"),2,median)
  sf<-sf/exp(mean(log(sf)))
  sf*exp(mean(log(colSums(y))))
}
gmt<-list()
for(f in list.files(file.path(root,"02-reanalysis","genesets"),pattern="\\.gmt$",full.names=TRUE)) {
  lines<-strsplit(readLines(f),"\t",fixed=TRUE)
  for(a in lines) {
    if(length(a)<3) next
    genes<-unique(toupper(trimws(sub(",.*$","",a[-c(1,2)]))))
    genes<-genes[nzchar(genes)]
    gmt[[paste0(sub("\\.gmt$","",basename(f)),"::",trimws(a[1]))]]<-genes
  }
}
panel<-read.csv("/data32TB/shared/ppmi-foundin/analyse/RNA-processed/shaohua-gene-list-2-10-25.csv")$gene
panel[panel=="PARK2"]<-"PRKN"

say("Reading cached donor counts and raw-verified culture counts")
cached<-read_counts(file.path(root,"02-reanalysis","pseudobulk_counts.tsv.gz"))
raw<-read_counts(file.path(here,"inputs","sample_pseudobulk_counts.tsv.gz"))
sm<-read.delim(file.path(here,"inputs","sample_pseudobulk_metadata.tsv"),check.names=FALSE)
stopifnot(identical(colnames(raw),as.character(sm$group)),identical(rownames(raw),rownames(cached)))
key<-paste(sm$PPMI_ID,sm$CellType,sep="|")
summed<-t(rowsum(t(raw),key,reorder=FALSE))
stopifnot(setequal(colnames(summed),colnames(cached)))
maxdiff<-max(abs(summed[,colnames(cached)]-cached))
stopifnot(maxdiff==0)
writeLines(c("All 34,960 genes and 879 donor/cell-type groups matched current raw-object extraction.",
             paste("Maximum absolute count difference",maxdiff)),file.path(out,"pseudobulk_source_validation.txt"))
rm(cached,summed);gc()

all_de<-list();all_camera<-list();all_roast<-list();summaries<-list();norm_rows<-list()
contexts<-c("iDA_pooled",sort(unique(sm$CellType)))
for(context in contexts) {
  selected<-if(context=="iDA_pooled") grepl("^iDA",sm$CellType) else sm$CellType==context
  m0<-sm[selected & sm$subtype %in% c("na_HC","na_iPD"),]
  y0<-raw[,m0$group,drop=FALSE]
  # Collapse clusters within each culture, then cultures within independent donors.
  culture_counts<-t(rowsum(t(y0),m0$SampleID,reorder=FALSE))
  culture_meta<-m0[match(colnames(culture_counts),m0$SampleID),]
  culture_meta$n_cells<-as.numeric(tapply(m0$n_cells,m0$SampleID,sum)[culture_meta$SampleID])
  ok<-culture_meta$n_cells>=20
  culture_counts<-culture_counts[,ok,drop=FALSE];culture_meta<-culture_meta[ok,]
  for(scenario in c("primary","exclude_reference_donor","original_30_iPD_grouping","sex_only")) {
    if(context!="iDA_pooled" && scenario!="primary") next
    eligible<-if(scenario=="original_30_iPD_grouping") rep(TRUE,nrow(culture_meta)) else !culture_meta$PPMI_ID %in% c(3954,4106)
    if(scenario=="exclude_reference_donor") eligible<-eligible & culture_meta$PPMI_ID!=3966
    cm<-culture_meta[eligible,]
    cy<-culture_counts[,eligible,drop=FALSE]
    donor_y<-t(rowsum(t(cy),as.character(cm$PPMI_ID),reorder=FALSE))
    md<-cm[match(colnames(donor_y),as.character(cm$PPMI_ID)),]
    md$group<-ifelse(md$subtype=="na_iPD","iPD","HC")
    # The covariates represent each batch's library contribution to a donor's
    # count pseudobulk. Reference-line cultures remain one independent donor.
    libs<-colSums(cy)
    for(b in 1:5) {
      contribution<-tapply(libs*as.integer(cm$BATCH==b),cm$PPMI_ID,sum)
      total<-tapply(libs,cm$PPMI_ID,sum)
      md[[paste0("batch_fraction_",b)]]<-as.numeric((contribution/total)[as.character(md$PPMI_ID)])
    }
    nh<-sum(md$group=="HC");np<-sum(md$group=="iPD")
    if(min(nh,np)<5) next
    x<-make_design(md,adjust_batch=scenario!="sex_only")
    lib<-median_ratio_lib(donor_y)
    cpm<-sweep(donor_y,2,lib,"/")*1e6
    keep<-rowSums(cpm>=1)>=min(nh,np)
    y<-donor_y[keep,,drop=FALSE]
    v<-voom(y,x,lib.size=lib,plot=FALSE)
    fit<-eBayes(lmFit(v,x))
    tt<-topTable(fit,coef="iPD",number=Inf,sort.by="none")
    tt$gene<-rownames(tt);tt$context<-context;tt$scenario<-scenario
    tt$n_HC<-nh;tt$n_iPD<-np;tt$model<-paste(colnames(x),collapse=" + ")
    all_de[[paste(context,scenario)]]<-tt
    nr<-data.frame(context=context,scenario=scenario,PPMI_ID=md$PPMI_ID,
                   raw_library_size=colSums(donor_y),effective_library_size=lib)
    norm_rows[[paste(context,scenario)]]<-nr
    idx<-intersect(panel,rownames(y))
    roast<-mroast(v,list(original_autophagy_panel=match(idx,rownames(y))),x,contrast=ncol(x),nrot=9999)
    roast$context<-context;roast$scenario<-scenario;roast$n_HC<-nh;roast$n_iPD<-np
    all_roast[[paste(context,scenario)]]<-roast
    summaries[[paste(context,scenario)]]<-data.frame(context=context,scenario=scenario,n_HC=nh,n_iPD=np,
        genes_tested=nrow(tt),DEG_FDR05=sum(tt$adj.P.Val<.05),panel_genes_tested=length(idx),
        panel_directional_P=roast$PValue,panel_mixed_P=roast$PValue.Mixed)
    say(context," ",scenario,": DEG=",sum(tt$adj.P.Val<.05),"; autophagy mixed P=",signif(roast$PValue.Mixed,3))
    if(scenario=="primary") {
      indices<-lapply(gmt,function(g) which(toupper(rownames(y)) %in% g))
      lengths<-lengths(indices)
      indices<-indices[lengths>=5 & lengths<=500]
      ca<-camera(v,indices,x,contrast=ncol(x),inter.gene.cor=NA,allow.neg.cor=FALSE,sort=FALSE)
      ca$term<-rownames(ca);ca$context<-context
      # Explicitly include every eligible set in the family.
      ca$q_within_context<-p.adjust(ca$PValue,method="BH",n=length(indices))
      all_camera[[context]]<-ca
      say(context," CAMERA complete: ",nrow(ca)," sets")
    }
  }
}
de<-do.call(rbind,all_de)
de$q_global_primary<-NA_real_
idx<-de$scenario=="primary"
de$q_global_primary[idx]<-p.adjust(de$P.Value[idx],"BH")
con<-gzfile(file.path(out,"pseudobulk_DE_all_results.csv.gz"),open="wt")
write.csv(de,con,row.names=FALSE);close(con)
write.csv(do.call(rbind,summaries),file.path(out,"pseudobulk_DE_summary.csv"),row.names=FALSE)
write.csv(do.call(rbind,norm_rows),file.path(out,"pseudobulk_normalization.csv"),row.names=FALSE)
ca<-do.call(rbind,all_camera);ca$q_global_contexts<-p.adjust(ca$PValue,"BH")
con<-gzfile(file.path(out,"camera_all_pathways.csv.gz"),open="wt")
write.csv(ca,con,row.names=FALSE);close(con)
ro<-do.call(rbind,all_roast)
ro$q_across_primary_contexts_directional<-NA_real_;ro$q_across_primary_contexts_mixed<-NA_real_
idx<-ro$scenario=="primary"
ro$q_across_primary_contexts_directional[idx]<-p.adjust(ro$PValue[idx],"BH")
ro$q_across_primary_contexts_mixed[idx]<-p.adjust(ro$PValue.Mixed[idx],"BH")
write.csv(ro,file.path(out,"autophagy_panel_roast.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"R_session_info.txt"))
say("Expression and pathway analyses complete")
