# libraries
library(stringr)
library(alluvial)
library(scales)

# input
ref_fn = "refOGs.csv"

ort_sets = list(
  list(id="tight", ort_fo="orthobench_trees/tight/"),
  list(id="raw", ort_fo="orthobench_trees/raw/"),
  list(id="own", ort_fo="results_orthology/")
)


# load data
ref = read.table(ref_fn, sep="\t", header = F, col.names = c("refOG", "gene"), stringsAsFactors=F)
ref$species = stringr::str_split(ref$gene, pattern="_", simplify = T)[,1]
fam_list = unique(ref$refOG)
sps_list = unique(ref$species)

# store diagnostics
dia = data.frame()

#### per-family evaluation ####


for (ort_set in ort_sets) {
  
  
  ort_fo = ort_set$ort_fo
  id=ort_set$id
  
  pdf(sprintf("eval_%s_classification_alluvial.pdf",id),height = 8, width = 7)
  for (fam in fam_list) {
    
    ort_fn = sprintf("%s/%s.possom.ortholog_groups.csv",ort_fo, fam)
    print(fam)
    
    if( file.exists(ort_fn)) {
      # read in possvm classification
      ort = read.table(ort_fn, sep="\t", header = T, stringsAsFactors = F)
      ort$species = stringr::str_split(ort$gene, pattern = "_", simplify = T)[,1]
      ort = ort[ort$species %in% sps_list,]
      
      # subset ref to refOG of interest
      rei = ref[ref$refOG == fam,c("gene","refOG")]
      rei$refOG_bool = T
      
      # add ref annot to possvm classification
      ort = merge(ort[,c("gene","orthogroup")],rei, by.x = "gene",by.y = "gene", all.y = T,all.x = T)
      ort[is.na(ort$refOG),"refOG"] = "other"
      ort[is.na(ort$refOG_bool),"refOG_bool"] = F
      
      # identify equivalent possvm orthogroup
      ort_xtab = as.data.frame(table(ort$refOG,ort$orthogroup), stringsAsFactors = F)
      colnames(ort_xtab) = c("refOG","Possvm","Freq")
      ort_xtab_pos = ort_xtab[ort_xtab$refOG == fam,]
      # this line will retrieve all OGs that have a shared gene with refOG (more inclusive, implies better recall)
      ort_xtab_max = which(ort_xtab_pos$Freq>1)
      # this line will the BEST OGs with shared genes with refOG (less inclusive, aims to identify a single hit, therefore lower recall)
      # ort_xtab_max = which.max(ort_xtab_pos$Freq)
      ort_equi = ort_xtab_pos[ort_xtab_max,"Possvm"]
      ort_equi_string = paste(ort_equi,collapse = ",")
      
      # assign bool to most likely Possvm OG
      ort$orthogroup_bool = (ort$orthogroup %in% ort_equi)
      ort[is.na(ort$orthogroup_bool),"orthogroup_bool"] = F
      
      # evaluate
      ev_TP = sum(ort$orthogroup_bool & ort$refOG_bool)
      ev_TN = sum(!ort$orthogroup_bool & !ort$refOG_bool)
      ev_FP = sum(ort$orthogroup_bool & !ort$refOG_bool)
      ev_FN = sum(!ort$orthogroup_bool & ort$refOG_bool)
      # evaluate precision = TP / (TP + FP)
      ev_precision = ev_TP / ( ev_TP + ev_FP )
      # recall = TP / (TP + FN)
      ev_recall = ev_TP / ( ev_TP + ev_FN )
      # F score = (2*Precision*Recall) / sum(Precision, Recall)
      ev_Fscore = (2*ev_precision*ev_recall) / sum(ev_recall, ev_precision)
      
      # alluvial plot
      if (nrow(ort_xtab)>1){
        alluvial(
          ort_xtab[,1:2], 
          freq = ort_xtab[,3], 
          col = c("lightblue"), gap.width = 0.1, blocks = T, 
          border=NA, cex=0.8)
      } else {
        plot(0,0, xlab = "", ylab = "", frame.plot = F, axes = F, col="blue", pch=19)
      }
      title(
        main=sprintf("%s\nn=%i",fam, length(rei$gene)),
        sub=sprintf("%s is %s\nPrecision = %.2f | Recall = %.2f | F-score = %.2f", fam, ort_equi_string,ev_precision, ev_recall, ev_Fscore),
        cex.main=0.8,
        cex.sub=0.7)
      
      # store diagnostics
      dia = rbind(
        dia, 
        data.frame(
          refOG = fam, 
          best_hit=ort_equi_string, 
          precision=ev_precision, 
          recall=ev_recall, 
          Fscore=ev_Fscore,
          TP=ev_TP,
          TN=ev_TN,
          FP=ev_FP,
          FN=ev_FN))
    }
    
  }
  dev.off()
  
  
  #### Summary plots ####
  
  # compare precsion and recall
  
  pdf(sprintf("eval_%s_summary.pdf",id),height = 6, width = 8)
  layout(matrix(1:6, nrow = 2))
  plot(dia$precision, dia$recall, xlim = c(0,1), ylim=c(0,1), xlab = "Precision", ylab="Recall",
       col=alpha("blue", 0.6), main="Precision & recall")
  text(dia$precision, dia$recall, labels = dia$refOG, col=alpha("lightblue",0.8))
  
  # dists
  hist(dia$Fscore, breaks = 10, xlim = c(0,1),main="F-score", col="blue", border = "white", xlab = "F-score")
  hist(dia$precision, breaks = 10, xlim = c(0,1),main="Precision", col="blue", border = "white", xlab = "Precision")
  hist(dia$recall, breaks = 10, xlim = c(0,1),main="Recall", col="blue", border = "white", xlab = "Recall")
  plot(sort(dia$precision), col="blue", ylab = "Precision", ylim = c(0,1), main="Precision")
  plot(sort(dia$recall), col="blue", ylab = "Recall", ylim = c(0,1), main="Recall")
  
  dev.off()
  
  # save table
  write.table(dia, file=sprintf("eval_%s_summary.csv",id), quote = F, sep="\t",row.names = F)
  
}