# libraries
library(stringr)
library(alluvial)
library(scales)

# input
ort_fo = "results_trees/"
out_fo = "results_evaluation/"

# list of datasets
set_list = list(
  list(id = "ANTPm", ref = "results_searches/ANTPm.seed.diamond.csv", eval_sps = c("Hsap","Mmus","Galgal","Xtro","Drer","Bflo","Dmel","Apimel","Tcas")),
  list(id = "PRDm", ref = "results_searches/PRDm.seed.diamond.csv", eval_sps = c("Hsap","Mmus","Galgal","Xtro","Drer","Bflo","Dmel","Apimel","Tcas")),
  list(id = "TALEm", ref = "results_searches/TALEm.seed.diamond.csv", eval_sps = c("Hsap","Mmus","Galgal","Xtro","Drer","Bflo","Dmel","Apimel","Tcas")),
  list(id = "ANTP", ref = "results_searches/ANTP.seed.diamond.csv", eval_sps = c("Hsap","Mmus","Ggal","Xtro","Drer","Bflo","Dmel","Tcas","Cele")),
  list(id = "PRD", ref = "results_searches/PRD.seed.diamond.csv", eval_sps = c("Hsap","Mmus","Ggal","Xtro","Drer","Bflo","Dmel","Tcas","Cele")),
  list(id = "TALE", ref = "results_searches/TALE.seed.diamond.csv", eval_sps = c("Hsap","Mmus","Ggal","Xtro","Drer","Bflo","Dmel","Tcas","Cele"))
)

for (set in set_list) {
  
  ref_fn = set$ref
  id = set$id
  
  # load data
  ref = read.table(ref_fn, sep="\t", header = F, col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"), stringsAsFactors=F)
  ref = ref[ref$pident>=99,]
  ref$refOG = stringr::str_split(ref$qseqid, "_", simplify = T)[,3]
  ref$species = stringr::str_split(ref$sseqid, pattern="_", simplify = T)[,1]
  ref = ref[,c("sseqid","refOG","species")]
  ref = ref[!duplicated(ref$sseqid),]
  write.table(ref, sprintf("reference_%s.csv", id), quote = F, sep="\t", row.names = F)
  
  fam_list = unique(ref$refOG)
  sps_list = set$eval_sps
  
  # store diagnostics
  
  #### per-family evaluation ####
  
  dia = data.frame()
  
  pdf(sprintf("%s/eval_%s_classification_alluvial.pdf",out_fo, id),height = 8, width = 7)
  for (fam in fam_list) {
    
    ort_fn = sprintf("%s/%s.possom.ortholog_groups.csv",ort_fo, id)
    print(paste(id,fam))
    
    # read in possvm classification
    ort = read.table(ort_fn, sep="\t", header = T, stringsAsFactors = F)
    ort$species = stringr::str_split(ort$gene, pattern = "_", simplify = T)[,1]
    
    # subset ref to refOG of interest
    rei = ref[ref$refOG == fam,c("sseqid","refOG")]
    rei$refOG_bool = T
    rei = rei[order(rei$refOG,rei$sseqid),]
    ort = ort[ort$species %in% sps_list,]
    
    
    # add ref annot to possvm classification
    ort = merge(ort[,c("gene","orthogroup")],rei, by.x = "gene",by.y = "sseqid", all.y = T,all.x = T)
    ort[is.na(ort$refOG),"refOG"] = "other"
    ort[is.na(ort$refOG_bool),"refOG_bool"] = F
    
    # identify equivalent possvm orthogroup
    ort_xtab = as.data.frame(table(ort$refOG,ort$orthogroup), stringsAsFactors = F)
    colnames(ort_xtab) = c("refOG","Possvm","Freq")
    ort_xtab_pos = ort_xtab[ort_xtab$refOG == fam,]
    # this line will retrieve all OGs that have a shared gene with refOG (more inclusive, implies better recall)
    ort_xtab_hits = which(ort_xtab_pos$Freq>0)
    ort_hits = ort_xtab_pos[ort_xtab_hits,"Possvm"]
    ort_hits_string = paste(ort_hits,collapse = ",")
    ort_hits_count = sum(ort_xtab_pos$Freq[ort_xtab_hits])
    
    # this line will the BEST OGs with shared genes with refOG (less inclusive, aims to identify a single hit, therefore lower recall)
    ort_xtab_max = which.max(ort_xtab_pos$Freq)
    ort_hits_best = ort_xtab_pos[ort_xtab_max,"Possvm"]
    ort_hits_best_count = sum(ort_xtab_pos$Freq[ort_xtab_max])
    
    
    # assign bool Possvm OGs
    ort$orthogroup_bool = (ort$orthogroup %in% ort_hits)
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
      main=sprintf("%s\nn=%i",fam, length(rei$sseqid)),
      sub=sprintf("%s is %s\nPrecision = %.2f | Recall = %.2f | F-score = %.2f", fam, ort_hits_string,ev_precision, ev_recall, ev_Fscore),
      cex.main=0.8,
      cex.sub=0.7)
    
    # store diagnostics
    dia = rbind(
      dia, 
      data.frame(
        refOG = fam, 
        hits=ort_hits_string, 
        num_hits=ort_hits_count,
        best_hit=ort_hits_best,
        num_best_hits=ort_hits_best_count,
        precision=ev_precision, 
        recall=ev_recall, 
        Fscore=ev_Fscore,
        TP=ev_TP,
        TN=ev_TN,
        FP=ev_FP,
        FN=ev_FN))
    
  }
  dev.off()
  
  #### Summary plots ####
  
  # compare precsion and recall
  
  pdf(sprintf("%s/eval_%s_summary.pdf",out_fo,id),height = 5, width = 7)
  layout(matrix(1:6, nrow = 2))
  plot(dia$precision, dia$recall, xlim = c(0,1), ylim=c(0,1), xlab = "Precision", ylab="Recall",
       col=alpha("blue", 0.6), main="Precision & recall", cex.axis=0.9, cex.lab=0.9)
  text(dia$precision, dia$recall, labels = dia$refOG, col=alpha("lightblue",0.8),cex=0.8)
  
  # dists of fscore, precision and recall
  hist(dia$Fscore, breaks = 10, xlim = c(0,1),main="F-score", col="blue", ylim=c(0,50),border = "white", xlab = "F-score", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i | mean = %.3f", sum(dia$Fscore >= 0.95), mean(dia$Fscore)))
  abline(v=0.95, lty=2, col="grey")
  
  hist(dia$precision, breaks = 10, xlim = c(0,1),main="Precision", col="blue",ylim=c(0,50), border = "white", xlab = "Precision", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i | mean = %.3f", sum(dia$precision >= 0.95), mean(dia$precision)))
  abline(v=0.95, lty=2, col="grey")
  
  hist(dia$recall, breaks = 10, xlim = c(0,1),main="Recall", col="blue",ylim=c(0,50), border = "white", xlab = "Recall", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i | mean = %.3f", sum(dia$recall >= 0.95), mean(dia$recall)))
  abline(v=0.95, lty=2, col="grey")
  
  plot(sort(dia$precision), col="blue", ylab = "Precision", ylim = c(0,1), main="Precision", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i", sum(dia$precision >= 0.95)))
  abline(h=0.95, lty=2, col="grey")
  
  plot(sort(dia$recall), col="blue", ylab = "Recall", ylim = c(0,1), main="Recall", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i", sum(dia$recall >= 0.95)))
  abline(h=0.95, lty=2, col="grey")
  
  # number of ogs per reference group
  plot(dia$num_hits, dia$num_best_hits, col="blue", xlab = "# genes in ref OG", ylab = "# genes in best OG", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(best >=0.9 all) = %i (%.3f | r=%.3f)\nn(best = all) = %i (p=%.3f | r=%.3f)\nn(best < all) = %i (p=%.3f | r=%.3f)", 
                    sum(dia$num_best_hits / dia$num_hits >= 0.9), 
                    mean(dia[dia$num_best_hits / dia$num_hits >= 0.9,"precision"]), 
                    mean(dia[dia$num_best_hits / dia$num_hits >= 0.9,"recall"]),
                    sum(dia$num_best_hits == dia$num_hits), 
                    mean(dia[dia$num_best_hits == dia$num_hits,"precision"]), 
                    mean(dia[dia$num_best_hits == dia$num_hits,"recall"]),
                    sum(dia$num_best_hits < dia$num_hits), 
                    mean(dia[dia$num_best_hits < dia$num_hits,"precision"]), 
                    mean(dia[dia$num_best_hits < dia$num_hits,"recall"])
  ),
  cex.sub=0.6)
  abline(a=0, b=1, lty=2, col="grey")
  
  dev.off()
  
  # save table
  write.table(dia, file=sprintf("%s/eval_%s_summary.csv",out_fo,id), quote = F, sep="\t",row.names = F)
}