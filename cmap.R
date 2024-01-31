#!/usr/bin/Rscript
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--twas_sumstats_dir", action="store", default=NA, type='character',
		help="Path to per chromosome target PLINK files [required]"),
make_option("--cmap_data", action="store", default=NA, type='character',
		help="Path to keep file for target [optional]"),
make_option("--output_dir", action="store", default=NA, type='character',
		help="Name of output [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# For testing set opt manually
opt$twas_sumstats_dir<-'/users/k1928822/cmap_phase2'
opt$cmap_data<-'/users/k1928822/cmap_phase2/AD_ctrln_MCF7'
opt$output_dir<-'/users/k1928822/cmap_phase2'

library(mygene)
library(dplyr)
library(Hmisc)
setwd(opt$twas_sumstats_dir)
noperm = 100
filenames <- list.files(".", pattern="*.txt", full.names=FALSE)
cmap = read.table(opt$cmap_data, header=T,row.names=1)
thres.N.vector = c(50,100,250,500)
no.tissues = length(filenames)
entrez_id = rownames(cmap)
cmap = cbind(entrez_id, cmap)
colnames(cmap)[1] = "entrezgene"
for (q in 1:no.tissues){

gwasres = read.table(filenames[q], header=TRUE)
gwasres <- na.omit(gwasres)
gwasres <- subset(gwasres, TWAS.Z!=-Inf & TWAS.Z!=Inf)
idres = queryMany(gwasres$ID, scopes="symbol", fields="entrezgene", species="human",returnall=TRUE, return.as="DataFrame")
result.mygene = idres$response  
ind.query = which(colnames(result.mygene)=="query")  
colnames(result.mygene)[ind.query] <- "symbol"
colnames(gwasres)[1] <- "symbol"
result.mygene.df = data.frame(result.mygene)  #convert to conventional dataframe instead of S4vector
gwasres2 = merge(gwasres, result.mygene.df ,by = "symbol" )
finalres = merge(gwasres2, cmap, by="entrezgene")
head(finalres)
nodrugs = ncol(finalres)-6  ##need to check, here we assume the first 12 columns of "finalres" are NOT drug-expression data; ie drug transcriptome data starts from 13th column
  
no.thres.N = length(thres.N.vector)
  
  cor.spearman=NULL
  cor.pearson=NULL
  p.spearman = NULL
  p.pearson = NULL
  ks.ES = NULL
  ks.p = NULL 
  
  ks.signed = matrix(ncol=no.thres.N,nrow=nodrugs)
  extreme.cor.spearman = matrix(ncol=no.thres.N,nrow=nodrugs)
  extreme.cor.pearson = matrix(ncol=no.thres.N,nrow=nodrugs)
  connect.score = vector()
  finalres = data.frame(finalres)
  
  #******************************************************
  #         comparing expression in drugs vs diseases
  #******************************************************
  for (i in 1:nodrugs) {
    disease.zscore = finalres$TWAS.Z
    spearman.obj = cor.test(finalres[,i+6],disease.zscore, method="spearman",use="na.or.complete")
    cor.spearman[i] = spearman.obj$estimate
    pearson.obj  = cor.test(finalres[,i+6],disease.zscore, method = "pearson",use="na.or.complete")
    cor.pearson[i]= pearson.obj$estimate
    p.spearman[i] = spearman.obj$p.value
    p.pearson[i]= pearson.obj$p.value

    nogenes = nrow(finalres)
    pos.zscore = disease.zscore[disease.zscore>=0]
    no.pos.genes = length(pos.zscore)
    neg.zscore = disease.zscore[disease.zscore<0]
    no.neg.genes = length(neg.zscore)
    
    
    for (j in 1:no.thres.N) {
      thres.N = thres.N.vector[j]
      
      #****************************************
      #   KS method
      #****************************************
      rank.pos = no.pos.genes+1-rank(pos.zscore)  
      up.thres = pos.zscore[rank.pos==thres.N] 
      rank.neg = rank(neg.zscore)
      down.thres = neg.zscore[rank.neg==thres.N]
      
      ind.upreg = order(disease.zscore,decreasing=T)[1:thres.N]
      
      geneset2 = rank(-finalres[,i+6])[ind.upreg] 
      geneset2 = sort(geneset2)
      a.up = max(   (1:thres.N)/thres.N - geneset2/nogenes      )
      b.up = max(   geneset2/nogenes -  (1:thres.N-1)/thres.N )                                           
      ks.test.obj.up = ifelse(a.up>b.up, a.up , -b.up)
      
      ind.downreg = order(disease.zscore,decreasing=F)[1:thres.N]
      
      geneset2 = rank(-finalres[,i+6])[ind.downreg]
      geneset2 = sort(geneset2)
      
      a.down = max(   (1:thres.N)/thres.N - geneset2/nogenes      )
      b.down = max(geneset2/nogenes - (1:thres.N-1)/thres.N )                                           
      ks.test.obj.down = ifelse(a.down>b.down, a.down , -b.down)
      
      
      ks.signed[i,j]=0
      if (  sign(ks.test.obj.up)!=sign(ks.test.obj.down)  ) 
      {
        ks.signed[i,j] = ks.test.obj.up - ks.test.obj.down
      }
      
      #****************************************
      #   spearman and Pearson correlations (either use all observations or only compare the most extreme observations)
      #****************************************
      spearman.obj.extr = cor.test( finalres[,i+6][c(ind.upreg,ind.downreg)], disease.zscore[c(ind.upreg,ind.downreg)] 
                                    ,method="spearman",use="na.or.complete"   )
      extreme.cor.spearman[i,j] =  spearman.obj.extr$estimate 
      pearson.obj.extr  = cor.test(finalres[,i+6][c(ind.upreg,ind.downreg)],disease.zscore[c(ind.upreg,ind.downreg)], method = "pearson",use="na.or.complete")
      extreme.cor.pearson[i,j]= pearson.obj.extr$estimate
      
    }
  }
  
  
  drugmat = cbind( colnames(finalres)[7:ncol(finalres)], 
                   cor.pearson ,cor.spearman, 
                   ks.signed, 
                   extreme.cor.spearman,
                   extreme.cor.pearson)
  drugmat = data.frame(drugmat)
  
  #*************************************************
  # averaging ranks over the thresholds and then averaging over all methods
  #*************************************************
  rank.ks = matrix(nrow=nodrugs, ncol=no.thres.N)
  rank.spearman = matrix(nrow=nodrugs, ncol=no.thres.N)
  rank.pearson = matrix(nrow=nodrugs, ncol=no.thres.N)
  
  for (j in 1:no.thres.N){
    rank.ks[,j]=rank(ks.signed[,j])
    rank.spearman[,j]=rank(extreme.cor.spearman[,j])
    rank.pearson[,j]=rank(extreme.cor.pearson[,j])
  }
  
  mean.rank.ks = rowMeans(rank.ks)
  mean.rank.spearman = rowMeans(rank.spearman)
  mean.rank.pearson = rowMeans(rank.pearson)
  
  five.method.rank = cbind(rank(cor.pearson),
                           rank(cor.spearman), 
                           rank(mean.rank.ks), 
                           rank(mean.rank.spearman), 
                           rank(mean.rank.pearson) )
  
  res.5method.avg = data.frame( colnames(finalres)[7:ncol(finalres)], rowMeans(five.method.rank) ) 
  colnames(res.5method.avg)   <- c("Drug","AvgRank")    
  
   
  #*********************************************************************
  # permutation to determine the significance of the avg. rank obtained (repeats the above code)
  #*********************************************************************
  cor.spearman=NULL
  cor.pearson=NULL
  
  ks.signed = matrix(ncol=no.thres.N,nrow=nodrugs)
  extreme.cor.spearman = matrix(ncol=no.thres.N,nrow=nodrugs)
  extreme.cor.pearson = matrix(ncol=no.thres.N,nrow=nodrugs)
  avgrank.perm <- matrix(nrow=nodrugs, ncol=noperm)
  
  for (r in 1:noperm) {
    disease.zscore = sample(finalres$TWAS.Z)
    
    for (i in 1:nodrugs) {
    spearman.obj = cor.test(finalres[,i+6],disease.zscore, method="spearman",use="na.or.complete")
    cor.spearman[i] = spearman.obj$estimate
    pearson.obj  = cor.test(finalres[,i+6],disease.zscore, method = "pearson",use="na.or.complete")
    cor.pearson[i]= pearson.obj$estimate

    nogenes = nrow(finalres)
    
    pos.zscore = disease.zscore[disease.zscore>=0]
    no.pos.genes = length(pos.zscore)
    neg.zscore = disease.zscore[disease.zscore<0]
    no.neg.genes = length(neg.zscore)
    
    
    for (j in 1:no.thres.N) {
      thres.N = thres.N.vector[j]
      ind.upreg = order(disease.zscore,decreasing=T)[1:thres.N]
      
      #****************************************
      #   KS method
      #******************************************
      geneset2 = rank(-finalres[,i+6])[ind.upreg] 
      geneset2 = sort(geneset2) 

      a.up = max(   (1:thres.N)/thres.N - geneset2/nogenes      )
      b.up = max(   geneset2/nogenes -  (1:thres.N-1)/thres.N )                                           
      ks.test.obj.up = ifelse(a.up>b.up, a.up , -b.up)
      
      
      ind.downreg = order(disease.zscore,decreasing=F)[1:thres.N]
      
      geneset2 = rank(-finalres[,i+6])[ind.downreg]
      geneset2 = sort(geneset2)
      a.down = max(   (1:thres.N)/thres.N - geneset2/nogenes      )
      b.down = max(geneset2/nogenes - (1:thres.N-1)/thres.N )                                           
      ks.test.obj.down = ifelse(a.down>b.down, a.down , -b.down)
      
      
      ks.signed[i,j]=0
      if (  sign(ks.test.obj.up)!=sign(ks.test.obj.down)  ) 
      {
        ks.signed[i,j] = ks.test.obj.up - ks.test.obj.down
      
      }
      
      #****************************************
      #   spearman and Pearson correlations (either use all observations or only compare the most extreme observations)
      #****************************************
      spearman.obj.extr = cor.test( finalres[,i+6][c(ind.upreg,ind.downreg)], disease.zscore[c(ind.upreg,ind.downreg)] 
                                    ,method="spearman",use="na.or.complete"   )
      extreme.cor.spearman[i,j] =  spearman.obj.extr$estimate 
      pearson.obj.extr  = cor.test(finalres[,i+6][c(ind.upreg,ind.downreg)],disease.zscore[c(ind.upreg,ind.downreg)], method = "pearson",use="na.or.complete")
      extreme.cor.pearson[i,j]= pearson.obj.extr$estimate
      
    }
  
  } #end of looping over all drugs
  
  drugmat = cbind( colnames(finalres)[7:ncol(finalres)], 
                   cor.pearson ,cor.spearman, 
                   ks.signed, 
                   extreme.cor.spearman,
                   extreme.cor.pearson)
  

  drugmat = data.frame(drugmat)
  
  #*************************************************
  # averaging ranks over the thresholds and then averaging over all methods
  #*************************************************
  rank.ks = matrix(nrow=nodrugs, ncol=no.thres.N)
  rank.spearman = matrix(nrow=nodrugs, ncol=no.thres.N)
  rank.pearson = matrix(nrow=nodrugs, ncol=no.thres.N)
  
  for (j in 1:no.thres.N){
    rank.ks[,j]=rank(ks.signed[,j])
    rank.spearman[,j]=rank(extreme.cor.spearman[,j])
    rank.pearson[,j]=rank(extreme.cor.pearson[,j])
  }
  
  mean.rank.ks = rowMeans(rank.ks)
  mean.rank.spearman = rowMeans(rank.spearman)
  mean.rank.pearson = rowMeans(rank.pearson)
  
  five.method.rank = cbind(rank(cor.pearson),
                           rank(cor.spearman), 
                           rank(mean.rank.ks), 
                           rank(mean.rank.spearman), 
                           rank(mean.rank.pearson) )
  
  avgrank.perm[,r] <- rowMeans(five.method.rank) 
    
    
  }  # end of permutation loop
  
  
  
  perm.finalres = cbind( as.numeric(res.5method.avg$AvgRank), avgrank.perm)
  #*************************************************
  # calculate permutation p-values (pooling over all drugs to obtain the permutation null distribution)
  #*************************************************
  perm.p =NULL
  allPermRank = c(avgrank.perm)
  realnoperm = length(allPermRank)
  for (s in 1:nodrugs){
    perm.p[s] = sum( allPermRank <= perm.finalres[s,1])/ (realnoperm)
  }
  
  perm.finalres = cbind(colnames(finalres)[7:ncol(finalres)], 
                        as.numeric(res.5method.avg$AvgRank), 
                        perm.p, 
                        perm.finalres[,-1])
  colnames(perm.finalres)[c(1,2)]<- c("Drug","AvgRank")
  #********************************************************
  # ranking the table by the original (un-permuted) ranking 
  #*********************************************************
  perm.finalres = data.frame(perm.finalres)
  perm.finalres$AvgRank <- as.numeric(as.character(perm.finalres$AvgRank))
  perm.finalres.arr = arrange(perm.finalres, AvgRank)
  
##writing the results##
write.table(avgrank.perm, paste(opt$output,"/avgrank.perm", gsub ("*.txt", "", x=filenames[q]), sep="_"))
write.table(cor.spearman, paste(opt$output,"/cor.spearman", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(extreme.cor.spearman, paste(opt$output,"/extreme.cor.spearman", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(cor.pearson, paste (opt$output,"/cor.pearson", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(five.method.rank, paste (opt$output,"/five.method.rank", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(ks.signed, paste (opt$output,"/ks.signed", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(mean.rank.ks, paste (opt$output,"/mean.rank.ks", gsub("*.txt", "", x=filenames[q]), sep="_")) 
write.table(p.spearman, paste (opt$output,"/p.spearman", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(rank.ks, paste (opt$output,"/rank.ks", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(res.5method.avg, paste (opt$output,"/res.5method.avg", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(mean.rank.pearson, paste (opt$output,"/mean.rank.pearson", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(rank.spearman, paste (opt$output,"/rank.spearman", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(extreme.cor.pearson, paste (opt$output,"/extreme.cor.pearson", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(ks.p, paste (opt$output,"/ks.p", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(mean.rank.spearman, paste (opt$output,"/mean.rank.spearman", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(p.pearson, paste (opt$output,"/p.pearson", gsub("*.txt", "", x=filenames[q]), sep="_"))
write.table(rank.pearson, paste (opt$output,"/rank.pearson", gsub("*.txt", "", x=filenames[q]), sep="_"))                      
write.table(perm.finalres.arr, paste(opt$output,"/AvgRank_5methods_w_permutation", gsub("*.txt", "", x=filenames[q]), sep="_") )
}
