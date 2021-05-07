#' Generate a kmer table from representative sequences
#'
#' Convert DNA sequences to tetranucleotide frequency table
#' @param repseq The fasta file path
#' @param file logical, default is TRUE, loading DNA sequences from file. if 'file=FALSE', the 'repseq' is loaded from object.
#' @return A 256 tetranucleotide combinations frequency table
#' @export
tetra.freq <- function(repseq,file=TRUE){
  rev_comp <- function(x){
    x.r <- paste0(rev(strsplit(x,"")[[1]]),collapse = "")
    x.r <- gsub("A","1",x.r);x.r <- gsub("T","2",x.r);x.r <- gsub("C","3",x.r);x.r <- gsub("G","4",x.r)
    x.r <- gsub("1","T",x.r);x.r <- gsub("2","A",x.r);x.r <- gsub("3","G",x.r);x.r <- gsub("4","C",x.r)
    return(x.r)
  }
  readseq <- function(seq){
    x <- seq
    x <- paste0(x,collapse = "")
    x.r <- rev_comp(x)
    return(list(x,x.r))
  }
  if(isTRUE(file)){
    scaffolds = Biostrings::readDNAStringSet(filepath = repseq, use.names = T)
  } else scaffolds = Biostrings::DNAStringSet(repseq)
  asv.id <- scaffolds@ranges@NAMES
  species <- lapply(scaffolds, function(x) readseq(as.character(x)))

  DNA <- c("A","T","C","G")
  tetra.mer <- expand.grid(DNA,DNA,DNA,DNA)
  tetra.mer <- do.call(paste0,tetra.mer)

  tetra.table <- matrix(nrow=256,ncol=length(species))
  rownames(tetra.table) <- tetra.mer

  #full length ASV
  for(i in 1:256){
    for(j in 1:length(species)){
      single.forward <- ifelse(length(grep(tetra.mer[i],species[[j]][[1]]))>0,grep(tetra.mer[i],species[[j]][[1]]),0)
      single.reverse <- ifelse(length(grep(tetra.mer[i],species[[j]][[2]]))>0,grep(tetra.mer[i],species[[j]][[2]]),0)
      tetra.table[i,j] <- (single.forward+single.reverse)
    }
  }
  colnames(tetra.table) <- asv.id
  tetra.table <- prop.table(tetra.table,2)
  return(tetra.table)
}


#' KTU clustering
#'
#' K-mer taxonomic unit clustering
#' @param repseq The fasta file path
#' @param feature.table (optional) 'data.frame' formated ASV/OTU table
#' @param write.fasta (optional) write out a representative KTU sequence fasta file
#' @param step Split searching range for optimal K. Two step searching process: large scale searching in the first round and smaller scale searching in the second round. Default is 'step=c(5,10)'.
#' @param search.min The minimum K for searching, default is 'NULL' = tip numbers at 0.03 height of cosine hierarchical clustering tree
#' @param search.max The maximum K for searching, default is 'NULL' = tip numbers at 0.015 height of cosine hierarchical clustering tree
#' @param cores Numbers of CPUs
#' @return KTU.table aggregated KTU table
#' @return ReqSeq Representative KTU sequences
#' @return kmer.table Tetranucleotide frequency table
#' @return clusters K-clusters of input features
#' @export
klustering <- function(repseq,feature.table=NULL,write.fasta=TRUE,step=c(5,10),search.min=NULL,search.max=NULL,cores=1){
  require(foreach,quietly = T)
  rev_comp <- function(x){
    x.r <- paste0(rev(strsplit(x,"")[[1]]),collapse = "")
    x.r <- gsub("A","1",x.r);x.r <- gsub("T","2",x.r);x.r <- gsub("C","3",x.r);x.r <- gsub("G","4",x.r)
    x.r <- gsub("1","T",x.r);x.r <- gsub("2","A",x.r);x.r <- gsub("3","G",x.r);x.r <- gsub("4","C",x.r)
    return(x.r)
  }
  readseq <- function(seq){
    x <- seq
    x <- paste0(x,collapse = "")
    x.r <- rev_comp(x)
    return(list(x,x.r))
  }
  aggregate2df <- function(data,groups,FUN){
    agg.data <- aggregate(data,list(groups),FUN)
    rownames(agg.data) <- agg.data[,1]
    agg.data <- as.data.frame(t(agg.data[,-1]))
    return(agg.data)
  }

  scaffolds = Biostrings::readDNAStringSet(filepath = repseq, use.names = T)
  asv.id <- scaffolds@ranges@NAMES
  species <- lapply(scaffolds, function(x) readseq(as.character(x)))

  DNA <- c("A","T","C","G")
  tetra.mer <- expand.grid(DNA,DNA,DNA,DNA)
  tetra.mer <- do.call(paste0,tetra.mer)

  tetra.table <- matrix(nrow=256,ncol=length(species))
  rownames(tetra.table) <- tetra.mer

  #full length ASV
  for(i in 1:256){
    for(j in 1:length(species)){
      single.forward <- ifelse(length(grep(tetra.mer[i],species[[j]][[1]]))>0,grep(tetra.mer[i],species[[j]][[1]]),0)
      single.reverse <- ifelse(length(grep(tetra.mer[i],species[[j]][[2]]))>0,grep(tetra.mer[i],species[[j]][[2]]),0)
      tetra.table[i,j] <- (single.forward+single.reverse)
    }
  }
  colnames(tetra.table) <- asv.id
  tetra.table <- prop.table(tetra.table,2)
  cos <- as.dist(1-coop::cosine(tetra.table))


  cl <- parallel::makeCluster(cores) #not to overload your computer
  doParallel::registerDoParallel(cl)

  tree <- hclust(cos) #### 20210428 update ####
  if(is.null(search.max)) search.max <- max(cutree(tree,h=0.015)) else search.max <- search.max #### 20210428 update ####
  if(is.null(search.min)) search.min <- max(cutree(tree,h=0.03)) else search.min <- search.min #### 20210428 update ####

  steps <- ceiling(seq(search.min,search.max,length.out = step[1]))
  repeat{
    steps.update <- steps
    asw <- foreach::foreach(k=steps.update, .combine=c) %dopar% {
      temp.asw = cluster::pam(cos,k,do.swap = F,pamonce = 5)$silinfo$avg.width #calling a function
      temp.asw #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    }
    k.best <- steps.update[order(asw,decreasing = T)][1:2]
    print(k.best)
    steps <- ceiling(seq(k.best[2],k.best[1],length.out = step[1]))
    if(all(steps==steps.update)) break
  }
  steps <- unique(ceiling(seq(k.best[2],k.best[1],length.out = step[2]))) #### 20210428 correct typo ####
  rm(k.best)
  repeat{
    steps.update <- steps
    asw <- foreach::foreach(k=steps.update, .combine=c) %dopar% {
      temp.asw = cluster::pam(cos,k,do.swap = F,pamonce = 5)$silinfo$avg.width #calling a function
      temp.asw #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    }
    if(length(steps.update)>1){
      k.best <- steps.update[order(asw,decreasing = T)][1:2]
    } else k.best <- c(steps.update[1],steps.update[1])
    print(k.best)
    steps <- unique(ceiling(seq(k.best[2],k.best[1],length.out = step[2])))
    if(all(steps==steps.update)) break
  }
  k.best <- steps.update[which.max(asw)]
  print(paste("Best KTU #:",k.best))
  parallel::stopCluster(cl) #stop cluster


  kms <- cluster::pam(cos,k.best,do.swap = F,pamonce = 5)
  #centoid rep seq
  centroid <- kms$id.med
  crepseq <- sapply(species[centroid], '[[',1)
  #names(crepseq) <- paste0("KTU",sprintf("%0004i",1:length(centroid)))
  for(i in 1:length(centroid)) names(crepseq)[i] <- digest::digest(crepseq[i],algo = "md5",serialize = F)

  #kmer table
  kmer.table <- data.frame(tetra.table[,centroid])
  colnames(kmer.table) <- names(crepseq)

  if(isTRUE(write.fasta)){
    repseq <- matrix(rbind(paste0(">",names(crepseq)),crepseq),ncol=1)
    write(repseq,file = "ktu-sequence.fasta")
  }

  if(!is.null(feature.table)){
    feature.table <- feature.table[match(asv.id,feature.table[,1],nomatch = 0),]
    otu <- feature.table[,-1]
    ktu <- as.data.frame(t(aggregate2df(otu,kms$clustering,sum)))
    rownames(ktu) <- names(crepseq)
    colnames(ktu) <- colnames(otu)
    return(list(KTU.table=ktu, ReqSeq=crepseq, kmer.table=kmer.table, clusters=kms$clustering))
  } else return(list(ReqSeq=crepseq, kmer.table=kmer.table, clusters=kms$clustering))
}


#' Prepare KTU taxonomy assignment database by trimming primers
#'
#' Trim primer sequences of full length 16S/18S/target gene database by PCR primers
#' @param fasta input database fasta file
#' @param forseq forward primer sequence (5' -> 3')
#' @param revseq reverse primer sequence (5' -> 3')
#' @param output.name (optional) the name for output file
#' @param progress show the processing progress
#' @export
trim.primer <- function(fasta,forseq,revseq,output.name=NULL,progress=TRUE){
  #read fasta file
  fastasets <- Biostrings::readDNAStringSet(fasta)
  seqNames <- fastasets@ranges@NAMES
  #format primer seq & reverse complement
  Lpattern=Biostrings::DNAString(forseq)
  Rpattern=Biostrings::reverseComplement(Biostrings::DNAString(revseq))
  #search primer position
  Forward <- Biostrings::vmatchPattern(pattern = Lpattern,subject = fastasets,fixed = F)
  Reverse <- Biostrings::vmatchPattern(pattern = Rpattern,subject = fastasets,fixed = F)

  trimmed <- c()
  for(i in 1:length(fastasets)){
    if(length(Forward[[i]]) > 0 & length(Reverse[[i]]) >0 ){
      Fstart <- ifelse(length(Forward[[i]])>1,Forward[[i]][1]@start,Forward[[i]]@start)
      Fwidth <- ifelse(length(Forward[[i]])>1,Forward[[i]][1]@width,Forward[[i]]@width)
      Rend   <- ifelse(length(Reverse[[i]])>1,Reverse[[i]][length(Reverse[[i]])]@start,Reverse[[i]]@start)
      if(Rend-(Fstart+Fwidth)>0){
      trimmed[i] <- as.character(fastasets[[i]][(Fstart+Fwidth):(Rend-1)])
      } else trimmed[i] <- NA
    } else trimmed[i] <- NA

    if(progress==TRUE){
      print(paste0(i,"/",length(fastasets)," sequences from DB are processed"))
      flush.console()
      Sys.sleep(0.01)
    }
  }

  names(trimmed) <- seqNames
  trimmed <- trimmed[!is.na(trimmed)]
  output <- matrix(rbind(paste0(">",names(trimmed)),trimmed),ncol=1)
  if(is.null(output.name)){
    write(output,file = "trimmed-sequence.fasta")
  } else if(!is.null(output.name)) write(output,file = paste0(output.name,"_trimmed-sequence.fasta"))
}

#' Build KTU taxonomy assignment database by tetranucleotide frequency table
#'
#' make KTU db
#' @param input.fasta input database fasta file (trimmed)
#' @param input.taxa input database taxonomy file
#' @param output.file (optional) the directory for output RDS format file
#' @export
makektudb <- function(input.fasta,input.taxa,output.file=NULL){
  rev_comp <- function(x){
    x.r <- paste0(rev(strsplit(x,"")[[1]]),collapse = "")
    x.r <- gsub("A","1",x.r);x.r <- gsub("T","2",x.r);x.r <- gsub("C","3",x.r);x.r <- gsub("G","4",x.r)
    x.r <- gsub("1","T",x.r);x.r <- gsub("2","A",x.r);x.r <- gsub("3","G",x.r);x.r <- gsub("4","C",x.r)
    return(x.r)
  }
  readseq <- function(seq){
    x <- seq
    x <- paste0(x,collapse = "")
    x.r <- rev_comp(x)
    return(list(x,x.r))
  }
  scaffolds = Biostrings::readDNAStringSet(filepath = input.fasta, use.names = T)

  DNA <- c("A","T","C","G")
  tetra.mer <- expand.grid(DNA,DNA,DNA,DNA)
  tetra.mer <- do.call(paste0,tetra.mer)

  species <- lapply(scaffolds, readseq)

  tetra.table <- matrix(nrow=256,ncol=length(species))
  rownames(tetra.table) <- tetra.mer

  for(i in 1:256){
    for(j in 1:length(species)){
      single.forward <- ifelse(length(grep(tetra.mer[i],species[[j]][[1]]))>0,grep(tetra.mer[i],species[[j]][[1]]),0)
      single.reverse <- ifelse(length(grep(tetra.mer[i],species[[j]][[2]]))>0,grep(tetra.mer[i],species[[j]][[2]]),0)
      tetra.table[i,j] <- (single.forward+single.reverse)
    }
  }

  colnames(tetra.table) <- scaffolds@ranges@NAMES
  tetra.table <- prop.table(tetra.table,2)

  taxa <- read.delim(input.taxa,header = F)
  taxa <- taxa[match(scaffolds@ranges@NAMES,taxa[,1]),]

  if(!is.null(output.file)){
    saveRDS(tetra.table,file = paste0(output.file,"_DB.RDS"))
    saveRDS(taxa,file = paste0(output.file,"_TX.RDS"))
  }
  return(list(tetra.table,taxa))
}

#' KTU annotation
#'
#' KTU-based alignment-free taxonomy assignment
#' @param dbRDS Kmer formated database RDS file
#' @param taxaRDS taxonomy ID RDS file
#' @param kmer.table tetranucleotide frequency table for annotation
#' @param cos.cutff cosine similarity cutoff, default=0.95
#' @param consensus taxonomy assignment by consensus [0-1]
#' @param candidate how many candidates (≥ cos.cutff) for taxonomy assignment
#' @param cores Numbers of CPUs
#' @export
kaxonomy <- function(dbRDS,taxaRDS,kmer.table,cos.cutff=0.95,consensus=0.8,candidate=10,cores=1){
  require(foreach,quietly = T)
  cl <- parallel::makeCluster(cores) #not to overload your computer
  doParallel::registerDoParallel(cl)

  db.tetra <- readRDS(dbRDS)
  db.taxa <- readRDS(taxaRDS)
  {
    taxa <- as.character(db.taxa[,2])
    sep <- ifelse(grepl("; ",taxa[1],fixed = T),"; ",";")
    taxa <- strsplit(taxa,sep)
    taxa <- data.frame(do.call(rbind,taxa),row.names = db.taxa[,1])
  }
  cos.table <- data.frame(matrix(data = NA,nrow = nrow(db.taxa),ncol=ncol(kmer.table)))
  cos.table <- foreach::foreach(cln=1:ncol(kmer.table), .combine=cbind) %dopar% {
    temp.cos.table = apply(db.tetra,2,function(x) coop::cosine(kmer.table[,cln],x)) #calling a function
    temp.cos.table #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
  }
  parallel::stopCluster(cl) #stop cluster

  machine.core <- function(x){
    xxx <- matrix(ncol = 8)
    xx <- cbind(taxa[which(x>=cos.cutff),],cos.score=x[which(x>=cos.cutff)])
    if(nrow(xx)>0){
      #xx <- xx[xx$cos.score>=quantile(xx$cos.score,top.pct,na.rm=TRUE),]
      ranks <- rank(1/(xx$cos.score))
      if(min(ranks)>candidate){
        xx <- xx[which(ranks==min(ranks)),]
      } else xx <- xx[which(ranks<=candidate),]
      xx[,-8] <- apply(xx[,-8], 2, as.character)
      xy <- apply(xx,2,function(x) names(table(x))[which(prop.table(table(x))>consensus)])
      xxx[,1:7] <- c(xy[[1]][1],xy[[2]][1],xy[[3]][1],xy[[4]][1],xy[[5]][1],xy[[6]][1],xy[[7]][1])
      #xxx[,8] <- mean(unlist(apply(xx,2,function(x) prop.table(table(x))[which(prop.table(table(x))>consensus)])))
      xxx[,8] <- mean(mean(xx[,8]) * unlist(apply(xx,2,function(x) prop.table(table(x))[which(prop.table(table(x))>consensus)])))
      xxx[,7] <- ifelse(!is.na(xxx[,7]) & is.na(xxx[,5]),NA,xxx[,7])
      for(i in 7:2) xxx[,i] <- ifelse(!is.na(xxx[,i]) & is.na(xxx[,i-1]),NA,xxx[,i])
    } else if(nrow(xx)==0) xxx[,c(1,8)] <- c("Unassigned",1)
    return(xxx)
  }

  taxa.consensus <- lapply(X = as.list(as.data.frame(cos.table)),
                                       FUN = machine.core)
  taxa.consensus <- data.frame(do.call(rbind,taxa.consensus))
  colnames(taxa.consensus) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Cos.score")
  taxa.consensus <- apply(taxa.consensus, 2, as.character)
  taxa.consensus[is.na(taxa.consensus)] <- "Unassigned"
  taxa.consensus <- as.data.frame(taxa.consensus)
  taxa.consensus[,8] <- as.numeric(as.character(taxa.consensus[,8]))
  rownames(taxa.consensus) <- colnames(kmer.table)
  return(taxa.consensus)
}


#' KTU evaluation
#'
#' Evaluate sequence similarity within KTUs
#' @param klusterRDS klustering output RDS file
#' @param ASVfasta representative ASV sequence fasta file
#' @export
KTUsim.eval <- function(klusterRDS,ASVfasta){
  kluster <- readRDS(klusterRDS)
  dna <- Biostrings::readDNAStringSet(ASVfasta,use.names = T)
  dna <- lapply(dna, as.character)
  dna <- sapply(dna, "[[", 1)
  dna <- dna[match(names(kluster$clusters),names(dna))]

  mean.simi <- c()
  Ks <- c()
  kn.stats <- table(kluster$clusters)
  for(k in 1:length(kn.stats)){
    if(kn.stats[k]>1){
      dna.cluster.n <- dna[which(kluster$clusters==k)]
      dna.cluster.cb <- combn(x = 1:length(dna.cluster.n),m = 2)
      dna.seq.cb <- apply(dna.cluster.cb,2,function(x) dna.cluster.n[x])
      within.simi <- c()
      for(i in 1:ncol(dna.seq.cb)) within.simi[i] <- Biostrings::pid(Biostrings::pairwiseAlignment(dna.seq.cb[1,i],dna.seq.cb[2,i]))
      mean.simi[k] <- mean(within.simi)

      tetra <- tetra.freq(dna[which(kluster$clusters==k)],file = F)
      Ks[k] <- mean(as.dist(1-coop::cosine(tetra)))
      Ks[k] <- ifelse(Ks[k]>0,Ks[k],0) # check no negative values
    } else if(kn.stats[k]==1){
      mean.simi[k] <- 100
      Ks[k] <- NA
    }
    print(paste("mean similarity within Kluster",k,"=",mean.simi[k],"% from",kn.stats[k],"seq features; divergence =",Ks[k]))
  }
  g.mean.simi <- mean(mean.simi)
  g.Ks <- mean(Ks[!is.na(Ks)])
  print(paste("Mean similarity of 'within KTU' of", length(kn.stats), "KTUs =",g.mean.simi,"; Mean divergence within KTU =",g.Ks))
  #hist(mean.simi)
  return(list(eachmean=data.frame(kluster=1:length(kn.stats),similarity.pct=mean.simi,divergence=Ks,n.feature=as.data.frame(kn.stats)$Freq),globalmean=g.mean.simi,globaldivergence=g.Ks))
}



