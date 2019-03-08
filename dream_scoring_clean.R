#Sets up the enviroment for scoring
#Defines the scoring metrics

library(DistMap)
library(purrr)
library(dplyr)
library(mccr)
library(caret)
library(synapser)
library(stringr)
library(readr)

#Use Attila's code to donwload, load data and initialize the environment for scoring <- syn16782361
initialize <- function(){
  if(!file.exists("init.RData")){
    if(!all(file.exists(c("dge_raw.txt.gz","dge_normalized.txt.gz","binarized_bdtnp.csv.gz","bdtnp.txt.gz","geometry.txt.gz")))){
      download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/dge_raw.txt.gz",destfile = "dge_raw.txt.gz")
      download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/dge_normalized.txt.gz",destfile = "dge_normalized.txt.gz")
      download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/binarized_bdtnp.csv.gz",destfile = "binarized_bdtnp.csv.gz")
      download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/bdtnp.txt.gz",destfile = "bdtnp.txt.gz")
      download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/geometry.txt.gz",destfile = "geometry.txt.gz")
    }
    
    raw.data = read.table(gzfile("dge_raw.txt.gz",'rt'),
                          sep = "\t",
                          row.names = NULL,
                          stringsAsFactors = F,
                          quote = "")
    raw.data.genes = raw.data$V1
    raw.data$V1 = NULL
    
    raw.data.genes = gsub("'","",raw.data.genes,fixed = T)
    
    raw.data = as.matrix(raw.data)
    rownames(raw.data) = raw.data.genes
    
    normalized.data = read.table(gzfile("dge_normalized.txt.gz",'rt'),
                                 sep = "\t",
                                 row.names = NULL,
                                 stringsAsFactors = F,
                                 quote = "")
    
    normalized.data.genes = normalized.data$row.names
    normalized.data$row.names = NULL
    
    normalized.data.genes = gsub("'","",normalized.data.genes,fixed = T)
    
    normalized.data = as.matrix(normalized.data)
    rownames(normalized.data) = normalized.data.genes
    
    stopifnot(all(normalized.data.genes == raw.data.genes))
    
    insitu.matrix = read.table(gzfile("binarized_bdtnp.csv.gz",'rt'), sep = ",",header = T)
    
    insitu.genes_orig <- colnames(insitu.matrix)
    
    #this is not needed for the normalized data
    insitu.genes = gsub(".","-",insitu.genes_orig,fixed = T)
    insitu.genes = gsub("-spl-","(spl)",insitu.genes,fixed = T)
    
    insitu.matrix = as.matrix(insitu.matrix)
    colnames(insitu.matrix) = insitu.genes
    
    stopifnot(all(insitu.genes %in% raw.data.genes))
    
    geometry = read.csv(gzfile("geometry.txt.gz",'rt'),sep = " ")
    
    colnames(geometry) = c("x","y","z")
    
    #close gz properly
    #closeAllConnections()
    
    dm <<- new("DistMap",
             raw.data=raw.data,
             data=normalized.data,
             insitu.matrix=insitu.matrix,
             geometry=as.matrix(geometry))
    
    dm <<- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
    
    dm <<- mapCells(dm)
    #Thank you Attila!
    
    #GROUND TRUTH
    
    ground.truth <<- t(apply(dm@mcc.scores,2,order,decreasing=TRUE))[,1:10]
    ambig.locations <<- t(apply(dm@mcc.scores,2,sort,decreasing=TRUE))[,1:2]
    ambig.locations <<- which(ambig.locations[,1] == ambig.locations[,2])
    
    
    #map every cell to its d84 value
    d84 <<- seq(nrow(ground.truth)) %>% map_dbl(function(j){
      #map every position to the norm of the difference in the geometry and calculate the mean
      ground.truth[j,] %>% map_dbl(~sqrt(sum((dm@geometry[.x,] - dm@geometry[ground.truth[j,1],])^2))) %>% mean
    })
    
    save(dm,ground.truth,d84,ambig.locations,file="init.RData")
    
  } else {
    load("init.RData",envir = .GlobalEnv)
  }
}


#Scoring function. 
#Input: path to the results .csv file (character), number of subchallenge (integer)
#Output: vector of scores (s1,s2,s3)
score <- function(path,sub){
  
  if (!exists("dm")) initialize()
  
  submission <- read.csv(path,header=FALSE,stringsAsFactors = FALSE)
  
  #separate the gene names from the location predictions
  gene.lines <- (4-sub)*2
  genes <- submission %>% slice(1:gene.lines)
  locations <- submission %>% slice(-1:-gene.lines)
  
  #preprocess genes and locations, remove NAs, sort locations by cellid
  genes <- genes %>% select(-1) %>% unlist %>% as.character
  locations <- locations[order(locations[,1]),] %>% select(-1) %>% apply(2,as.numeric)
  
  #do the same mapping for the submission as for d48
  dsub <- seq(nrow(locations)) %>% map_dbl(function(j){
    vals <- locations[j,] %>% map_dbl(~sqrt(sum((dm@geometry[.x,] - dm@geometry[ground.truth[j,1],])^2))) %>% mean
  })
  
  #calculate relative precision
  pk <- d84/dsub
  
  #s1
  
  #select fluorescence data only for the submitted subset of genes
  reduced.insitu <- data.frame(dm@insitu.matrix) %>% select(genes)
  #get binarized data from distmap
  ts <- data.frame(t(dm@binarized.data))
  #select binarized data only for the submitted subset of genes
  reduced.ts <- ts %>% select(genes)
  
  #map every cell location prediction to the MCC between the ground truth location and the predicted most likely position, using the submitted subset of genes 
  mccrs <- seq(nrow(locations)) %>% map_dbl(~mccr(reduced.insitu[ground.truth[.x,1],],reduced.insitu[locations[.x,1],]))
  
  #do not take into account the cells with ambiguous locations
  s1 <- sum(((pk/sum(pk)) * mccrs)[-ambig.locations])
  
  #s2
  #do not take into account the cells with ambiguous locations
  s2<- mean(pk[-ambig.locations])
  
  #s3
  
  #comparing rnaseq and fluorescence data using true locations
  true.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[ground.truth[-ambig.locations,1],.x],reduced.ts[-ambig.locations,.x]))
  #.. using submitted locations
  competitor.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[locations[-ambig.locations,1],.x],reduced.ts[-ambig.locations,.x]))
  
  #do not take into account the cells with ambiguous locations
  s3 <- sum(((true.mccs/sum(true.mccs)) * competitor.mccs))
  
  return(c(s1,s2,s3))
}

#Scoring function with bootstraping
#Input: path to the results .csv file (character), number of subchallenge (integer), number of bootstraps (integer, optional)
#Output: data frame with scores
score.bootstrapped <- function(path,sub,nboot=1000){
  
  if (!exists("dm")) initialize()
  
  submission <- read.csv(path,header=FALSE,stringsAsFactors = FALSE)
  
  #separate the gene names from the location predictions
  gene.lines <- (4-sub)*2
  genes <- submission %>% slice(1:gene.lines)
  locations <- submission %>% slice(-1:-gene.lines)
  
  #preprocess genes and locations, remove NAs, sort locations by cellid
  genes <- genes %>% select(-1) %>% unlist %>% as.character
  locations <- locations[order(locations[,1]),] %>% select(-1) %>% apply(2,as.numeric)
  
  #fix incompatibility
  genes = gsub("-",".",genes,fixed = T)
  genes = gsub("(spl)",".spl.",genes,fixed = T)
  
  #remove ambiguous locations
  locations.n <- locations[-ambig.locations,]
  ground.truth.n <- ground.truth[-ambig.locations,]
  
  #do the same mapping for the submission as for d48
  dsub <- seq(nrow(locations.n)) %>% map_dbl(function(j){
    vals <- locations.n[j,] %>% map_dbl(~sqrt(sum((dm@geometry[.x,] - dm@geometry[ground.truth.n[j,1],])^2))) %>% mean
  })
  
  #calculate relative precision
  pk <- d84[-ambig.locations]/dsub
  
  #s1
  
  #select fluorescence data only for the submitted subset of genes
  reduced.insitu <- data.frame(dm@insitu.matrix) %>% select(genes)
  #get binarized data from distmap and remove abiguous locations
  ts <- data.frame(t(dm@binarized.data))[-ambig.locations,]
  #select binarized data only for the submitted subset of genes
  reduced.ts <- ts %>% select(genes)
  
  #map every cell location prediction to the MCC between the ground truth location and the predicted most likely position, using the submitted subset of genes 
  mccrs <- seq(nrow(locations.n)) %>% map_dbl(~mccr(reduced.insitu[ground.truth.n[.x,1],],reduced.insitu[locations.n[.x,1],]))
  
  #bootstrapping
  samples <- seq(nboot) %>% map_dfr(function(seed){
    set.seed(seed)
    bootstrap <- sample.int(nrow(locations.n),replace=TRUE)
    #s1
    s1.b <- sum((pk[bootstrap]/sum(pk[bootstrap])) * mccrs[bootstrap])
    #s2
    s2.b <- mean(pk[bootstrap])
    
    #s3
    #since we bootstrap by locations, these must be recalculated
    true.mccs.b <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[ground.truth.n[bootstrap,1],.x],reduced.ts[bootstrap,.x]))
    #submitted locations
    competitor.mccs.b <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[locations.n[bootstrap,1],.x],reduced.ts[bootstrap,.x]))
    
    #here i assumed that the denominator is the sum of true mccs
    s3.b <- sum((true.mccs.b/sum(true.mccs.b)) * competitor.mccs.b)
    
    
    data.frame(s1=s1.b,s2=s2.b,s3=s3.b)
  })
  
  return(samples)
}

#wrapper for summarising the bootstrapped scores
score.bootstrapped.summary <- function(path,sub,nboot=1000){
  score.bootstrapped(path,sub,nboot) %>% summarise(mean(s1), sd(s1), mean(s2), sd(s2), mean(s3), sd(s3)) %>% as.numeric
}

#computes the Bayes factor between two results using bootstrapped scores
bayes.bootstrap <- function(path1,path2,sub,nboot=1000){
  
  samples1 <- score.bootstrapped(path1,sub,nboot)
  samples2 <- score.bootstrapped(path2,sub,nboot)
  
  wins <- colSums(samples1 >= samples2)
  B <- wins/(nboot-wins)
  
  return(B)
}

#input: a csv file with a sid column and a team column, requires synapse login
bootstraped.ranks <- function(submissions,sub){
  
  s <- read.csv(submissions,stringsAsFactors = F)
  s$sid <- as.character(s$sid)
  
  #if not logged in synapse, load in the variable files the paths to the submissions
  files <- s$sid %>% map_chr(~synGetSubmission(.x)$filePath)
  
  #evaluate
  eval.boot <-files %>% map(~score.bootstrapped(.x,sub))
  
  #rank on each score separately and reduce to sum
  #need to make scores negative in order to properly rank them
  ranks <- seq(3) %>% map(function(score){
    ranks<- eval.boot %>% map_dfc(~-.x[,score]) %>% apply(1,rank) %>% t
    colnames(ranks) <- s$team  
    return(ranks)
  }) %>% reduce(`+`)
  
  ranks <- (ranks/3) %>% apply(1,rank) %>% t
  save(ranks,file=paste0("sc",sub,"ranks.Rdata"))
  
  #draw the boxplot
  avg.ranks <- ranks %>% colMeans %>% rank  
  ordering <- order(avg.ranks)
  pdf(paste0("sc",sub,"_final_boxplot.pdf"), width=11, height=8)
  par(mar=c(5,10.5,4,2) + 0.1)
  boxplot(ranks[,ordering],horizontal=T,las=2,at=rev(1:ncol(ranks)), xlab="Rank")
  
  factors <- map2_dfr(ordering[-length(ordering)],ordering[-1],function(c1,c2){
    win <- sum(ranks[,c1] < ranks[,c2])
    lose <- sum(ranks[,c2] < ranks[,c1])
    BF <- win/lose
    data.frame(c1=c1,c2=c2,BF=BF)
  })
  
  abh <- ncol(ranks) - factors$c1[which(factors$BF>=3)[1]] + 0.5
  abline(h=abh,lwd=2)
  
  dev.off()
  
  #report final ranking
  result <- mutate(s, rank=avg.ranks) %>% arrange(rank)
  write.csv(result,file=paste0("sc",sub,"_final_table.csv"),row.names=F)
  
  return(mutate(s, rank=avg.ranks) %>% arrange(rank))
}


folds <- function(k=10){
  if (!exists("dm")) initialize()
  
  set.seed(42)
  folds.train <- createFolds(ground.truth[,1],k = k,returnTrain = T)
  set.seed(42)
  folds.test <- createFolds(ground.truth[,1],k = k, returnTrain = F)
  
  walk(folds.train, ~cat(.x,"\n", file="folds_train.csv", sep=",", append=T))
  walk(folds.test, ~cat(.x,"\n", file="folds_test.csv", sep=",", append=T))
  
  ids <- colnames(dm@data)
  walk(folds.train, ~cat(.x,"\n", file="folds_train.csv", sep=",", append=T))
  
  write_csv(data.frame(id=seq(1297),cell=ids),path="cell_ids.csv")
  
  #the union of the test locations is the set of all locations
  #submit 10 files (one per fold) per subchallenge
}


#score function for the 10 fold crossvalidation
#the pattern should be provided in glue format, where the fold number is represented by {.} or {.x}
score.post.folds <- function(pattern, sub){
  if (!exists("dm")) initialize()
  
  #read all folds
  all.submissions <- seq(10) %>% map(~read.csv(str_glue(pattern),header=FALSE,stringsAsFactors = FALSE))
  
  #for each fold get scores
  all.submissions %>% map_dfr(function(submission){
    #separate the gene names from the location predictions
    gene.lines <- (4-sub)*2
    genes <- submission %>% slice(1:gene.lines)
    locations <- submission %>% slice(-1:-gene.lines)
    #BCBU hack
    #locations <- locations[which(locations[,2] != ""),]
    
    #preprocess genes and locations, remove NAs, sort locations by cellid
    #mutations cause cancer, but anyway ...
    genes <- genes %>% select(-1) %>% unlist %>% as.character
    
    #fix incompatibility
    genes = gsub("-",".",genes,fixed = T)
    genes = gsub("(spl)",".spl.",genes,fixed = T)
    
    locations <- locations[order(locations[,1]),] 
    cell.ids <- locations[,1]
    locations <- locations %>% select(-1) %>% apply(2,as.numeric)
    
    #select from d84 the falues from the current test fold 
    d84.fold <- d84[cell.ids]
    
    #do the same mapping for the submission as for d84
    dsub <- cell.ids %>% imap_dbl(function(j,iter){
      vals <- locations[iter,] %>% map_dbl(~sqrt(sum((dm@geometry[.x,] - dm@geometry[ground.truth[j,1],])^2))) %>% mean
    })
    
    #calculate relative precision
    pk <- d84.fold/dsub
    
    #s1
    
    #select fluorescence data only for the submitted subset of genes
    reduced.insitu <- data.frame(dm@insitu.matrix) %>% select(genes)
    #get binarized data from distmap
    ts <- data.frame(t(dm@binarized.data))
    #select binarized data only for the submitted subset of genes
    reduced.ts <- ts %>% select(genes)
    
    
    #map every cell location prediction to the MCC between the ground truth location and the predicted most likely position, using the submitted subset of genes 
    mccrs <- cell.ids %>% imap_dbl(~mccr(reduced.insitu[ground.truth[.x,1],],reduced.insitu[locations[.y,1],]))
    
    ambig.indexes <- which(cell.ids %in% ambig.locations)
    
    #do not take into account the cells with ambiguous locations
    s1 <- sum(((pk/sum(pk)) * mccrs)[-ambig.indexes])
    
    #s2
    #do not take into account the cells with ambiguous locations
    s2 <- mean(pk[-ambig.indexes])
    
    #s3
    
    #comparing rnaseq and fluorescence data using true locations
    true.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[ground.truth[cell.ids[-ambig.indexes],1],.x],reduced.ts[cell.ids[-ambig.indexes],.x]))
    #.. using submitted locations
    competitor.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[locations[-ambig.indexes,1],.x],reduced.ts[cell.ids[-ambig.indexes],.x]))
    
    #do not take into account the cells with ambiguous locations
    s3 <- sum(((true.mccs/sum(true.mccs)) * competitor.mccs))
    
    data.frame(s1=s1,s2=s2,s3=s3)
    
  })
  
}


#random competitor is random
#selects a random subset of genes and assigns 10 random locations for each cell
random.competitor <- function(seed=42){

  set.seed(seed)
  
  isgenes <- colnames(dm@insitu.matrix)
  
  selected.genes <- sample(isgenes, 60)
  result.locations <- t(replicate(ncol(dm@data), sample(nrow(dm@insitu.matrix),10))) %>% data.frame
  
  result60 <- list(genes = selected.genes, locations = result.locations)
  
  selected.genes <- sample(isgenes, 40)
  result.locations <- t(replicate(ncol(dm@data), sample(nrow(dm@insitu.matrix),10))) %>% data.frame
  
  
  result40 <- list(genes = selected.genes, locations = result.locations)
  
  selected.genes <- sample(isgenes, 20)
  result.locations <- t(replicate(ncol(dm@data), sample(nrow(dm@insitu.matrix),10))) %>% data.frame
  
  result20 <- list(genes = selected.genes, locations = result.locations)
  
  return(list(r60=result60, r40=result40, r20=result20))
  
}


#the path points to the folder with csv results
rank.post.folds <- function(path){
  #read calculated scores
  scores <- list.files(path,".csv$",full.names = T) %>% map_dfr(~read_delim(.,delim = " ",col_types=cols()) %>% mutate(fold=seq(10), team=str_extract(.x,"([a-zA-Z0-9 _])+(?=.csv)")))
  #rank
  ranks <- scores %>% group_by(fold) %>% mutate_if(is.numeric,"-") %>% mutate_if(is.numeric,rank,na.last="keep") %>% mutate(sub1 = (s1.1+s2.1+s3.1)/3,sub2 = (s1.2+s2.2+s3.2)/3, sub3 = (s1.3+s2.3+s3.3)/3)
  #summarize sub
  sub.summary <- ranks %>% ungroup %>% group_by(team) %>% select(contains("sub"))
  
  #final ranking + draw boxplots
  avranks <- sub.summary %>% summarise(ar1 = mean(sub1),ar2=mean(sub2),ar3=mean(sub3))
  
  sub.summary %<>% ungroup %>% mutate(ar1 = rep(avranks$ar1,each=10), ar2 = rep(avranks$ar2,each=10), ar3 = rep(avranks$ar3,each=10))
  
  #avoid the same average rank
  jitter <- rep(rnorm(length(unique(sub.summary$team)),0,1e-8),each=10)
  
  sub.summary$ar1 <- sub.summary$ar1 + jitter
  sub.summary$ar2 <- sub.summary$ar2 + jitter
  sub.summary$ar3 <- sub.summary$ar3 + jitter
  
  par(mar=c(5,10.5,4,2) + 0.1)
  boxplot(sub1~ar1,sub.summary,horizontal=T,las=2,xlim=c(13,0),xlab="Rank",names=avranks$team[order(avranks$ar1)], main="Subchallenge 1")
  boxplot(sub2~ar2,sub.summary,horizontal=T,las=2,xlim=c(13,0),xlab="Rank",names=avranks$team[order(avranks$ar2,na.last=NA)], main="Subchallenge 2")
  boxplot(sub3~ar3,sub.summary,horizontal=T,las=2,xlim=c(13,0),xlab="Rank",names=avranks$team[order(avranks$ar3,na.last=NA)], main="Subchallenge 3")
  
}


jaccard.genes <- function(pattern, sub){
  all.submissions <- seq(10) %>% map(~read.csv(str_glue(pattern),header=FALSE,stringsAsFactors = FALSE))
  
  #for each fold get genes
  selected.genes <- all.submissions %>% map_dfc(function(submission){
    #separate the gene names from the location predictions
    gene.lines <- (4-sub)*2
    genes <- submission %>% slice(1:gene.lines) %>% select(-1) %>% unlist %>% make.names
  
    data.frame(f = genes %>% as.character)
  })
  
  combn(seq(ncol(selected.genes)), 2, function(x){
    length(intersect(selected.genes[,x[1]],selected.genes[,x[2]]))/length(union(selected.genes[,x[1]],selected.genes[,x[2]]))
  })
}


team_all <- c("BCBU", "Challengers18", "Challengers18 Official", "Christoph Hafemeister", "DeepCMC", "MLB", "NAD", "OmicsEngineering", "random", "Thin Nguyen", "WhatATeam", "Zho")

pattern1_all <- c("post/BCBU_bin/BCBU_60_{.}.csv", "post/challengers18_CV/60genes_{.}.csv", "post/challengers18_official_submisison_CV/60genes_{.}.csv", "post/cv_results_christoph_hafemeister/60genes_{.}.csv", "post/DeepCMC_SubmissionFiles_CV_DistMapTrainTest/60genes_a_{.}_test.csv", 
                  "post/Team_MLB_10foldCV_predictions_v2/60genes_{.}.csv", "post/NAD_Result/60genes_{.}.csv", "post/OmicsEngineering_new/subchallenge1_submission/gene60_testfold{.}_index.csv", "post/random/random_{.}_.r60.genes.csv", "post/thinnguyen_ten_fold_binarized/60genes_{.}.csv",
                  "post/WhatATeam/60genes_{.}.csv", "post/post_challenge_zho_team/results/60genes_{.}_.csv")
pattern2_all <- c("post/BCBU_bin/BCBU_40_{.}.csv", "post/challengers18_CV/40genes_{.}.csv", "post/challengers18_official_submisison_CV/40genes_{.}.csv", "post/cv_results_christoph_hafemeister/40genes_{.}.csv", "post/DeepCMC_SubmissionFiles_CV_DistMapTrainTest/40genes_a_{.}_test.csv", 
                  "post/Team_MLB_10foldCV_predictions_v2/40genes_{.}.csv", "post/NAD_Result/40genes_{.}.csv", "post/OmicsEngineering_new/subchallenge2_submission/gene40_testfold{.}_index.csv", "post/random/random_{.}_.r40.genes.csv", "post/thinnguyen_ten_fold_binarized/40genes_{.}.csv",
                  "post/WhatATeam/40genes_{.}.csv", "post/post_challenge_zho_team/results/40genes_{.}_.csv")
pattern3_all <- c("post/BCBU_bin/BCBU_20_{.}.csv", "post/challengers18_CV/20genes_{.}.csv", "post/challengers18_official_submisison_CV/20genes_{.}.csv", "post/cv_results_christoph_hafemeister/20genes_{.}.csv", "post/DeepCMC_SubmissionFiles_CV_DistMapTrainTest/20genes_a_{.}_test.csv", 
                  "post/Team_MLB_10foldCV_predictions_v2/20genes_{.}.csv", "post/NAD_Result/20genes_{.}.csv", "post/OmicsEngineering_new/subchallenge3_submission/gene20_testfold{.}_index.csv", "post/random/random_{.}_.r20.genes.csv", "post/thinnguyen_ten_fold_binarized/20genes_{.}.csv",
                  "post/WhatATeam/20genes_{.}.csv", "post/post_challenge_zho_team/results/20genes_{.}_.csv")

patterns <- data.frame(Team = team_all, pattern1 = pattern1_all, pattern2 = pattern2_all, pattern3 = pattern3_all)

team <- team[1]
pattern1 <- pattern1_all[1]
pattern2 <- pattern2_all[1]
pattern3 <- pattern3_all[1]


result <- c(pattern1,pattern2,pattern3) %>% imap(~score.post.folds(.x,.y))

write.table(result,paste0(team,".csv"),row.names = F,col.names = apply(expand.grid("s",seq(3),".",seq(3)),1,paste0,collapse=""))

stats <- c(team, result %>% map(~apply(.x,2,function(r) c(mean(r),sd(r)))) %>% unlist)

write(stats,"stats.csv",ncolumns = length(stats), append=TRUE, sep=",")

rank.post.folds("./scores")

#expected jaccard similarity of 2 sets of m element chosen randomly from a set of n elements
expected.jaccard <- function(n,m){
  seq(m) %>%  map_dbl(function(k){
    choose(m,k)*choose(n-m,m-k)/choose(n,m) * k/(2*m - k)
  }) %>% sum
}

.pardefault <- par()
par(mfrow=c(1,3), mar = c(5,2,4,2) + 0.1, oma = c(1,8,1,1))

sub1_genes <- map2_dfc(team_all, pattern1_all, function(x,y){
  team.stats <- data.frame(jaccard.genes(y,1))
  colnames(team.stats)[1] <- x
  team.stats
})

sub1_genes <- sub1_genes[,rev(colnames(sub1_genes))]

boxplot(sub1_genes, horizontal = T, ylim = c(0.1,1),las=2, xlab= "Jaccard similarity", main="Subchallenge 1")
abline(v=expected.jaccard(84,60), lty=3)

sub2_genes <- map2_dfc(team_all, pattern2_all, function(x,y){
  team.stats <- data.frame(jaccard.genes(y,2))
  colnames(team.stats)[1] <- x
  team.stats
})
sub2_genes <- sub2_genes[,rev(colnames(sub2_genes))]

boxplot(sub2_genes, horizontal = T, ylim = c(0.1,1),las=2, xlab= "Jaccard similarity", main="Subchallenge 2", names = F)
abline(v=expected.jaccard(84,40), lty=3)

sub3_genes <- map2_dfc(team_all, pattern3_all, function(x,y){
  team.stats <- data.frame(jaccard.genes(y,3))
  colnames(team.stats)[1] <- x
  team.stats
})

sub3_genes <- sub3_genes[,rev(colnames(sub3_genes))]

boxplot(sub3_genes, horizontal = T, ylim = c(0.1,1),las=2, xlab= "Jaccard similarity", main="Subchallenge 3", names = F)
abline(v=expected.jaccard(84,20), lty=3)

raw.selected.genes <- function(pattern, sub){
  all.submissions <- seq(10) %>% map(~read.csv(str_glue(pattern),header=FALSE,stringsAsFactors = FALSE))
  
  #for each fold get genes
  selected.genes <- all.submissions %>% map_dfc(function(submission){
    #separate the gene names from the location predictions
    gene.lines <- (4-sub)*2
    genes <- submission %>% slice(1:gene.lines) %>% select(-1) %>% unlist %>% make.names
    
    data.frame(f = genes %>% as.character)
  })
  
  return(selected.genes)
}

par(.pardefault)

frequencies1 <- sort(table(pattern1_all %>% map_dfc(~raw.selected.genes(.x,1)) %>% apply(1,as.character))/130,decreasing=T)
plot(frequencies1, ylab="Frequency", xlab="", las=2, main= "Subchallenge 1")

frequencies2 <- sort(table(pattern2_all %>% map_dfc(~raw.selected.genes(.x,2)) %>% apply(1,as.character))/130,decreasing=T)
plot(frequencies2, ylab="Frequency", xlab="", las=2, main= "Subchallenge 2")

frequencies3 <- sort(table(pattern3_all %>% map_dfc(~raw.selected.genes(.x,3)) %>% apply(1,as.character))/130,decreasing=T)
plot(frequencies3, ylab="Frequency", xlab="", las=2, main= "Subchallenge 3")


#top-k plot

overlap12 <- seq(min(length(frequencies1),length(frequencies2))) %>% map_dbl(function(c){
  s1 <- names(frequencies1)[1:c]
  s2 <- names(frequencies2)[1:c]
  s1
  length(intersect(s1,s2))/length(union(s1,s2))
})

overlap23 <- seq(min(length(frequencies2),length(frequencies3)))  %>% map_dbl(function(c){
  s1 <- names(frequencies2)[1:c]
  s2 <- names(frequencies3)[1:c]
  length(intersect(s1,s2))/length(union(s1,s2))
})

overlap13 <- seq(min(length(frequencies1),length(frequencies3)))  %>% map_dbl(function(c){
  s1 <- names(frequencies1)[1:c]
  s2 <- names(frequencies3)[1:c]
  length(intersect(s1,s2))/length(union(s1,s2))
})

plot(overlap12, type="l", col="green", ylab="Jaccard", xlab="Top K selected genes by frequency", ylim=c(0,1))
lines(overlap23, col="blue")
lines(overlap13, col="red")

lines(seq(84) %>% map_dbl(~expected.jaccard(84,.x)), lty=3)


legend(60,0.35,c("Subchallenge 1-2","Subchallenge 2-3","Subchallenge 1-3", "random"), col=c("green","blue","red", "black"), lty=c(1,1,1,3))
