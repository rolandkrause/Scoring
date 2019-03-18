source("dream_scoring_clean.R")
library(FNN)

if (!exists("dm")) initialize()

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

saveRDS(patterns,"patterns.rds")

#team <- team[1]
#pattern1 <- pattern1_all[1]
#pattern2 <- pattern2_all[1]
#pattern3 <- pattern3_all[1]


#result <- c(pattern1,pattern2,pattern3) %>% imap(~score.post.folds(.x,.y))

#write.table(result,paste0(team,".csv"),row.names = F,col.names = apply(expand.grid("s",seq(3),".",seq(3)),1,paste0,collapse=""))

#stats <- c(team, result %>% map(~apply(.x,2,function(r) c(mean(r),sd(r)))) %>% unlist)

#write(stats,"stats.csv",ncolumns = length(stats), append=TRUE, sep=",")

#rank.post.folds("./scores")

#expected jaccard similarity of 2 sets of m element chosen randomly from a set of n elements
expected.jaccard <- function(n,m){
  seq(m) %>%  map_dbl(function(k){
    choose(m,k)*choose(n-m,m-k)/choose(n,m) * k/(2*m - k)
  }) %>% sum
}

pattern1_all <- as.character(patterns$pattern1[-9])
pattern2_all <- as.character(patterns$pattern2[-9])
pattern3_all <- as.character(patterns$pattern3[-9])
team_all <- as.character(patterns$Team[-9])

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

frequencies1 <- sort(table(pattern1_all %>% map_dfc(~raw.selected.genes(.x,1)) %>% apply(1,as.character))/(length(pattern1_all)*10),decreasing=T)
plot(frequencies1, ylab="Frequency", xlab="", las=2, main= "Subchallenge 1", ylim=c(0,1))

frequencies2 <- sort(table(pattern2_all %>% map_dfc(~raw.selected.genes(.x,2)) %>% apply(1,as.character))/(length(pattern2_all)*10),decreasing=T)
plot(frequencies2, ylab="Frequency", xlab="", las=2, main= "Subchallenge 2", ylim=c(0,1))

frequencies3 <- sort(table(pattern3_all %>% map_dfc(~raw.selected.genes(.x,3)) %>% apply(1,as.character))/(length(pattern3_all)*10),decreasing=T)
plot(frequencies3, ylab="Frequency", xlab="", las=2, main= "Subchallenge 3", ylim=c(0,1))



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


legend(50,0.2,c("Subchallenge 1-2","Subchallenge 2-3","Subchallenge 1-3", "random"), col=c("green","blue","red", "black"), lty=c(1,1,1,3))


#Features in subchallenge intersection
write("Top k\tNumber of features in intersection\tFeatures",file = "intersect_features_cv.txt")

seq(10,60) %>% walk(function(k){
  feat <- intersect(intersect(names(frequencies1)[1:k],names(frequencies2)[1:k]),names(frequencies3)[1:k])
  write_delim(tibble(k,length(feat),paste0(feat,collapse = "\t")),path = "intersect_features_cv.txt", col_names = F, append = T)
})


insitus <- dm@insitu.matrix
geom <- data.frame(cbind(dm@geometry,apply(insitus,2,rep,2)))


neighbors <- get.knnx(dm@geometry, dm@geometry, k=11)$nn.index[,-1]

spatial.stats <- colnames(insitus) %>% map_dfr(function(gene){
  probs <- table(insitus[,gene])/nrow(insitus)
  entropy <- sum(-probs * log2(probs))
  
  xbar <- mean(insitus[,gene])
  invvar <- 1/var(insitus[,gene])
  
  #spatial weight 1 for the 10 nearest neighbors
  MoranI <- (seq(nrow(dm@geometry)) %>% map_dbl(function(id1){
    sum((geom[id1,gene] - xbar) * (geom[neighbors[id1,],gene] - xbar))
  }) %>% sum) * invvar*nrow(dm@geometry)/((nrow(dm@geometry)*10)*(nrow(dm@geometry) - 1))
  
  GearyC <- (seq(nrow(dm@geometry)) %>% map_dbl(function(id1){
    sum((geom[id1,gene] - geom[neighbors[id1,],gene])^2)
  }) %>% sum) * invvar/(2*nrow(dm@geometry)*10)
  
  print(gene)
  data.frame(gene = gene, entropy = entropy, MoranI=MoranI, GearyC=GearyC)
})
