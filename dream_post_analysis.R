source("dream_scoring_clean.R")
library(FNN)
library(gridExtra)

if (!exists("dm")) initialize()

jaccard.genes <- function(pattern, sub) {
  all.submissions <- seq(10) %>% map(~ read.csv(str_glue(pattern), header = FALSE, stringsAsFactors = FALSE))

  # for each fold get genes
  selected.genes <- all.submissions %>% map_dfc(function(submission) {
    # separate the gene names from the location predictions
    gene.lines <- (4 - sub) * 2
    genes <- submission %>% slice(1:gene.lines) %>% select(-1) %>% unlist() %>% make.names()

    data.frame(f = genes %>% as.character())
  })

  combn(seq(ncol(selected.genes)), 2, function(x) {
    length(intersect(selected.genes[, x[1]], selected.genes[, x[2]])) / length(union(selected.genes[, x[1]], selected.genes[, x[2]]))
  })
}


team_all <- c("BCBU", "Challengers18", "Challengers18 Official", "Christoph Hafemeister", "DeepCMC", "MLB", "NAD", "OmicsEngineering", "random", "Thin Nguyen", "WhatATeam", "Zho")

pattern1_all <- c(
  "post/BCBU_bin/BCBU_60_{.}.csv", "post/challengers18_CV/60genes_{.}.csv", "post/challengers18_official_submisison_CV/60genes_{.}.csv", "post/cv_results_christoph_hafemeister/60genes_{.}.csv", "post/DeepCMC_SubmissionFiles_CV_DistMapTrainTest/60genes_a_{.}_test.csv",
  "post/Team_MLB_10foldCV_predictions_v2/60genes_{.}.csv", "post/NAD_Result/60genes_{.}.csv", "post/OmicsEngineering_new/subchallenge1_submission/gene60_testfold{.}_index.csv", "post/random/random_{.}_.r60.genes.csv", "post/thinnguyen_ten_fold_binarized/60genes_{.}.csv",
  "post/WhatATeam/60genes_{.}.csv", "post/post_challenge_zho_team/results/60genes_{.}_.csv"
)
pattern2_all <- c(
  "post/BCBU_bin/BCBU_40_{.}.csv", "post/challengers18_CV/40genes_{.}.csv", "post/challengers18_official_submisison_CV/40genes_{.}.csv", "post/cv_results_christoph_hafemeister/40genes_{.}.csv", "post/DeepCMC_SubmissionFiles_CV_DistMapTrainTest/40genes_a_{.}_test.csv",
  "post/Team_MLB_10foldCV_predictions_v2/40genes_{.}.csv", "post/NAD_Result/40genes_{.}.csv", "post/OmicsEngineering_new/subchallenge2_submission/gene40_testfold{.}_index.csv", "post/random/random_{.}_.r40.genes.csv", "post/thinnguyen_ten_fold_binarized/40genes_{.}.csv",
  "post/WhatATeam/40genes_{.}.csv", "post/post_challenge_zho_team/results/40genes_{.}_.csv"
)
pattern3_all <- c(
  "post/BCBU_bin/BCBU_20_{.}.csv", "post/challengers18_CV/20genes_{.}.csv", "post/challengers18_official_submisison_CV/20genes_{.}.csv", "post/cv_results_christoph_hafemeister/20genes_{.}.csv", "post/DeepCMC_SubmissionFiles_CV_DistMapTrainTest/20genes_a_{.}_test.csv",
  "post/Team_MLB_10foldCV_predictions_v2/20genes_{.}.csv", "post/NAD_Result/20genes_{.}.csv", "post/OmicsEngineering_new/subchallenge3_submission/gene20_testfold{.}_index.csv", "post/random/random_{.}_.r20.genes.csv", "post/thinnguyen_ten_fold_binarized/20genes_{.}.csv",
  "post/WhatATeam/20genes_{.}.csv", "post/post_challenge_zho_team/results/20genes_{.}_.csv"
)

patterns <- data.frame(Team = team_all, pattern1 = pattern1_all, pattern2 = pattern2_all, pattern3 = pattern3_all)

saveRDS(patterns, "patterns.rds")

# team <- team[1]
# pattern1 <- pattern1_all[1]
# pattern2 <- pattern2_all[1]
# pattern3 <- pattern3_all[1]


# result <- c(pattern1,pattern2,pattern3) %>% imap(~score.post.folds(.x,.y))

# write.table(result,paste0(team,".csv"),row.names = F,col.names = apply(expand.grid("s",seq(3),".",seq(3)),1,paste0,collapse=""))

# stats <- c(team, result %>% map(~apply(.x,2,function(r) c(mean(r),sd(r)))) %>% unlist)

# write(stats,"stats.csv",ncolumns = length(stats), append=TRUE, sep=",")

# rank.post.folds("./scores")

# expected jaccard similarity of 2 sets of m element chosen randomly from a set of n elements
expected.jaccard <- function(n, m) {
  seq(m) %>% map_dbl(function(k) {
    choose(m, k) * choose(n - m, m - k) / choose(n, m) * k / (2 * m - k)
  }) %>% sum()
}

sub1_genes <- map2_dfc(team_all, pattern1_all, function(x, y) {
  team.stats <- data.frame(jaccard.genes(y, 1))
  colnames(team.stats)[1] <- x
  team.stats
})

sub1_genes <- sub1_genes[, rev(colnames(sub1_genes))]
s1g <- ggplot(sub1_genes %>% select(-random) %>% gather(key = "team", value = "jac", factor_key = TRUE), aes(x = team, y = jac)) + geom_boxplot() +
  geom_hline(yintercept = expected.jaccard(84, 60), linetype = "dashed") +
  ylim(0.1, 1) + labs(title = "Subchallenge 1", y = "Jaccard similarity", x = element_blank()) + theme_classic() + coord_flip()


sub2_genes <- map2_dfc(team_all, pattern2_all, function(x, y) {
  team.stats <- data.frame(jaccard.genes(y, 2))
  colnames(team.stats)[1] <- x
  team.stats
})

sub2_genes <- sub2_genes[, rev(colnames(sub2_genes))]
s2g <- ggplot(sub2_genes %>% select(-random) %>% gather(key = "team", value = "jac", factor_key = TRUE), aes(x = team, y = jac)) + geom_boxplot() +
  geom_hline(yintercept = expected.jaccard(84, 40), linetype = "dashed") +
  ylim(0.1, 1) + labs(title = "Subchallenge 2", y = "Jaccard similarity", x = element_blank()) + theme_classic() + coord_flip()


sub3_genes <- map2_dfc(team_all, pattern3_all, function(x, y) {
  team.stats <- data.frame(jaccard.genes(y, 3))
  colnames(team.stats)[1] <- x
  team.stats
})

sub3_genes <- sub3_genes[, rev(colnames(sub3_genes))]
s3g <- ggplot(sub3_genes %>% select(-random) %>% gather(key = "team", value = "jac", factor_key = TRUE), aes(x = team, y = jac)) + geom_boxplot() +
  geom_hline(yintercept = expected.jaccard(84, 20), linetype = "dashed") +
  ylim(0.1, 1) + labs(title = "Subchallenge 3", y = "Jaccard similarity", x = element_blank()) + theme_classic() + coord_flip()

grid.arrange(s1g, s2g, s3g, ncol = 3)

raw.selected.genes <- function(pattern, sub) {
  all.submissions <- seq(10) %>% map(~ read.csv(str_glue(pattern), header = FALSE, stringsAsFactors = FALSE))

  # for each fold get genes
  selected.genes <- all.submissions %>% map_dfc(function(submission) {
    # separate the gene names from the location predictions
    gene.lines <- (4 - sub) * 2
    genes <- submission %>% slice(1:gene.lines) %>% select(-1) %>% unlist() %>% make.names()

    data.frame(f = genes %>% as.character())
  })

  return(selected.genes)
}

par(.pardefault)

frequencies1 <- sort(table(pattern1_all %>% map_dfc(~ raw.selected.genes(.x, 1)) %>% apply(1, as.character)) / (length(pattern1_all) * 10), decreasing = T)


ggplot(tibble(genes = factor(names(frequencies1), levels = names(frequencies1)), freq = as.numeric(frequencies1)), aes(x = genes, y = freq)) +
  geom_point(stat = "identity", shape = 21, size = 3, fill = alpha("blue", 0.4)) + geom_segment(aes(x = genes, xend = genes, y = 0, yend = freq - 0.005)) +
  ylim(0, 1) + labs(title = "Subchallenge 1", x = element_blank(), y = "Frequency") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

write(names(frequencies1)[1:60], "most.freq.txt", ncolumns = 60, append = T)

frequencies2 <- sort(table(pattern2_all %>% map_dfc(~ raw.selected.genes(.x, 2)) %>% apply(1, as.character)) / (length(pattern2_all) * 10), decreasing = T)

ggplot(tibble(genes = factor(names(frequencies2), levels = names(frequencies2)), freq = as.numeric(frequencies2)), aes(x = genes, y = freq)) +
  geom_point(stat = "identity", shape = 21, size = 3, fill = alpha("green", 0.4)) + geom_segment(aes(x = genes, xend = genes, y = 0, yend = freq - 0.005)) +
  ylim(0, 1) + labs(title = "Subchallenge 2", x = element_blank(), y = "Frequency") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

write(names(frequencies2)[1:40], "most.freq.txt", ncolumns = 40, append = T)

frequencies3 <- sort(table(pattern3_all %>% map_dfc(~ raw.selected.genes(.x, 3)) %>% apply(1, as.character)) / (length(pattern3_all) * 10), decreasing = T)

ggplot(tibble(genes = factor(names(frequencies3), levels = names(frequencies3)), freq = as.numeric(frequencies3)), aes(x = genes, y = freq)) +
  geom_point(stat = "identity", shape = 21, size = 3, fill = alpha("red", 0.4)) + geom_segment(aes(x = genes, xend = genes, y = 0, yend = freq - 0.005)) +
  ylim(0, 1) + labs(title = "Subchallenge 3", x = element_blank(), y = "Frequency") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

write(names(frequencies3)[1:20], "most.freq.txt", ncolumns = 20, append = T)


# top-k plot

overlap12 <- seq(min(length(frequencies1), length(frequencies2))) %>% map_dbl(function(c) {
  s1 <- names(frequencies1)[1:c]
  s2 <- names(frequencies2)[1:c]
  length(intersect(s1, s2)) / length(union(s1, s2))
})

overlap23 <- seq(min(length(frequencies2), length(frequencies3))) %>% map_dbl(function(c) {
  s1 <- names(frequencies2)[1:c]
  s2 <- names(frequencies3)[1:c]
  length(intersect(s1, s2)) / length(union(s1, s2))
})

overlap13 <- seq(min(length(frequencies1), length(frequencies3))) %>% map_dbl(function(c) {
  s1 <- names(frequencies1)[1:c]
  s2 <- names(frequencies3)[1:c]
  length(intersect(s1, s2)) / length(union(s1, s2))
})

overlap12 <- c(overlap12, rep(NA, 84 - (length(overlap12) %% 84)))
overlap23 <- c(overlap23, rep(NA, 84 - (length(overlap23) %% 84)))
overlap13 <- c(overlap13, rep(NA, 84 - (length(overlap13) %% 84)))

ggplot(tibble(topk = seq(84), overlap12, overlap23, overlap13), aes(x = topk)) +
  geom_line(aes(y = overlap12, color = "a")) + geom_line(aes(y = overlap23, color = "b")) + geom_line(aes(y = overlap13, color = "c")) +
  geom_line(aes(y = seq(84) %>% map_dbl(~ expected.jaccard(84, .x)), color = "d"), linetype = "dashed") +
  labs(x = "Top k genes per subchallenge", y = "Jaccard similarity") +
  scale_colour_manual(name = "", values = c("a" = alpha("blue", 0.6), "b" = alpha("green", 0.6), "c" = alpha("red", 0.6), d = "gray"), labels = c("Subchallenge 1-2", "Subchallenge 2-3", "Subchallenge 1-3", "Random")) +
  theme_classic() #+ theme(legend.position = c(50,0.2))


# Features in subchallenge intersection
write("Top k\tNumber of features in intersection\tFeatures", file = "intersect_features_cv.txt")

seq(10, 60) %>% walk(function(k) {
  feat <- intersect(intersect(names(frequencies1)[1:k], names(frequencies2)[1:k]), names(frequencies3)[1:k])
  write_delim(tibble(k, length(feat), paste0(feat, collapse = "\t")), path = "intersect_features_cv.txt", col_names = F, append = T)
})


insitus <- dm@insitu.matrix
geom <- data.frame(cbind(dm@geometry, apply(insitus, 2, rep, 2)))


neighbors <- get.knnx(dm@geometry, dm@geometry, k = 11)$nn.index[, -1]

spatial.stats <- colnames(insitus) %>% map_dfr(function(gene) {
  probs <- table(insitus[, gene]) / nrow(insitus)
  entropy <- sum(-probs * log2(probs))

  xbar <- mean(insitus[, gene])
  invvar <- 1 / var(insitus[, gene])

  # spatial weight 1 for the 10 nearest neighbors
  MoranI <- (seq(nrow(dm@geometry)) %>% map_dbl(function(id1) {
    sum((geom[id1, gene] - xbar) * (geom[neighbors[id1, ], gene] - xbar))
  }) %>% sum()) * invvar * nrow(dm@geometry) / ((nrow(dm@geometry) * 10) * (nrow(dm@geometry) - 1))

  GearyC <- (seq(nrow(dm@geometry)) %>% map_dbl(function(id1) {
    sum((geom[id1, gene] - geom[neighbors[id1, ], gene])^2)
  }) %>% sum()) * invvar / (2 * nrow(dm@geometry) * 10)

  print(gene)
  data.frame(gene = gene, entropy = entropy, MoranI = MoranI, GearyC = GearyC)
})


challenge.frequencies1 <- list.files("challenge/final_submissions/", "60genes", full.names = TRUE) %>% 
  map(~read_csv(.x, col_names = F, n_max = 6, col_types = cols()) %>% select(-1) %>% unlist %>% make.names()) %>%
  unlist %>% table %>% sort(decreasing = T)


ggplot(tibble(genes = factor(names(challenge.frequencies1), levels = names(challenge.frequencies1)), freq = as.numeric(challenge.frequencies1/34)), aes(x = genes, y = freq)) +
  geom_point(stat = "identity", shape = 21, size = 3, fill = alpha("blue", 0.4)) + geom_segment(aes(x = genes, xend = genes, y = 0, yend = freq - 0.005)) +
  ylim(0, 1) + labs(title = "Subchallenge 1", x = element_blank(), y = "Frequency") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

challenge.frequencies2 <- list.files("challenge/final_submissions/", "40genes", full.names = TRUE) %>% 
  map(~read_csv(.x, col_names = F, n_max = 4, col_types = cols()) %>% select(-1) %>% unlist %>% make.names()) %>%
  unlist %>% table %>% sort(decreasing = T)


ggplot(tibble(genes = factor(names(challenge.frequencies2), levels = names(challenge.frequencies2)), freq = as.numeric(challenge.frequencies2/34)), aes(x = genes, y = freq)) +
  geom_point(stat = "identity", shape = 21, size = 3, fill = alpha("green", 0.4)) + geom_segment(aes(x = genes, xend = genes, y = 0, yend = freq - 0.005)) +
  ylim(0, 1) + labs(title = "Subchallenge 2", x = element_blank(), y = "Frequency") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

challenge.frequencies3 <- list.files("challenge/final_submissions/", "20genes", full.names = TRUE) %>% 
  map(~read_csv(.x, col_names = F, n_max = 2, col_types = cols()) %>% select(-1) %>% unlist %>% make.names()) %>%
  unlist %>% table %>% sort(decreasing = T)


write(names(challenge.frequencies1)[1:60], "post analysis/chal.most.freq.txt", ncolumns = 60, append = T)
write(names(challenge.frequencies2)[1:40], "post analysis/chal.most.freq.txt", ncolumns = 40, append = T)
write(names(challenge.frequencies3)[1:20], "post analysis/chal.most.freq.txt", ncolumns = 20, append = T)

ggplot(tibble(genes = factor(names(challenge.frequencies3), levels = names(challenge.frequencies3)), freq = as.numeric(challenge.frequencies3/34)), aes(x = genes, y = freq)) +
  geom_point(stat = "identity", shape = 21, size = 3, fill = alpha("red", 0.4)) + geom_segment(aes(x = genes, xend = genes, y = 0, yend = freq - 0.005)) +
  ylim(0, 1) + labs(title = "Subchallenge 3", x = element_blank(), y = "Frequency") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))


prepost1 <- seq(min(length(frequencies1), length(challenge.frequencies1))) %>% map_dbl(function(c) {
  s1 <- names(frequencies1)[1:c]
  s2 <- names(challenge.frequencies1)[1:c]
  length(intersect(s1, s2)) / length(union(s1, s2))
})

prepost2 <- seq(min(length(frequencies2), length(challenge.frequencies2))) %>% map_dbl(function(c) {
  s1 <- names(frequencies2)[1:c]
  s2 <- names(challenge.frequencies2)[1:c]
  length(intersect(s1, s2)) / length(union(s1, s2))
})

prepost3 <- seq(min(length(frequencies3), length(challenge.frequencies3))) %>% map_dbl(function(c) {
  s1 <- names(frequencies3)[1:c]
  s2 <- names(challenge.frequencies3)[1:c]
  length(intersect(s1, s2)) / length(union(s1, s2))
})


prepost2 <- c(prepost2, rep(NA, 84 - (length(prepost2) %% 84)))
prepost3 <- c(prepost3, rep(NA, 84 - (length(prepost3) %% 84)))


ggplot(tibble(topk = seq(84), prepost1, prepost2, prepost3), aes(x = topk)) +
  geom_line(aes(y = prepost1, color = "a")) + geom_line(aes(y = prepost2, color = "b")) + geom_line(aes(y = prepost3, color = "c")) +
  geom_line(aes(y = seq(84) %>% map_dbl(~ expected.jaccard(84, .x)), color = "d"), linetype = "dashed") +
  labs(x = "Top k genes per subchallenge", y = "Jaccard similarity") +
  scale_colour_manual(name = "", values = c("a" = alpha("blue", 0.6), "b" = alpha("green", 0.6), "c" = alpha("red", 0.6), d = "gray"), labels = c("Subchallenge 1", "Subchallenge 2", "Subchallenge 3", "Random")) +
  theme_classic() #+ theme(legend.position = c(50,0.2))


trans.data <- t(dm@data)
colnames(trans.data) <- make.names(colnames(trans.data))
tsne1 <- Rtsne(as_tibble(trans.data) %>% select(names(frequencies1)[1:60]), theta = 0.01)
tsne2 <- Rtsne(as_tibble(trans.data) %>% select(names(frequencies2)[1:40]), theta = 0.01)
tsne3 <- Rtsne(as_tibble(trans.data) %>% select(names(frequencies3)[1:20]), theta = 0.01)

tsneClus <- function(tsne, eps, sub){
  assignment <- as.factor(dbscan(tsne$Y,eps,minPts = 10)$cluster)
  
  tsneplot <- ggplot(data.frame(x=tsne$Y[,1], y=tsne$Y[,2], Cluster=assignment), aes(x=x, y=y, color=Cluster))  + 
    geom_point() + labs(x = element_blank(), y= element_blank()) +
    theme_classic()
  
  print(tsneplot)
  
  locations <- read_csv(paste0("valid_topk",sub,".txt"),col_names = F, col_types = cols(), skip = (4-sub)*2) %>% pull(2)
  
  embryo <- sort(unique(assignment)) %>% map(function(c){
    geom <- cbind(dm@geometry[locations,c(1,3)],color=(assignment==c))
    ggplot(as_tibble(geom), aes(x=x,y=z,color=as.factor(color))) + 
      geom_point() + 
      labs(title = paste0("Cluster ",c)) + 
      scale_color_brewer(palette="Set1", direction=-1) + 
      theme_classic() + theme(legend.position = "none")
  })
  
  grid.arrange(grobs=embryo)
}

tsneClus(tsne1,5,1)
tsneClus(tsne2,6,2)
tsneClus(tsne3,5,3)


