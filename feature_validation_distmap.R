source("dream_scoring_clean.R")


reduced.DistMap <- function(gene.names, output.file){

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
  raw.data.genes <- make.names(raw.data.genes)
  
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
  normalized.data.genes <- make.names(normalized.data.genes)
  
  rownames(normalized.data) = normalized.data.genes
  
  stopifnot(all(normalized.data.genes == raw.data.genes))
  
  insitu.matrix = read.table(gzfile("binarized_bdtnp.csv.gz",'rt'), sep = ",",header = T)
  
  insitu.genes_orig <- colnames(insitu.matrix)
  
  #insitu.genes <- gsub(".","-",insitu.genes_orig,fixed = T)
  #insitu.genes <- gsub("-spl-","(spl)",insitu.genes,fixed = T)
  
  insitu.genes <- make.names(insitu.genes_orig)
  
  
  stopifnot(all(insitu.genes %in% raw.data.genes))
  
  colnames(insitu.matrix) <- insitu.genes
  
  #Subset here!
  insitu.matrix <- insitu.matrix %>% select(gene.names)
  insitu.matrix <- as.matrix(insitu.matrix)
  
  geometry = read.csv(gzfile("geometry.txt.gz",'rt'),sep = " ")
  
  closeAllConnections()
  
  colnames(geometry) = c("x","y","z")
  
  dm <- new("DistMap",
             raw.data=raw.data,
             data=normalized.data,
             insitu.matrix=insitu.matrix,
             geometry=as.matrix(geometry))
  
  dm <- binarizeSingleCellData(dm, seq(0.15, 0.75, 0.01))
  
  dm <- mapCells(dm)
  
  predictions <- t(apply(dm@mcc.scores,2,order,decreasing=TRUE))[,1:10]
  
  
  in.format <- c(gene.names, rep(NA, ifelse(length(gene.names) %% 20 == 0, 0, 20 - (length(gene.names) %% 20))))
  
  dim(in.format) <- c(length(in.format)/10,10)
  in.format <- cbind(rep(NA,nrow(in.format)),in.format)
  #sub <- (4 - dim(in.format)[1]/2)
  write.table(in.format, file=output.file, na = "",sep=",", row.names=FALSE, col.names = FALSE)                   
  in.format <- cbind(c(seq(nrow(predictions))),predictions)
  write.table(in.format,file=output.file, na ="",sep=",",row.names=FALSE, col.names = FALSE, append = TRUE) 
  
}

#path to file with intersection of genes across subchallenges
lines <- read_lines("post analysis/intersect_features_cv.txt", skip = 1)

lines %>% walk(function(line){

  splitted <- str_split(line, "\\s")[[1]]
  out.filename <- paste0("validation/", paste(splitted[1:2],collapse="."), ".txt")
  reduced.DistMap(splitted[-(1:2)],out.filename)

  if(as.numeric(splitted[2]) <= 60) sub <- 1
  if(as.numeric(splitted[2]) <= 40) sub <- 2
  if(as.numeric(splitted[2]) <= 20) sub <- 3

  scores <- score(out.filename,sub)

  write(c(out.filename,scores), file = "validation_results.txt", ncolumns = 4,  append = T)
})


insitu.matrix = read.table(gzfile("binarized_bdtnp.csv.gz",'rt'), sep = ",",header = T)

colnames(insitu.matrix) <- make.names(colnames(insitu.matrix))


#generate and score random samples
seq(100) %>% walk(function(seed){
  set.seed(seed)
  gene.sample <- sample(colnames(insitu.matrix),60)
  
  
  reduced.DistMap(gene.sample,paste0("temp",seed))
  scores <- score(paste0("temp",seed),1)
  
  write_lines(paste(c(gene.sample,scores),collapse = ","),path = "validation_distmap_sub1.csv",append = T)
  
  reduced.DistMap(gene.sample[1:40],paste0("temp",seed))
  scores <- score(paste0("temp",seed),2)
  
  write_lines(paste(c(gene.sample[1:40],scores),collapse = ","),path = "validation_distmap_sub2.csv",append = T)
  
  reduced.DistMap(gene.sample[1:20], paste0("temp",seed))
  scores <- score(paste0("temp",seed),3)
  
  write_lines(paste(c(gene.sample[1:20],scores),collapse = ","),path = "validation_distmap_sub3.csv",append = T)
  
  file.remove(paste0("temp",seed))
})


#collect and plot
scores <- seq(3) %>% map_dfr(~read_csv(paste0("post analysis/validation_distmap_sub",.x,".csv"), col_names = FALSE, col_types = cols()) %>%
                  select_if(is.numeric) %>% mutate(subc=.x) %>% setNames(c("s1","s2","s3","subc")) %>% gather("score","value", -subc))

top.sub <- read_lines("post analysis/most.freq.txt") %>% iwalk(~reduced.DistMap(str_split(.x," ",simplify = T),output.file=paste0("valid_topk",.y,".txt"))) %>% imap(~score(paste0("valid_topk",.y,".txt"),.y))

ggplot(scores, aes(as.factor(subc),value, fill=score)) + geom_violin(draw_quantiles = c(0.25,0.5,0.75), position = position_dodge(0.6)) +
  scale_fill_brewer(palette="Blues") + labs(title="DistMap scores of top-k selected genes vs. null distributions", x="Subchallenge", y="Value", fill="Score") + 
  geom_point(data = tibble(subc = rep(seq(3), each=3) + rep(c(-0.2,0,0.2),3), vals = top.sub %>% unlist), shape=21, size=2, stroke=1, inherit.aes=FALSE, aes(x=subc, y=vals), fill="red") + theme_classic()

#percentile scores
ecdf(scores %>% filter(subc==1, score=="s1") %>% pull(value))(top.sub[[1]][1])
ecdf(scores %>% filter(subc==1, score=="s2") %>% pull(value))(top.sub[[1]][2])
ecdf(scores %>% filter(subc==1, score=="s3") %>% pull(value))(top.sub[[1]][3])

ecdf(scores %>% filter(subc==2, score=="s1") %>% pull(value))(top.sub[[2]][1])
ecdf(scores %>% filter(subc==2, score=="s2") %>% pull(value))(top.sub[[2]][2])
ecdf(scores %>% filter(subc==2, score=="s3") %>% pull(value))(top.sub[[2]][3])

ecdf(scores %>% filter(subc==3, score=="s1") %>% pull(value))(top.sub[[3]][1])
ecdf(scores %>% filter(subc==3, score=="s2") %>% pull(value))(top.sub[[3]][2])
ecdf(scores %>% filter(subc==3, score=="s3") %>% pull(value))(top.sub[[3]][3])


ss <- read_csv("spatial.stats.csv")
tidyss <- ss %>% select(-4) %>% gather("measure","value",-gene)

cv_genes <- str_split(lines[[37]],"\\s")[[1]][-(1:2)]
cv_genes_ss <- tidyss %>% filter(gene %in% cv_genes)

ggplot(tidyss, aes(measure,value, fill = measure)) + geom_violin(trim=F, draw_quantiles = c(0.25,0.5,0.75)) + labs(y="Value", x=element_blank(), fill = "Measure") + 
  geom_jitter(data=cv_genes_ss, aes(x=measure, y=value), colour="black", width=0.1)  + 
  coord_flip() + theme_classic() + theme(axis.text.y = element_blank())


shapiro.test(tidyss %>%  filter(measure=="entropy") %>% pull(value))
shapiro.test(tidyss %>%  filter(measure=="MoranI") %>% pull(value))

#Mann-Whitney U test since the distributions are not normal Shapiro-Wilk rejects normality
wilcox.test(tidyss %>%  filter(measure=="entropy") %>% pull(value), cv_genes_ss %>% filter(measure=="entropy") %>% pull(value), alternative = "less")
wilcox.test(tidyss %>%  filter(measure=="MoranI") %>% pull(value), cv_genes_ss %>% filter(measure=="MoranI") %>% pull(value), alternative = "less")

topk.genes <- read_lines("post analysis/most.freq.txt") %>% map(~str_split(.x, " ",simplify=T))

spate1 <- ggplot(tidyss, aes(measure,value, fill = measure)) + geom_violin(draw_quantiles = c(0.25,0.5,0.75)) + 
   labs(title="Subchallenge 1", y="Value", x=element_blank(), fill = "Measure") + scale_fill_brewer(palette="Greens") +
   geom_jitter(data = tidyss %>% filter(gene %in% topk.genes[[1]]), aes(x=measure, y=value), colour="black", shape=21, fill="red", size=2, stroke=1, width=0.1) +
   theme_classic() + theme(axis.text.x = element_blank())

spate2 <- ggplot(tidyss, aes(measure,value, fill = measure)) + geom_violin(draw_quantiles = c(0.25,0.5,0.75)) + 
  labs(title="Subchallenge 2", y="Value", x=element_blank(), fill = "Measure") + scale_fill_brewer(palette="Greens") +
  geom_jitter(data = tidyss %>% filter(gene %in% topk.genes[[2]]), aes(x=measure, y=value), colour="black", shape=21, fill="red", size=2, stroke=1, width=0.1) +
  theme_classic() + theme(axis.text.x = element_blank())

spate3 <- ggplot(tidyss, aes(measure,value, fill = measure)) + geom_violin(draw_quantiles = c(0.25,0.5,0.75)) + 
  labs(title="Subchallenge 3", y="Value", x=element_blank(), fill = "Measure") + scale_fill_brewer(palette="Greens") +
  geom_jitter(data = tidyss %>% filter(gene %in% topk.genes[[3]]), aes(x=measure, y=value), colour="black", shape=21, fill="red", size=2, stroke=1, width=0.1) +
  theme_classic() + theme(axis.text.x = element_blank())

grid.arrange(spate1,spate2,spate3,ncol=3)

wilcox.test(tidyss %>%  filter(measure=="entropy") %>% pull(value), tidyss %>% filter(gene %in% topk.genes[[1]]) %>% filter(measure=="entropy") %>% pull(value), alternative = "less")
wilcox.test(tidyss %>%  filter(measure=="MoranI") %>% pull(value), tidyss %>% filter(gene %in% topk.genes[[1]]) %>% filter(measure=="MoranI") %>% pull(value), alternative = "less")

wilcox.test(tidyss %>%  filter(measure=="entropy") %>% pull(value), tidyss %>% filter(gene %in% topk.genes[[2]]) %>% filter(measure=="entropy") %>% pull(value), alternative = "less")
wilcox.test(tidyss %>%  filter(measure=="MoranI") %>% pull(value), tidyss %>% filter(gene %in% topk.genes[[2]]) %>% filter(measure=="MoranI") %>% pull(value), alternative = "less")

wilcox.test(tidyss %>%  filter(measure=="entropy") %>% pull(value), tidyss %>% filter(gene %in% topk.genes[[3]]) %>% filter(measure=="entropy") %>% pull(value), alternative = "less")
wilcox.test(tidyss %>%  filter(measure=="MoranI") %>% pull(value), tidyss %>% filter(gene %in% topk.genes[[3]]) %>% filter(measure=="MoranI") %>% pull(value), alternative = "less")
