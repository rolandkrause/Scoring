library(FNN)
library(igraph)
library(gridExtra)
library(VennDiagram)
library(Rtsne)
library(dbscan)
library(vioplot)

jaccard.genes <- function(pattern, sub) {
  all.submissions <- seq(10) %>%
    map(~ read.csv(str_glue(pattern), header = FALSE, stringsAsFactors = FALSE))

  # for each fold get genes
  selected.genes <- all.submissions %>% map_dfc(function(submission) {
    # separate the gene names from the location predictions
    gene.lines <- (4 - sub) * 2
    genes <- submission %>%
      slice(1:gene.lines) %>%
      select(-1) %>%
      unlist() %>%
      make.names()

    data.frame(f = genes %>% as.character())
  })

  combn(seq(ncol(selected.genes)), 2, function(x) {
    length(intersect(selected.genes[, x[1]], selected.genes[, x[2]])) /
      length(union(selected.genes[, x[1]], selected.genes[, x[2]]))
  })
}

expected.jaccard <- function(n, m) {
  seq(m) %>%
    map_dbl(function(k) {
      choose(m, k) * choose(n - m, m - k) / choose(n, m) * k / (2 * m - k)
    }) %>%
    sum()
}

raw.selected.genes <- function(pattern, sub) {
  all.submissions <- seq(10) %>%
    map(~ read.csv(str_glue(pattern), header = FALSE, stringsAsFactors = FALSE))

  # for each fold get genes
  selected.genes <- all.submissions %>% map_dfc(function(submission) {
    # separate the gene names from the location predictions
    gene.lines <- (4 - sub) * 2
    genes <- submission %>%
      slice(1:gene.lines) %>%
      select(-1) %>%
      unlist() %>%
      make.names()

    data.frame(f = genes %>% as.character())
  })

  return(selected.genes)
}

plot_insitu <- function(gene.name, color = "cyan") {
  pattern <- tibble(
    x = dm@geometry[seq(1, 6078, by = 2), 1], z = dm@geometry[seq(1, 6078, by = 2), 3],
    gene = as.factor(rep(dm@insitu.matrix[, gene.name], 2))[seq(1, 6078, by = 2)]
  )
  ggplot(pattern, aes(x = x, y = z, fill = gene, color = gene)) + geom_point(shape = 21) +
    scale_fill_manual(values = c("white", color)) + scale_color_manual(values = c("gray80", "gray30")) + theme_void() +
    theme(legend.position = "none")
}


tsneClus <- function(tsne, eps, sub) {
  assignment <- dbscan(tsne$Y, eps, minPts = 10)$cluster

  set3.ext <- c(RColorBrewer::brewer.pal(length(unique(assignment)), "Set3"), "#E5E5E5")

  tsneplot <- ggplot(data.frame(x = tsne$Y[, 1], y = tsne$Y[, 2], Cluster = as.factor(assignment)), aes(x = x, y = y, color = Cluster)) +
    geom_point() + labs(x = element_blank(), y = element_blank()) + scale_color_manual(values = set3.ext) +
    theme_void()

  print(tsneplot)

  locations <- read_csv(paste0("valid_topk", sub, ".txt"), col_names = F, col_types = cols(), skip = (4 - sub) * 2) %>% pull(2)
  remaining <- which(!(seq(3039) %in% locations))

  geom <- tibble(x = dm@geometry[locations, 1], z = dm@geometry[locations, 3], Cluster = assignment)
  geom <- rbind(geom, tibble(x = dm@geometry[remaining, 1], z = dm@geometry[remaining, 3], Cluster = 13))



  all.clusters <- ggplot(geom, aes(x = x, y = z, color = as.factor(Cluster))) + 
    geom_point() + scale_color_manual(values = set3.ext) +
    theme_void() + theme(legend.position = "none")
  print(all.clusters)



  embryo <- sort(unique(assignment)) %>% map(function(c) {
    geom <- cbind(dm@geometry[locations, c(1, 3)], color = (assignment == c))
    ggplot(as_tibble(geom), aes(x = x, y = z, color = as.factor(color))) +
      geom_point() +
      labs(title = paste0("Cluster ", c)) +
      scale_color_manual(values = c("gray90", "red")) +
      theme_void() + theme(legend.position = "none")
  })

  grid.arrange(grobs = embryo)
}

# Run DistMap with custom gene list
reduced.DistMap <- function(gene.names, output.file) {
  if (!all(file.exists(c("dge_raw.txt.gz", "dge_normalized.txt.gz", "binarized_bdtnp.csv.gz", "bdtnp.txt.gz", "geometry.txt.gz")))) {
    download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/dge_raw.txt.gz", destfile = "dge_raw.txt.gz")
    download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/dge_normalized.txt.gz", destfile = "dge_normalized.txt.gz")
    download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/binarized_bdtnp.csv.gz", destfile = "binarized_bdtnp.csv.gz")
    download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/bdtnp.txt.gz", destfile = "bdtnp.txt.gz")
    download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/geometry.txt.gz", destfile = "geometry.txt.gz")
  }

  raw.data <- read.table(gzfile("dge_raw.txt.gz", "rt"),
    sep = "\t",
    row.names = NULL,
    stringsAsFactors = F,
    quote = ""
  )
  raw.data.genes <- raw.data$V1
  raw.data$V1 <- NULL

  raw.data.genes <- gsub("'", "", raw.data.genes, fixed = T)
  raw.data.genes <- make.names(raw.data.genes)

  raw.data <- as.matrix(raw.data)
  rownames(raw.data) <- raw.data.genes

  normalized.data <- read.table(gzfile("dge_normalized.txt.gz", "rt"),
    sep = "\t",
    row.names = NULL,
    stringsAsFactors = F,
    quote = ""
  )

  normalized.data.genes <- normalized.data$row.names
  normalized.data$row.names <- NULL

  normalized.data.genes <- gsub("'", "", normalized.data.genes, fixed = T)

  normalized.data <- as.matrix(normalized.data)
  normalized.data.genes <- make.names(normalized.data.genes)

  rownames(normalized.data) <- normalized.data.genes

  stopifnot(all(normalized.data.genes == raw.data.genes))

  insitu.matrix <- read.table(gzfile("binarized_bdtnp.csv.gz", "rt"), sep = ",", header = T)

  insitu.genes_orig <- colnames(insitu.matrix)

  # insitu.genes <- gsub(".","-",insitu.genes_orig,fixed = T)
  # insitu.genes <- gsub("-spl-","(spl)",insitu.genes,fixed = T)

  insitu.genes <- make.names(insitu.genes_orig)


  stopifnot(all(insitu.genes %in% raw.data.genes))

  colnames(insitu.matrix) <- insitu.genes

  # Subset here!
  insitu.matrix <- insitu.matrix %>% select(gene.names)
  insitu.matrix <- as.matrix(insitu.matrix)

  geometry <- read.csv(gzfile("geometry.txt.gz", "rt"), sep = " ")

  closeAllConnections()

  colnames(geometry) <- c("x", "y", "z")

  dm <- new("DistMap",
    raw.data = raw.data,
    data = normalized.data,
    insitu.matrix = insitu.matrix,
    geometry = as.matrix(geometry)
  )

  dm <- binarizeSingleCellData(dm, seq(0.15, 0.75, 0.01))

  dm <- mapCells(dm)

  predictions <- t(apply(dm@mcc.scores, 2, order, decreasing = TRUE))[, 1:10]


  in.format <- c(gene.names, rep(NA, ifelse(length(gene.names) %% 20 == 0, 0, 20 - (length(gene.names) %% 20))))

  dim(in.format) <- c(length(in.format) / 10, 10)
  in.format <- cbind(rep(NA, nrow(in.format)), in.format)
  
  write.table(in.format, file = output.file, na = "", sep = ",", row.names = FALSE, col.names = FALSE)
  in.format <- cbind(c(seq(nrow(predictions))), predictions)
  write.table(in.format, file = output.file, na = "", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
}

plot_insitu <- function(gene.name, color="cyan"){
  pattern <- tibble(x=dm@geometry[seq(1,6078,by=2),1], z=dm@geometry[seq(1,6078,by=2),3], 
                    gene = as.factor(rep(dm@insitu.matrix[,gene.name], 2))[seq(1,6078,by=2)])
  ggplot(pattern, aes(x=x,y=z,fill=gene,color=gene)) + geom_point(shape=21) + 
    scale_fill_manual(values=c("white",color)) + scale_color_manual(values=c("gray80","gray30")) + theme_void() + 
    theme(legend.position = "none")
}
