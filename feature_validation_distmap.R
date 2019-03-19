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
lines <- read_lines("post analysis/intersect_features.txt", skip = 1)

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


