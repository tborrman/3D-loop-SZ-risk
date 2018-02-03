
# Read data tables
NPC_intxns <- read.table("S12.NPC-specific_PGC_intrxns.txt", sep="\t", header=TRUE, comment.char="")
chromSizes <- read.table("chrom.sizes", sep="\t", header=FALSE)
# Remove M chrom
chromSizes <- chromSizes[1:nrow(chromSizes) -1, ]
colnames(chromSizes) <- c("chrom", "bp")
# generate cumulative genome table
chrom <- as.character(chromSizes$chrom)
bp <- cumsum(as.numeric(chromSizes$bp))
cmChromSizes <- data.frame(chrom, bp)



is_cis <- function(row) {
  # Determine if significant interaction is in
  # cis or trans
  # Returns: boolean; True if in cis"
  chr1 <- sub(":.*", "", row$anchor.bin.coord)
  chr2 <- sub(":.*", "", row$target.bin.coord)
  if(chr1 == chr2){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

get_distance <- function(row) {
  # Returns: bp distance between significant interaction
  start <- as.numeric(sub(".*:(.*)-.*", "\\1", row$anchor.bin.coord))
  end <- as.numeric(sub(".*:(.*)-.*", "\\1", row$target.bin.coord))
  d <- abs(start - end)
  return(d)
}

rand_genomic_pos <- function(d, c, cm) {
  # Args: '
  # d = bp distance
  # c = chromosome size table
  # cm = cumulative chromosome size table 
  # Returns: list coords
  #   chrom = chromosome
  #   g1 = random genomic coordinate
  #   g2 = paired genomic coordinate d bp from g1
  total <- cm$bp[length(cm$bp)]
  chrStarts <- c(1,cm$bp[1:length(cm$bp) -1] + 1)
  rand_pos <- sample(total, 1)
  chr_idx <- max(which(chrStarts <= rand_pos))
  chrom <- as.character(c$chrom[chr_idx])
  g1 <- rand_pos - chrStarts[chr_idx] + 1
  g2 <- g1 + d
  coords <- list(chrom, g1, g2)
  if (g2 > c$bp[chr_idx]) {
    print(coords)
    coords<- rand_genomic_pos(d, c, cm)
  }
  return(coords)
}

# Get random PGC interacion
r <- sample(1:nrow(NPC_intxns), 1)
randRow = NPC_intxns[r,]
if (is_cis(randRow)) {
  d <- get_distance(randRow)
  x <- rand_genomic_pos(d, chromSizes, cmChromSizes)

} else {
  print("ERROR: trans interaction")
  stop()
}




