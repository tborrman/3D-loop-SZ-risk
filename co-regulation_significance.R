NPC_intxns <- read.table("S12.NPC-specific_PGC_intrxns.txt", sep="\t", header=TRUE, comment.char="")


is_cis <- function(row) {
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
  start <- as.numeric(sub(".*:(.*)-.*", "\\1", row$anchor.bin.coord))
  end <- as.numeric(sub(".*:(.*)-.*", "\\1", row$target.bin.coord))
  d <- abs(start - end)
  return(d)
}

# Get random PGC interacion
r <- sample(1:nrow(NPC_intxns), 1)
randRow = NPC_intxns[r,]
if (is_cis(randRow)) {
  d <- get_distance(randRow)
  print(randRow)
  print(d)

} else {
  print("ERROR: trans interaction")
  stop()
}






