library(VennDiagram)

signif <- function(d, cells) {
  for (i in 1:length(cells)) {
    d <- subset(d, d[cells[i]] < -1)
  }
  return(nrow(d))
}

clean_qvals <- function(d) {
  print(nrow(d))
  # Remove NA cases
  d <- na.omit(d)
  print(nrow(d))
  print(head(d))
  # Fix zeros
  d[d$AstrofdrDonut == 0, "AstrofdrDonut"] = 1.0e-42
  d[d$GMfdrDonut == 0, "GMfdrDonut"] = 1.0e-42
  d[d$NeufdrDonut == 0, "NeufdrDonut"] = 1.0e-42
  d[d$NPCfdrDonut == 0, "NPCfdrDonut"] = 1.0e-42
  print(head(d))
  # Change q values to log10(qvals)
  d$AstrofdrDonut <- log10(d$AstrofdrDonut)
  d$GMfdrDonut <- log10(d$GMfdrDonut)
  d$NeufdrDonut <- log10(d$NeufdrDonut)
  d$NPCfdrDonut <- log10(d$NPCfdrDonut)
  print(head(d))
  
  # Remove duplicates
  d <- d[!duplicated(d[,c('chr1','x1','x2','chr2','y1','y2')]),]
  print(nrow(d))
  return(d)
  
  
}



df <- read.table('../master_loops_Schahram/master_requested_loops', header=TRUE, sep="\t")

df <- clean_qvals(df)

cell = c('AstrofdrDonut', 'GMfdrDonut', 'NeufdrDonut', 'NPCfdrDonut')
png('4wayVenn.png', width=3000, height=2700, res=300)
draw.quad.venn(signif(df, cell[1]), signif(df, cell[2]), signif(df, cell[3]), signif(df, cell[4]),
               signif(df, cell[1:2]), signif(df, cell[c(1,3)]), signif(df, cell[c(1,4)]), 
               signif(df, cell[2:3]), signif(df, cell[c(2,4)]), signif(df, cell[3:4]),
               signif(df, cell[1:3]), signif(df, cell[c(1,2,4)]), signif(df, cell[c(1,3,4)]),
               signif(df, cell[2:4]), signif(df, cell),
               category=c('Astro', 'GM', 'Neu', 'NPC'),
               fill = c('seagreen2', 'lightskyblue', 'orangered', 'mediumorchid4')
)
dev.off()

# Get brain specific loops
bs <- df
for (cell in c('AstrofdrDonut', 'NeufdrDonut', 'NPCfdrDonut')) {
  bs <- subset(bs, bs[cell] < -1)
}

bs <- subset(bs, bs['GMfdrDonut'] > -1)

write.table(bs, 'brain_specific_loops.txt', row.names=FALSE, sep='\t', quote=FALSE)

gm <- df
for (cell in c('AstrofdrDonut', 'NeufdrDonut', 'NPCfdrDonut')) {
  gm <- subset(gm, gm[cell] > -1)
}
gm <- subset(gm, gm['GMfdrDonut'] < -1)
write.table(gm, 'GM_specific_loops.txt', row.names=FALSE, sep='\t', quote=FALSE)