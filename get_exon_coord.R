###################################################################################
###
### Get UCSC downloaded RefSeq (curated) exon coordinates
### Go to UCSC table browser   https://genome.ucsc.edu/cgi-bin/hgTables
###  select: Group = Genes and Gene Predictions, 
###          Assembly = GRCh37/hg19
###          Track = RefSeq Genes
###          Table = refSeq Genes
###          Upload 591 onco actionable gene list
###  Download retrieved data and save as txt file
###
#################################################################################### 
### Load UCSC cooddintes data.
d = read.delim("ucsc_exon_coord_032818.txt", header=T, sep="\t", stringsAsFactors=F)

### Split Exon coordinates and calculate exon length
d = t(apply(d, 1, function(x) {
  x3 = x[3]
  x4 = x[4]
  x3 = as.numeric(strsplit(x3, ",")[[1]])
  x4 = as.numeric(strsplit(x4, ",")[[1]])
  L = sum(x4-x3+1, na.rm=T)
  x = c(x, L)
  names(x) =  c("chrom", "exonCount", "exonStarts", "exonEnds", "Gene", "Length")
  x
}))

### Format data
d = data.frame(d, stringsAsFactors=F)
d$exonCount = as.numeric(d$exonCount)
d$Length = as.numeric(d$Length)

### Remove any zro length rows
d = subset(d, subset=Length>0)


################################################################################
###
### Illumina TruSeq Coordinates
### https://support.illumina.com/downloads/truseq-exome-product-files.html
###   unzip bed file 
###
################################################################################

exCoord.truseq = read.delim("truseq-exome-targeted-regions-manifest-v1-2.bed", header=F, sep="\t", stringsAsFactors=F)
names(exCoord.truseq) = c("chrom", "exonStart", "exonEnd", "Info")

### Parse trueseq data and calculate exon length
exCoord.truseq = t(apply(exCoord.truseq, 1, function(x) {
  x4 = x[4]
  x4 = strsplit(x4, "-")[[1]][1]
  x = c(x, Gene=x4)
  x
}))

exCoord.truseq = data.frame(exCoord.truseq, stringsAsFactors=F)
exCoord.truseq$exonStart = as.numeric(exCoord.truseq$exonStart)
exCoord.truseq$exonEnd = as.numeric(exCoord.truseq$exonEnd)
exCoord.truseq$Length = exCoord.truseq$exonEnd - exCoord.truseq$exonStart + 1



################################################################################
###
### Illumina TrieSight 170 tumor set BED Coordinates
### https://support.illumina.com/downloads/trusight_cancer_product_files.html
###   unzip bed file 
###
################################################################################
exCoord = read.delim("trusight_cancer_manifest_a.bed", header=F, sep="\t", stringsAsFactors=F)
names(exCoord) = c("chrom", "exonStart", "exonEnd", "Info")

exCoord = t(apply(exCoord, 1, function(x) {
  x4 = x[4]
  x4 = strsplit(x4, "\\.")[[1]][1]
  x = c(x, Gene=x4)
  x
}))

exCoord = data.frame(exCoord, stringsAsFactors=F)
exCoord$exonStart = as.numeric(exCoord$exonStart)
exCoord$exonEnd = as.numeric(exCoord$exonEnd)
exCoord$Length = exCoord$exonEnd - exCoord$exonStart + 1


save.image("C:\\Data\\Misc\\Gastric Ca\\20180326_GDX\\exome_coords.RData")
