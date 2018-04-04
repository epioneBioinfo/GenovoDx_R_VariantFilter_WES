##########################################################################################
#
#  Import GenovoDx variant calling data
#  Receive data in Excel format 
#      Take Excel sheet = "All Variants"
#      Rename header names to a single row and easy to header names
#      Import data in to R
#  Import other supporting files
#      Selected Consequnce File (Moderate and High groups (SO_Impact.txt)
#      Actionable 591 gene panel list (gene_comb_list.txt)
#
########################################################################################
setwd("C:/Data/2018-03 CTC Gastric Ca/2018-03-26_GenovoDx/20180326_variants_analysis")

### Import variant data
raw = read.delim("data_variants_03262018.txt", header=T, sep="\t", stringsAsFactors = F)

### Import Sample information
SampleInfo = read.delim("SampleInfo.txt", header=T, sep="\t", stringsAsFactors = F)

### Sort variant data by Chromosome and Position
chr = as.character(raw$CHROM)
chr = gsub("chr", "", chr)
chr[chr=="X"] = "24"
chr[chr=="Y"] = "25"
chr[chr=="M"] = "26"
chr = as.integer(chr)
oo = order(chr, raw$POS)
raw = raw[oo, ]

remove(list=c("chr", "oo"))

### Import 591 actionable onco-gene list
gcomb = read.delim("gene_comb_list.txt", header=F, sep="")
gcomb = as.character(gcomb[,1])

### Import Consequence information
so.file = read.delim("SO_Impact.txt", header=T, sep="\t", stringsAsFactors=F)
so.list = as.character(so.file$SO.term)


save.image("C:/Data/2018-03 CTC Gastric Ca/2018-03-26_GenovoDx/20180326_variants_analysis/20180326_variants.RData")

### Function to split imported variant data for Tumor/Control pair if available, and combined for identity columns
split.data.by.sample = function(dat, sample.names, indexLeadCols=1:10) {
  info = dat[, indexLeadCols]
  d = lapply(sample.names, function(x) {
               j = grep(x, colnames(dat))
               dd = cbind(info, dat[, j])
               dd
            })
  d
}


### Get unique sample names
sample.names=as.character(unique(sapply(SampleInfo$SampleID, function(x)strsplit(x, "-")[[1]][1])))

### perform sample pair data splitting
dat = split.data.by.sample(raw, sample.names = sample.names, indexLeadCols=1:10)
ii = sapply(sample.names, function(x)grep(x, SampleInfo$SampleID)[1])
names(dat) = paste(sample.names, SampleInfo$Alias[ii], sep="_")
                   

### Function that filters Variants
filter.variants = function(inputDat, i.case, i.ctrl, AF.thres=0.2, DP.thres=30, useConsequence=T) {
           
  outputDat = inputDat
      
  d.case = inputDat[, i.case]
  d.ctrl = inputDat[, i.ctrl]
  af1 = as.numeric(d.case[, 1])
  af0 = as.numeric(d.ctrl[, 1])
  dp1 = as.numeric(d.case[, 2])
  dp0 = as.numeric(d.ctrl[, 2])
      
  ii1 = (af1 >= AF.thres & !is.na(af1)) & dp1 >= DP.thres
  ii0 = (af0 < AF.thres | is.na(af0)) & dp0 >= DP.thres
  f = ii1 & ii0
  
  jj = grep("TumorPresence", colnames(outputDat))
  outputDat[, jj] = ifelse(f, "Yes", "No")
  
  if (useConsequence) 
    f = f & outputDat$Gene!="" & outputDat$Consequence%in%so.list
  else
    f = f & outputDat$Gene!=""

  outputDat = subset(outputDat, subset=f)    
  outputDat      
   
}

#### Perform variant filtering
### Filter FFPE samples for mutations that meet the filtering criteria -> Called Variants
df.fe = lapply(dat[1:3], function(x) {
             y = filter.variants(x, i.case=c(11, 12), i.ctrl=c(13, 14), AF.thres=0.05, DP.thres=30, useConsequence=T)
          })

### Filter CTC samples for mutations that meet the filtering criteria -> Called Variants
df.ctc = lapply(dat[4:5], function(x) {
                y = filter.variants(x, i.case=c(11, 12), i.ctrl=c(13, 14), AF.thres=0.2, DP.thres=30, useConsequence=T)
               })


### 
### Match called variants to actionable mutations
### FFPE
df.fe.act = lapply(df.fe, function(d_i) {
                   ii = unique(unlist(lapply(gcomb, function(x)grep(x, d_i$Gene))))
                   ii = ii[order(ii)]
                   dgcomb = d_i[ii, ]
                 })

### CTC
df.ctc.act = lapply(df.ctc, function(d_i) {
                    ii = unique(unlist(lapply(gcomb, function(x)grep(x, d_i$Gene))))
                    ii = ii[order(ii)]
                    dgcomb = d_i[ii, ]
                 })


### Load Stomach Ca Good and Poor Prognosis targets
gene.good = read.delim("stomach_good_prognosis.txt", header=T, sep="\t", stringsAsFactors = F)
gene.poor = read.delim("stomach_poor_prognosis.txt", header=T, sep="\t", stringsAsFactors = F)
g.good = gene.good$Gene
g.poor = gene.poor$Gene

### Function to build the unique key for each row to match mutations.
gene.key = function(inputDat) { paste(inputDat$CHROM, inputDat$Gene, inputDat$POS, inputDat$REF, inputDat$ALT, sep=":") }

### Match to good prognosis gene list
### FFPE
df.fe.good = lapply(df.fe, function(di) {
  ii = lapply(paste(":", g.good, ":", sep=""), function(x)grep(x, gene.key(di)))
  ii = unique(unlist(ii))
  ii= ii[order(ii)]  
  di[ii, ]
})

### Match to poor prognosis gene list
### FFPE
df.fe.poor = lapply(df.fe, function(di) {
  ii = lapply(paste(":", g.poor, ":", sep=""), function(x)grep(x, gene.key(di)))
  ii = unique(unlist(ii))
  ii= ii[order(ii)]  
  di[ii, ]
})

### Match to good prognosis gene list
### CTC
df.ctc.good = lapply(df.ctc, function(di) {
  ii = lapply(paste(":", g.good, ":", sep=""), function(x)grep(x, gene.key(di)))
  ii = unique(unlist(ii))
  ii= ii[order(ii)]  
  di[ii, ]
})

### Match to poor prognosis gene list
### CTC
df.ctc.poor = lapply(df.ctc, function(di) {
  ii = lapply(paste(":", g.poor, ":", sep=""), function(x)grep(x, gene.key(di)))
  ii = unique(unlist(ii))
  ii= ii[order(ii)]  
  di[ii, ]
})



### Export Results into tab-delimited files
#### Called Variants
### FFPE
for (i in 1:length(df.fe)) {
   write.table(df.fe[[i]], paste(names(df.fe)[i], "ffpe_called_variants.txt", sep="_"), row.names=F, col.names=T, se="\t") 
}

### CTC
for (i in 1:length(df.ctc)) {
  write.table(df.ctc[[i]], paste(names(df.ctc)[i], "ctc_called_variants.txt", sep="_"), row.names=F, col.names=T, se="\t") 
}


#### Actionable Variants - FFPE
for (i in 1:length(df.fe.act)) {
  write.table(df.fe.act[[i]], paste(names(df.fe.act)[i], "ffpe_action_variants.txt", sep="_"), row.names=F, col.names=T, se="\t") 
}

#### Actionable Variants - CTC
for (i in 1:length(df.ctc.act)) {
  write.table(df.ctc.act[[i]], paste(names(df.ctc.act)[i], "ctc_action_variants.txt", sep="_"), row.names=F, col.names=T, se="\t") 
}


#### Good Prognosis Variants - FFPE
for (i in 1:length(df.fe.good)) {
  write.table(df.fe.good[[i]], paste(names(df.fe.good)[i], "ffpe_good_variants.txt", sep="_"), row.names=F, col.names=T, se="\t") 
}

#### Good Prognosis Variants - CTC
for (i in 1:length(df.ctc.good)) {
  write.table(df.ctc.good[[i]], paste(names(df.ctc.good)[i], "ctc_good_variants.txt", sep="_"), row.names=F, col.names=T, se="\t") 
}


#### Poor Prognosis Variants - FFPE
for (i in 1:length(df.fe.poor)) {
  write.table(df.fe.poor[[i]], paste(names(df.fe.poor)[i], "ffpe_poor_variants.txt", sep="_"), row.names=F, col.names=T, se="\t") 
}

#### Poor Prognosis Variants - CTC
for (i in 1:length(df.ctc.poor)) {
  write.table(df.ctc.poor[[i]], paste(names(df.ctc.poor)[i], "ctc_poor_variants.txt", sep="_"), row.names=F, col.names=T, se="\t") 
}



###################################################################################
#
#        TMB calculation and MMR extraction
#
###################################################################################
### Calculate TMB
### Load exon coordinates data - prestored in RData (TruSeq coordinates, TruSight-170 coordinates,
### and UCSC RefSeq exon coordinates.)
load("C:/Data/Misc/Gastric Ca/20180326_GDX/exon_coords.RData")

### Define TMB calculation with TruSeq coordinates
truseq.calc = function(di, ex) {
  chroms = unique(di$CHROM)
  cs = 0
  for (ci in chroms) {
    dci = subset(di, subset=CHROM==ci)
    exi = ex[[ci]]
    idx = apply(dci,1,function(r){
                ii = which(exi$exonStart<=r[[2]] & exi$exonEnd>=r[[2]])
                })
    cs = cs + sum(sapply(idx, function(xi)length(xi)>0))
  }
  cs
}

### Split TruSeq coordinates by Chromodome
ex = split(exCoord.truseq, exCoord.truseq$chrom)

### FFPE - with TruSeq coordinates
TMB.fe.truseq = lapply(df.fe, function(x) {
  countMuts = truseq.calc(x, ex)
  L = sum(exCoord.truseq$Length, na.rm=T)
  r = countMuts * 1e6 / L
  r
})

### FFPE - with UCSC RefSeq coordinates
TMB.fe.ucsc = lapply(df.fe, function(x) {
  ii = lapply(x$Gene, function(x)grep(x, exCoord.ucsc$Gene))
  countMuts = sum(sapply(ii, function(x)length(x)>0), na.rm=T)
  ii = unique(unlist(ii))
  L = sum(exCoord.ucsc$Length, na.rm=T)
  r = countMuts * 10^6 / L
  r
})


### CTC - with TruSeq coordinates
TMB.ctc.truseq = lapply(df.ctc, function(x) {
  countMuts = truseq.calc(x, ex)
  L = sum(exCoord.truseq$Length, na.rm=T)
  r = countMuts * 1e6 / L
  r
})

### CTC - with UCSC coordinates
TMB.ctc.ucsc = lapply(df.ctc, function(x) {
  ii = lapply(x$Gene, function(x)grep(x, exCoord.ucsc$Gene))
  countMuts = sum(sapply(ii, function(x)length(x)>0), na.rm=T)
  ii = unique(unlist(ii))
  L = sum(exCoord.ucsc$Length, na.rm=T)
  r = countMuts * 10^6 / L
  r
})

TMB.fe.truseq1 = do.call("rbind", TMB.fe.truseq)
TMB.fe.ucsc1 = do.call("rbind", TMB.fe.ucsc)
TMB.ctc.truseq1 = do.call("rbind", TMB.ctc.truseq)
TMB.ctc.ucsc1 = do.call("rbind", TMB.ctc.ucsc)

plot(TMB.fe.truseq1[,1], TMB.fe.ucsc1[,1])
plot(c(TMB.fe.truseq1[,1], TMB.ctc.truseq1[2,1]), c(TMB.fe.ucsc1[,1], TMB.ctc.ucsc1[2,1]))


### Calculate MMR
mmr = c("MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2")

### FFPE MMR
df.fe.mmr = lapply(df.fe, function(di) {
  ii = lapply(paste(":", mmr, ":", sep=""), function(x)grep(x, gene.key(di)))
  ii = unique(unlist(ii))
  ii= ii[order(ii)]  
  di[ii, ]
})

### CTC MMR
df.ctc.mmr = lapply(df.ctc, function(di) {
  ii = lapply(paste(":", mmr, ":", sep=""), function(x)grep(x, gene.key(di)))
  ii = unique(unlist(ii))
  ii= ii[order(ii)]  
  di[ii, ]
})


#### FFPE-MMR Variants
for (i in 1:length(df.fe.mmr)) {
  write.table(df.fe.mmr[[i]], paste(names(df.fe.mmr)[i], "ffpe_mmr_hits.txt", sep="_"), row.names=F, col.names=T, se="\t") 
}

#### CTC-MMR Variants
for (i in 1:length(df.ctc.mmr)) {
  write.table(df.ctc.mmr[[i]], paste(names(df.ctc.mmr)[i], "ctc_mmr_hits.txt", sep="_"), row.names=F, col.names=T, se="\t") 
}





