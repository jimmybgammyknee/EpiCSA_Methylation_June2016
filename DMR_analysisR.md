## Differentially methylated regions



```
library(methylKit)
library(magrittr)

# reading multiple files
bigList <- list(CpG = list("colWT_CpG.bedGraph", "met1_CpG.bedGraph"),
              CHG = list("colWT_CHG.bedGraph", "met1_CHG.bedGraph"),
              CHH = list("colWT_CHH.bedGraph", "met1_CHH.bedGraph"))

 methTypes <- names(bigList)                 

# read in annotation now -> for later
gene.obj <- read.transcript.features("VV_annotation.bed")

# get new objects of each (characterise each column using the pipeline option)
Objects <- methTypes %>%
  lapply(function(x){
    read(bigList[[x]],
         sample.id=list("Control", "Test"),
         assembly = "Vv145",
         header = FALSE,
         context = x,
         resolution = "base",
         pipeline = list(fraction = FALSE, chr.col = 1, start.col = 2, end.col = 3,
                         coverage.col = 5, strand.col = 4, freqC.col = 6),
                         treatment = c(0,1))
  })
names(bhObjects) <- methTypes

# Make plots/pdfs
methTypes %>%
  lapply(function(x){
    pdf(paste("Rplot_old", x, "bh.pdf", sep = "_"))
    getMethylationStats(bhObjects[[x]][[1]], plot = TRUE, both.strands = FALSE)
    dev.off()

    pdf(paste("Rplot_young", x, "bh.pdf", sep = "_"))
    getMethylationStats(bhObjects[[x]][[2]], plot = TRUE, both.strands = FALSE)
    dev.off()
  })


methObjects <-  bhObjects %>%
  lapply(unite)
names(methObjects) <- methTypes

diffObjects <- methObjects %>%
  lapply(calculateDiffMeth)
names(diffObjects) <- methTypes

# Write all
methTypes %>%
  lapply(function(x){
    write.table(diffObjects[[x]], paste("BH", x, "diffMethylation_all.tsv", sep = "_"), sep="\t")
  })

                   # cluster all samples using correlation distance and plot hiarachical clustering
                   # JB: not needed with two samples -> could change if I use all 3 replicates instead of pooled

                   #clusterSamples(meth.cpg, dist="correlation", method="ward", plot=FALSE)
                   #clusterSamples(meth.chg, dist="correlation", method="ward", plot=FALSE)
                   #clusterSamples(meth.chh, dist="correlation", method="ward", plot=FALSE)

diff25p <- methTypes %>%
  lapply(function(x){
    get.methylDiff(diffObjects[[x]], difference = 25, qvalue = 0.01)
  })
names(diff25p) <- methTypes

# Export all the tables
methTypes %>%
  lapply(function(x){

    # Export Significant DMRs
    write.table(diff25p[[x]],
                paste("BH", x, "diffMethylation_25p_q01.tsv", sep="_"),
                sep = "\t")

    # Significant DMRs on each chromosome
    pdf(paste("Rplot", x, "bh_chrom.pdf", sep= "_"))
    diffMethPerChr(diff25p[[x]], plot = TRUE, qvalue.cutoff=0.01, meth.cutoff=25)
    dev.off()

    # Plot gene annotation summary
    gene.ann <- annotate.WithGenicParts(diff25p[[x]], gene.obj)
    pdf(paste("Rplot", x, "bh_geneSummary.pdf", sep="_"))
    plotTargetAnnotation(gene.ann)
    dev.off()

    tss.ann <- getAssociationWithTSS(gene.ann[[x]])
    pdf(paste("Rplot", x, "bh_TSS.pdf", sep="_"))
    hist(tss.ann$dist.to.feature[abs(tss.ann$dist.to.feature)<=100000],
         main="distance to nearest TSS",
         xlab="distance in bp",
         breaks=50,
         col="brown4")
    dev.off()

  })
```
