library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(RColorBrewer)

getPalette.1 <- colorRampPalette(brewer.pal(9, "Set1"))


# 1. 2 figures: Replicates consistency Venn diagram for Donor1 vs Donor2


DC552_NI18 <- toGRanges("meth_CpG_hmr/hg38_16/SRR7449826.hmr", format="BED", header=FALSE)
DC552_TB18 <- toGRanges("meth_CpG_hmr/hg38_16/SRR7449830.hmr", format="BED", header=FALSE)
DC555_NI18 <- toGRanges("meth_CpG_hmr/hg38_16/SRR7449838.hmr", format="BED", header=FALSE)
DC555_TB18 <- toGRanges("meth_CpG_hmr/hg38_16/SRR7449842.hmr", format="BED", header=FALSE)

hyper_DC552_NI18 <- toGRanges("meth_CpG_hypermr/hg38_16/SRR7449826.hypermr", format="BED", header=FALSE)
hyper_DC552_TB18 <- toGRanges("meth_CpG_hypermr/hg38_16/SRR7449830.hypermr", format="BED", header=FALSE)
hyper_DC555_NI18 <- toGRanges("meth_CpG_hypermr/hg38_16/SRR7449838.hypermr", format="BED", header=FALSE)
hyper_DC555_TB18 <- toGRanges("meth_CpG_hypermr/hg38_16/SRR7449842.hypermr", format="BED", header=FALSE)


ol_hypo_ni <- findOverlapsOfPeaks(DC552_NI18, DC555_NI18)
ol_hyper_ni <- findOverlapsOfPeaks(hyper_DC552_NI18, hyper_DC555_NI18)

ol_hypo_tb <- findOverlapsOfPeaks(DC552_TB18, DC555_TB18)
ol_hyper_tb <- findOverlapsOfPeaks(hyper_DC552_TB18, hyper_DC555_TB18)

hypo <- findOverlapsOfPeaks(DC552_NI18, DC555_NI18, DC552_TB18, DC555_TB18)
hyper <- findOverlapsOfPeaks(hyper_DC552_NI18, hyper_DC555_NI18, hyper_DC552_TB18, hyper_DC555_TB18)

#   1 figure: Hypo methylated regions, NI

makeVennDiagram(ol_hypo_ni, fill=getPalette.1(2), # circle fill color
                col=getPalette.1(2), # circle border color
                cat.col=getPalette.1(2), main='hypomethylated regions of NI18 samples') # label color
makeVennDiagram(ol_hypo_tb, fill=getPalette.1(2), # circle fill color
                col=getPalette.1(2), # circle border color
                cat.col=getPalette.1(2), main='hypomethylated regions of TB18 samples') # label color


#   1 figure: Hyper methylated regions

makeVennDiagram(ol_hyper_ni, fill=getPalette.1(2), # circle fill color
                col=getPalette.1(2), # circle border color
                cat.col=getPalette.1(2), main='hypermethylated regions of NI18 samples') # label color
makeVennDiagram(ol_hyper_tb, fill=getPalette.1(2), # circle fill color
                col=getPalette.1(2), # circle border color
                cat.col=getPalette.1(2), main='hypermethylated regions of TB18 samples') # label color


makeVennDiagram(hypo, fill=getPalette.1(4), # circle fill color
                col=getPalette.1(4), # circle border color
                cat.col=getPalette.1(4), main='hypomethylated regions: all samples') # label color
makeVennDiagram(hyper, fill=getPalette.1(4), # circle fill color
                col=getPalette.1(4), # circle border color
                cat.col=getPalette.1(4), main='hypermethylated regions: all samples') # label color

# 2. 2 figures: Distribution with respect to genes for 4 files ( 2x Donors * 2x Conditions) using

#   1 figure:Hypo methylated regions

hypo_peaks <- GRangesList(DC552_NI18 = DC552_NI18, DC552_TB18 = DC552_TB18,
                          DC555_NI18 = DC555_NI18, DC555_TB18 = DC555_TB18)

genomicElementDistribution(hypo_peaks,
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))

#   1 figure:Hyper methylated regions

hyper_peaks <- GRangesList(DC552_NI18 = hyper_DC552_NI18, DC552_TB18 = hyper_DC552_TB18,
                          DC555_NI18 = hyper_DC555_NI18, DC555_TB18 = hyper_DC555_TB18)
genomicElementDistribution(hyper_peaks,
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))


# 3. 4 figures: ‘Distribution vs genomic features’. For selected donor:

#   2 figures: Hypo methylated regions x 2 conditions

# non-infected
genomicElementDistribution(DC552_NI18,
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000),
                           promoterLevel=list(
                             # from 5' -> 3', fixed precedence 3' -> 5'
                             breaks = c(-2000, -1000, -500, 0, 500),
                             labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", "upstream <500b", "TSS - 500b"),
                             colors = c("#FFE5CC", "#FFCA99", "#FFAD65", "#FF8E32")))

# infected
genomicElementDistribution(DC552_TB18,
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000),
                           promoterLevel=list(
                             # from 5' -> 3', fixed precedence 3' -> 5'
                             breaks = c(-2000, -1000, -500, 0, 500),
                             labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", "upstream <500b", "TSS - 500b"),
                             colors = c("#FFE5CC", "#FFCA99", "#FFAD65", "#FF8E32")))


#   2 figures: Hyper methylated regions x 2 conditions

# non-infected
genomicElementDistribution(hyper_DC552_NI18,
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000),
                           promoterLevel=list(
                             breaks = c(-2000, -1000, -500, 0, 500),
                             labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", "upstream <500b", "TSS - 500b"),
                             colors = c("#FFE5CC", "#FFCA99", "#FFAD65", "#FF8E32")))

# infected
genomicElementDistribution(hyper_DC552_TB18,
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000),
                           promoterLevel=list(
                             breaks = c(-2000, -1000, -500, 0, 500),
                             labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", "upstream <500b", "TSS - 500b"),
                             colors = c("#FFE5CC", "#FFCA99", "#FFAD65", "#FF8E32")))


# 4. 3 figures: distribution vs genomic features for  DMRs (Up / Down / Up+Down)

dmr_down <- toGRanges('meth_CpG_dmrs/dmrs_avgdiff_0.025_ncyto_3.down.bed', format = 'BED', header = FALSE)

genomicElementDistribution(dmr_down,
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000),
                           promoterLevel=list(
                             breaks = c(-2000, -1000, -500, 0, 500),
                             labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", "upstream <500b", "TSS - 500b"),
                             colors = c("#FFE5CC", "#FFCA99", "#FFAD65", "#FF8E32")))



