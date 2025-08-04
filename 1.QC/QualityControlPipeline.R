library(plinkQC)

package.dir <- find.package('plinkQC')
indir <- 'rawdata'
qcdir <- 'results'
name <- 'GSAdata_raw'
path2plink <- "/usr/bin/plink"


args <-commandArgs(trailingOnly=True)


prefixMergedDataset <- paste("GSAdata_raw.merge.HapMapIII_CGRCh37", sep="")


pdf(file="individualQC.pdf")
fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
                                    refSamplesFile=paste(indir, "/HapMap_ID2Pop.txt",
                                                         sep=""), 
                                    refColorsFile=paste(indir, "/HapMap_PopColors.txt",
                                                        sep=""),
                                    prefixMergedDataset=prefixMergedDataset,
                                    path2plink=path2plink, dont.check_sex=TRUE,
                                    interactive=TRUE, verbose=TRUE)
dev.off()


pdf(file="individualQC_Overview.pdf")
overview_individuals <- overviewPerIndividualQC(fail_individuals,
                                                interactive=TRUE)
dev.off()

pdf(file="markerQC.pdf")
fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
                            path2plink=path2plink,
                            verbose=TRUE, interactive=TRUE,
                            showPlinkOutput=FALSE)
dev.off()

pdf(file="markerQC_Overview.pdf")
overview_marker <- overviewPerMarkerQC(fail_markers, interactive=TRUE)
dev.off()

Ids  <- cleanData(indir=indir, qcdir=qcdir, name=name, path2plink=path2plink,
                  verbose=TRUE, showPlinkOutput=FALSE)
