GFFtoBED<- function(gtf, assembly_report) {

### SELECT AND MANE ### GTF

MANE_RefSeqSelect <- read.delim(gtf, header=FALSE, comment.char="#")
MANE_RefSeqSelect = MANE_RefSeqSelect[MANE_RefSeqSelect$V3=="exon",]
Refseq_Mane_select = MANE_RefSeqSelect[grepl("RefSeq Select|MANE Select",MANE_RefSeqSelect$V9,),]
rm(MANE_RefSeqSelect)

library(stringr)

information = as.data.frame(str_split_fixed(Refseq_Mane_select$V9, ";", 10))
###tag= as.data.frame(str_split_fixed(information, "tag ", 2))
gene= as.data.frame(str_split_fixed(information$V1, "gene_id ", 2))
trascript= as.data.frame(str_split_fixed(information$V2, "transcript_id", 2))
exon= as.data.frame(str_split_fixed(Refseq_Mane_select$V9, "exon_number ", 2))
exon$V2= str_replace(exon$V2, "; ", "")
name=str_c(gene$V2,exon$V2,sep = "_")

#hg38.p13.chromAlias <- read.delim("~/Downloads/Panel_Roche/hg38.p13.chromAlias.txt", header=FALSE, row.names=NULL, comment.char="#")
GCF_000001405_39_GRCh38_p13_assembly_report <- read_delim(assembly_report, 
                                                          delim = "\t", escape_double = FALSE, 
                                                          col_names = FALSE, comment = "#", trim_ws = TRUE)



BED_MANE = data.frame(chr=Refseq_Mane_select$V1, start=Refseq_Mane_select$V4, 
                      end=Refseq_Mane_select$V5, name=name,strand=Refseq_Mane_select$V7 ,gen=gene$V2, exon=exon$V2)

chr= GCF_000001405_39_GRCh38_p13_assembly_report$X10[match(BED_MANE$chr, GCF_000001405_39_GRCh38_p13_assembly_report$X7)]
BED_MANE$chr = chr
return(BED_MANE)
}

#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################

panel_gens = c("ATM","BRCA1")

Exons_to_introns<- function(BEDMANE, panel_gens) {
Intrones = data.frame(matrix(ncol = 7))
names(Intrones)<- c("chr","start","end","name","strand","gen","intron")
for (gen in panel_gens) {
  prueba=BEDMANE[ which(BEDMANE$gen==gen),]
  if (prueba$strand[1]=="+") {
    for (num in c(1:length(prueba$exon))) {
      Intrones=rbind(Intrones,(c(prueba$chr[num], prueba$end[num], prueba$start[num+1],prueba$name[num],
                                 prueba$strand[num], prueba$gen[num], prueba$exon[num] )))
    }
  }
  else{
    for (num in c(1:length(prueba$exon))) {
      Intrones=rbind(Intrones,(c(prueba$chr[num],prueba$end[num+1],prueba$start[num],prueba$name[num],
                                 prueba$strand[num], prueba$gen[num], prueba$exon[num] )))
    }
  }
}
Intrones = na.omit(Intrones)
name=paste(Intrones$gen,Intrones$intron,"Intron",sep = "_")
Introns_bed= data.frame(chr=Intrones$chr,start=Intrones$start,end=Intrones$end, name=name)
return(Introns_bed)
}


