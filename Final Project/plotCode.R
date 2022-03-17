library(vcfR)
#library(pinfsc50)

# Find the files.
#vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50") #Zaire_Ebola_mapped_snarl_genotypes.vcf.gz
#dna_file <- system.file("extdata", "pinf_sc50.fasta", package = "pinfsc50") #Zaire_Ebola_genomic.fna
#gff_file <- system.file("extdata", "pinf_sc50.gff", package = "pinfsc50") #Zaire_Ebola_GFF.gff

#vcf_file <- system.file(file.path(), mustWork=TRUE, p_ackage="base") #Zaire_Ebola_mapped_snarl_genotypes.vcf.gz
#vcf_file <- "/Users/shwetajones/Desktop/Zaire_Ebola_mapped_snarl_genotypes.vcf.gz"
vcf_file <- "/Users/shwetajones/Desktop/Ebola_Vars_Merged.vcf.gz"
dna_file <- "/Users/shwetajones/Desktop/Zaire_Ebola_genomic.fna" #Zaire_Ebola_genomic.fna
gff_file <- "/Users/shwetajones/Desktop/Zaire_Ebola_GFF.gff" #Zaire_Ebola_GFF.gff

# Input the files.
vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")

x=2
while(x<1000){
  x = x^2
  print(x)
}

# Create a chromR object.
#chrom <- create.chromR(name="contig", vcf=vcf, seq=dna, ann=gff, verbose=TRUE)

#chrom <- masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
#chrom <- proc.chromR(chrom, verbose = TRUE)

#chromoqc(chrom, dp.alpha = 22)

#chrom <- proc.chromR(chrom, verbose=FALSE, win.size=1e4)
#chromoqc(chrom, dp.alpha = 22)

#chrom <- proc.chromR(chrom, verbose=FALSE, win.size=1e3)

#head(chrom)

#plot(chrom)
### Visualization 2
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chrom <- masker(chrom, min_DP = 300, max_DP = 700)
chrom <- proc.chromR(chrom, verbose = FALSE)

head(chrom)

dp <- extract.gt(chrom, element="DP", as.numeric=TRUE)
rownames(dp) <- 1:nrow(dp)
head(dp)

#heatmap.bp(dp[1:55,], col.ramp = colorRampPalette(c("red", "yellow", "#008000"))(100))
heatmap.bp(dp[1:55,])
