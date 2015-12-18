# Convert shapeit file to ihs 
# Released 
#
# shapeit (http://www.shapeit.fr/pages/m02_formats/hapssample.html)
# ihs: (http://hgdp.uchicago.edu/Software/)
# Uses xor operation to speed up conversion
#
# Input:
# 1. shapeit haplotype file (http://www.shapeit.fr/pages/m02_formats/hapssample.html)
# 2. genetic map file (cols: bp, recomb rate, map)
# 3. ancestry allele file (cols: snp, ancestral allele)
# 4. output filename
#
# Output:
# ihsmap file and ihshapfile
#
# Example Usage: 
# Rscript convert-shapeit-to-ihs.R \ 
# chr1.haps genetic_map_chr1_b37.txt \
# chr1.derived_allele.txt outname_chr1

args <- commandArgs(TRUE)
file_hap <- args[1]
file_map <- args[2]
file_anc <- args[3]

hap <- read.table(file_hap)
names(hap)[2] <- "snp"

ancestry <- read.table(file_anc)
names(ancestry) <- c("snp","aa")
hapm <- merge(hap,ancestry,by="snp",sort=FALSE)
names(hapm)[4:5] <- c("A1","A2")

# Merge in ancestral allele
hapm$aa <- toupper(hapm$aa)
# remove rows with mismatching ancestral allele
hapm <- hapm[hapm$A1==hapm$aa & hapm$A2==hapm$aa,]

names(hapm)[3] <- "bp"
map <- read.table(file_map,header=T)
names(map) <- c("bp","comb","map")
map <- map[c("bp","map")]
mapfile <- merge(hapm,map,by="bp",sort=FALSE)
nn <- ncol(mapfile)
mapout <- mapfile[,c(2,1,nn,4,5,nn-2)]

# WRITE MAP FILE
# Convert A1 to other allele (0), A2 to ancestral allele (1)
# (ie. if A2 isn't the ancestral allele, swap them - name only not genotype values)
ndx <- mapout$aa!=mapout$A2
mapout[ndx,]$A1 <- mapout[ndx,]$A2
mapout[ndx,]$A2 <- mapout[ndx,]$aa
mapout <- mapout[,1:5]
write.table(mapout,paste(args[4],'.ihsmapfile',sep=""),quote=F,row.names=F,col.names=F)

# to use the XOR operation to quickly match alleles (instead of NOT XOR)
# the which values are inversed from above representation (in map output file)
# hence 1 is A1, 0 is A2
mapfile[mapfile$A1==mapfile$aa,]$which <- 1
mapfile[mapfile$A2==mapfile$aa,]$which <- 0
nn  <- ncol(mapfile)
out <- cbind(mapfile[c("which")],mapfile[,6:(nn-3)])

# first column is 'which'
# Using xor, write 0 for not ancestral, 1 for ancestral allele
nn <- ncol(out)-1
for (ind in 2:nn) {
    tmp <- bitwXor(as.numeric(out[,ind]),out[,1])
    out[,ind] <- tmp
}

# merging columns into one before transposing 
ihshap <- c()
for (ind in seq(2,nn,2)){
    out[,ind] <- paste(out[,ind],out[,ind+1])
}

# only take second rows - which is merged combination of both, before turning cols to rows
out_b <- out[,seq(2,nn,2)]
# cols to rows
output <- t(out_b)
# WRITE HAP FILE
write.table(output,paste(args[4],'.ihshap',sep=""),quote=F,row.names=F,col.names=F)
