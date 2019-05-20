source("http://bioconductor.org/biocLite.R")
biocLite("SNPRelate")
library(SNPRelate)
install.packages('Rcpp')
library(ggplot2)
library(tidyverse)
library(dplyr)
library(qqman)

eigenvector <- read.table("rd_eye_pca.eigenvec")
eigenvalue <- read.table("rd_eye_pca.eigenval")
tekst <- read.table("eye_color.txt")


#PCA plot
plot(eigenvector$V3, eigenvector$V4, xlab=("PC1"), ylab=("PC2"))

#Creating categories for 4 eye colors.
pheno_only_brown <- filter(tekst, tekst$V2=="brown"|  tekst$V2=="dark_brown"|  tekst$V2=="olive-brown_ringing_burnt_umber-brown"|  tekst$V2=="brown/black"| tekst$V2=="brown_-_brown_and_green_in_bright_sunlight"| tekst$V2=="grey_brown"| tekst$V2=="amber_-_(yellow/ocre__brown)")
pheno_hazel <- filter(tekst,tekst$V2=="hazel"| tekst$V2=="hazel/light_brown"|tekst$V2=="hazel_(brown/green)"| tekst$V2=="hazel_(light_brown,_dark_green,_dark_blue)"|tekst$V2=="hazel/yellow")
pheno_green <- filter(tekst, tekst$V2=="blue-green"| tekst$V2=="brown-green"| tekst$V2=="green-hazel" |tekst$V2=="green-brown" | tekst$V2=="green"|tekst$V2=="light-mixed_green"| tekst$V2=="green-gray"|tekst$V2=="amber-green"| tekst$V2=="blue-green-grey"| tekst$V2=="blue-green;_amber_collarette,_and_blue-grey_ringing")
pheno_only_blue<- filter(tekst, tekst$V2=="blue"|  tekst$V2=="blue-grey"| tekst$V2=="grey"| tekst$V2=="blue-grey;_broken_amber_collarette"|  tekst$V2=="blue_with_yellow_parts")

pheno_4_done <- mutate(tekst, Phenotype=if_else(tekst$V1 %in% pheno_only_brown$V1, 1, if_else(tekst$V1 %in% pheno_only_blue$V1, 2, if_else(tekst$V1 %in% pheno_green$V1, 3, 4))))
#Making it output correctly
pheno_4_FID <- pheno_4_done 
pheno_4_FID$V2 <- pheno_4_done$V1
names(pheno_4_FID) <- c("FID", "IID", "Phenotype")

write.table(pheno_4_FID, "eye_color_all_colors.txt", sep="\t", row.names = FALSE, quote = FALSE)

blue_brown <- filter(pheno_4_FID, Phenotype <= 2)
green_hazel <- filter(pheno_4_FID, Phenotype >= 3)

write.table(blue_brown, "eye_color_blue_brown.txt", sep="\t", row.names = FALSE, quote = FALSE)
write.table(green_hazel, "eye_color_green_hazel.txt", sep="\t", row.names = FALSE, quote = FALSE)

#Manhattan plot based on the fishers exact test
fisher <- read.table("fisher.assoc.fisher", head=T)

min(fisher$P, na.rm = TRUE)
manhattan(fisher)
qq(fisher$P)
head(fisher)
summary(fisher)
#significant SNPs
sum(fisher$P <= 0.05)
#Total SNPs
sum(fisher$P <= 1)
#Significant SNPs with bonferroni
sum(fisher$P <= 0.05/nrow(fisher))
#finding the most significant SNP.
which.min(fisher$P)
fisher$BP[749376]
#inflation
fish_mutate <- mutate(fisher, chi_s = qchisq(P, df=1, lower.tail = F))
median(fish_mutate$chi_s, na.rm = TRUE)/0.5


#logistic regression
logistic_3 <- read.table("logistic_3.assoc.logistic", head=T)

logistic_3_filter <- na.omit(filter(logistic_3, TEST == "ADD"))
which.min(logistic_3_filter$P)
logistic_3_filter$SNP[324534]
sum(logistic_3_filter$P <= 0.05/nrow(logistic_3_filter))
manhattan(logistic_3_filter)

qq(logistic_3_filter$P)
sum(logistic_3_filter$P <= 0.05/nrow(logistic_3_filter), na.rm = TRUE)
log_mutate <- mutate(logistic_3_filter, chi_s = qchisq(P, df=1, lower.tail = F))
median(log_mutate$chi_s, na.rm = TRUE)/0.5
which.min(logistic_3_filter$P)
logistic_3_filter$SNP[324534]
logistic_3_filter$BP[324534]


logistic_10 <- read.table("logistic_10.assoc.logistic", head=T)
logistic_10_filter <- na.omit(filter(logistic_10, TEST == "ADD"))
qq(logistic_10_filter$P)
manhattan(logistic_10_filter)
log_10_mutate <- mutate(logistic_10_filter, chi_s = qchisq(P, df=1, lower.tail = F))
median(log_10_mutate$chi_s, na.rm = TRUE)/0.5
sum(logistic_10_filter$P <= 0.05/nrow(logistic_10_filter))
which.min(logistic_10_filter$P)
logistic_10_filter$SNP[691572]
logistic_10_filter$BP[691572]

#additive, dominant or recessive?
recode <- read.table("snp_recode.raw", head=T)
recode_filter <- filter(recode, PHENOTYPE != -9, is.na(rs1667394_C)== FALSE)
Genotypes <- matrix(c(sum(recode_filter$PHENOTYPE == 1 & recode_filter$rs1667394_C==0),sum(recode_filter$PHENOTYPE == 1 & recode_filter$rs1667394_C==1),sum(recode_filter$PHENOTYPE == 1 & recode_filter$rs1667394_C==2),
                      sum(recode_filter$PHENOTYPE == 2 & recode_filter$rs1667394_C==0),sum(recode_filter$PHENOTYPE == 2 & recode_filter$rs1667394_C==1),sum(recode_filter$PHENOTYPE == 2 & recode_filter$rs1667394_C==2)), ncol=3, byrow = TRUE)
sum(recode_filter$PHENOTYPE == 1 & recode_filter$rs1667394_C==0)
sum(recode_filter$PHENOTYPE == 2 & recode_filter$rs1667394_C==0)
colnames(Genotypes) <- c(0,1,2)
rownames(Genotypes) <- c("Brown Eyes", "Blue Eyes")
barplot(Genotypes, ylab = "Individuals", xlab = "Genotypes")
percentage <- c(Genotypes[1]/sum(Genotypes[1:2]), Genotypes[3]/sum(Genotypes[3:4]),Genotypes[5]/sum(Genotypes[5:6]))
names(percentage) <- c(0,1,2)

#conditioning
condition_10_additive <- read.table("logistic_condition.assoc.logistic", head=T)
condition_10_filter <- na.omit(filter(condition_10_additive, TEST == "ADD"))
manhattan(condition_10_filter)
qq(condition_10_filter$P)
sum(p.adjust(condition_10_filter$P))
sum(condition_10_filter$P <= 0.05/nrow(condition_10_filter))
sum(condition_10_filter$P <= 0.05/nrow(condition_10_filter))
which.min(condition_10_filter$P)
condition_10_filter$SNP[713540]
#conditioning for dominant or recessive does not work out.
condition_10_dominant <- read.table("logistic_condition.assoc.logistic", head=T)
condition_10_filter_dominant <- na.omit(filter(condition_10_dominant, TEST == "ADD"))
manhattan(condition_10_filter_dominant)

condition_10_recessive <- read.table("logistic_condition_recessive.assoc.logistic", head=T)
condition_10_filter_recessive <- na.omit(filter(condition_10_recessive, TEST == "ADD"))
manhattan(condition_10_filter_recessive)

#green vs hazel
gh_10 <- read.table("logistic_10_green_hazel.Phenotype.assoc.linear", head=T)
gh_10_filter <- na.omit(filter(gh_10, TEST == "ADD"))
manhattan(gh_10_filter)
sum(gh_10_filter$P <= 0.05/nrow(gh_10_filter))
(gh_10_filter$P <= 0.05/nrow(gh_10_filter))
gh_10_filter$BP[709528]

gh_10_chr15 <- read.table("logistic_10_green_hazel_chr15.Phenotype.assoc.linear", head=T)
gh_10_chr15_filter <- na.omit(filter(gh_10_chr15, TEST == "ADD"))
manhattan(gh_10_chr15_filter)
qq(gh_10_chr15_filter$P)
sum(gh_10_chr15_filter$P <= 0.05/nrow(gh_10_chr15_filter))

#gh conditioning.
gh_cond <- read.table("logistic_10_green_hazel.Phenotype.assoc.linear", head=T)
gh_cond_filter <- na.omit(filter(gh_cond, TEST == "ADD"))
sum(gh_cond_filter$P <= 0.05/nrow(gh_cond_filter))


#additional analysis, 4 phenotypes and epistasis.
pheno_4 <- read.table("logistic_10_all.Phenotype.assoc.linear", head=T)
pheno_4_filter <- na.omit(filter(pheno_4, TEST == "ADD"))
sum(pheno_4_filter$P <= 0.05/nrow(pheno_4_filter))
manhattan(pheno_4_filter)

#epistasis
epistasis_chr15 <- read.table("first_set_epistasis.epi.qt", head=T)

#All significant SNP epistasis
epistasis_chr15_significant <- read.table("set_epistasis_significant.epi.qt", head=T)
sum(epistasis_chr15_significant$P <= 0.25/189956)

#All significant, only against each other
epistasis_set_by_set <- read.table("set_by_set_epistasis.epi.qt", head=T)
sum(p.adjust(epistasis_set_by_set$P, method="none") <= 0.05)

#The two most significant, against all in the range defined by the found significant SNPs.
epistasis_restricted <- read.table("set_restricted_epistasis.epi.qt", head=T)
sum(epistasis_restricted$P <= 0.10/39)
#finding phenotypes
snp_epistasis <- na.omit(read.table("snp_4_pheno.raw", head=T))

x=1
eye_color_brown <- (c(sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==0),
                      sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==1),
                      sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==2),
                      sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==0),
                      sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==1),
                      sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==2),
                      sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==0),
                      sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==1),
                      sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==2)))
x=2
eye_color_blue <- (c(sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==0),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==1),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==2),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==0),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==1),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==2),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==0),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==1),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==2)))
x=3
eye_color_green <- (c(sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==0),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==1),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==2),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==0),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==1),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==2),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==0),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==1),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==2)))
x=4
eye_color_hazel <- (c(sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==0),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==1),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 0 & snp_epistasis$rs11636232_T ==2),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==0),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==1),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 1 & snp_epistasis$rs11636232_T ==2),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==0),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==1),
                            sum(snp_epistasis$PHENOTYPE == x & snp_epistasis$rs1129038_C == 2 & snp_epistasis$rs11636232_T ==2)))
total_matrix <- rbind(eye_color_brown, eye_color_blue, eye_color_green, eye_color_hazel)

colnames(total_matrix) <- c("00", "01", "02", "10", "11", "12", "20", "21", "22")

barplot(total_matrix, main = "Genotype and Epistasis", xlab = "Phenotype", col = c("chocolate4", "skyblue1", "seagreen3", "tan2"))


#Not significant
snp_epistasis_ns <- na.omit(read.table("snp_4_pheno_not_significant.raw", head=T))
x=1
eye_color_brown_ns <- (c(sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==2),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==2),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==2)))
x=2
eye_color_blue_ns <- (c(sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==2),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==2),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==2)))
x=3
eye_color_green_ns <- (c(sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==2),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==2),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==2)))
x=4
eye_color_hazel_ns <- (c(sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 0 & snp_epistasis_ns$rs8043146_T ==2),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 1 & snp_epistasis_ns$rs8043146_T ==2),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==0),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==1),
                         sum(snp_epistasis_ns$PHENOTYPE == x & snp_epistasis_ns$rs916977_T == 2 & snp_epistasis_ns$rs8043146_T ==2)))
total_matrix_ns <- rbind(eye_color_brown_ns, eye_color_blue_ns, eye_color_green_ns, eye_color_hazel_ns)

colnames(total_matrix_ns) <- c("00", "01", "02", "10", "11", "12", "20", "21", "22")

barplot(total_matrix_ns, main = "Genotype and Epistasis", xlab = "Phenotype", col = c("chocolate4", "skyblue1", "seagreen3", "tan2"), beside = T)
