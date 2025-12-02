library(LEA)

geno_file = "structure_Ivom384.txt.geno"
snmf_nuc3 <- snmf(geno_file,
               K = 1:32,
               entropy = TRUE,
               repetitions = 10,
               CPU=8,
               project = "new")


save(snmf_nuc3, "snmf_Ivom384.RData")