library(LEA)

structure_file="structure_Ivom384.txt" #nuclear loci only

struct2geno(input.file=structure_file, 
            ploidy=2, #for diploids
            FORMAT=2, #2 rows of data per individual
            extra.row = 1, #rows before genotype data (locusID)
            extra.col=2) #columns before genotype data (sampleID and pop)