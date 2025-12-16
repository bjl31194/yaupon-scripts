library(readxl)
negative <- read_excel("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/writing-grants/Yaupon Metabolomics/S5_HRMS_compound_metadata.xlsx", sheet = "negative")

negative <- negative[, -c(6,7,8,9,10,11)]

syn_neg <- read_excel("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/writing-grants/Yaupon Metabolomics/S5_HRMS_compound_metadata.xlsx", sheet = "synthetic_neg", col_names=FALSE)

neg_filtered <- filter(negative, !name %in% syn_neg$...1)
neg_synthetic <- filter(negative, name %in% syn_neg$...1)
write.csv(neg_filtered, "neg_filtered.csv")
write.csv(neg_synthetic, "neg_synthetic.csv")


positive <- read_excel("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/writing-grants/Yaupon Metabolomics/S5_HRMS_compound_metadata.xlsx", sheet = "positive")

positive <- positive[, -c(6)]

syn_pos <- read_excel("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/writing-grants/Yaupon Metabolomics/S5_HRMS_compound_metadata.xlsx", sheet = "synthetic_pos", col_names=FALSE)
pos_filtered <- filter(positive, !name %in% syn_pos$...1)
pos_synthetic <- filter(positive, name %in% syn_pos$...1)

write.csv(pos_filtered, "pos_filtered.csv")
write.csv(pos_synthetic, "pos_synthetic.csv")

S1 <- read_excel("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/writing-grants/Yaupon Metabolomics/S1_neg_sig_assn_age.xlsx")
S1_filtered <- filter(S1,!Compound %in% syn_neg$...1)
S1_syn <- filter(S1, Compound %in% syn_neg$...1)

S2 <- read_excel("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/writing-grants/Yaupon Metabolomics/S2_pos_sig_assn_age.xlsx")
S2_filtered <- filter(S2,!Compound %in% syn_pos$...1)
S2_syn <- filter(S2, Compound %in% syn_pos$...1)

S3 <- read_excel("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/writing-grants/Yaupon Metabolomics/S3_neg_sig_roasting.xlsx")
S3_filtered <- filter(S3,!Compound %in% syn_neg$...1)
S3_syn <- filter(S3, Compound %in% syn_neg$...1)
write.csv(S3_filtered, "S3_filtered.csv")

S4 <- read_excel("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/writing-grants/Yaupon Metabolomics/S4_pos_sig_roasting.xlsx")
S4_filtered <- filter(S4,!Compound %in% syn_pos$...1)
S4_syn <- filter(S4, Compound %in% syn_pos$...1)
write.csv(S4_filtered, "S4_filtered.csv")
