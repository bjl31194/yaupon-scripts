library(LEA)

geno_file = "structure_Ivom384.txt.geno"
snmf_nuc3 <- snmf(geno_file,
               K = 1:32,
               entropy = TRUE,
               repetitions = 10,
               CPU=8,
               project = "new")


save(snmf_nuc3, file = "snmf_Ivom384.RData")

load("/Users/ben/Library/CloudStorage/OneDrive-UniversityofGeorgia/yaupon/vcf/plates1234/snmf_data_Ivom384.RData")

names <- read.table("names.txt")

pop_order <- c("AR","TX","LA","MS","AL","FL","GA","SC","NC","VA","MC")

K2_Q <- as.data.frame(snmf_nuc_Q_K2) %>% #convert to a dataframe
  mutate(sample = names$V1) %>% #change the name to sample
  mutate(site = states$V1) # make a site column by shortening sample names - this will be different for you if you do it
K2_Q$site <- factor(K2_Q$site, levels=pop_order)
K2_Q <- K2_Q[order(K2_Q$site),]
K2_Q$sample <- factor(K2_Q$sample, levels = K2_Q$sample) #set the order in which "sample" variable is plotted
K2_Q$group <- max.col(K2_Q[,1:2]) #For my map, I assigned each individual to it's maximum identity group. 
K2Q <- melt(K2_Q, id.vars = c("sample", "site"),
            measure.vars = c("V1", "V2"),
            variable.name = "Cluster",
            value.name = "Ancestry")

cols <- c("red","blue")
K2_plot <- ggplot(K2Q, aes(x = sample, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Individuals", y = "Ancestry proportion")
K2_plot

## K=3 ##
K3_Q <- as.data.frame(snmf_nuc_Q_K3) %>% #convert to a dataframe
  mutate(sample = names$V1) %>% #change the name to sample
  mutate(site = states$V1) # make a site column by shortening sample names - this will be different for you if you do it
K3_Q$site <- factor(K3_Q$site, levels=pop_order)
K3_Q <- K3_Q[order(K3_Q$site),]
K3_Q$sample <- factor(K3_Q$sample, levels = K3_Q$sample) #set the order in which "sample" variable is plotted
K3_Q$group <- max.col(K3_Q[,1:3]) #For my map, I assigned each individual to it's maximum identity group. 
K3Q <- melt(K3_Q, id.vars = c("sample", "site"),
            measure.vars = c("V1", "V2","V3"),
            variable.name = "Cluster",
            value.name = "Ancestry")

cols <- c("red","blue","darkgreen")
K3_plot <- ggplot(K3Q, aes(x = sample, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Individuals", y = "Ancestry proportion")
K3_plot

## K=4 ##
K4_Q <- as.data.frame(snmf_nuc_Q_K4) %>% #convert to a dataframe
  mutate(sample = names$V1) %>% #change the name to sample
  mutate(site = states$V1) # make a site column by shortening sample names - this will be different for you if you do it
K4_Q$site <- factor(K4_Q$site, levels=pop_order)
K4_Q <- K4_Q[order(K4_Q$site),]
K4_Q$sample <- factor(K4_Q$sample, levels = K4_Q$sample) #set the order in which "sample" variable is plotted
K4_Q$group <- max.col(K4_Q[,1:4]) #For my map, I assigned each individual to it's maximum identity group. 
K4Q <- melt(K4_Q, id.vars = c("sample", "site"),
            measure.vars = c("V1", "V2","V3", "V4"),
            variable.name = "Cluster",
            value.name = "Ancestry")

cols <- c("red","blue","darkgreen", "orange")
K4_plot <- ggplot(K4Q, aes(x = sample, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Individuals", y = "Ancestry proportion")
K4_plot

## K=5 ##

K5_Q <- as.data.frame(snmf_nuc_Q_K5) %>% #convert to a dataframe
  mutate(sample = names$V1) %>% #change the name to sample
  mutate(site = states$V1) # make a site column by shortening sample names - this will be different for you if you do it
K5_Q$site <- factor(K5_Q$site, levels=pop_order)
K5_Q <- K5_Q[order(K5_Q$site),]
K5_Q$sample <- factor(K5_Q$sample, levels = K5_Q$sample) #set the order in which "sample" variable is plotted
K5_Q$group <- max.col(K5_Q[,1:5]) #For my map, I assigned each individual to it's maximum identity group. 
K5Q <- melt(K5_Q, id.vars = c("sample", "site"),
            measure.vars = c("V1", "V2","V3","V4","V5"),
            variable.name = "Cluster",
            value.name = "Ancestry")

cols <- c("red","blue","darkgreen","orange","darkviolet")
K5_plot <- ggplot(K5Q, aes(x = sample, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = cols) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Individuals", y = "Ancestry proportion")
K5_plot






