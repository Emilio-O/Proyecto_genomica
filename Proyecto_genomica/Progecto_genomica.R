setwd("C:/Users/52442/Desktop/Proyecto_genomica/")
memory.size(max = T)

library(Biostrings)
library(msa)

hCoV19_Mexico <- readDNAStringSet(c("Wuhan_refseq.fasta","1_Estado_de_Mexico.fasta", "2_CDMX.fasta", 
                                    "3_Jalisco.fasta", "4_Veracruz.fasta", "5_Puebla.fasta", 
                                    "6_Guanajuato.fasta", "7_Nuevo_Leon.fasta", "8_Chiapas.fasta", 
                                    "9_Michoacan.fasta", "10_Oaxaca.fasta"))


# install.packages("remotes")
# remotes::install_github("mhahsler/rMSA", force = T)
library("rMSA")
Mafft_aln <- mafft(x = hCoV19_Mexico)
Mafft_aln

# save(Mafft_aln, file = "Mafft_aln.RData")
load("Mafft_aln.RData")




#install.packages("ips")
library(ips)

PI_1 <- as.DNAbin(Mafft_aln)

PI_1_frac <- pis(PI_1, what = "frac", use.ambiguities = FALSE)
PI_1_abs <- pis(PI_1, what = "abs", use.ambiguities = FALSE)
PI_1_ind <- pis(PI_1, what = "ind", use.ambiguities = FALSE)


PI_1_frac
PI_1_abs
PI_1_ind

library(ape)
Mafft_aln_nogaps <- deleteGaps(PI_1, gap.max = nrow(PI_1)-4)

dim(Mafft_aln)
dim(Mafft_aln_nogaps)


PI_2_frac <- pis(Mafft_aln_nogaps, what = "frac", use.ambiguities = FALSE)
PI_2_abs <- pis(Mafft_aln_nogaps, what = "abs", use.ambiguities = FALSE)
PI_2_ind <- pis(Mafft_aln_nogaps, what = "ind", use.ambiguities = FALSE)

PI_2_frac
PI_2_abs
PI_2_ind


PI_hCoV19 <- Mafft_aln_nogaps[ , c(PI_2_ind)]

Estados <- c("Wuh", paste0("E_Mex_", seq("1", "20")), paste0("CDMX_", seq("1", "20")), 
             paste0("Jal_", seq("1", "20")), paste0("Ver_", seq("1", "20")), 
             paste0("Pue_", seq("1", "20")), paste0("Gua_", seq("1", "20")), 
             paste0("NL_", seq("1", "20")), paste0("Chi_", seq("1", "20")), 
             paste0("Mich_", seq("1", "20")), paste0("Oax_", seq("1", "20")))

rownames(PI_hCoV19) <- Estados
rownames(PI_hCoV19)
# save(PI_hCoV19, file = "PI_hCoV19.RData")
# write.nexus.data(PI_hCoV19, file = "PI_hCoV19.nex")
# write.FASTA(PI_hCoV19, file = "PI_hCoV19.fa")









library(seqinr)
library(msa)

PI_hCoV19 <- readDNAMultipleAlignment("PI_hCoV19.fa")

PI_hCoV19 <- msaConvert(PI_hCoV19, type="seqinr::alignment")

D <- dist.alignment(PI_hCoV19, "similarity")
as.matrix(D)



Dendograma <- hclust(D, "average")
plot(Dendograma)


Dendo <- as.phylo(Dendograma)
# write.tree(phy = Dendo, file="Dendo.newick")


# install.packages("png")
library(png)
knitr::include_graphics("Figtree_IPS.png")





# install.packages("bios2mds")
library(bios2mds)

PI_hCoV19 <- import.fasta("PI_hCoV19.fa")


Nombres <- Estados

Estados <- c("Wuh", rep("E_Mex", 20), rep("CDMX", 20), 
             rep("Jal", 20), rep("Ver", 20), rep("Pue", 20), 
             rep("Gua", 20), rep("NL", 20), rep("Chi", 20), 
             rep("Mich", 20), rep("Oax", 20))

Colores <- c("blue", rep("turquoise", 20), rep("purple3", 20), 
             rep("mediumspringgreen", 20), rep("lightpink1", 20), rep("sienna1", 20), 
             rep("magenta", 20), rep("darkred", 20), rep("navyblue", 20), 
             rep("darkgoldenrod1", 20), rep("lightcoral", 20))

Nombres <- sapply(Nombres, function(x) gsub("\"", "", x))
Estados <- sapply(Estados, function(x) gsub("\"", "", x))
Colores <- sapply(Colores, function(x) gsub("\"", "", x))


Grupos <- cbind(Nombres, Estados, Colores)
Grupos <- Grupos[ ,-1]
# write.csv(Grupos, file = "Grupos.csv")


Dist_mat <- mat.dif(PI_hCoV19, PI_hCoV19)
Dist_mat


mmds_hCoV19 <- mmds(Dist_mat, group.file = "C:/Users/52442/Desktop/Proyecto_genomica/Grupos.csv")

scree.plot(mmds_hCoV19$eigen.perc, lab = TRUE, title = "IPS de hCoV19")
mmds.2D <- mmds.2D.plot(mmds_hCoV19, title = "IPS de hCoV19")

mmds.3D <- mmds.3D.plot(mmds_hCoV19, title = "IPS de hCoV19")










library(igraph)
Ady_list <- read.table("FinalytalTABLE.txt")[ ,c(-2)]
Ady_list <- as.matrix(Ady_list)

colnames(Ady_list) <- c("V1", "V2")



Red_popart <- graph_from_edgelist(Ady_list, directed = F)
layout <- layout.fruchterman.reingold(Red_popart)
plot(Red_popart, layout = layout, vertex.size = 4,
     vertex.label = NA, edge.arrow.size = .1, vertex.color="gray50")

plot(degree.distribution(Red_popart))

eb <- edge.betweenness.community(Red_popart)
lb <- label.propagation.community(Red_popart)
cs <- cluster_spinglass(Red_popart)

plot(eb, Red_popart, layout = layout, vertex.size = 4,
     vertex.label = NA, edge.arrow.size = .1, vertex.color="gray50")

plot(lb, Red_popart, layout = layout, vertex.size = 4,
     vertex.label = NA, edge.arrow.size = .1, vertex.color="gray50")

plot(cs, Red_popart, layout = layout, vertex.size = 4,
     vertex.label = NA, edge.arrow.size = .1, vertex.color="gray50")


eccentricity(Red_popart)[1]

betweenness(Red_popart)[1]

closeness(Red_popart)[1]






