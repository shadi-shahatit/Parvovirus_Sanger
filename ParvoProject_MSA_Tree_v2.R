#### Sanger sequence analysis MSA - v2
#### Shadi Shahatit, MBV Lab, JUST, 2026
# libraries ---------------------------------------------------------------

sys_dir <- "/home/shadi/Desktop/ParvoProject/"

library(ggbreak)
library(systemPipeR)
library(Biostrings)
library(msa)
library(ape)
library(phangorn)
library(ggplot2)
library(ggrepel)
library(ggalt)
library(ggmsa)
library(pegas)
library(haplotypes)
library(ggtree)
library(viridis)
library(patchwork)
library(tidyverse)    # dplyr, tidyr, ggplot2, readr
library(stringr)
library(ggnewscale)
library(readxl)
library(writexl)

# MSA, Trees, and others - OLD ---------------------------------------------------------------------

## load unaligned sequences

fasta_dir <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\Samples\\"
fa_files <- list.files(fasta_dir, pattern = ".fa$", full.names = TRUE)
seqs_list <- lapply(fa_files, function(f) {
  s <- readDNAStringSet(f)
  clean_name <- gsub("\\.fa$", "", basename(f))
  names(s) <- clean_name
  s
})

seqs <- do.call(c, seqs_list)

stats <- do.call(rbind, lapply(fa_files, function(f) {
  s <- readDNAStringSet(f)
  seq <- s[[1]]
  data.frame(
    sample    = gsub("\\.fa$", "", basename(f)),
    length    = Biostrings::width(seq),
    A_base         = letterFrequency(seq, "A"),
    T_base         = letterFrequency(seq, "T"),
    G_base         = letterFrequency(seq, "G"),
    C_base         = letterFrequency(seq, "C"),
    N_base         = letterFrequency(seq, "N"),
    GC_pct    = round(letterFrequency(seq, "GC", as.prob = TRUE) * 100, 2))}))

## translate 

aa <- translate(seqs)

orfs <- predORF(seqs,
                n = "all",
                type = "gr",
                mode = "ORF",
                strand = "both",
                longest_disjoint = TRUE)

orfs <- predORF(seqs,
                n = "all",
                type = "df",
                mode = "ORF",
                strand = "both",
                longest_disjoint = TRUE)

## multiple sequence alignment via ClustalW

aln <- msa(seqs, method = "ClustalW", type = "dna")
aln_dna <- as(aln, "DNAStringSet")

ggmsa(aln_dna, start = 1, end = 500,
      color = "Chemistry_NT", font = "DroidSansMono", seq_name = TRUE)

## PCA on distance matrix

pcoa  <- cmdscale(dist_mat, k = 2, eig = TRUE)
var_explained <- round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)

pca_df <- data.frame(
  PC1    = pcoa$points[, 1],
  PC2    = pcoa$points[, 2],
  sample = rownames(pcoa$points),
  group  = ifelse(grepl("feline", rownames(pcoa$points), ignore.case = TRUE), "Feline",
                  ifelse(grepl("canine", rownames(pcoa$points), ignore.case = TRUE), "Canine",
                         "Lion")))

ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_manual(values = c("Feline" = "steelblue", "Canine" = "tomato", "Lion" = "forestgreen")) +
  labs(
    title = "PCoA of pairwise DNA distances",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = "Group") +
  theme_bw()

## tree for our samples

## convert to phyDat for distance/tree

aln_phyDat <- as.phyDat(as.DNAbin(aln_dna))

## NJ tree

dist_mat <- dist.dna(as.DNAbin(aln_dna), model = "K80")   # Kimura 2-param
nj_tree  <- NJ(dist_mat)
nj_tree  <- ladderize(nj_tree)

# plot tree

plot(nj_tree, main = "Neighbor-Joining Tree (K80)",
     cex = 0.8, tip.color = "black")
add.scale.bar()
dev.off()

## tree for our samples and other regional seq and ref seq

fasta_dir_all <- "C:\\Shadi's File\\Shadi_Transcriptomics_RA\\ParvoProject\\ParvoGenomeRefs\\parvo_allseq.fasta"
ref_seqs <- readDNAStringSet(fasta_dir_all)

seqs_all <- c(seqs, ref_seqs)

aln      <- msa(seqs_all, method = "ClustalW", type = "dna")
aln_dna  <- as(aln, "DNAStringSet")

aln_dnabin_all <- as.matrix(as.DNAbin(aln_dna))
dist_mat_all   <- dist.dna(aln_dnabin_all, model = "K80")
nj_tree_all    <- NJ(dist_mat_all)
nj_tree_all    <- ladderize(nj_tree_all)

tip_colors <- ifelse(grepl("feline", nj_tree_all$tip.label, ignore.case = TRUE), "steelblue",
                     ifelse(grepl("canine", nj_tree_all$tip.label, ignore.case = TRUE), "tomato",
                            ifelse(grepl("lion",   nj_tree_all$tip.label, ignore.case = TRUE), "forestgreen",
                                   ifelse(grepl("FJ222821|M24003|M38245|M38246|M74849", nj_tree_all$tip.label), "purple",
                                          "gray40"))))

plot(nj_tree_all, main = "Neighbor-Joining Tree (K80)",
     cex = 0.6, tip.color = tip_colors)
add.scale.bar()
legend("bottomleft",
       legend = c("Feline", "Canine", "Lion", "Reference (known)", "Other"),
       text.color = c("steelblue", "tomato", "forestgreen", "purple", "gray40"),
       bty = "n")

# tip_df <- data.frame(
#   label  = nj_tree_all$tip.label,
#   group  = ifelse(grepl("sample", nj_tree_all$tip.label, ignore.case = TRUE) &
#                     grepl("feline", nj_tree_all$tip.label, ignore.case = TRUE), "My Feline",
#                   ifelse(grepl("sample", nj_tree_all$tip.label, ignore.case = TRUE) &
#                            grepl("canine", nj_tree_all$tip.label, ignore.case = TRUE), "My Canine",
#                          ifelse(grepl("sample", nj_tree_all$tip.label, ignore.case = TRUE) &
#                                   grepl("lion",   nj_tree_all$tip.label, ignore.case = TRUE), "My Lion",
#                                 ifelse(grepl("FJ222821|M24003|M38245|M38246|M74849", nj_tree_all$tip.label), "Reference",
#                                        "Database"))))
# )

tip_df <- data.frame(
  label = nj_tree_all$tip.label,
  group = ifelse(grepl("sample", nj_tree_all$tip.label, ignore.case = TRUE) & grepl("feline", nj_tree_all$tip.label, ignore.case = TRUE), "My_Feline",
                 ifelse(grepl("sample", nj_tree_all$tip.label, ignore.case = TRUE) & grepl("canine", nj_tree_all$tip.label, ignore.case = TRUE), "My_Canine",
                        ifelse(grepl("sample", nj_tree_all$tip.label, ignore.case = TRUE) & grepl("lion",   nj_tree_all$tip.label, ignore.case = TRUE), "My_Lion",
                               ifelse(grepl("FJ222821|M24003|M38245|M38246|M74849", nj_tree_all$tip.label), "Reference",
                                      ifelse(grepl("feline", nj_tree_all$tip.label, ignore.case = TRUE), "Database_Feline",
                                             ifelse(grepl("canine", nj_tree_all$tip.label, ignore.case = TRUE), "Database_Canine",
                                                    "Database"))))))
)

ggtree(nj_tree_all) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 2, fontface = "bold", align = FALSE)+
  scale_color_manual(values = c(
    "My_Feline" = "tomato",
    "My_Canine" = "steelblue",
    "My_Lion"   = "purple",
    "Reference" = "forestgreen",
    "Database_Feline" = "gray80",
    "Database_Canine" = "gray40",
    "Database"  = "gray20")) +
  theme_tree2() +
  labs(title = "Neighbor-Joining Tree (K80)", color = "Group") +
  theme(legend.position = "bottom")

## pop gene stats (viral, haploid)

## prep: convert alignment to DNAbin and define groups

aln_dnabin <- as.DNAbin(aln_dna)
group <- ifelse(grepl("feline", names(aln_dna), ignore.case = TRUE), "Feline",
                ifelse(grepl("canine", names(aln_dna), ignore.case = TRUE), "Canine", "Lion"))
aln_dnabin <- as.matrix(as.DNAbin(aln_dna))

## π and Hd per group

groups <- c("Feline", "Canine", "Lion")

div_stats <- do.call(rbind, lapply(groups, function(g) {
  sub_aln <- aln_dnabin[group == g, ]
  nuc     <- nuc.div(sub_aln)
  hap     <- pegas::haplotype(sub_aln)
  hd      <- pegas::hap.div(hap)
  data.frame(group = g, n = sum(group == g),
             pi = round(nuc, 6), Hd = round(hd, 4))}))

## Dxy between Canine and Feline

feline_idx <- which(group == "Feline")
canine_idx <- which(group == "Canine")
both_aln   <- aln_dnabin[c(feline_idx, canine_idx), ]
d_mat      <- dist.dna(both_aln, model = "raw", as.matrix = TRUE)
dxy        <- mean(d_mat[seq_along(feline_idx), (length(feline_idx) + 1):ncol(d_mat)])
da         <- dxy - (nuc.div(aln_dnabin[feline_idx, ]) + nuc.div(aln_dnabin[canine_idx, ])) / 2
cat("Dxy (Canine vs Feline):", round(dxy, 6), "\n")
cat("Da  (Canine vs Feline):", round(da,  6), "\n")

## haplotype network (median joining via pegas)

hap_all <- pegas::haplotype(aln_dnabin)
net     <- pegas::haploNet(hap_all)
plot(net, size = attr(net, "freq"), show.mutation = 1,
     main = "Haplotype Network")

## dN/dS per pairwise comparison

dn_ds <- dnds(aln_dna)

"red""red""red""red""red""red""red""red""red""red""red""red""red""red""red""red"

aa_dir <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\Samples\\aa_translate\\SS_aa_translate_alin_with3refs_with4vac.fa"

aa_all <- readAAStringSet(aa_dir)

library(Biostrings)
library(msa)
library(seqinr)
library(ape)
library(ggtree)

# read
aa_dir <- "C:\\Users\\Shadi Shahatit\\OneDrive\\Desktop\\Samples\\aa_translate\\SS_aa_translate_with3refs_with4vac.fa"
aa_all <- readAAStringSet(aa_dir)

# dedupe identical seqs (handles duplicate REF_46 if truly identical), else just uniquify names
if (any(duplicated(as.character(aa_all)))) {
  aa_all <- aa_all[!duplicated(as.character(aa_all))]
}
names(aa_all) <- make.unique(names(aa_all))

# align
aln <- msa(aa_all, method = "Muscle")
aln_mat <- as.matrix(as(aln, "AAStringSet"))   # rows=seqs, cols=alignment positions
aln_seqinr <- msaConvert(aln, type = "seqinr::alignment")

# tree
d <- dist.alignment(aln_seqinr, matrix = "identity")
tree_ClustalW <- nj(d)
tree_ClustalO <- nj(d)
tree_MUSCLE <- nj(d)

tree_ClustalW
tree_ClustalO
tree_MUSCLE 

tree_rooted_ClustalW <- root(tree_ClustalW, outgroup = ref_name, resolve.root = TRUE)
tree_rooted_ClustalO <- root(tree_ClustalO, outgroup = ref_name, resolve.root = TRUE)
tree_rooted_MUSCLE <- root(tree_MUSCLE, outgroup = ref_name, resolve.root = TRUE)

ggtree(tree_rooted_ClustalW) + geom_tiplab() + geom_treescale()
ggtree(tree_rooted_ClustalO) + geom_tiplab() + geom_treescale()
ggtree(tree_rooted_MUSCLE) + geom_tiplab() + geom_treescale()

ref_name <- "REF_46"     # set reference

names(aa_all)

ref_name <- "vaccineFeline_Zoetis"     # set reference
ref_name <- "vaccine_MSD"     # set reference
ref_name <- "vaccineCanine_MSD"     # set reference
ref_name <- "vaccineCanine_Zoetis"     # set reference

ref_name <- "vaccineFeline_Zoetis"     # set reference

tree_rooted <- root(tree, outgroup = ref_name, resolve.root = TRUE)
ggtree(tree_rooted) + geom_tiplab() + geom_treescale()

# variant calling vs reference
ref_seq <- aln_mat[ref_name, ]
variants <- lapply(rownames(aln_mat), function(sn) {
  q <- aln_mat[sn, ]
  diffs <- which(q != ref_seq & q != "-" & ref_seq != "-")
  if (length(diffs) == 0) return(NULL)
  data.frame(seq = sn, pos = diffs, ref_aa = ref_seq[diffs], alt_aa = q[diffs],
             mutation = paste0(ref_seq[diffs], diffs, q[diffs]))
})
variants_df <- do.call(rbind, variants)

# MSA - DNA & Protein analysis --------------------------------------------------------

#### DNA and VP2 protein analysis

## load the files

dna_file <- file.path(sys_dir, "Parvo_sanger_full_seq","combined_sample_vac_ref_ncbi_out.fasta")
protein_file <- file.path(sys_dir,"Parvo_sanger_full_seq/Expasy_translate/Expasy_output/best_orf","combined_best_orf.fasta")

dna_sequences <- readDNAStringSet(dna_file)
protein_sequences <- readAAStringSet(protein_file)

dna_sequences                     
length(dna_sequences)              
width(dna_sequences)        
names(dna_sequences)       
dna_sequences[[1]]               
as.character(dna_sequences[1])     
letterFrequency(dna_sequences, letters = "GC", as.prob = TRUE)

## naming issues

## clean names of FASTA HEADERS
## spaces/commas/slashes -> underscore
## accession_species_country for NCBI viruses

dna_names <- names(dna_sequences)
pipe_hdrs <- dna_names[grepl("\\|", dna_names)]
pipe_parts <- strsplit(pipe_hdrs, "\\|")
acc     <- trimws(sapply(pipe_parts, `[`, 1))
species <- trimws(sapply(pipe_parts, `[`, 3))
country <- trimws(sapply(pipe_parts, `[`, 4))
length(unique(acc))
unique(species)
unique(country)
ncbi_viruses_map <- setNames(paste0(gsub("_+", "_", gsub("[ ,/]+", "_", trimws(species))), "_",
                                    gsub("_+", "_", gsub("[ ,/]+", "_", trimws(country)))),
                             acc)

outgroup_map <- c(
  "NC_001510.1" = "Protoparvovirus_rodent1", 
  "NC_001718.1" = "Protoparvovirus_ungulate1",
  "NC_029797.1" = "Megabat_bufavirus",
  "NC_038544.1" = "Primate_protoparvovirus1"
  )

rename_fasta_headers <- function(hdrs) {
  first_tok <- trimws(sub("\\s*\\|.*", "", hdrs))   # strip pipe metadata
  first_tok <- sub("\\s.*", "", first_tok)          # strip anything after first space
  out <- first_tok
  in_ncbi <- first_tok %in% names(ncbi_viruses_map)
  in_outgroup <- first_tok %in% names(outgroup_map)
  out[in_ncbi]     <- paste0(first_tok[in_ncbi], "_", ncbi_viruses_map[first_tok[in_ncbi]])
  out[in_outgroup] <- paste0(first_tok[in_outgroup], "_", outgroup_map[first_tok[in_outgroup]])
  out
}

dna_sequences_mod <- dna_sequences
protein_sequences_mod <- protein_sequences

names(dna_sequences_mod) <- rename_fasta_headers(dna_names)
names(protein_sequences_mod) <- rename_fasta_headers(names(protein_sequences))

names(dna_sequences)
names(dna_sequences_mod) 
names(protein_sequences) 
names(protein_sequences_mod) 
dna_sequences_mod                     
length(dna_sequences_mod)              
width(dna_sequences_mod)        
names(dna_sequences_mod)   
protein_sequences_mod                     
length(protein_sequences_mod)              
width(protein_sequences_mod)        
names(protein_sequences_mod) 

## Alignment

## DNA

dna_alignment_ClustalO <- msa(dna_sequences_mod, method = "ClustalOmega")
dna_phy_ClustalO <- phyDat(as.matrix(dna_alignment_ClustalO), type = "DNA")

dna_alignment_Muscle <- msa(dna_sequences_mod, method = "Muscle")
dna_phy_Muscle <- phyDat(as.matrix(dna_alignment_Muscle), type = "DNA")

# writeXStringSet(DNAStringSet(dna_alignment_ClustalO), file.path(sys_dir, "Parvo_sanger_full_seq/aln/Parvovirus_DNA_alignment_ClustalO.fasta"))
# writeXStringSet(DNAStringSet(dna_alignment_Muscle), file.path(sys_dir, "Parvo_sanger_full_seq/aln/Parvovirus_DNA_alignment_Muscle.fasta"))

## we will stick with MUSCLE; medium-sized datasets, highly similar sequences, and protein-coding DNA

## trim MSA_Muscle

# M38246.1_ref_FPV VP2 DNA cord = 2655 - 4585

dna_alignment_Muscle_trimmed <- DNAStringSet(dna_alignment_Muscle)
ref_FPV_name <- grep("M38246.1_ref_FPV", names(dna_alignment_Muscle_trimmed), value = TRUE)
ref_FPV_seq <- dna_alignment_Muscle_trimmed[[ref_FPV_name]]
ref_FPV_chars <- strsplit(as.character(ref_FPV_seq), "")[[1]]
non_gap_positions <- which(ref_FPV_chars != "-")

start_col <- non_gap_positions[2655]
end_col <- non_gap_positions[4585]
# ref_FPV_seq[2655]
# ref_FPV_seq[4585]
# ref_FPV_seq[start_col]
# ref_FPV_seq[end_col]

dna_alignment_Muscle_trimmed <- subseq(dna_alignment_Muscle_trimmed, start = start_col, end = end_col)
width(dna_alignment_Muscle_trimmed)[1]
dna_phy_Muscle_trimmed <- phyDat(as.matrix(dna_alignment_Muscle_trimmed), type = "DNA")

# writeXStringSet(DNAStringSet(dna_alignment_Muscle_trimmed), file.path(sys_dir, "Parvo_sanger_full_seq/aln/Parvovirus_DNA_alignment_Muscle_trimmed.fasta"))

## Alignment

## Protein 

protein_alignment_ClustalO <- msa(protein_sequences_mod, method = "ClustalOmega")
protein_phy_ClustalO <- phyDat(as.matrix(protein_alignment_ClustalO), type = "AA")

protein_alignment_Muscle <- msa(protein_sequences_mod, method = "Muscle")
protein_phy_Muscle <- phyDat(as.matrix(protein_alignment_Muscle), type = "AA")

# writeXStringSet(AAStringSet(protein_alignment_ClustalO), file.path(sys_dir, "Parvo_sanger_full_seq/aln/Parvovirus_Protein_alignment_ClustalO.fasta"))
# writeXStringSet(AAStringSet(protein_alignment_Muscle), file.path(sys_dir, "Parvo_sanger_full_seq/aln/Parvovirus_Protein_alignment_Muscle.fasta"))

# Neighbor Joining trees - DNA & Protein analysis --------------------------------------------------

## Neighbor Joining trees

## DNA

set.seed(32)
dna_NJ_Muscle_trimmed_JC <- NJ(dist.ml(dna_phy_Muscle_trimmed))
dna_NJ_Muscle_trimmed_F81 <- NJ(dist.ml(dna_phy_Muscle_trimmed,model = "F81"))
dna_NJ_Muscle_trimmed_K80 <- NJ(dist.dna(as.DNAbin(dna_phy_Muscle_trimmed), model = "K80"))

model_test <- modelTest(dna_phy_Muscle_trimmed, model = c("JC","F81","K80"), G = FALSE, I = FALSE)
# Model        df  logLik   AIC      BIC
# JC 279 -13280.04 27118.09 28700.01 
# F81 282 -13214.39 26992.78 28591.71 
# K80 280 -13240.38 27040.76 28628.35 

parsimony(dna_NJ_Muscle_trimmed_JC, dna_phy_Muscle_trimmed)
parsimony(dna_NJ_Muscle_trimmed_F81, dna_phy_Muscle_trimmed)
parsimony(dna_NJ_Muscle_trimmed_K80, dna_phy_Muscle_trimmed)

## F81 is best
dna_NJ_Muscle_trimmed <- dna_NJ_Muscle_trimmed_F81

## plot the trees

tip_df <- data.frame(
  label = dna_NJ_Muscle_trimmed$tip.label,
  group = ifelse(grepl("sample", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE) & grepl("feline", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "My_Feline",
                 ifelse(grepl("sample", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE) & grepl("canine", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "My_Canine",
                        ifelse(grepl("sample", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE) & grepl("lion",   dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "My_Lion",
                               ifelse(grepl("NC_001510|NC_001718|NC_029797|NC_038544", dna_NJ_Muscle_trimmed$tip.label), "Outgroup",
                                      ifelse(grepl("_ref_", dna_NJ_Muscle_trimmed$tip.label), "Reference",
                                             ifelse(grepl("vaccine", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "Vaccine",
                                                    ifelse(grepl("feline", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "Regional_Feline",
                                                           ifelse(grepl("canine", dna_NJ_Muscle_trimmed$tip.label, ignore.case = TRUE), "Regional_Canine",
                                                                  "Other"))))))))
  )

## filtering and rooting

## remove zero branch
dna_NJ_Muscle_trimmed_collapse <- di2multi(dna_NJ_Muscle_trimmed, tol = 1e-8)
length(dna_NJ_Muscle_trimmed$tip.label)
length(dna_NJ_Muscle_trimmed_collapse$tip.label)
## none were removed; use manual filtering

outgroup_tips <- tip_df$label[tip_df$group == "Outgroup"]
dna_NJ_Muscle_trimmed_rooted <- root(dna_NJ_Muscle_trimmed, outgroup = outgroup_tips, resolve.root = TRUE)
dna_NJ_Muscle_trimmed_midpoint <- midpoint(dna_NJ_Muscle_trimmed)

dna_NJ_Muscle_trimmed_NOoutgroup <- drop.tip(dna_NJ_Muscle_trimmed, outgroup_tips)
dna_NJ_Muscle_trimmed_NOoutgroup_midpoint <- midpoint(dna_NJ_Muscle_trimmed_NOoutgroup)

parsimony(dna_NJ_Muscle_trimmed, dna_phy_Muscle_trimmed)
parsimony(dna_NJ_Muscle_trimmed_rooted, dna_phy_Muscle_trimmed)
parsimony(dna_NJ_Muscle_trimmed_midpoint, dna_phy_Muscle_trimmed)
parsimony(dna_NJ_Muscle_trimmed_NOoutgroup, dna_phy_Muscle_trimmed)
parsimony(dna_NJ_Muscle_trimmed_NOoutgroup_midpoint, dna_phy_Muscle_trimmed)

ggtree(dna_NJ_Muscle_trimmed_rooted, layout = "rectangular",
       # branch.length = "none",
       open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick", #E63946
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#f66",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - NJ - ", deparse(substitute(dna_NJ_Muscle_trimmed_NOoutgroup))), color = "Group") +
  theme(legend.position = "below",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

## more filtering

plot_data <- ggtree_plot$data
plot_data <- plot_data[plot_data$isTip, ]
tip_order <- plot_data$label[order(plot_data$y, decreasing = TRUE)]

tips_to_keep <- c(
  # "NC_029797.1_Megabat_bufavirus",
  "MZ056887.1_Canine_parvovirus_2_Egypt",
  "OM721656.1_Canine_parvovirus_2b_Turkey",
  "OM100702.1_Canine_parvovirus_2_Egypt",
  "MZ056892.1_Canine_parvovirus_2_Egypt",
  "OM721655.1_Canine_parvovirus_2b_Turkey",
  "OL330979.1_Canine_parvovirus_2_Iran",
  "MW653250.1_Canine_parvovirus_2a_Iran",
  "KX268117.1_Canine_parvovirus_2_Turkey",
  "MW653251.1_Canine_parvovirus_2b_Iran",
  "KX268118.1_Canine_parvovirus_2_Turkey",
  "MW653252.1_Canine_parvovirus_2b_Iran",
  "MW653256.1_Canine_parvovirus_2a_Iran",
  "M24003.1_ref_CPV2A",
  "sample8_canine_79_rev",
  "sample18_canine", "sample68_canine",  "Sample47_canine",
  "OM100699.1_Canine_parvovirus_2_Egypt",
  "OM100696.1_Canine_parvovirus_2_Egypt",
  "Sample77_canine", "sample76_canine", "sample57_canine", "sample51_canine",  "sample72_canine",
  "OM100701.1_Canine_parvovirus_2_Egypt",
  "OM100698.1_Canine_parvovirus_2_Egypt",
  "OM100697.1_Canine_parvovirus_2_Egypt",
  "OR451707.1_Canine_parvovirus_2_Iraq",
  "OR667802.1_Canine_parvovirus_2_Iraq",
  "MW653253.1_Canine_parvovirus_2c_Iran",
  "KX268109.1_Canine_parvovirus_2_Turkey",
  "AY380577.1_ref_CPV2C", "FJ222821.1_ref_CPV2C",
  "KX268114.1_Canine_parvovirus_2_Turkey",
  "M74849.1_ref_CPV2B",
  "vaccineCanine_Zoetis",
  "M38245.1_ref_CPV2",
  "vaccineCanine_MSD",
  "M38246.1_ref_FPV",
  "KP081409.1_Feline_parvovirus_Iran",
  "vaccineFeline_Zoetis",
  "vaccineFeline_MSD",
  "PV521955.1_Feline_parvovirus_Egypt",
  "PP663056.2_Feline_parvovirus_Egypt",
  "sample67_feline",
  "PV521947.1_Feline_parvovirus_Egypt",
  "PP663048.2_Feline_parvovirus_Egypt",
  "Sample69_Feline", "FPV_JO24_Sample_11",
  "sample_lion","Sample25_Feline", "Sample70_Feline", "sample27_feline", "sample38_canine", 
  "Sample_44_feline","sample7_feline", "Sample12_Feline"
  )                   

tips_to_drop <- setdiff(dna_NJ_Muscle_trimmed$tip.label, tips_to_keep)
dna_NJ_Muscle_trimmed_filtered <- drop.tip(dna_NJ_Muscle_trimmed, tips_to_drop)
dna_NJ_Muscle_trimmed_filtered_midpoint <- midpoint(dna_NJ_Muscle_trimmed_filtered)

length(dna_NJ_Muscle_trimmed$tip.label)
length(dna_NJ_Muscle_trimmed_filtered_midpoint$tip.label)

ggtree(dna_NJ_Muscle_trimmed_filtered_midpoint, layout = "fan",
       open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 3, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick", #E63946
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#F4978E",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - NJ - dna_MS_Mustrim_fil_mid", deparse(substitute())), color = "Group") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) 

ggtree(dna_NJ_Muscle_trimmed_filtered_midpoint, layout = "rectangular",
       # branch.length = "none",
       open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick", #E63946
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#f66",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - NJ - ", deparse(substitute(dna_NJ_Muscle_trimmed_filtered_midpoint))), color = "Group") +
  theme(legend.position = "below",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) 

## Protein

## Neighbor Joining trees

set.seed(32)
protein_NJ_Muscle <- NJ(dist.ml(protein_phy_Muscle))
protein_NJ_Muscle_WAG <- NJ(dist.ml(protein_phy_Muscle,model = "WAG"))
protein_NJ_Muscle_JTT <- NJ(dist.ml(protein_phy_Muscle,model = "JTT"))
protein_NJ_Muscle_LG <- NJ(dist.ml(protein_phy_Muscle,model = "LG"))

model_test <- modelTest(protein_phy_Muscle, model = c("WAG","JTT","LG"), G = FALSE, I = FALSE)
# Model        df  logLik   AIC      BIC
# WAG 279 -8078.763 16715.53 18008.96 
# JTT 279 -8128.357 16814.71 18108.14 
# LG 279 -8125.194 16808.39 18101.82 

parsimony(protein_NJ_Muscle_WAG, protein_phy_Muscle)
parsimony(protein_NJ_Muscle_JTT, protein_phy_Muscle)
parsimony(protein_NJ_Muscle_LG,  protein_phy_Muscle)

## WAG is best
protein_NJ_Muscle <- protein_NJ_Muscle_WAG

## plot the trees

## filtering and rooting

## remove zero branch
protein_NJ_Muscle_collapse <- di2multi(protein_NJ_Muscle, tol = 1e-8)
length(protein_NJ_Muscle$tip.label)
length(protein_NJ_Muscle_collapse$tip.label)
## none were removed; use manual filtering

protein_NJ_Muscle_rooted <- root(protein_NJ_Muscle, outgroup = outgroup_tips, resolve.root = TRUE)
protein_NJ_Muscle_midpoint <- midpoint(protein_NJ_Muscle)

protein_NJ_Muscle_NOoutgroup <- drop.tip(protein_NJ_Muscle, outgroup_tips)
protein_NJ_Muscle_NOoutgroup_midpoint <- midpoint(protein_NJ_Muscle_NOoutgroup)

ggtree(protein_NJ_Muscle_rooted, layout = "rectangular",
       # branch.length = "none",
       open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick", #E63946
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#F4978E",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - NJ - ", deparse(substitute(protein_NJ_Muscle_rooted))), color = "Group") +
  theme(legend.position = "below",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +                  
  hexpand(0.5, direction = 1)

ggtree(protein_NJ_Muscle_rooted, layout = "fan",
       open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 3, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick", #E63946
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#F4978E",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - NJ - ", deparse(substitute(protein_NJ_Muscle_rooted))), color = "Group") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) 

## more filtering

protein_NJ_Muscle_filtered <- drop.tip(protein_NJ_Muscle, tips_to_drop)
protein_NJ_Muscle_filtered_midpoint <- midpoint(protein_NJ_Muscle_filtered)

length(protein_NJ_Muscle$tip.label)
length(protein_NJ_Muscle_filtered_midpoint$tip.label)

ggtree(protein_NJ_Muscle_filtered_midpoint, layout = "rectangular",
       open.angle = 10, size = 0.3) %<+% tip_df +
  geom_tiplab(aes(color = group), size = 3, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick", #E63946
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#F4978E",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - NJ - ", deparse(substitute(protein_NJ_Muscle_filtered_midpoint))), color = "Group") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) 

# Maximum Likelihood trees - DNA & Protein analysis -----------------------

## Maximum Likelihood trees from iqtree 

## DNA

dna_ML_Muscle_trimmed_tree <- read.tree(file.path(sys_dir,"Parvo_sanger_full_seq/iqtree_output/Maxll_DNA_aln_Muscle_trimmed.contree"))

tip_df_dna <- data.frame(
  label = dna_ML_Muscle_trimmed_tree$tip.label,
  group = ifelse(grepl("sample", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE) & grepl("feline", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "My_Feline",
                 ifelse(grepl("sample", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE) & grepl("canine", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "My_Canine",
                        ifelse(grepl("sample", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE) & grepl("lion",   dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "My_Lion",
                               ifelse(grepl("NC_001510|NC_001718|NC_029797|NC_038544", dna_ML_Muscle_trimmed_tree$tip.label), "Outgroup",
                                      ifelse(grepl("_ref_", dna_ML_Muscle_trimmed_tree$tip.label), "Reference",
                                             ifelse(grepl("vaccine", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "Vaccine",
                                                    ifelse(grepl("feline", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "Regional_Feline",
                                                           ifelse(grepl("canine", dna_ML_Muscle_trimmed_tree$tip.label, ignore.case = TRUE), "Regional_Canine",
                                                                  "Other"))))))))
)

## rooting

dna_ML_Muscle_trimmed_tree_rooted <- root(dna_ML_Muscle_trimmed_tree, outgroup = outgroup_tips, resolve.root = TRUE)
dna_ML_Muscle_trimmed_tree_midpoint <- midpoint(dna_ML_Muscle_trimmed_tree)

dna_ML_Muscle_trimmed_tree_NOoutgroup <- drop.tip(dna_ML_Muscle_trimmed_tree, outgroup_tips)
dna_ML_Muscle_trimmed_tree_NOoutgroup_midpoint <- midpoint(dna_ML_Muscle_trimmed_tree_NOoutgroup)

parsimony(dna_ML_Muscle_trimmed_tree, dna_phy_Muscle_trimmed)
parsimony(dna_ML_Muscle_trimmed_tree_rooted, dna_phy_Muscle_trimmed)
parsimony(dna_ML_Muscle_trimmed_tree_midpoint, dna_phy_Muscle_trimmed)
parsimony(dna_ML_Muscle_trimmed_tree_NOoutgroup, dna_phy_Muscle_trimmed)
parsimony(dna_ML_Muscle_trimmed_tree_NOoutgroup_midpoint, dna_phy_Muscle_trimmed)

parsimony(dna_ML_Muscle_trimmed_tree_NOoutgroup_midpoint, dna_phy_Muscle_trimmed)

ggtree(dna_ML_Muscle_trimmed_tree_NOoutgroup_midpoint, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick", #E63946
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#f66",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_trimmed_tree_NOoutgroup_midpoint))), color = "Group") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

## more filtering

tips_to_drop <- setdiff(dna_ML_Muscle_trimmed_tree$tip.label, tips_to_keep)
dna_ML_Muscle_trimmed_tree_filtered <- drop.tip(dna_ML_Muscle_trimmed_tree, tips_to_drop)
dna_ML_Muscle_trimmed_tree_filtered_midpoint <- midpoint(dna_ML_Muscle_trimmed_tree_filtered)

length(dna_ML_Muscle_trimmed_tree$tip.label)
length(dna_ML_Muscle_trimmed_tree_filtered_midpoint$tip.label)

ggtree(dna_ML_Muscle_trimmed_tree_filtered_midpoint, layout = "rectangular", size = 0.3) %<+% tip_df_dna +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick", #E63946
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#f66",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  theme_tree2() +
  labs(title = paste0("Parvovirus DNA - ML", deparse(substitute(dna_ML_Muscle_trimmed_tree_filtered_midpoint))), color = "Group") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8))

## Protein

protein_ML_Muscle_tree <- read.tree(file.path(sys_dir,"Parvo_sanger_full_seq/iqtree_output/Maxll_Protein_aln_Muscle.contree"))

tip_df_protein <- data.frame(
  label = protein_ML_Muscle_tree$tip.label,
  group = ifelse(grepl("sample", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE) & grepl("feline", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "My_Feline",
                 ifelse(grepl("sample", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE) & grepl("canine", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "My_Canine",
                        ifelse(grepl("sample", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE) & grepl("lion",   protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "My_Lion",
                               ifelse(grepl("NC_001510|NC_001718|NC_029797|NC_038544", protein_ML_Muscle_tree$tip.label), "Outgroup",
                                      ifelse(grepl("_ref_", protein_ML_Muscle_tree$tip.label), "Reference",
                                             ifelse(grepl("vaccine", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "Vaccine",
                                                    ifelse(grepl("feline", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "Regional_Feline",
                                                           ifelse(grepl("canine", protein_ML_Muscle_tree$tip.label, ignore.case = TRUE), "Regional_Canine",
                                                                  "Other"))))))))
)

## rooting

protein_ML_Muscle_tree_rooted <- root(protein_ML_Muscle_tree, outgroup = outgroup_tips, resolve.root = TRUE)
protein_ML_Muscle_tree_midpoint <- midpoint(protein_ML_Muscle_tree)

protein_ML_Muscle_tree_NOoutgroup <- drop.tip(protein_ML_Muscle_tree, outgroup_tips)
protein_ML_Muscle_tree_NOoutgroup_midpoint <- midpoint(protein_ML_Muscle_tree_NOoutgroup)

parsimony(protein_ML_Muscle_tree, protein_phy_Muscle)
parsimony(protein_ML_Muscle_tree_rooted, protein_phy_Muscle)
parsimony(protein_ML_Muscle_tree_midpoint,  protein_phy_Muscle)
parsimony(protein_ML_Muscle_tree_NOoutgroup, protein_phy_Muscle)
parsimony(protein_ML_Muscle_tree_NOoutgroup_midpoint, protein_phy_Muscle)

parsimony(protein_NJ_Muscle_NOoutgroup_midpoint, protein_phy_Muscle)

ggtree(protein_ML_Muscle_tree, layout = "rectangular", size = 0.3) %<+% tip_df_protein +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick", #E63946
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#f66",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - ML", deparse(substitute(protein_ML_Muscle_tree))), color = "Group") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8)) +
  coord_cartesian(clip = "off") +
  hexpand(0.5, direction = 1)

## more filtering

protein_ML_Muscle_tree_filtered <- drop.tip(protein_ML_Muscle_tree, tips_to_drop)
protein_ML_Muscle_tree_filtered_midpoint <- midpoint(protein_ML_Muscle_tree_filtered)

length(protein_ML_Muscle_tree$tip.label)
length(protein_ML_Muscle_tree_filtered_midpoint$tip.label)

ggtree(protein_ML_Muscle_tree_filtered_midpoint, layout = "rectangular", size = 0.3) %<+% tip_df_protein +
  geom_tiplab(aes(color = group), size = 2, offset = 0.002, align = TRUE, linesize = 0.1) +
  geom_tippoint(aes(color = group), size = 1.2) +
  geom_text2(aes(label = label, subset = !isTip), size = 1.8, hjust = 1.2, vjust = -0.4) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick", #E63946
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#f66",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  theme_tree2() +
  labs(title = paste0("Parvovirus Protein - ML", deparse(substitute(protein_ML_Muscle_tree_filtered_midpoint))), color = "Group") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 8))

# PCA + π - DNA & Protein analysis --------------------------------------------

## PCA

## DNA

dna_alignment_Muscle_trimmed_NOoutgroup <- dna_alignment_Muscle_trimmed[!(names(dna_alignment_Muscle_trimmed) %in% outgroup_tips)]
length(dna_alignment_Muscle_trimmed)
length(dna_alignment_Muscle_trimmed_NOoutgroup)

dna_distance <- dist.dna(as.DNAbin(dna_alignment_Muscle_trimmed_NOoutgroup), model = "F81")
dna_pcoa <- cmdscale(dna_distance, k = 2, eig = TRUE)
dna_var_explained <- round(dna_pcoa$eig / sum(dna_pcoa$eig[dna_pcoa$eig > 0]) * 100, 1)
dna_pcoa_df <- data.frame(
  PC1 = dna_pcoa$points[,1],
  PC2 = dna_pcoa$points[,2],
  Sample = rownames(dna_pcoa$points))
dna_pcoa_df <- merge(dna_pcoa_df, tip_df, by.x = "Sample", by.y = "label", all.x = TRUE)

ggplot(dna_pcoa_df, aes(PC1, PC2, color = group, label = Sample)) +
  geom_point(size = 3) +
  # geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick",
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#F4978E",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  labs(
    title = "DNA PCoA",
    x = paste0("PC1 (", dna_var_explained[1], "%)"),
    y = paste0("PC2 (", dna_var_explained[2], "%)"),
    color = "Group") +
  theme_bw()

## Protein

protein_alignment_Muscle_NOoutgroup <- AAStringSet(protein_alignment_Muscle)
protein_alignment_Muscle_NOoutgroup <- protein_alignment_Muscle_NOoutgroup[!(names(protein_alignment_Muscle_NOoutgroup) %in% outgroup_tips)]
length(protein_alignment_Muscle_NOoutgroup)

protein_phy_Muscle_NOoutgroup <- phyDat(as.matrix(protein_alignment_Muscle_NOoutgroup), type = "AA")
protein_distance <- dist.ml(protein_phy_Muscle_NOoutgroup,model = "WAG")
protein_pcoa <- cmdscale(protein_distance, k = 2, eig = TRUE)
protein_var_explained <- round(protein_pcoa$eig / sum(protein_pcoa$eig[protein_pcoa$eig > 0]) * 100, 1)
protein_pcoa_df <- data.frame(
  PC1 = protein_pcoa$points[,1],
  PC2 = protein_pcoa$points[,2],
  Sample = rownames(protein_pcoa$points))
protein_pcoa_df <- merge(protein_pcoa_df, tip_df, by.x = "Sample", by.y = "label", all.x = TRUE)

ggplot(protein_pcoa_df, aes(PC1, PC2, color = group, label = Sample)) +
  geom_point(size = 3) +
  # geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_manual(values = c(
    "My_Feline"        = "firebrick",
    "My_Canine"        = "#1D4E89",
    "My_Lion"          = "#6A0DAD",
    "Outgroup"         = "black",
    "Reference"        = "forestgreen",
    "Vaccine"          = "#F4A261",
    "Regional_Feline"  = "#F4978E",
    "Regional_Canine"  = "#A9D6E5",
    "Other"            = "gray20")) +
  labs(
    title = "VP2 Protein PCoA",
    x = paste0("PC1 (", protein_var_explained[1], "%)"),
    y = paste0("PC2 (", protein_var_explained[2], "%)"),
    color = "Group") +
  theme_bw()

## Nucleotide diversity (π)

## π is the avg # of nucleotide differences per site between all pairs of sequences per group

dna_alignment_Muscle_trimmed_NOoutgroup

dna_aln_bins <- as.DNAbin(dna_alignment_Muscle_trimmed_NOoutgroup)
group_df <- tip_df[match(names(dna_aln_bins), tip_df$label), ]
group <- group_df$group

dna_aln_bins <- as.matrix(as.DNAbin(dna_alignment_Muscle_trimmed_NOoutgroup))
is.matrix(dna_aln_bins) 
dim(dna_aln_bins)

nuc.div(dna_aln_bins)

nucleotide_diversity <- data.frame(
  Group = unique(group),
  N = as.numeric(table(group)))
nucleotide_diversity$Pi <- sapply(nucleotide_diversity$Group,
                                  function(x) { round(nuc.div(dna_aln_bins[group == x, , drop = FALSE]), 6)})
nucleotide_diversity

# Variant calling - Protein analysis --------------------------------------

## calling variants

## using 
# protein_sequences_mod
# protein_alignment_Muscle

dna_names_mod <- names(protein_sequences_mod)
SamplesVacRef_seq <- dna_names_mod[grepl("sample|vaccine|ref", dna_names_mod, ignore.case = TRUE)]

SamplesVacRef_seq_order <- c(
  "M38246.1_ref_FPV", "M38245.1_ref_CPV2",
  "M24003.1_ref_CPV2A", "M74849.1_ref_CPV2B", "AY380577.1_ref_CPV2C" , "FJ222821.1_ref_CPV2C",
  "sample38_canine", "sample18_canine", "Sample47_canine","sample51_canine", "sample57_canine" , "sample68_canine", 
  "sample72_canine"   ,    "sample76_canine"     ,  "Sample77_canine",   "sample8_canine_79_rev",
  "Sample12_Feline","Sample25_Feline", "sample27_feline","Sample_44_feline","sample67_feline","Sample70_Feline",
  "Sample69_Feline", "sample7_feline", "sample_lion", "FPV_JO24_Sample_11",
  "vaccineCanine_MSD" , "vaccineCanine_Zoetis" , "vaccineFeline_MSD"  ,  "vaccineFeline_Zoetis"
)

protein_alignment_Muscle_ss <- AAStringSet(protein_alignment_Muscle)
protein_phy_Muscle_ss <- phyDat(as.matrix(protein_alignment_Muscle_ss), type = "AA")

protein_sequences_SamplesVacRef <- protein_sequences_mod[names(protein_sequences_mod) %in% SamplesVacRef_seq]
protein_alignment_Muscle_ss_SamplesVacRef <- protein_alignment_Muscle_ss[
  names(protein_alignment_Muscle_ss) %in% SamplesVacRef_seq,]

ref_FPV <- "M38246.1_ref_FPV"     # set reference

names(protein_alignment_Muscle_ss_SamplesVacRef)
protein_alignment_Muscle_ss_SamplesVacRef_mat <- as.matrix(as(protein_alignment_Muscle_ss_SamplesVacRef, "AAStringSet"))
# rows = seqs, cols = alignment positions

## map alignment columns -> ungapped reference residue numbers for publication

ref_seq <- protein_alignment_Muscle_ss_SamplesVacRef_mat[ref_FPV, ]
ref_nongap_cols <- which(ref_seq != "-")
col_to_refpos <- setNames(seq_along(ref_nongap_cols), ref_nongap_cols)

## variant calling vs a reference (aka, FPV)

ProMSA_SamplesVacRef_Var_ref_FPV <- lapply(rownames(protein_alignment_Muscle_ss_SamplesVacRef_mat), function(sn) {
  q <- protein_alignment_Muscle_ss_SamplesVacRef_mat[sn, ]
  diffs <- which(q != ref_seq & q != "-" & ref_seq != "-")
  if (length(diffs) == 0) return(NULL)
  refpos <- col_to_refpos[as.character(diffs)]
  data.frame(seq = sn, pos = refpos, ref_aa = ref_seq[diffs], alt_aa = q[diffs],
             mutation = paste0(ref_seq[diffs], refpos, q[diffs]))
})

ProMSA_SamplesVacRef_Var_ref_FPV_df <- do.call(rbind, ProMSA_SamplesVacRef_Var_ref_FPV)

## loop through all 6 refs and get all divergence in a variant table

ref_names <- grep("_ref_", rownames(protein_alignment_Muscle_ss_SamplesVacRef_mat), value = TRUE)
# SamplesVac_seq <- setdiff(rownames(protein_alignment_Muscle_ss_SamplesVacRef_mat), ref_names)
SamplesVacRef_seq

all_cols <- as.numeric(names(col_to_refpos))

ProMSA_SamplesVacRef_Var_allref_list <- lapply(SamplesVacRef_seq, function(sn) {
  sample_row <- protein_alignment_Muscle_ss_SamplesVacRef_mat[sn, all_cols]
  ref_matrix <- protein_alignment_Muscle_ss_SamplesVacRef_mat[ref_names, all_cols, drop = FALSE]
  
  disagrees <- sapply(seq_along(all_cols), function(i) {
    sv <- sample_row[i]
    rv <- ref_matrix[, i]
    sv != "-" && any(sv != rv & rv != "-")
  })
  
  if (!any(disagrees)) return(NULL)
  
  hit_cols <- all_cols[disagrees]
  hit_refpos <- col_to_refpos[as.character(hit_cols)]
  
  out <- data.frame(seq = sn, pos = hit_refpos, sample_aa = sample_row[disagrees])
  ref_vals <- t(ref_matrix[, disagrees, drop = FALSE])
  colnames(ref_vals) <- paste0("aa_", ref_names)
  cbind(out, ref_vals)
})

ProMSA_SamplesVacRef_Var_allref_df <- do.call(rbind, ProMSA_SamplesVacRef_Var_allref_list)
ProMSA_SamplesVacRef_Var_allref_df <- ProMSA_SamplesVacRef_Var_allref_df[order(ProMSA_SamplesVacRef_Var_allref_df$seq, ProMSA_SamplesVacRef_Var_allref_df$pos), ]
rownames(ProMSA_SamplesVacRef_Var_allref_df) <- NULL

ProMSA_SamplesVacRef_Var_allref_df

## pivoting to a table

all_cols <- as.numeric(names(col_to_refpos))
all_refpos <- col_to_refpos[as.character(all_cols)]

var_cols_final <- sort(unique(ProMSA_SamplesVacRef_Var_allref_df$pos))
var_cols_alncol <- all_cols[match(var_cols_final, all_refpos)]

var_table_mat <- sapply(var_cols_alncol, function(c) {
  fpv_aa <- ref_seq[c]
  col_vals <- protein_alignment_Muscle_ss_SamplesVacRef_mat[, c]
  out <- ifelse(col_vals == "-", "-", ifelse(col_vals == fpv_aa, ".", col_vals))
  out[ref_FPV] <- fpv_aa
  out
})

colnames(var_table_mat) <- var_cols_final
var_table_df <- data.frame(Strain = rownames(protein_alignment_Muscle_ss_SamplesVacRef_mat), var_table_mat, check.names = FALSE)
rownames(var_table_df) <- NULL

var_table_df <- var_table_df[match(SamplesVacRef_seq_order, var_table_df$Strain), ]
rownames(var_table_df) <- NULL
var_table_df

# "." = matches FPV at this position
# "-" = no sequence here (aka, missing data)

## AA mapping 

aa_3to1 <- c(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Glu="E", Gln="Q", Gly="G",
             His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F", Pro="P", Ser="S",
             Thr="T", Trp="W", Tyr="Y", Val="V")
aa_1to3 <- setNames(names(aa_3to1), aa_3to1)

ProMSA_SamplesVacRef_Var_ref_FPV_df$ref_aa_3 <- aa_1to3[ProMSA_SamplesVacRef_Var_ref_FPV_df$ref_aa]
ProMSA_SamplesVacRef_Var_ref_FPV_df$alt_aa_3 <- aa_1to3[ProMSA_SamplesVacRef_Var_ref_FPV_df$alt_aa]
ProMSA_SamplesVacRef_Var_ref_FPV_df$mutation_3 <- paste0(
  ProMSA_SamplesVacRef_Var_ref_FPV_df$ref_aa_3,
  ProMSA_SamplesVacRef_Var_ref_FPV_df$pos,
  ProMSA_SamplesVacRef_Var_ref_FPV_df$alt_aa_3)

var_table_df_3letter <- var_table_df
cols_to_convert <- setdiff(names(var_table_df_3letter), "Strain")

var_table_df_3letter[cols_to_convert] <- lapply(var_table_df_3letter[cols_to_convert], function(col) {
  ifelse(col %in% c(".", "-"), col, aa_1to3[col])
})

# write.xlsx(var_table_df_3letter, file.path(sys_dir, "Parvo_sanger_full_seq/ProMSA_SamplesVacRef_Var_ref_FPV.xlsx"))
# write.xlsx(var_table_df_3letter, file.path(sys_dir, "Parvo_sanger_full_seq/VP2_aa_var_table_3letter.xlsx"))


