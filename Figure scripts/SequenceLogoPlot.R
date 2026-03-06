# Load libraries
library(stringr)
library(seqinr)
library(Biostrings)
library(ggseqlogo)
library(ggplot2)

# Read aligned sequences
seqs <- read.fasta("MatPep_noPseudo.fasta", seqtype = "AA", as.string = TRUE)
seqs <- read.fasta("Aug4_Drs_all_incl_Katydid.fasta", seqtype = "AA", as.string = TRUE)

seqs_gapped<-seqs

# replace gaps with X to make alphanumeric (ggseqlogo unhappy with "-")
seqs_gapped_X <- gsub("-", "X", seqs_gapped)
aa_with_X_gap <- c(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]], "X")

# Mark's custom colour palette (mimics Geneious' colour palette, which I like)
aa_colors <- c(
  A = "lightgrey", #lightgrey
  C = "gold", #yellow
  D = "salmon", #salmon
  E = "blue", #
  F = "khaki",
  G = "darkolivegreen3",
  H = "orange",
  I = "springgreen3",
  K = "lightblue",
  L = "hotpink",
  M = "magenta",
  N = "#00FFFF",
  P = "wheat3",
  Q = "lightskyblue",
  R = "lightblue1",
  S = "magenta",
  T = "wheat1",
  V = "green",
  W = "pink",
  Y = "palegreen",
  X = "white"  # make X invisible
)

# Create a custom color scheme object so ggseqlogo is happy to use it
### ggseqlogo doesn't take lists etc... apparently
custom_scheme <- make_col_scheme(chars = names(aa_colors), cols = aa_colors)

# indicate which residues asterisks should be drawn over
asterisk_positions <- c(6, 8, 18, 19, 23, 25, 32, 39)  # positions where no 75% majority present

# indicate which residues asterisks should be drawn over
serine_positions <- c(4, 33, 37, 40)  # fixed serines

# plot everything
### various workarounds in the plot text below to customize text and ticks placement etc...
p <- ggseqlogo(seqs_gapped_X, method = "prob", seq_type = "aa",
               namespace = aa_with_X_gap, col_scheme = custom_scheme) +
  # labs(caption = "Asterisks denote residues with <75% agreement, and \"X\" indicates an alignment gap.", size = 32) +  # move title to caption
  labs(x = "Residue position", y = "Frequency") +
  theme(
    axis.title.x = element_text(size = 32, margin = margin(t = 15)),
    axis.title.y = element_text(size = 32, margin = margin(l = 10), angle = 90, vjust = 5, hjust = 0.5),
    axis.text.x  = element_text(size = 32),
    axis.text.y  = element_text(size = 32, margin = margin(r = -150)),
    plot.title   = element_text(size = 32, face = "bold", vjust = 1.5),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks.y = element_blank(),  # hide default y-axis ticks
    plot.caption = element_text(size = 32, face = "bold", hjust = 0, vjust = 1.5, margin = margin(t = 20)),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 40)  # add left padding
  ) +
  annotate(
    "text", x = asterisk_positions, y = 1.01, label = "*", size = 24, vjust = 0
  ) +
  annotate(
    "text", x = serine_positions, y = 1.01, label = "S", size = 14, vjust = -.8
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25)
  ) +
  coord_cartesian(ylim = c(0, 1.3)) +
  geom_segment(data = data.frame(y = seq(0, 1, by = 0.25)),
               aes(x = -1, xend = -.8, y = y, yend = y), # draw lines from 0 to 0.5
               inherit.aes = FALSE,
               size = 1, color = "black")

# Save as PNG
png("Aug4MatPep_logo_noPseudo.png", width = 15000, height = 1500, res = 300)  # Add resolution for clarity
p  # actually draw the plot to the file
dev.off()  # close the device and write the file
pdf("Fig3A_CA.pdf", width = 60, height = 6)
p
dev.off()
