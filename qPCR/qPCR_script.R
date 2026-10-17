library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

ge <- read_excel("Dbus_HGT_qPCR_analysis.xlsx") # comment
str(ge) # check

ge <- ge %>% 
  select(Well,`Well Position`,spp,treatment,gene,replicate,CT)

ge<- ge %>% 
  filter(!is.na(spp), spp != "NA")

dotcolor = c("grey40","blue","firebrick","purple3")

str(ge) # note: tibbles round to 1 decimal, but full data retained
unique(ge$spp)
unique(ge$treatment)
unique(ge$gene)
unique(ge$replicate)

ge <- ge %>%
  mutate(spp = recode(spp,
                      Dbus = "D. busckii",
                      Dmel = "D. melanogaster")) %>%
  mutate(treatment = recode(treatment,
                      UC = "Unchallenged",
                      Ec = "HK E. coli",
                      Ml = "HK M. luteus",
                      EcMl = "HK cocktail"))

head(ge)

ge$replicate<-factor(ge$replicate)

# get average of technical replicates and delta ct
ge_dct <- ge %>%
  group_by(spp, treatment, replicate, gene) %>%
  summarise(
    CT_avg = mean(CT, na.rm = TRUE),
    CT_sd = sd(CT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # next get reference gene CT, here RpL32
  group_by(spp, treatment, replicate) %>%
  mutate(
    ref_ct = mean(CT_avg[gene == "RpL32"], na.rm = TRUE),
    dCT = CT_avg - ref_ct
  ) %>%
  ungroup()

head(ge_dct)

# Standardize BomBc naming
ge_plot <- ge_dct %>%
  filter(gene != "RpL32") %>%       # omit reference gene
  mutate(
    gene_plot = case_when(
      gene %in% c("BomBc", "BomBc3") ~ "BomBc",
      TRUE ~ gene
    )
  )

# set order for levels
ge_plot <- ge_plot %>%
  mutate(
    treatment = factor(treatment, levels = c("Unchallenged", "HK E. coli", "HK M. luteus", "HK cocktail")),
    gene_plot = factor(gene, levels = c("DptA","DptC", "BomBc", "BomBc3", "Dso1", "Drs","Drsl1",
                                   "Drsl2", "Drsl3", "Drsl4", "Drsl5"))
  )

str(ge_plot)
head(ge_plot)

ge_plot <- ge_plot %>%
  group_by(spp, treatment, gene_plot) %>%
  mutate(dCT_sd = sd(dCT, na.rm = TRUE)) %>%
  ungroup()

print(ge_plot[ge_plot$gene=="DptC",],n=40) # seems likely that Dbus samples are reverse order between 1 and 2, and I must reverse datasheet
print(ge_plot[ge_plot$gene=="DptA",],n=40)
print(ge_plot[ge_plot$gene=="Dso1",],n=40)

# for error bars
ge_summary <- ge_plot %>%
  group_by(spp, gene_plot, treatment) %>%
  summarise(
    dCT_mean = mean(dCT, na.rm = TRUE),
    dCT_sd   = sd(dCT, na.rm = TRUE),
    .groups = "drop"
  )

# Plot for Dbus
p_dbus <- ge_plot %>%
  filter(spp == "D. busckii") %>%
  ggplot(aes(x = treatment, y = dCT, ), aes(fill = treatment)) +
  geom_col(data = ge_summary %>% filter(spp == "D. busckii"),
    aes(y = dCT_mean, fill = treatment),
    position = position_dodge(width = 0.9),
    alpha = 0.6) +
  scale_fill_manual(values = dotcolor) +
  geom_errorbar(
    data = ge_summary %>% filter(spp == "D. busckii"),
    aes(y = dCT_mean, ymin = dCT_mean - dCT_sd, ymax = dCT_mean + dCT_sd),
    position = position_dodge(width = 0.9),
    width = 0.2) +
  geom_jitter(aes(color = treatment, shape = replicate),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 1, show.legend = FALSE) +
  scale_color_manual(values = dotcolor) +
  facet_wrap(~ gene_plot, scales = "free_y", ncol = 4) +
  labs(title = bquote(bolditalic("D. busckii")), x = "", y = expression("raw "*Delta*CT)) +
  theme_bw() +
  theme(strip.background = element_blank(),
  strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(color = "black"), plot.title = element_text(size = 24))

# Plot for Dmel
p_dmel <- ge_plot %>%
  filter(spp == "D. melanogaster") %>%
  ggplot(aes(x = treatment, y = dCT, ), aes(fill = treatment)) +
  geom_col(data = ge_summary %>% filter(spp == "D. melanogaster"),
           aes(y = dCT_mean, fill = treatment),
           position = position_dodge(width = 0.9),
           alpha = 0.6) +
  scale_fill_manual(values = dotcolor) +
  geom_errorbar(
    data = ge_summary %>% filter(spp == "D. melanogaster"),
    aes(y = dCT_mean, ymin = dCT_mean - dCT_sd, ymax = dCT_mean + dCT_sd),
    position = position_dodge(width = 0.9),
    width = 0.2) +
  geom_jitter(aes(color = treatment, shape = replicate),
                               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                               size = 1, show.legend = FALSE) +
  scale_color_manual(values = dotcolor) +
  facet_wrap(~ gene_plot, scales = "free_y", ncol = 4) +
  labs(title = bquote(bolditalic("D. melanogaster")), x = "", y = expression("raw "*Delta*CT)) +
  theme_bw() +
  theme(strip.background = element_blank(),
  strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(color = "black"), plot.title = element_text(size = 24))

p_combined <- p_dmel / p_dbus +
  plot_layout(heights = c(1, 2))  # stack vertically, gives two rows of plot height to Dbus and one to Dmel
p_combined

ggsave("Drs_HGT_dCT_qPCR.png", p_combined, dpi = 300, width = 11, height = 9)


# reminder
print(ge_plot, n = 150)
print(ge_summary, n = 50)

## transform to ddCT
ge_ddct <- ge_plot %>%
  # UC rows to use as baseline
  filter(treatment == "Unchallenged") %>%
  select(spp, gene, replicate, dCT_uc = dCT) %>%
  # join back onto full dataset
  right_join(ge_dct, by = c("spp", "gene", "replicate")) %>%
  # calculate ddCT
  mutate(ddCT = dCT - dCT_uc)

# First compute mean dCT of the reference gene within each species
ref_means <- ge_ddct %>%
  filter(gene == "RpL32") %>%      
  filter(!is.na(dCT)) %>%          # just in case
  group_by(spp) %>%
  summarise(ref_mean_dCT = mean(dCT))

# Now join this back and compute ddCT relative to that mean
ge_ddct <- ge_ddct %>%
  left_join(ref_means, by = "spp") %>%
  mutate(
    ddCT = dCT - ref_mean_dCT,
    fold = 2^(-ddCT)
  )
# ge_ddct <- ge_ddct %>%
#   mutate(fold = 2^(-ddCT))

print(ge_ddct, n = 180)


ge_ddct <- ge_ddct %>%
  filter(gene != "RpL32"
         #gene != "Dso1"
         ) %>%       # omit reference gene, and Dso1 for simplicity (Dso1 shows same pattern as BomBc3, but more noise)
  mutate(
    gene_plot = case_when(
      gene %in% c("BomBc", "BomBc3") ~ "BomBc",
      TRUE ~ gene
    )
  )

# set order for levels
ge_ddct <- ge_ddct %>%
  mutate(
    treatment = factor(treatment, levels = c("Unchallenged", "HK E. coli", "HK M. luteus", "HK cocktail")),
    gene_plot = factor(gene, levels = c("DptA","DptC", "BomBc", "BomBc3", "Dso1", "Drs","Drsl1",
                                        "Drsl2", "Drsl3", "Drsl4", "Drsl5"))
  )

ge_ddct <- ge_ddct %>%
  group_by(spp, treatment, gene_plot) %>%
  mutate(fold_sd = sd(fold, na.rm = TRUE)) %>%
  ungroup()

# for error bars
ge_summary <- ge_ddct %>%
  group_by(spp, gene_plot, treatment) %>%
  summarise(
    fold_mean = mean(fold, na.rm = TRUE),
    fold_sd   = sd(fold, na.rm = TRUE),
    .groups = "drop"
  )

print(ge_ddct, n = 150)
print(ge_summary, n = 50)

print(ge_ddct[ge_ddct$gene=="DptC",], n = 150)
print(ge_summary[ge_summary$gene_plot=="DptC",], n = 50)
head(ge_ddct)
head(ge_summary)

### normalise so mean(Unchallenged = 1) per gene
# Compute mean fold of Unchallenged for each species+gene
ref_fold_means <- ge_ddct %>%
  filter(treatment == "Unchallenged") %>%
  group_by(spp, gene_plot) %>%
  summarise(ref_fold_mean = mean(fold), .groups = "drop")

# Join and compute normalized values
ge_ddct <- ge_ddct %>%
  left_join(ref_fold_means, by = c("spp", "gene_plot")) %>%
  mutate(
    fold_normalised = fold / ref_fold_mean,
    fold_sd_normalised = fold_sd / ref_fold_mean
  )

# Compute reference mean fold from ge_summary (Unchallenged only)
ref_summary_means <- ge_summary %>%
  filter(treatment == "Unchallenged") %>%
  group_by(spp, gene_plot) %>%
  summarise(ref_fold_mean = mean(fold_mean), .groups = "drop")

# Join and normalize all fold_means + fold_sd
ge_summary <- ge_summary %>%
  left_join(ref_summary_means, by = c("spp", "gene_plot")) %>%
  mutate(
    fold_mean_normalised = fold_mean / ref_fold_mean,
    fold_sd_normalised   = fold_sd / ref_fold_mean
  )


# Plot for Dbus
p_dbus <- ge_ddct %>%
  filter(spp == "D. busckii") %>%
  ggplot(aes(x = treatment, y = fold_normalised, ), aes(fill = treatment)) +
  geom_col(data = ge_summary %>% filter(spp == "D. busckii"),
           aes(y = fold_mean_normalised, fill = treatment),
           position = position_dodge(width = 0.9),
           alpha = 0.6) +
  scale_fill_manual(values = dotcolor) +
  geom_errorbar(
    data = ge_summary %>% filter(spp == "D. busckii"),
    aes(y = fold_mean_normalised, ymin = fold_mean_normalised - fold_sd_normalised, ymax = fold_mean_normalised + fold_sd_normalised),
    position = position_dodge(width = 0.9),
    width = 0.2) +
  geom_jitter(aes(color = treatment, shape = replicate),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 1, show.legend = FALSE) +
  scale_color_manual(values = dotcolor) +
  facet_wrap(~ gene_plot, scales = "free_y", ncol = 3) +
  labs(title = bquote(bolditalic("D. busckii")), x = "", y = expression(`2-`*Delta*Delta*CT)) +
  theme_bw() +
  theme(strip.background = element_blank(),
  strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(color = "black"), plot.title = element_text(size = 24))


# Plot for Dmel
p_dmel <- ge_ddct %>%
  filter(spp == "D. melanogaster") %>%
  ggplot(aes(x = treatment, y = fold_normalised, ), aes(fill = treatment)) +
  geom_col(data = ge_summary %>% filter(spp == "D. melanogaster"),
           aes(y = fold_mean_normalised, fill = treatment),
           position = position_dodge(width = 0.9),
           alpha = 0.6) +
  scale_fill_manual(values = dotcolor) +
  geom_errorbar(
    data = ge_summary %>% filter(spp == "D. melanogaster"),
    aes(y = fold_mean_normalised, ymin = fold_mean_normalised - fold_sd_normalised, ymax = fold_mean_normalised + fold_sd_normalised),
    position = position_dodge(width = 0.9),
    width = 0.2) +
  geom_jitter(aes(color = treatment, shape = replicate),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 1, show.legend = FALSE) +
  scale_color_manual(values = dotcolor) +
  facet_wrap(~ gene_plot, scales = "free_y", ncol = 3) +
  labs(title = bquote(bolditalic("D. melanogaster")), x = "", y = expression(`2-`*Delta*Delta*CT)) +
  theme_bw() +
  theme(strip.background = element_blank(),
  strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(color = "black"), plot.title = element_text(size = 24))

p_combined2 <-  p_dmel / p_dbus +
  plot_layout(heights = c(1, 3))  # stack vertically, gives two rows of plot height to Dbus and one to Dmel
p_combined2

ggsave("Drs_HGT_ddCT_qPCR.png", p_combined2, dpi = 300, width = 11, height = 9)


#####################
#### stats
#####################

str(ge_ddct)
ge_ddct <- ge_ddct %>%
  mutate(treatment = relevel(treatment, ref = "Unchallenged"))

# load additional libraries
library(purrr)
library(broom)

tukey_results <- ge_ddct %>%
  group_by(spp, gene_plot) %>%
  group_modify(~ {
    # Fit model
    fit <- aov(ddCT ~ treatment, data = .x)
    
    # Run Tukey HSD
    tk <- TukeyHSD(fit)
    
    # Extract treatment table
    tk_df <- tk$treatment %>%
      as.data.frame() %>%
      mutate(contrast = rownames(.)) %>%
      select(contrast, diff, lwr, upr, p_adj = `p adj`)
    
    # Keep only comparisons vs Unchallenged
    tk_df <- tk_df %>%
      filter(grepl("Unchallenged", contrast))
    
    tk_df
  }) %>%
  ungroup()

print(tukey_results, n = 30)

##########################
## plot with stats in
##########################
# generate statistical groups
tukey_results <- tukey_results %>%
  mutate(sig = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01  ~ "**",
    p_adj < 0.05  ~ "*",
    # p_adj < 0.1 ~ "°",
    p_adj < .10 ~ "˙",
    TRUE ~ ""
  ))

# extract comparisons vs Unchallenged
tukey_sig <- tukey_results %>%
  separate(contrast, into = c("treatment", "ref"), sep = "-") %>%
  filter(ref == "Unchallenged", sig != "") %>%
  select(spp, gene_plot, treatment, diff, lwr, upr, p_adj, sig)


######
# let's plot
######

# build stat annotations df
sig_positions <- ge_summary %>%
  select(spp, gene_plot, treatment, fold_mean_normalised) %>%
  left_join(tukey_sig, by = c("spp", "gene_plot", "treatment")) %>%
  group_by(spp, gene_plot) %>%
  mutate(y_pos = max(fold_mean_normalised) * 1.25 +
           row_number() * max(fold_mean_normalised) * 0.15)

# Filter only significant comparisons for each species
sig_positions_db <- sig_positions %>% filter(spp == "D. busckii", !is.na(sig))
sig_positions_dm <- sig_positions %>% filter(spp == "D. melanogaster", !is.na(sig))

# Plot D. busckii
p_dbus <- ge_ddct %>%
  filter(spp == "D. busckii") %>%
  ggplot(aes(x = treatment, y = fold_normalised, fill = treatment)) +
  geom_col(data = ge_summary %>% filter(spp == "D. busckii"),
           aes(y = fold_mean_normalised, fill = treatment),
           position = position_dodge(width = 0.9), alpha = 0.6) +
  geom_errorbar(data = ge_summary %>% filter(spp == "D. busckii"),
                aes(y = fold_mean_normalised,
                    ymin = fold_mean_normalised - fold_sd_normalised,
                    ymax = fold_mean_normalised + fold_sd_normalised),
                position = position_dodge(width = 0.9), width = 0.2) +
  geom_jitter(aes(color = treatment, shape = replicate),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 2, show.legend = FALSE) +
  scale_fill_manual(values = dotcolor) +
  scale_color_manual(values = dotcolor) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  facet_wrap(~ gene_plot, scales = "free_y", ncol = 4) +
  labs(title = bquote(bolditalic("D. busckii")), x = "", y = "Fold-change") + # expression(`2-`*Delta*Delta*CT)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 24))

# Plot D. melanogaster
p_dmel <- ge_ddct %>%
  filter(spp == "D. melanogaster") %>%
  ggplot(aes(x = treatment, y = fold_normalised, fill = treatment)) +
  geom_col(data = ge_summary %>% filter(spp == "D. melanogaster"),
           aes(y = fold_mean_normalised, fill = treatment),
           position = position_dodge(width = 0.9), alpha = 0.6) +
  geom_errorbar(data = ge_summary %>% filter(spp == "D. melanogaster"),
                aes(y = fold_mean_normalised,
                    ymin = fold_mean_normalised - fold_sd_normalised,
                    ymax = fold_mean_normalised + fold_sd_normalised),
                position = position_dodge(width = 0.9), width = 0.2) +
  geom_jitter(aes(color = treatment, shape = replicate),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 2, show.legend = FALSE) +
  scale_fill_manual(values = dotcolor) +
  scale_color_manual(values = dotcolor) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  facet_wrap(~ gene_plot, scales = "free_y", ncol = 4) +
  labs(title = bquote(bolditalic("D. melanogaster")), x = "", y = "Fold-change") + # expression(`2-`*Delta*Delta*CT)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 24))

### Increase the y-axis range so the asterisks later don't get cut by the margins
# Compute max y per facet for D. busckii
max_y_db <- ge_summary %>%
  filter(spp == "D. busckii") %>%
  group_by(gene_plot) %>%
  summarise(max_y = max(fold_mean_normalised, na.rm = TRUE)) %>%
  mutate(max_y = max_y * 1.1)

# Compute max y per facet for D. melanogaster
max_y_dm <- ge_summary %>%
  filter(spp == "D. melanogaster") %>%
  group_by(gene_plot) %>%
  summarise(max_y = max(fold_mean_normalised, na.rm = TRUE)) %>%
  mutate(max_y = max_y * 1.1)

# Add invisible points to extend y-axis for D. busckii
p_dbus <- p_dbus +
  geom_blank(data = max_y_db %>%
               mutate(treatment = factor(levels(ge_ddct$treatment)[1],
                                         levels = levels(ge_ddct$treatment))),
             aes(x = treatment, y = max_y))

# Add invisible points to extend y-axis for D. melanogaster
p_dmel <- p_dmel +
  geom_blank(data = max_y_dm %>%
               mutate(treatment = factor(levels(ge_ddct$treatment)[1],
                                         levels = levels(ge_ddct$treatment))),
             aes(x = treatment, y = max_y))

# Fix the asterisks
# Convert treatment factors to numeric positions for plotting annotations
sig_positions_db <- sig_positions_db %>%
  mutate(
    x_num = as.numeric(factor(treatment, levels = levels(ge_ddct$treatment))),
    xend_num = as.numeric(factor("Unchallenged", levels = levels(ge_ddct$treatment))),
    x_mid = (x_num + xend_num) / 2
  )

sig_positions_dm <- sig_positions_dm %>%
  mutate(
    x_num = as.numeric(factor(treatment, levels = levels(ge_ddct$treatment))),
    xend_num = as.numeric(factor("Unchallenged", levels = levels(ge_ddct$treatment))),
    x_mid = (x_num + xend_num) / 2
  )

# Add segments and asterisks to ggplot
p_dbus <- p_dbus +
  geom_segment(data = sig_positions_db,
               aes(x = x_num, xend = xend_num, y = y_pos, yend = y_pos),
               inherit.aes = FALSE, linewidth = 0.15) +
  geom_text(data = sig_positions_db,
            aes(x = x_mid, y = y_pos, label = sig),
            vjust = .5, size = 5, inherit.aes = FALSE)

p_dmel <- p_dmel +
  geom_segment(data = sig_positions_dm,
               aes(x = x_num, xend = xend_num, y = y_pos, yend = y_pos),
               inherit.aes = FALSE, linewidth = 0.15) +
  geom_text(data = sig_positions_dm,
            aes(x = x_mid, y = y_pos, label = sig),
            vjust = .5, size = 5, inherit.aes = FALSE)

### if using ncol = 3
# Combine plots vertically
p_combined3 <- p_dmel / p_dbus + plot_layout(heights = c(1, 2))
p_combined3

ggsave("Drs_HGT_ddCT_qPCR_stats.png", p_combined3, dpi = 300, width = 11, height = 11)




########
# main text plot
########

# generate statistical groups
tukey_results <- tukey_results %>%
  mutate(sig = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01  ~ "**",
    p_adj < 0.05  ~ "*",
    TRUE ~ ""
  ))

# extract comparisons vs Unchallenged
tukey_sig <- tukey_results %>%
  separate(contrast, into = c("treatment", "ref"), sep = "-") %>%
  filter(ref == "Unchallenged", sig != "") %>%
  select(spp, gene_plot, treatment, diff, lwr, upr, p_adj, sig)

# Remove unnecessary genes
ge_ddct <- ge_ddct %>%
  filter(gene != "RpL32",
         gene != "Dso1",
         gene != "Drsl1",
         gene != "Drsl2",
         gene != "Drsl3",
         gene != "Drsl5")
# Standardize BomBc naming
ge_summary <- ge_summary %>%
  filter(gene_plot != "Dso1",
         gene_plot != "Drsl1",
         gene_plot != "Drsl2",
         gene_plot != "Drsl3",
         gene_plot != "Drsl5")

######
# let's plot
######

# build stat annotations df
sig_positions <- ge_summary %>%
  select(spp, gene_plot, treatment, fold_mean_normalised) %>%
  left_join(tukey_sig, by = c("spp", "gene_plot", "treatment")) %>%
  group_by(spp, gene_plot) %>%
  mutate(y_pos = max(fold_mean_normalised) * 1.25 +
           row_number() * max(fold_mean_normalised) * 0.15)

# Filter only significant comparisons for each species
sig_positions_db <- sig_positions %>% filter(spp == "D. busckii", !is.na(sig))
sig_positions_dm <- sig_positions %>% filter(spp == "D. melanogaster", !is.na(sig))

# Plot D. busckii
p_dbus <- ge_ddct %>%
  filter(spp == "D. busckii") %>%
  ggplot(aes(x = treatment, y = fold_normalised, fill = treatment)) +
  geom_col(data = ge_summary %>% filter(spp == "D. busckii"),
           aes(y = fold_mean_normalised, fill = treatment),
           position = position_dodge(width = 0.9), alpha = 0.6) +
  geom_errorbar(data = ge_summary %>% filter(spp == "D. busckii"),
                aes(y = fold_mean_normalised,
                    ymin = fold_mean_normalised - fold_sd_normalised,
                    ymax = fold_mean_normalised + fold_sd_normalised),
                position = position_dodge(width = 0.9), width = 0.2) +
  geom_jitter(aes(color = treatment, shape = replicate),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 2, show.legend = FALSE) +
  scale_fill_manual(values = dotcolor) +
  scale_color_manual(values = dotcolor) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  facet_wrap(~ gene_plot, scales = "free_y", ncol = 4) +
  labs(title = bquote(bolditalic("D. busckii")), x = "", y = "Fold-change") + # expression(`2-`*Delta*Delta*CT)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16))

# Plot D. melanogaster
p_dmel <- ge_ddct %>%
  filter(spp == "D. melanogaster") %>%
  ggplot(aes(x = treatment, y = fold_normalised, fill = treatment)) +
  geom_col(data = ge_summary %>% filter(spp == "D. melanogaster"),
           aes(y = fold_mean_normalised, fill = treatment),
           position = position_dodge(width = 0.9), alpha = 0.6) +
  geom_errorbar(data = ge_summary %>% filter(spp == "D. melanogaster"),
                aes(y = fold_mean_normalised,
                    ymin = fold_mean_normalised - fold_sd_normalised,
                    ymax = fold_mean_normalised + fold_sd_normalised),
                position = position_dodge(width = 0.9), width = 0.2) +
  geom_jitter(aes(color = treatment, shape = replicate),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 2, show.legend = FALSE) +
  scale_fill_manual(values = dotcolor) +
  scale_color_manual(values = dotcolor) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  facet_wrap(~ gene_plot, scales = "free_y", ncol = 4) +
  labs(title = bquote(bolditalic("D. melanogaster")), x = "", y = "Fold-change") + # expression(`2-`*Delta*Delta*CT)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16))

### Increase the y-axis range so the asterisks later don't get cut by the margins
# Compute max y per facet for D. busckii
max_y_db <- ge_summary %>%
  filter(spp == "D. busckii") %>%
  group_by(gene_plot) %>%
  summarise(max_y = max(fold_mean_normalised, na.rm = TRUE)) %>%
  mutate(max_y = max_y * 1.1)

# Compute max y per facet for D. melanogaster
max_y_dm <- ge_summary %>%
  filter(spp == "D. melanogaster") %>%
  group_by(gene_plot) %>%
  summarise(max_y = max(fold_mean_normalised, na.rm = TRUE)) %>%
  mutate(max_y = max_y * 1.1)

# Add invisible points to extend y-axis for D. busckii
p_dbus <- p_dbus +
  geom_blank(data = max_y_db %>%
               mutate(treatment = factor(levels(ge_ddct$treatment)[1],
                                         levels = levels(ge_ddct$treatment))),
             aes(x = treatment, y = max_y))

# Add invisible points to extend y-axis for D. melanogaster
p_dmel <- p_dmel +
  geom_blank(data = max_y_dm %>%
               mutate(treatment = factor(levels(ge_ddct$treatment)[1],
                                         levels = levels(ge_ddct$treatment))),
             aes(x = treatment, y = max_y))

# Fix the asterisks
# Convert treatment factors to numeric positions for plotting annotations
sig_positions_db <- sig_positions_db %>%
  mutate(
    x_num = as.numeric(factor(treatment, levels = levels(ge_ddct$treatment))),
    xend_num = as.numeric(factor("Unchallenged", levels = levels(ge_ddct$treatment))),
    x_mid = (x_num + xend_num) / 2
  )

sig_positions_dm <- sig_positions_dm %>%
  mutate(
    x_num = as.numeric(factor(treatment, levels = levels(ge_ddct$treatment))),
    xend_num = as.numeric(factor("Unchallenged", levels = levels(ge_ddct$treatment))),
    x_mid = (x_num + xend_num) / 2
  )

# Add segments and asterisks to ggplot
p_dbus <- p_dbus +
  geom_segment(data = sig_positions_db,
               aes(x = x_num, xend = xend_num, y = y_pos, yend = y_pos),
               inherit.aes = FALSE, linewidth = 0.15) +
  geom_text(data = sig_positions_db,
            aes(x = x_mid, y = y_pos, label = sig),
            vjust = .5, size = 4, inherit.aes = FALSE)

p_dmel <- p_dmel +
  geom_segment(data = sig_positions_dm,
               aes(x = x_num, xend = xend_num, y = y_pos, yend = y_pos),
               inherit.aes = FALSE, linewidth = 0.15) +
  geom_text(data = sig_positions_dm,
            aes(x = x_mid, y = y_pos, label = sig),
            vjust = .5, size = 4, inherit.aes = FALSE)

### if using ncol = 3
# Combine plots vertically
p_combined4 <- p_dmel / p_dbus + plot_layout(heights = c(1, 1))
p_combined4

ggsave("Drs_HGT_ddCT_qPCR_main.png", p_combined4, dpi = 300, width = 4.5, height = 5.5)
