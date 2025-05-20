PTM_matcher <- function(delta_masses, tolerance = 0.025) {
  ptm_list <- c(
    Methylation = 14.01565,
    Dimethylation = 28.0313,
    Trimethylation = 42.04695,
    Acetylation = 42.0106,
    Propionylation = 56.0262,
    Butyrylation = 70.0419,
    Succinylation = 100.0160,
    Malonylation = 86.0004,
    Formylation = 27.9949,
    Oxidation = 15.9949,
    Dioxidation = 31.9898,
    Nitrosylation = 29.9979,
    Nitration = 44.9851,
    Hydroxylation = 15.9949,
    Phosphorylation = 79.9663,
    Sulfation = 79.9568,
    DisulfideLoss = -2.0157,
    Carbamidomethyl = 57.0215,
    Carboxymethylation = 58.0055,
    Carbamylation = 43.0058,
    Iodoacetamide = 57.0215,
    Hexose = 162.0528,
    HexNAc = 203.0794,
    Deoxyhexose = 146.0579,
    SialicAcid = 291.0954,
    O_Glycosylation = 203.08,
    N_GlycanCore = 1202.42,
    Glycation = 162.0528,
    Palmitoylation = 238.2297,
    Myristoylation = 210.1984,
    Farnesylation = 204.1878,
    Geranylgeranylation = 272.2504,
    Ubiquitination = 114.0429,
    NEDDylation = 114.0429,
    SUMOylation = 383.2281
  )
  
  ptm_combos <- unlist(lapply(1:3, function(n) combn(names(ptm_list), n, simplify = FALSE)), recursive = FALSE)
  combo_masses <- sapply(ptm_combos, function(combo) sum(ptm_list[combo]))
  names(combo_masses) <- sapply(ptm_combos, paste, collapse = " + ")
  
  best_match_for_mass <- function(delta_mass) {
    trials <- 0
    current_tol <- tolerance
    repeat {
      trials <- trials + 1
      diffs <- abs(combo_masses - delta_mass)
      best_idx <- which(diffs <= current_tol)
      if (length(best_idx) > 0 || trials > 5) break
      current_tol <- current_tol * 2
    }
    if (length(best_idx) == 0) {
      return(data.frame(
        delta_mass = delta_mass,
        match = NA,
        mass = NA,
        delta = NA,
        used_tolerance = current_tol,
        trials = trials
      ))
    }
    best <- best_idx[which.min(diffs[best_idx])]
    data.frame(
      delta_mass = delta_mass,
      match = names(combo_masses)[best],
      mass = combo_masses[best],
      delta = combo_masses[best] - delta_mass,
      used_tolerance = current_tol,
      trials = trials
    )
  }
  
  do.call(rbind, lapply(delta_masses, best_match_for_mass))
}

# XML integration and delta mass computation
library(xml2)
library(dplyr)
library(ggplot2)

file_path <- "C:/Users/chanm/Downloads/01CPTAC_COprospective_W_VU_20150901_01CO006_f01_uncalibrated.pepXML"
doc <- read_xml(file_path)
ns <- xml_ns(doc)

spectrum_queries <- xml_find_all(doc, ".//d1:spectrum_query", ns)
total_psms <- length(spectrum_queries)

chunk_size <- floor(total_psms / 5)
sample_size <- 25

sample_indices <- sort(unlist(lapply(0:4, function(i) {
  start_idx <- i * chunk_size + 1
  end_idx <- min((i + 1) * chunk_size, total_psms)
  sample(start_idx:end_idx, sample_size)
})))

print(sample_indices)

psm_list <- list()
delta_masses <- c()

for (i in seq_along(spectrum_queries)) {
  if (!(i %in% sample_indices)) next
  
  spectrum <- spectrum_queries[[i]]
  precursor_mass <- as.numeric(xml_attr(spectrum, "precursor_neutral_mass"))
  hits <- xml_find_all(spectrum, ".//d1:search_hit", ns)
  
  for (hit in hits) {
    pep_mass <- as.numeric(xml_attr(hit, "calc_neutral_pep_mass"))
    delta <- precursor_mass - pep_mass
    delta_masses <- c(delta_masses, delta)
    
    peptide <- xml_attr(hit, "peptide")
    protein <- xml_attr(hit, "protein")
    spectrum_id <- xml_attr(spectrum, "spectrum")
    
    psm <- data.frame(
      spectrum = spectrum_id,
      peptide = peptide,
      protein = protein,
      precursor_mass = precursor_mass,
      pep_mass = pep_mass,
      delta_mass = delta,
      stringsAsFactors = FALSE
    )
    
    psm_list[[length(psm_list) + 1]] <- psm
  }
}

psms_df <- bind_rows(psm_list)

ptm_matches <- PTM_matcher(psms_df$precursor_mass - psms_df$pep_mass, tolerance = 0.025)
final_results <- bind_cols(psms_df, ptm_matches[, c("match", "mass", "delta", "used_tolerance", "trials")])

# sort by used_tolerance (ascending)
final_results <- final_results[order(final_results$used_tolerance), ]

#  column names
colnames(final_results) <- c("spectrum", "peptide", "protein", "precursor_mass", "pepfrag_mass", "delta_mass", "match", "PTM_combo_mass", "delta_fit", "used_tolerance", "#trials")

# output columns
print(final_results[, c("peptide", "protein", "precursor_mass", "pepfrag_mass", "delta_mass", "PTM_combo_mass", "delta_fit", "match", "used_tolerance", "#trials")], row.names = FALSE)

cat("\nTotal sampled PSMs:", nrow(final_results), "out of", total_psms, "\n")



############################################################ PLOTTING #######################################################

par(mfrow = c(1, 1))

# plot 1: PTM frequency 
all_ptms <- unlist(strsplit(na.omit(final_results$match), " \\+ "))
all_ptms <- trimws(all_ptms)
all_ptms <- tolower(all_ptms)
all_ptms <- tools::toTitleCase(all_ptms)

ptm_freq <- as.data.frame(table(all_ptms))
colnames(ptm_freq) <- c("PTM", "Count")
ptm_freq <- ptm_freq[order(-ptm_freq$Count), ]

barplot(ptm_freq$Count, names.arg = ptm_freq$PTM, las = 2, col = "lightblue",
        main = "PTM Frequencies", ylab = "Count")

# plot 2: histogram of PTM combo masses 
hist(final_results$PTM_combo_mass, breaks = seq(min(final_results$PTM_combo_mass, na.rm = TRUE), 
                                                max(final_results$PTM_combo_mass, na.rm = TRUE) + 15, by = 15),
     col = "lightgreen", main = "PTM Combo Mass frequency",
     xlab = "Mass Range", ylab = "Count")

# plot 3: tolerance usage frequency
tol_freq <- table(final_results$used_tolerance)
barplot(tol_freq, main = "Tolerance Usage Frequency",
        xlab = "Tolerance", ylab = "Count", col = "salmon")

hist(final_results$delta_mass[final_results$match == "DisulfideLoss"],
     breaks = 30, main = "Delta mass deviation ~ Disulfide loss artefact ",
     xlab = "Delta Mass (Da)", col = "orange")
abline(v = -2.0157, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = "True Disulfide loss -2.0157 Da", col = "blue", lwd = 2, lty = 2)



################_Co-occurrence analysis_##########

library(tidyr)
library(stats)
library(RColorBrewer)

# prepare PTM matches
ptm_split <- strsplit(as.character(na.omit(final_results$match)), " \\+ ")
unique_ptms <- sort(unique(unlist(ptm_split)))

# create binary PTM matrix
ptm_matrix <- sapply(unique_ptms, function(ptm) {
  sapply(ptm_split, function(psm_ptms) ptm %in% psm_ptms)
})

# compute jaccard-distance matrix
ptm_dist <- dist(t(ptm_matrix), method = "binary")

# perform classical MDS
mds_coords <- cmdscale(ptm_dist, k = 2)

# assign colors to each PTM
ptm_names <- colnames(ptm_matrix)
num_ptms <- length(ptm_names)
colors <- brewer.pal(min(num_ptms, 8), "Set1")
if (num_ptms > 8) {
  colors <- colorRampPalette(colors)(num_ptms)
}

#K-means clustering on MDS coordinates
set.seed(666)  
k <- 6  
km <- kmeans(mds_coords, centers = k)
kmeans_clusters <- km$cluster








# assign cluster-based colors
kmeans_colors <- brewer.pal(k, "Dark2")[kmeans_clusters]


par(mfrow = c(1, 1), mar = c(5, 5, 4, 10), xpd = TRUE)

# first plott:text labels
plot(mds_coords, type = "n", main = "PTM Co-occurrence",
     xlab = "Dimension 1", ylab = "Dimension 2")
points(mds_coords, pch = 19, col = colors, cex = 2)
text(mds_coords, labels = ptm_names, pos = 3, cex = 1.2)

# second plot: colored dots 
plot(mds_coords, type = "n", main = "PTM Co-occurrence",
     xlab = "Dimension 1", ylab = "Dimension 2")
points(mds_coords, pch = 19, col = colors, cex = 2.5)
legend("topright", inset = c(-0.3, 0), legend = ptm_names, col = colors, pch = 19, cex = 1.8)

#third plot:K-means clusters
plot(mds_coords, type = "n", main = "MDS K-means Clusters",
     xlab = "Dimension 1", ylab = "Dimension 2")
points(mds_coords, pch = 19, col = kmeans_colors, cex = 2.5)
legend("topright", inset = c(-0.3, 0), legend = paste("Cluster", 1:k), col = brewer.pal(k, "Dark2"), pch = 19, cex = 1.5)


# Co-occurrence statistical significance (Fisher's Exact Test)
fisher_results <- list()

for (i in 1:(ncol(ptm_matrix) - 1)) {
  for (j in (i + 1):ncol(ptm_matrix)) {
    a <- ptm_matrix[, i] & ptm_matrix[, j]
    b <- ptm_matrix[, i] & !ptm_matrix[, j]
    c <- !ptm_matrix[, i] & ptm_matrix[, j]
    d <- !ptm_matrix[, i] & !ptm_matrix[, j]
    contingency <- matrix(c(sum(a), sum(b), sum(c), sum(d)), nrow = 2)
    test <- fisher.test(contingency)
    pair <- paste(colnames(ptm_matrix)[i], colnames(ptm_matrix)[j], sep = " & ")
    fisher_results[[pair]] <- c(p_value = test$p.value, odds_ratio = test$estimate)
  }
}

fisher_df <- as.data.frame(do.call(rbind, fisher_results))
fisher_df$pair <- rownames(fisher_df)
fisher_df <- fisher_df[order(fisher_df$p_value), ]

cat("\nTop significant PTM co-occurrences:\n")
print(head(fisher_df, 50))

