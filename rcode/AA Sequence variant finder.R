# load required packages
library(xml2)
library(stringr)

#read fasta 
read_fasta <- function(filepath) {
  lines <- readLines(filepath)
  ids <- character()
  seqs <- character()
  current_seq <- ""
  current_id <- ""
  
  for (line in lines) {
    line <- str_trim(line)
    if (line == "") next
    
    if (startsWith(line, ">")) {
      if (current_id != "") {
        ids <- c(ids, current_id)
        seqs <- c(seqs, current_seq)
      }
      current_id <- sub(">", "", line)
      current_seq <- ""
    } else {
      current_seq <- paste0(current_seq, line)
    }
  }
  if (current_id != "") {
    ids <- c(ids, current_id)
    seqs <- c(seqs, current_seq)
  }
  return(data.frame(id = ids, sequence = seqs, stringsAsFactors = FALSE))
}


# -------------main script ---------------

pepXML_path <- "C:/Users/chanm/Downloads/01CPTAC_COprospective_W_VU_20150901_01CO006_f01_uncalibrated.pepXML"
fasta_path <- "C:/Users/chanm/Downloads/UP000005640_9606.fasta"

# read data
doc <- read_xml(pepXML_path)
ns <- xml_ns(doc)
fasta_df <- read_fasta(fasta_path)

# extract all PSM's 
spectrum_queries <- xml_find_all(doc, ".//d1:spectrum_query", ns)
total_psms <- length(spectrum_queries)

# subsampling strategy
chunk_size <- max(1, floor(total_psms / 5)) 
sample_size <- 150
sample_indices <- sort(unlist(lapply(0:2, function(i) {
  start_idx <- i * chunk_size + 1
  end_idx <- min((i + 1) * chunk_size, total_psms)
  if(start_idx > end_idx) return(integer(0))
  sample(start_idx:end_idx, min(sample_size, end_idx - start_idx + 1))
})))
print(sample_indices) #optioneeel om te zien welke da ge hebt genomen

sampled_queries <- spectrum_queries[sample_indices]

# process sampld PSM's
psm_data <- lapply(sampled_queries, function(query) {
  hit <- xml_find_first(query, ".//d1:search_hit", ns)
  if (length(hit)) {
    data.frame(
      peptide = xml_attr(hit, "peptide"),
      protein = xml_attr(hit, "protein"),
      stringsAsFactors = FALSE
    )
  }
})

# remove non-hits/null entries 
psm_data <- Filter(Negate(is.null), psm_data)
psm_df <- do.call(rbind, psm_data)

# find variant sequences - CORE LOGIC 
results <- list()
if (nrow(psm_df) > 0) {
  for (i in 1:nrow(psm_df)) {
    peptide <- psm_df$peptide[i]
    protein <- psm_df$protein[i]
    
    protein_seq <- fasta_df$sequence[grepl(protein, fasta_df$id)]
    if (length(protein_seq) > 0) {
      protein_seq <- protein_seq[1]  # Take first match
      peptide_len <- nchar(peptide)
      variants <- list()
      
      # slide window thru sequence
      for (j in 1:(nchar(protein_seq) - peptide_len + 1)) {
        window <- substr(protein_seq, j, j + peptide_len - 1)
        
        # count differences
        diffs <- strsplit(peptide, "")[[1]] != strsplit(window, "")[[1]]
        diff_count <- sum(diffs)
        
        # restrict to only 1 difference 
        if (diff_count == 1) {
          changes <- which(diffs)
          subst <- paste0(
            substr(peptide, changes, changes),
            "→",
            substr(window, changes, changes)
          )
          
          variants[[length(variants) + 1]] <- data.frame(
            original = peptide,
            variant = window,
            position = j + changes - 1,
            substitutions = subst,
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (length(variants) > 0) {
        results[[length(results) + 1]] <- cbind(
          data.frame(protein = protein),
          do.call(rbind, variants)
        )
      }
    }
  }
}

# generate and print output
if (length(results) > 0) {
  final_results <- do.call(rbind, results)
  cat("\n==== SEQUENCE VARIANTS FOUND ====\n")
  print(final_results)
} else {
  cat("\nNo sequence variants found in sampled PSMs.\n")
}