# Load libraries
library(readr)
library(stringr)
library(dplyr)
library(openxlsx)
library(tibble)

#################### The following script is scripted for a single marker (Dloop), change the marker accordingly ###############

# Read genbank file
gb_text <- read_file("GenBank_Dloop.gb")

# Split by individual entries
entries <- str_split(gb_text, pattern = fixed("//"))[[1]]

# Main extraction function
extract_features_with_extra_info <- function(entry) {
  # ACCESSION
  accession <- str_match(entry, "ACCESSION\\s+(\\S+)")[,2]
  
  # ORGANISM
  org_match <- str_match(entry, "ORGANISM\\s+([^\n]+)\n((?:\\s{12,}[^\\n]+\\n?)+)")
  species <- org_match[,2]
  taxo_raw <- org_match[,3]
  
  taxo_cols <- rep(NA, 15)
  if (!is.na(taxo_raw)) {
    taxo_clean <- taxo_raw %>%
      str_replace_all("\n", " ") %>%
      str_squish() %>%
      str_split(";\\s*") %>%
      unlist()
    taxo_cols[seq_along(taxo_clean)] <- taxo_clean
  }
  names(taxo_cols) <- paste0("Tax_", seq_along(taxo_cols))
  
  # JOURNALs
  journal_lines <- str_match_all(entry, "JOURNAL\\s+((?:[^\\n]+(?:\\n\\s{12,}[^\\n]+)*))")[[1]][,2]
  journal_lines <- str_replace_all(journal_lines, "\n\\s{12,}", " ")  # joins broken lines
  journal1 <- ifelse(length(journal_lines) >= 1, str_squish(journal_lines[1]), NA)
  journal2 <- ifelse(length(journal_lines) >= 2, str_squish(journal_lines[2]), NA)
  
  # FEATURES
  features_block <- str_match(entry, "(?s)FEATURES(.+?)ORIGIN")[,2]
  if (is.na(features_block)) return(NULL)
  
  lines <- unlist(str_split(features_block, "\n"))
  feature_indexes <- grep("^ {5}\\S", lines)
  feature_indexes <- c(feature_indexes, length(lines) + 1)
  
  results <- list()
  source_qualifiers <- NULL
  
  for (i in seq_along(feature_indexes[-length(feature_indexes)])) {
    start <- feature_indexes[i]
    end <- feature_indexes[i + 1] - 1
    block <- lines[start:end]
    
    header <- str_match(block[1], "^ {5}([\\w\\-]+)\\s+(.+)$")
    type <- header[,2]
    location <- header[,3]
    
    qualifiers <- block[-1]
    qualifiers_text <- paste(trimws(qualifiers), collapse = " ")
    matches <- str_match_all(qualifiers_text, "/([\\w\\-]+)=\\\"(.*?)\\\"")[[1]]
    
    if (nrow(matches) > 0) {
      qual_df <- as_tibble(t(matches[,3]))
      colnames(qual_df) <- matches[,2]
    } else {
      qual_df <- tibble()
    }
    
    if (type == "source") {
      source_qualifiers <- qual_df
    } else {
      result <- tibble(
        Accession = accession,
        Feature_type = type,
        Location = location,
        Organism = species,
        Journal_1 = journal1,
        Journal_2 = journal2
      )
      
      if (ncol(qual_df) > 0) {
        result <- bind_cols(result, qual_df)
      }
      
      result <- bind_cols(result, as_tibble_row(taxo_cols))
      results[[length(results) + 1]] <- result
    }
  }
  
  if (!is.null(source_qualifiers) && length(results) > 0) {
    results <- lapply(results, function(df) {
      bind_cols(df, source_qualifiers)
    })
  }
  
  bind_rows(results)
}

# Applies the function
final_table <- lapply(entries, extract_features_with_extra_info) %>% bind_rows()


# Auxiliary function to remove line breaks in all columns
clean_breaks <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ str_replace_all(., "[\r\n]+", " ")))
}

# Cleans breaks before exporting
final_table_clean <- clean_breaks(final_table)

# Export csv
write.csv(final_table_clean, "GenBank_totalsearch.csv", row.names = FALSE)

# Defines the products we want to keep according to marker 
products_Dloop <- c("D-loop", "control region", "D-loop partial sequence", "mitochondrial control region", "mitochondrial DNA control region")
products_COI <- c("cytochrome c oxidase subunit I", "cytochrome c oxidase subunit 1",  "cytochrome oxidase subunit I",  "cytochrome oxidase subunit 1", "CO1", "COI", "COX1")
products_Cytb <- c("CYTB", "cytochrome b")
products_12S <- c("12S ribosomal RNA", "12S ribosonal RNA", "12S small subunit ribosomal RNA", "small subunit ribosomal RNA", "s-rRNA")
products_16s <- c("16S ribosomal RNA", "large subunit ribosomal RNA", "l-rRNA")

# Filter only the lines containing these products (in 'product', 'gene' or 'note')
table_Dloop <- final_table_clean %>%
  filter(if_any(c(Feature_type, product, gene, note), ~ .x %in% products_Dloop))

# Export csv
write.csv(table_Dloop, "GenBank_Dloop_filtered.csv", row.names = FALSE)

cat("✅ CSV filtered")
