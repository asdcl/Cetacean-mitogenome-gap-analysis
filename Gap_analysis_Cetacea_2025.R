# Load libraries
library(readxl)    # read Excel
library(writexl)   # export to Excel
library(bold)      # BOLD API
library(rentrez)   # NCBI API
library(dplyr)     # data wrangling
library(stringr)   # string operations
library(purrr)     # functional helpers
library(tibble)    # tibble data frame structure

# Species list retrieved from WoRMS: Excel must have at least scientificname, valid_name, AphiaID, valid_AphiaID
species_df <- read_excel("SpeciesList_Cetacea_WoRMS.xlsx", col_types = "text")
species_vec <- species_df$scientificname
length(species_vec) #adapt chunks in the next section accordingly


########## RETRIEVE BOLD DATA #########################

#Download bold data for our species list in chunks
listBOLD1 <- bold_seqspec(taxon = species_vec[1:300], format = "tsv")
listBOLD2 <- bold_seqspec(taxon = species_vec[301:600], format = "tsv")
listBOLD3 <- bold_seqspec(taxon = species_vec[601:900], format = "tsv")
listBOLD4 <- bold_seqspec(taxon = species_vec[901:1200], format = "tsv")
listBOLD5 <- bold_seqspec(taxon = species_vec[1201:1280], format = "tsv")
listBOLD <- bind_rows(listBOLD1, listBOLD2, listBOLD3, listBOLD4, listBOLD5)

######## Optional if previous function does not work ##############

# helper to fetch one chunk and normalise types if the previous function works
prep_chunk <- function(x) {
  df <- suppressWarnings(bold_seqspec(taxon = x, format = "tsv"))
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(tibble())
  df %>%
    #list-cols (if any) -> plain strings
    mutate(across(where(is.list), ~ vapply(.x, toString, character(1)))) %>%
    # make every column character so bind_rows can't complain
    mutate(across(everything(), as.character))
}

n <- length(species_vec)

# --- keep your original chunk layout, but guard the indices ---
listBOLD1 <- prep_chunk(species_vec[1:min(300, n)])
listBOLD2 <- if (n >= 301) prep_chunk(species_vec[301:min(600, n)])  else tibble()
listBOLD3 <- if (n >= 601) prep_chunk(species_vec[601:min(900, n)])  else tibble()
listBOLD4 <- if (n >= 901) prep_chunk(species_vec[901:min(1200, n)]) else tibble()
listBOLD5 <- if (n >= 1201) prep_chunk(species_vec[1201:min(1280, n)]) else tibble()
listBOLD <- bind_rows(listBOLD1, listBOLD2, listBOLD3, listBOLD4, listBOLD5)

########## End of optional #############################

#See which markers you have in your data
unique(listBOLD$markercode)

#Remove GenBank-mined records explicitly
listBOLD <- listBOLD %>%
  mutate(GenbankMined = institution_storing == "Mined from GenBank, NCBI") %>%
  filter(!GenbankMined)

#Rename collecting_date to collection_date if present
if("collecting_date" %in% names(listBOLD)) {
  listBOLD <- rename(listBOLD, collection_date = collecting_date)
}

#Check markers again after cleaning
unique(listBOLD$markercode)

#Creates a new column called number with the length of the DNA sequence
listBOLD$length <- str_count(listBOLD$nucleotides, "[A-Z]")

#Save total BOLD data
write_xlsx(listBOLD, "listBOLD_total.xlsx")

# Save COI BOLD data
write_xlsx(listBOLD[listBOLD$markercode == "COI-5P", ], "listBOLD_COI.xlsx")


###### GETTING BOLD STATS FOR EACH MOLECULAR MARKER IN CASE YOU USE MORE THAN ONE #####################

#Define safe taxon status function once - bold_tax_id2() fetches taxon-level stats from the BOLD database for a given species ID. The wrapper ensures that if any error occurs (e.g., no data for a species), it returns a dummy row with NA values instead of stopping the script.
#if there are warnings of content type it's okay to ignore as they aren't errors
safe_tax_id <- possibly(
  .f = function(id) bold_tax_id2(id = id, dataTypes = "stats"),
  otherwise = tibble(AphiaID = NA, total = NA, bin_count = NA)
)

#marker-specific processing function for stats for each molecular marker
process_marker <- function(marker_code, output_prefix, length_filter = NULL) {
  message("Processing ", marker_code)
  
  #Filter for marker
  marker_data <- filter(listBOLD, markercode == marker_code)
  
  #Optional sequence length filter
  if (!is.null(length_filter)) {
    marker_data <- filter(marker_data, length >= length_filter)
  }
  
  #Save raw marker data
  write_xlsx(marker_data, path = paste0("listBOLD_", output_prefix, ".xlsx"))
  
  #Summarize stats per species
  stats <- marker_data %>%
    group_by(phylum_name, class_name, species_name, species_taxID) %>%
    summarize(
      Public_BOLD = n(),
      countries = if("country" %in% names(marker_data)) paste(unique(na.omit(country)), collapse = "; ") else NA_character_,
      collectiondate_start = if("collectiondate_start" %in% names(marker_data)) {
        dates <- as.Date(na.omit(collectiondate_start), "%Y-%m-%d")
        if (length(dates) > 0) min(dates) else as.Date(NA)
      } else as.Date(NA),
      
      collectiondate_end = if("collectiondate_end" %in% names(marker_data)) {
        dates <- as.Date(na.omit(collectiondate_end), "%Y-%m-%d")
        if (length(dates) > 0) max(dates) else as.Date(NA)
      } else as.Date(NA),
      
      .groups = "drop"
    )
  
  # Fetch BOLD taxon-level stats
  stats2 <- map_dfr(
    stats$species_taxID,
    function(x) {
      res <- suppressWarnings(safe_tax_id(x))
      if (!"species_taxID" %in% colnames(res)) {
        res <- tibble(species_taxID = x, AphiaID = NA, total = NA, bin_count = NA)
      } else {
        res$species_taxID <- x
      }
      res
    }
  )
  
  #Join summarized stats with taxon-level stats
  stats_final <- left_join(stats, stats2, by = "species_taxID")
  
  #Makes dates character to match detailed_columns
  stats_final <- stats_final %>% 
    dplyr::mutate(across(c(collectiondate_start, collectiondate_end), as.character))
  
  #Selection of detailed columns of interest
  detailed_columns <- marker_data %>%
    select(processid, sampleid, recordID, bin_uri, phylum_name, class_name, order_name, family_name,
           collectiondate_start, collectiondate_end, lat, lon, country, province_state, region, 
           sequenceID, markercode, nucleotides, length, species_name)
  
  #Merge detailed data with stats_final (robust join)
  full_stats_final <- left_join(stats_final, detailed_columns, 
                                by = c("phylum_name", "class_name", "species_name", 
                                       "collectiondate_start", "collectiondate_end", "countries" = "country"))
  
  #Final export
  write_xlsx(full_stats_final, path = paste0("listBOLD_", output_prefix, "_statsfinal.xlsx"))
}

#Run stats
process_marker("COI-5P", "COI")


########## RETRIEVE GENBANK DATA #########################

#this is an API KEY which you can ask for in NCBI to increase your power of search
Sys.setenv(ENTREZ_KEY = "YOUR_KEY_HERE") 

#The queries are more extensive to make sure we don't miss out on any important data
genbank_queries <- list(
  COI   = 'Cetacea[Organism] AND (COI[Gene] OR CO1[Gene] OR COXI[Gene] OR COX1[Gene] OR "complete genome"[All Fields] OR "mitochondrial genome"[All Fields])',
  `12S` = 'Cetacea[Organism] AND ("12S"[All Fields] OR "12S ribosomal RNA"[Title] OR "12S rRNA"[Title] OR "complete genome"[All Fields] OR "mitochondrial genome"[All Fields])',
  `16S` = 'Cetacea[Organism] AND ("16S"[All Fields] OR "16S ribosomal RNA"[Title] OR "16S rRNA"[Title] OR "complete genome"[All Fields] OR "mitochondrial genome"[All Fields])',
  Dloop = 'Cetacea[Organism] AND ("D-loop"[All Fields] OR "control region"[All Fields] OR "complete genome"[All Fields] OR "mitochondrial genome"[All Fields])',
  CytB  = 'Cetacea[Organism] AND (CYTB[Gene] OR "cytochrome b"[Gene] OR "cyt b"[All Fields] OR "complete genome"[All Fields] OR "mitochondrial genome"[All Fields])'
)

#Function to retrieve GenBank data
fetch_genbank_records <- function(query, marker, retmax = 500) {
  message("🔍 Searching GenBank for: ", marker)
  search <- tryCatch(
    entrez_search(db = "nuccore", term = query, use_history = TRUE, retmax = 0),
    error = function(e) return(NULL)
  )
  if(is.null(search)||is.null(search$web_history)||search$count==0) {
    message("⚠️ No records for: ", marker);
    return(NULL)
  }
  total <- search$count; message("  📄 Found: ", total)
  results <- character()
  for(start in seq(0, total-1, by=retmax)){
    message("  ⏳ Fetching ", start+1, " to ", min(start+retmax, total))
    xml_batch <- tryCatch(
      entrez_fetch(db="nuccore", web_history=search$web_history, rettype="gb", retmode="text", retstart=start, retmax=retmax),
      error = function(e) NULL
    )
    if(!is.null(xml_batch)) results <- c(results, xml_batch)
    Sys.sleep(0.5)
  }
  file_out <- paste0("GenBank_", marker, ".gb")
  writeLines(results, file_out)
  message("✅ Saved to: ", file_out)
}

#Download GenBank records for all markers
for(marker in names(genbank_queries)) {
  fetch_genbank_records(genbank_queries[[marker]], marker)
}