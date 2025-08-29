# Set working directory
setwd("C:/Users/aslav/Desktop/CIIMAR/Script Luís Cetáceos")

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
length(species_vec) #adapt chuncks in the next section accordingly


########## RETRIEVE BOLD DATA #########################

#Download bold data for our species list in chuncks
listBOLD1 <- bold_seqspec(taxon = species_vec[1:300], format = "tsv")
listBOLD2 <- bold_seqspec(taxon = species_vec[301:600], format = "tsv")
listBOLD3 <- bold_seqspec(taxon = species_vec[601:900], format = "tsv")
listBOLD4 <- bold_seqspec(taxon = species_vec[901:1200], format = "tsv")
listBOLD5 <- bold_seqspec(taxon = species_vec[1201:1280], format = "tsv")
listBOLD <- bind_rows(listBOLD1, listBOLD2, listBOLD3, listBOLD4, listBOLD5)

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

# Save marker-specific BOLD data
write_xlsx(listBOLD[listBOLD$markercode == "COI-5P", ], "listBOLD_COI.xlsx")
write_xlsx(listBOLD[listBOLD$markercode == "12S", ], "listBOLD_12S.xlsx")
write_xlsx(listBOLD[listBOLD$markercode == "16S", ], "listBOLD_16S.xlsx")
write_xlsx(listBOLD[listBOLD$markercode == "D-loop", ], "listBOLD_Dloop.xlsx")
write_xlsx(listBOLD[listBOLD$markercode == "CYTB", ], "listBOLD_CytB.xlsx")


###### GETTING BOLD STATS FOR EACH MOLECULAR MARKER #####################

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

#Run stats for all markers
process_marker("COI-5P", "COI")
process_marker("12S",    "12S")
process_marker("16S",    "16S")
process_marker("D-loop", "Dloop")
process_marker("CYTB", "CYTB")


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

########### EXTRACT GENBANK DATA FROM GB FILES ###################
parse_genbank_flatfile <- function(file_path, marker_name) {
  lines <- readLines(file_path, warn = FALSE)
  records <- unlist(strsplit(paste(lines, collapse = "\n"), "//\\s*\n"))
  
  result <- map_dfr(records, function(rec) {
    acc   <- str_match(rec, "ACCESSION\\s+(\\S+)")[,2]
    org   <- str_match(rec, "  ORGANISM\\s+([^\n]+)")[,2]
    desc  <- str_trim(str_match(rec, "DEFINITION\\s+(.+?)(\\n[A-Z ]|$)")[,2])
    date  <- str_match(rec, "/collection_date=\"([^\"]+)\"")[,2]
    note  <- str_match(rec, "/note=\"([^\"]+)\"")[,2]
    
    # Get country or geo_loc_name
    geo_loc_matches <- c(
      str_match_all(rec, "/geo_loc_name=\"([^\"]+)\"")[[1]][,2],
      str_match_all(rec, "/country=\"([^\"]+)\"")[[1]][,2]
    )
    geo <- if (length(geo_loc_matches) > 0) geo_loc_matches[1] else NA_character_
    
    # Get /lat_lon as coordinates string
    coord_str <- str_match(rec, "/lat_lon=\"([^\"]+)\"")[,2]
    
    # Optionally also extract as numeric lat/lon
    if (!is.na(coord_str)) {
      parts <- unlist(str_split(coord_str, " "))
      lat <- as.numeric(parts[1]) * ifelse(parts[2] == "S", -1, 1)
      lon <- as.numeric(parts[3]) * ifelse(parts[4] == "W", -1, 1)
    } else {
      lat <- NA_real_
      lon <- NA_real_
    }
    
    # Extract all genes and products
    genes <- str_match_all(rec, "/gene=\"([^\"]+)\"")[[1]][,2]
    prods <- str_match_all(rec, "/product=\"([^\"]+)\"")[[1]][,2]
    
    # Handle mismatched lengths
    if (length(genes) == 0 && length(prods) > 0) {
      genes <- rep(NA, length(prods))
    }
    if (length(prods) == 0 && length(genes) > 0) {
      prods <- rep(NA, length(genes))
    }
    
    max_len <- max(length(genes), length(prods))
    genes <- rep_len(genes, max_len)
    prods <- rep_len(prods, max_len)
    
    tibble(
      marker          = marker_name,
      accession       = acc,
      species         = org,
      gene            = genes,
      product         = prods,
      country         = geo,
      coordinates     = coord_str,
      latitude        = lat,
      longitude       = lon,
      collection_date = date,
      note            = note,
      description     = desc
    )
  }) %>% filter(!is.na(accession))
  
  return(result)
}

filter_genbank_marker <- function(df, marker) {
  marker <- tolower(marker)
  
  if (marker == "coi") {
    df <- df %>%
      filter(tolower(gene) %in% c("coi", "co1", "cox1", "coxi") |
               str_detect(tolower(product), "coi|co1|cox1|coxi"))
    
  } else if (marker == "12s") {
    df <- df %>%
      filter(str_detect(tolower(product), "12s") |
               str_detect(tolower(description), "12s"))
    
  } else if (marker == "16s") {
    df <- df %>%
      filter(str_detect(tolower(product), "16s") |
               str_detect(tolower(description), "16s"))
    
  } else if (marker == "dloop" || marker == "d-loop") {
    df <- df %>%
      filter(str_detect(tolower(product), "d[- ]?loop|control region") |
               str_detect(tolower(description), "d[- ]?loop|control region"))
    
  } else if (marker == "cytb") {
    df <- df %>%
      filter(str_detect(tolower(gene), "cytb|cytochrome b") |
               str_detect(tolower(product), "cytb|cytochrome b"))
  }
  
  return(df)
}

# Run parsing and filtering for each marker
markers <- c("COI", "12S", "16S", "Dloop", "CytB")
all_data <- list()
metadata_all <- list()

for (m in markers) {
  file <- paste0("GenBank_", m, ".gb")
  if (file.exists(file)) {
    df <- parse_genbank_flatfile(file, m)
    metadata_all[[m]] <- df
    write_xlsx(df, paste0("GenBank_", m, "_Metadata.xlsx"))
    message("✅ Saved: GenBank_", m, "_Metadata.xlsx")
    
    # Apply filtering
    df_filtered <- filter_genbank_marker(df, m)
    all_data[[m]] <- df_filtered
    write_xlsx(df_filtered, paste0("GenBank_", m, "_Filtered.xlsx"))
    message("✅ Saved filtered: GenBank_", m, "_Filtered.xlsx")
    
  } else {
    message("⚠️ Missing file: ", file)
  }
}

# Export all metadata in one Excel (1 sheet per marker)
write_xlsx(metadata_all, "GenBank_Metadata_AllMarkers.xlsx")

# Export all filtered in one Excel (1 sheet per marker)
write_xlsx(all_data, "GenBank_Filtered_AllMarkers.xlsx")

# Create a summary per accession with the filtered data
all_combined <- bind_rows(all_data, .id = "marker")

summary_table <- all_combined %>%
  group_by(accession, species) %>%
  summarize(
    markers_present = paste(unique(marker), collapse = "; "),
    genes_present   = paste(unique(na.omit(gene)), collapse = "; "),
    products        = paste(unique(na.omit(product)), collapse = "; "),
    countries       = paste(unique(na.omit(country)), collapse = "; "),
    coordinates     = paste(unique(na.omit(coordinates)), collapse = "; "),
    .groups = "drop"
  )

write_xlsx(summary_table, "GenBank_Accessions_Summary.xlsx")
message("✅ Saved summary per accession.")


############OPTIONAL#######
# Create a summary per accession with the raw data
all_combined <- bind_rows(metadata_all, .id = "marker")

summary_table <- all_combined %>%
  group_by(accession, species) %>%
  summarize(
    markers_present = paste(unique(marker), collapse = "; "),
    genes_present   = paste(unique(na.omit(gene)), collapse = "; "),
    products        = paste(unique(na.omit(product)), collapse = "; "),
    countries       = paste(unique(na.omit(country)), collapse = "; "),
    coordinates     = paste(unique(na.omit(coordinates)), collapse = "; "),
    .groups = "drop"
  )

write_xlsx(summary_table, "GenBank_Accessions_Metadata_Summary.xlsx")
message("✅ Saved summary per accession.")