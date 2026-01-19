# Cetacean-mitogenome-gap-analysis

README — Cetacea Barcode Data Retrieval & Parsing Scripts 

Overview
These R scripts automates:
1. Retrieval of barcode metadata from BOLD Systems (for a WoRMS species list)
2. Fetching GenBank flatfiles via NCBI (rentrez)
3. Parsing GenBank features/ qualifiers and filtering target products
4. Cleaning, summary stats, and Excel/CSV exports for five mitochondrial markers: COI, 12S, 16S, D-loop, CytB

Use it to build local reference datasets or perform gap analyses on cetacean barcode coverage.

Requires a species list file to retrieve data from BOLD
•	Name: SpeciesList_Cetacea_WoRMS.xlsx
•	Must include at least the column: scientificname (valid names retrieved from WoRMS)

Dependencies
•	Install these R packages if not already installed:
install.packages(c("readxl", "writexl", "bold", "rentrez", "dplyr", "stringr", "purrr", "tibble", "openxlsx"))

•	NCBI API key (optional, but recommended for large downloads). Request one at https://www.ncbi.nlm.nih.gov/account/settings/. 
Add here:
> Sys.setenv(ENTREZ_KEY="YOUR_KEY_HERE") 

Citation:
If you use or modify this code, please cite our paper: (manuscript submitted - reference will be updated when published)
