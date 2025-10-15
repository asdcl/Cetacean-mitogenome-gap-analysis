# Cetacean-mitogenome-gap-analysis

README — Cetacea Barcode Data Retrieval Script 

Overview
This R script automates:
1. Retrieval of barcode sequences from BOLD Systems
2. Fetching and parsing GenBank records via NCBI
3. Filtering, cleaning, and compiling metadata
4. Export of summary statistics and Excel tables
for five mitochondrial markers: COI, 12S, 16S, D-loop, and CytB.

Use it to build local reference datasets or perform gap analyses on cetacean barcode coverage.

Requires a species list file to retrieve data from BOLD
•	Name: SpeciesList_Cetacea_WoRMS.xlsx
•	Must include at least the column: scientificname (valid names retrieved from WoRMS)

Dependencies
•	Install these R packages if not already installed:
install.packages(c("readxl", "writexl", "bold", "rentrez", "dplyr", "stringr", "purrr", "tibble"))

•	NCBI API key (optional, but recommended for large downloads). Request one at https://www.ncbi.nlm.nih.gov/account/settings/. 
Add here:
> Sys.setenv(ENTREZ_KEY="YOUR_KEY_HERE") 

Citation:
If you use or modify this code, please cite our paper:
