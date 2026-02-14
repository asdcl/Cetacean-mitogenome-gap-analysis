# Cetacean-mitogenome-gap-analysis

README — Cetacea Barcode Data Retrieval & Parsing Scripts 

Overview
These R scripts automate:
1. Retrieval of barcode metadata from BOLD Systems (for a WoRMS species list)
2. Fetching GenBank flatfiles via NCBI (rentrez)
3. Parsing GenBank features/ qualifiers and filtering target products
4. Cleaning, summary stats, and Excel/CSV exports for five mitochondrial markers: COI, 12S, 16S, D-loop, CytB

Use it to build local reference datasets or perform gap analyses on cetacean barcode coverage.

Note: Public databases evolve. Re-running retrieval at a later date may produce slightly different totals. For transparency, this repository includes frozen outputs used in the manuscript (see “Reproducibility data”).

Requires a species list file to retrieve data from BOLD
•	Name: SpeciesList_Cetacea_WoRMS.xlsx (included in data/raw/)
•	Must include at least the column: scientificname (valid names retrieved from WoRMS)

Dependencies
•	Install these R packages if not already installed:
install.packages(c("readxl", "writexl", "bold", "rentrez", "dplyr", "stringr", "purrr", "tibble", "openxlsx"))

•	NCBI API key (optional, but recommended for large downloads). Request one at https://www.ncbi.nlm.nih.gov/account/settings/. 
Add here:
> Sys.setenv(ENTREZ_KEY="YOUR_KEY_HERE") 

Repository structure:
•	R/ : contains R scripts
•	data/raw/ : contains required inputs
•	data/derived/ : contains frozen outputs and derived tables

Important: `Genbank_files_data_extraction.R` is configured for a single marker at a time (e.g., D-loop). To run other markers, update the input `.gb` filename and the marker keyword list used for filtering.

Reproducibility data:
Included in this repository is:
•	Input species list: data/raw/SpeciesList_Cetacea_WoRMS.xlsx
•	Frozen BOLD exports: data/derived/listBOLD_*.xlsx
•	Frozen GenBank derived tables:
  o	data/derived/GenBank_*_Metadata.xlsx
  o	data/derived/GenBank_*_Filtered.xlsx
  o	Large raw files (provided as GitHub Release assets)
  
GitHub blocks very large files in standard commits. The following raw inputs used in the manuscript are provided as GitHub Release assets (Release tag: data-v1):
•	GenBank_12S.gb, GenBank_16S.gb, GenBank_COI.gb, GenBank_CytB.gb, GenBank_Dloop.gb
•	GenBank_totalsearch.csv

After downloading, place:
•	GenBank flatfiles (*.gb) in data/raw/
•	GenBank_totalsearch.csv in data/derived/

Quick start: 
1.	Clone the repository.
2.	Confirm the small files exist in data/raw/ and data/derived/.
3.	Download large files from Release data-v1 and place them in data/raw/ and data/derived/ as described above.
4.	Run scripts in R/ from the repository root in R/RStudio.

Citation:
If you use or modify this code, please cite our paper: (manuscript under review - reference will be updated when published)
