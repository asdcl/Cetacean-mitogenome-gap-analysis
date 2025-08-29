# gap-analysis-cetacea

README — Cetacea Barcode Data Retrieval Script (Criar repositório no meu Github e depois alterar para o título do paper)

Overview
This R script automates:
1. Retrieval of barcode sequences from BOLD Systems
2. Fetching and parsing GenBank records via NCBI
3. Filtering, summary statistics, and Excel exports
for five mitochondrial markers: COI, 12S, 16S, D-loop, and CytB.

Use it to build local reference datasets or perform gap analyses on cetacean barcode coverage.

Requires a species list file to retrieve data from BOLD
•	Name: SpeciesList_Cetacea_WoRMS.xlsx
•	Must include at least the column: scientificname

Dependencies
•	Install these R packages if not already installed:
install.packages(c("readxl", "writexl", "bold", "rentrez", "dplyr", "stringr", "purrr", "tibble"))

•	NCBI API key (optional, but recommended for large downloads). Request one at https://www.ncbi.nlm.nih.gov/account/settings/. Add here:

> Sys.setenv(ENTREZ_KEY="YOUR_KEY_HERE") (retirar isto do script antes de colocar no github)

Citation:
If you use or modify this code, please cite our paper:

