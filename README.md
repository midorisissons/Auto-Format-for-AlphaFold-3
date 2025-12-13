# Auto-Format-for-AlphaFold-3
Automated formatting of JSON files for AlphaFold 3, with fragmenting of protein sequence

## Installation
### 1. Install Python
[Download and install Python](https://www.python.org/downloads/)
### 2. Download code from GitHub
[Download and extract zip from GitHub](https://github.com/midorisissons/Auto-Format-for-AlphaFold-3/archive/refs/heads/main.zip)
### 3. Extract the folder
Opening the zip file will make a folder called Auto-Format-for-AlphaFold-3-main. All AlphaFold protein sequences, fragments, and jobs will be stored here, so move the folder to where you wish on your computer.
### 4. Setup
For MacOS/Linux, run setup.sh 🤫 in the extracted folder by `sh setup.sh` directly in your terminal
For Windows run setup.bat 🦇 by double clicking the file 

### 5. Ready to run
For MacOS/Linux, run AutoFormatAF.sh directly in your terminal by `sh AutoFormatAF.sh`. Make sure you are currently in the directory of the extracted folder.
For Windows, run AutoFormatAF.bat by double clicking the file.

## Functionalities
Auto-Format-for-AlphaFold-3 has three main functionalities, described below.

### 1. Get a protein sequence
The program is able to retrieve and store protein sequences of interest into a storage CSV file "store_sequences.csv" within the extracted folder.
*Protein sequences can be saved via two methods:*
1. Manual input - manually input the protein ID and sequence
2. Download from UniProt - enter the UniProt ID (eg. Q8IUZ5) to automatically retrieve the sequence
Stored protein sequences can be retrieved via their stored ID

### 2. Fragment a protein sequence
The program can take a stored protein sequence and fragment it. The user can input the desired fragment size and overlap size, as well as section of the protein sequence of interest (start residue:end residue). If there is a remainder, some of the fragments' lengths are extended by one to use the full sequence.

Fragments are stored in the "fragments" folder in the extracted folder, under a CSV file named `ID_size[size]_overlap[overlapsize]_residues[startresidue]-[endresidue]`.

### 3. Automatically make AlphaFold files
AlphaFold Server takes JSON format files to make tasks. This program can automatically output JSON files with stores sequences/fragments, ions and ligands.
1. Full length protein sequences - retrieved from the CSV store
2. Fragment protein sequences - full CSV files are retrieved and read fragment-by-fragment from "fragments" folder. AlphaFold can only take up to 100 tasks in one JSON, if this is exceeded, multiple files may be outputted.

[How to specify ions and ligands by CCD codes](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)
