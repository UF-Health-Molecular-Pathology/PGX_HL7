## PGX TaqMan Genotyping Data Analysis and HL7 Generation

This repository contains Python code for analyzing PGX (Pharmacogenomics) TaqMan Genotyping data from Quantstudios. The code processes raw data files, translates results, and generates various reports for review and reporting purposes.

## Table of Contents

    Introduction
    Installation
    Project Structure    
    Description
    Usage
    Functionality
    Results
    Contributing
    License

## Introduction

    Pharmacogenomics (PGX) is the study of how an individual's genetic makeup affects their response to drugs. 
    PGX testing provides information about how a person's genes may impact their response to specific medications. 
    This repository focuses on the analysis of PGX TaqMan Genotyping data generated by Quantstudios and the generation 
    of an HL7 formated result that is compatible with EPIC beaker. 
    
    Important note: The tables provided here are labratory/assay specific and will need to be updated to refelct the 
    probes and 


## Installation

    Clone the repository to your local machine:

    git clone https://github.com/kja3/PGX_HL7.git
    cd PGX_HL7

    ## Install the required Python packages:

    Python 3
    NumPy
    pandas
    hl7
    Matplotlib

## Project Structure
    
    PGX_HL7/
    ├── HL7_Orders/
    │ ├── 18279-85bfc99a-0f1a-426e-add1-c6ecbaacab03.txt
    │ ├── 18281-f964c4cf-f41f-4fdf-a971-439acb9ec736.txt
    │ └── 18957-4d514a60-3de0-4aec-a55e-253ac239b7af.txt
    ├── HL7_Results/
    │ ├── hl7-GatorPGX-100067795-output.txt
    │ ├── hl7-GatorPGX-100067800-output.txt
    │ └── hl7-GatorPGX-100079914-output.txt
    ├── PGX_Tables/
    │ ├── GT_PT_Translator.csv
    │ └── QS_Translator.csv
    ├── QS_Data/
    │ ├── CN GPGX-QS-23-30b Results.csv
    │ ├── GPGX-QS-23-30b_20230328_112910_Export.csv
    │ └── PGX_TaqMan_Genotyping_Data_Review.pdf
    ├── PGX_HL7.py
    └── README.md

## Description

    - **HL7_Orders:** Contains HL7 order files for pharmacogenomics analysis.
    - **HL7_Results:** Stores the output HL7 files after PGx analysis. The files present here are for example output.
    - **PGX_Tables:** Holds translation tables for genotyping and other relevant data. This is labratory dependant and the hardest part! 
    - **QS_Data:** This is the Input folder for the TaqMan genotyping and copy number results.

## Usage

       1. Clone the repository.
       2. Navigate to the `PGX_HL7` directory.
       3. Place your raw PGX TaqMan Genotyping data files (CSV format) in the QS_Data folder. Example data is provided.
       3. Run the `PGX_HL7.py` script for PGx analysis.
       4. Review the generated reports in the output directories HL7_Results for the HL7 files and QS_Data for the 
          PGX_TaqMan_Genotyping_Data_Review PDF file.

   
## Functionality

    The code provides the following functionality:

    Data Processing: Reads raw PGX TaqMan Genotyping data, cleans and organizes the information.
    Translation: Translates the results of the TaqMan Genotyping assays to Genotype-Phenotype (GT-PT) information.
    Report Generation: Generates review reports, full data reports, and reports suitable for HL7 reporting.
    PDF Generation: Creates a PDF file containing tables for reviewing PGX TaqMan Genotyping results.
    HL7 Reporting: Prepares Genotype-Phenotype (GT-PT) results DataFrame for HL7 reporting.

## Results

    The results include review reports, full data reports, and HL7-compatible reports. The PDF file provides a visual 
    representation of the data for easy review.
    
## Contributing
    
    Contributions to this project are welcome. Feel free to open issues for bug reports, feature requests, or general 
    feedback. If you would like to contribute code, please fork the repository and submit a pull request.

## License

This project is licensed under the MIT License. #PGX_HL7
