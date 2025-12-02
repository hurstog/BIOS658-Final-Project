# BIOS658-Final-Project
#This project contains two scripts: "project_code.r" and "functional_analysis.r"
# project_code.r:
  #Description: Main code file of the analysis for this project. It downloads data from GEO and downloads the raw gene expression counts and conducts most of     the analytical steps of the project. It contains multiple sections of results and accompanying plots/figures. The main ouput from this file is a CSV file of    the differential expression analysis.
  #Inputs:
    #GEO data (via getGEO)
    #Raw data of gene counts (via downloaded CSV)
    #Results of functional enrichment analysis (GSEA2_proj_deg.xlsx): from "functional_analysis.r"
  #Outputs:
    #CSV of gene expression analysis results (DGE_results.csv): contains log2fold change, p-value, adjusted p-value, FDR, ensemble gene IDs, HGNC Gene Symbols,     etc.)

# functional_analysis.r:
  #Description: This file takes the results of the differential gene expression and conducts functional enrichment analysis. Results are generated fro KEGG       pathways and three gene ontologies using both GSEA and overrepresentation analysis. Results are saved as an excel file.
  #Inputs:
    #CSV of gene expression results (DGE_results.csv)
  #Outputs:
    #XLSX of the results of analysis (GSEA2_proj_deg.xlsx): One sheet for eahc possible combination of analysis method and domain/pathways (8 total sheets)
  
