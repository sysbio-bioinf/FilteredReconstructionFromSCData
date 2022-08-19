# FilteredReconstructionFromSCData

Source code to reconstruct Boolean networks from times series of data. The pipeline specified in this code combines the algorithm by Maucher et al. (Bioninformatics, 2011) with the best-fit extension algorithm by Lähdesmäki et al. (Bioinformatics, 2003). The algorithm by Maucher et al. is used as preprocessing step to reduce the number of inputs that will be tested for the succeeding best-fit extension algorithm.

Furthermore, the respository contains an R Markdown File (ReconstructingBooleanNetworksEnsemblesExample.Rmd) and RData objects from the analyses of a publicly available single-cell RNA sequencing data from Ratliff and colleagues (Immun Aging, 2020; NCBI Gene Expression Omnibus GSE138544). The PDF with a description of the analyses, the code, and its results can also be found in the repository (ReconstructingBooleanNetworksEnsemblesExample.pdf)

Source data files are available from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138544). Intermediate files can be created by rendering the rmarkdown file.
