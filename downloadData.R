#' Download, subset, and tidy pharmacogenomics datasets 
#'
#' (c) 2017 Alejandro Reyes, Keegan Korthauer
#' (c) 2019 Patrick Kimes, Kelly Street
#' (c) 2020 Patrick Kimes

## use earlier stable version of PharmacoGx package
remotes::install_github("bhklab/PharmacoGx", ref = "RELEASE_3_10")

## load libraries
library(Biobase)
library(PharmacoGx)
library(tidyverse)

## download parsed data for GDSC, CCLE studies from Haibe-Kains lab
GDSC <- downloadPSet("GDSC")
CCLE <- downloadPSet("CCLE")

## only use intersection of two studies
common <- intersectPSet(list('CCLE' = CCLE, 'GDSC' = GDSC),
                        intersectOn = c("cell.lines", "drugs"),
                        strictIntersect = TRUE)

## specify output directory
dir.create("data", showWarnings = FALSE)


## ##############################################################################
## Parse raw dose and viability measures merged by Haibe-Kains et al. (2014)
## ##############################################################################

## helper function to transform 'raw' data in PSet to data.frame 
rawToDataFrame <- function(dataset) {
    rawdat <- common[[dataset]]@sensitivity$raw
    names(dimnames(rawdat)) <- c("cellLine_drug", "doseID", "measure")
    rawdat <- as.data.frame.table(rawdat, responseName = "value",
                                  stringsAsFactors = FALSE)
    rawdat <- dplyr::mutate(rawdat, value = as.numeric(as.character(value)))
    rawdat <- tidyr::spread(rawdat, measure, value)
    rawdat <- tidyr::separate(rawdat, cellLine_drug,
                              c("cellLine", "drug"), sep = "_")
    rawdat <- dplyr::rename(rawdat, concentration = Dose, viability = Viability)
    rawdat <- dplyr::filter(rawdat, !is.na(concentration), !is.na(viability))
    rawdat
}

## transform PSets to data.frame (long)
rawData <- dplyr::bind_rows(CCLE = rawToDataFrame("CCLE"),
                            GDSC = rawToDataFrame("GDSC"),
                            .id = "study")

## reorder columns
colOrder <- c("cellLine", "drug", "doseID", "concentration", "viability", "study")
rawData <- rawData[, colOrder]

## determine set of cell lines in our data subset
allCellLines <- unique(rawData$cellLine)

## save as RDS file
saveRDS(rawData, file = file.path("data", "rawPharmacoData.rds"))

## ## alernatively, save as csv files
## write.table(rawData, sep = ",", quote = FALSE, col.names = TRUE,
##             row.names = FALSE, file = file.path("data", "rawPharmacoData.csv"))


## ##############################################################################
## Parse summarized sensitivity values from Haibe-Kains et al. (2014)
## -- will use "published" values (data also includes recomputed values)
## ##############################################################################

## helper function to transform summarized 'profile' data in PSet to data.frame 
profilesToDataFrame <- function(dataset) {
    sumdat <- common[[dataset]]@sensitivity$profiles
    keepCols <- c("ic50_published", "auc_published")
    sumdat <- sumdat[, keepCols]
    sumdat <- tibble::rownames_to_column(sumdat, "cellLine_drug")
    sumdat <- dplyr::rename_if(sumdat, is.numeric, ~gsub("_published", "", .))
    sumdat <- tidyr::separate(sumdat, cellLine_drug,
                              c("cellLine", "drug"), sep = "_")
    sumdat
}

## transform PSets to data.frame (wide)
sumData <- dplyr::inner_join(profilesToDataFrame("CCLE"),
                             profilesToDataFrame("GDSC"),
                             by = c("cellLine", "drug"),
                             suffix = c("_CCLE", "_GDSC"))

## save as RDS file
saveRDS(sumData, file = file.path("data", "summarizedPharmacoData.rds"))

## ## alernatively, save as csv files
## write.table(sumData, sep = ",", quote = FALSE, col.names = TRUE,
##             row.names = FALSE, file = file.path("data", "summarizedPharmacoData.csv"))


## ##############################################################################
## Parse gene expression measurements merged by Haibe-Kains et al. (2014)
## ##############################################################################

## to parse gene expression data set to TRUE
if (FALSE) { 
    ## helper function to transform gene expression data in PSet to data.frame 
    exprToDataFrame <- function(dataset) {
        rna_pdat <- pData(common[[dataset]]@molecularProfiles$rna)
        rna_vals <- exprs(common[[dataset]]@molecularProfiles$rna)
        ## only keep first instance of each cell line
        rna_keep <- !duplicated(rna_pdat$cellid)
        
        rna_vals <- rna_vals[, rna_keep]
        colnames(rna_vals) <- rna_pdat$cellid[rna_keep]
        
        rna_vals <- as.data.frame(rna_vals, optional = TRUE)
        rna_vals <- tibble::rownames_to_column(rna_vals, "gene")
        rna_vals <- tidyr::gather(rna_vals, cellLine, expression, -gene)
        rna_vals
    }

    ## transofmr PSets to data.frame (wide)
    exprData <- dplyr::inner_join(exprToDataFrame("CCLE"),
                                  exprToDataFrame("GDSC"),
                                  by = c("gene", "cellLine"),
                                  suffix = c("_CCLE", "_GDSC"))

    ## note: expression data includes more cell lines than the intersection
    ## subset on cell lines with sensitivity data
    exprData <- dplyr::filter(exprData, cellLine %in% !! allCellLines)

    ## save as RDS file
    saveRDS(exprData, file = file.path("data", "exprPharmacoData.rds"))
    
    ## ## alernatively, save as csv files
    ## write.table(exprData, sep = ",", quote = FALSE, col.names = TRUE,
    ##             row.names = FALSE, file = file.path("data", "exprPharmacoData.csv"))
}


## ##############################################################################
## Parse cell line metadata cleaned by Haibe-Kains et al. (2014)
## -- parse cell line metadata in CCLE, GDSC objects separately
## -- overlap is moderate, will keep data separate and only use subset of columns
## ##############################################################################

## to parse cell line metadata set to TRUE
if (FALSE) { 
    ## parse cell line metadata from CCLE
    clDataCCLE <- common[["CCLE"]]@cell
    clDataCCLE <- dplyr::select(clDataCCLE, cellid, tissueid, Gender,
                                Site.Primary, Histology, Hist.Subtype1)
    clDataCCLE <- dplyr::filter(clDataCCLE, cellid %in% !! allCellLines)
    clDataCCLE <- dplyr::rename(clDataCCLE,
                                cellLine = cellid,
                                tissue = tissueid,
                                gender = Gender, 
                                primarySite = Site.Primary,
                                histology = Histology,
                                subtype = Hist.Subtype1)

    ## parse cell line metadata from GDSC
    clDataGDSC <- common[["GDSC"]]@cell
    clDataGDSC <- dplyr::select(clDataGDSC, cellid, tissueid, Primary.site,
                                Primary.histology, Histology.subtype)
    clDataGDSC <- dplyr::filter(clDataGDSC, cellid %in% !! allCellLines)
    clDataGDSC <- dplyr::rename(clDataGDSC,
                                cellLine = cellid,
                                tissue = tissueid, 
                                primarySite = Primary.site,
                                histology = Primary.histology,
                                subtype = Histology.subtype)

    ## save as RDS file
    saveRDS(clDataCCLE, file = file.path("data", "cellLineDataCCLE.rds"))
    saveRDS(clDataGDSC, file = file.path("data", "cellLineDataGDSC.rds"))

    ## ## alernatively, save as csv files
    ## write.table(clDataCCLE, sep = ",", quote = FALSE, col.names = TRUE,
    ##             row.names = FALSE, file = file.path("data", "cellLineDataCCLE.csv"))
    ## write.table(clDataGDSC, sep = ",", quote = FALSE, col.names = TRUE,
    ##             row.names = FALSE, file = file.path("data", "cellLineDataGDSC.csv"))
}


## ##############################################################################
## Parse subset of drug metadata cleaned by Haibe-Kains et al. (2014)
## -- parse cell line metadata in CCLE, GDSC objects separately
## ##############################################################################

## to parse drug metadata set to TRUE
if (FALSE) { 
    ## parse drug metadata from CCLE
    drugMapCCLE <- common[["CCLE"]]@curation$drug
    drugMapCCLE <- dplyr::mutate_if(drugMapCCLE, is.factor, as.character)
    drugDataCCLE <- common[["CCLE"]]@drug
    drugDataCCLE <- dplyr::left_join(drugDataCCLE, drugMapCCLE, by = c("drug.name" = "CCLE.drugid"))
    drugDataCCLE <- dplyr::select(drugDataCCLE,
                                  drug = unique.drugid,
                                  class = Class)

    ## parse drug metadata from GDSC
    drugMapGDSC <- common[["GDSC"]]@curation$drug
    drugMapGDSC <- dplyr::mutate_if(drugMapGDSC, is.factor, as.character)
    drugDataGDSC <- common[["GDSC"]]@drug
    drugDataGDSC <- dplyr::left_join(drugDataGDSC, drugMapGDSC, by = c("drug.name" = "GDSC.drugid"))
    drugDataGDSC <- dplyr::select(drugDataGDSC,
                                  drug = unique.drugid,
                                  type = Drug.type,
                                  targetted = Drug.class.II)

    ## merge drug data
    drugData <- dplyr::full_join(drugDataGDSC, drugDataCCLE, by = "drug")

    ## save as RDS file
    saveRDS(drugData, file = file.path("data", "drugData.rds"))

    ## ## alernatively, save as csv files
    ## write.table(drugData, sep = ",", quote = FALSE, col.names = TRUE,
    ##             row.names = FALSE, file = file.path("data", "drugData.csv"))
}






