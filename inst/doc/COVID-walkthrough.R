## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----setup--------------------------------------------------------------------
#  library(MOCHA)
#  library(ArchR)
#  library(TxDb.Hsapiens.UCSC.hg38.refGene)
#  library(org.Hs.eg.db)
#  library(BSgenome.Hsapiens.UCSC.hg19)

## -----------------------------------------------------------------------------
#  # You should substitute this with your own ArchR project.
#  # You must have completed cell labeling with your ArchR project.
#  
#  ArchRProj <- ArchR::loadArchRProject("/home/jupyter/FullCovid")
#  
#  metadata <- data.table::as.data.table(ArchR::getCellColData(ArchRProj))
#  studySignal <- median(metadata$nFrags)
#  
#  # Get metadata information at the sample level
#  lookup_table <- unique(
#    metadata[, c(
#      "Sample",
#      "COVID_status",
#      "Visit",
#      "days_since_symptoms"
#    ),
#    with = FALSE
#    ]
#  )
#  
#  # Subset to visit 1 and extract samples
#  samplesToKeep <- lookup_table$Sample[
#    lookup_table$Visit == "FH3 COVID-19 Visit 1" &
#      lookup_table$days_since_symptoms <= 15 |
#      is.na(lookup_table$days_since_symptoms)
#  ]
#  
#  # subset ArchR Project
#  idxSample <- BiocGenerics::which(ArchRProj$Sample %in% samplesToKeep)
#  cellsSample <- ArchRProj$cellNames[idxSample]
#  ArchRProj <- ArchRProj[cellsSample, ]
#  

## -----------------------------------------------------------------------------
#  # Parameters for calling open tiles.
#  cellPopLabel <- "CellSubsets"
#  cellPopulations <- c("CD16 Mono")
#  numCores <- 20

## -----------------------------------------------------------------------------
#  tileResults <- MOCHA::callOpenTiles(
#    ArchRProj,
#    cellPopLabel = cellPopLabel,
#    cellPopulations = cellPopulations,
#    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
#    Org = "org.Hs.eg.db",
#    numCores = numCores,
#    studySignal = studySignal,
#    outDir = tempdir()
#  )

## -----------------------------------------------------------------------------
#  # Computing the TSAM can take into account groupings of
#  # samples when determining consensus tiles.
#  # Our samples can be grouped by the metadata column 'COVID_status'
#  # into 'Positive' and 'Negative' groups.
#  # Since these groupings may have unique biology and we expect differences
#  # in accessibility, we want to compute consensus tiles on each
#  # group independently and take the union of consensus tiles from each group.
#  groupColumn <- "COVID_status"
#  
#  # We set the threshold to require a tile must be open in at least
#  # (0.2 * the number of samples in each group) samples to be
#  # retained
#  threshold <- 0.2
#  
#  # Alternatively, you can set the threshold to 0 to keep the union of
#  # all samples' open tiles.
#  # This is equivalent to setting a threshold that would retain
#  # tiles that are open in at least one sample.
#  
#  SampleTileMatrices <- MOCHA::getSampleTileMatrix(
#    tileResults,
#    cellPopulations = "CD16 Mono",
#    groupColumn = groupColumn,
#    threshold = threshold,
#    verbose = FALSE
#  )
#  

## -----------------------------------------------------------------------------
#  # This function can also take any GRanges object
#  # and add annotations to its metadata.
#  SampleTileMatricesAnnotated <- MOCHA::annotateTiles(SampleTileMatrices)
#  
#  # Load a curated motif set from library(chromVARmotifs)
#  # included with ArchR installation
#  data(human_pwms_v2)
#  SampleTileMatricesAnnotated <- MOCHA::addMotifSet(
#    SampleTileMatricesAnnotated,
#    pwms = human_pwms_v2,
#    w = 7 # weight parameter for motifmatchr
#  )

## -----------------------------------------------------------------------------
#  regionToPlot = "chr4:XXX-XXXX"
#  
#  countSE <- MOCHA::extractRegion(
#    SampleTileObj = SampleTileMatrices,
#    cellPopulations = "CD16 Mono",
#    region = regionToPlot,
#    groupColumn = "COVID_status",
#    numCores = numCores,
#    sampleSpecific = FALSE
#  )
#  dev.off()
#  pdf("ExamplePlot.pdf")
#  # Note that to show specific genes with the option' whichGene'
#  # you must have the package RMariaDB installed
#  MOCHA::plotRegion(countSE = countSE, whichGene = "MYD88")
#  dev.off()
#  

## -----------------------------------------------------------------------------
#  cellPopulation <- "CD16 Mono"
#  groupColumn <- "COVID_status"
#  foreground <- "Positive"
#  background <- "Negative"
#  
#  # Choose to output a GRanges or data.frame.
#  # Default is TRUE
#  outputGRanges <- TRUE
#  
#  # Optional: Standard output will display the number of tiles found
#  # below a false-discovery rate threshold.
#  # This parameter does not filter results and only affects the
#  # afforementioned message.
#  fdrToDisplay <- 0.2
#  
#  differentials <- MOCHA::getDifferentialAccessibleTiles(
#    SampleTileObj = SampleTileMatrices,
#    cellPopulation = cellPopulation,
#    groupColumn = groupColumn,
#    foreground = foreground,
#    background = background,
#    fdrToDisplay = fdrToDisplay,
#    outputGRanges = outputGRanges,
#    numCores = numCores
#  )
#  
#  # The output contains a GRanges with all tiles and their differential
#  # test results. We can filter by FDR to get our set of
#  # differentially accessible tiles:
#  
#  res = head(plyranges::filter(differentials, seqnames =='chr4' & FDR < 0.2))
#  

## -----------------------------------------------------------------------------
#  regions = res$Tile
#  
#  # Alternatively, define regions as a character vector
#  # of region strings in the format "chr:start-end"
#  # regions <- c(
#  # "chr4:7326500-7326999",
#  # "chr4:7327000-7327499",
#  # "chr4:7339500-7339999",
#  # "chr4:7344500-7344999"
#  # )
#  
#  links <- MOCHA::getCoAccessibleLinks(
#    SampleTileObj = SampleTileMatrices,
#    cellPopulation = cellPopulation,
#    regions = regions,
#    windowSize = 1 * 10^6,
#    numCores = numCores,
#    verbose = TRUE
#  )
#  
#  # Optionally filter these links by their absolute
#  # correlation - this output also adds the chromosome,
#  # start, and end site of each link to the table.
#  
#  MOCHA::filterCoAccessibleLinks(links, threshold = 0.4)

