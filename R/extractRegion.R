#' @title \code{extractRegion}
#'
#' @description \code{extractRegion} will extract the coverage files created
#'   by callOpenTiles and return a specific region's coverage
#'
#' @param SampleTileObj The SummarizedExperiment object output from getSampleTileMatrix
#' @param region a GRanges object or vector or strings containing the regions on which to compute co-accessible links. Strings must be in the format "chr:start-end", e.g. "chr4:1300-2222".
#' @param cellPopulations vector of strings. Cell subsets for which to call
#'   peaks. This list of group names must be identical to names that appear in
#'   the SampleTileObj.  Optional, if cellPopulations='ALL', then peak
#'   calling is done on all cell populations. Default is 'ALL'.
#' @param groupColumn Optional, the column containing sample group labels for returning coverage within sample groups. Default is NULL, all samples will be used.
#' @param subGroups a list of subgroup(s) within the groupColumn from the metadata. Optional, default is NULL, all labels within groupColumn will be used.
#' @param sampleSpecific If TRUE, get a sample-specific count dataframe out. Default is FALSE, average across samples and get a dataframe out.
#' @param approxLimit Optional limit to region size, where if region is larger than approxLimit basepairs, binning will be used. Default is 100000.
#' @param binSize Optional, size of bins in basepairs when binning is used. Default is 250.
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'   multiple cell populations
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return countSE a SummarizedExperiment containing coverage for the given input cell populations.
#'
#' @examples
#' \dontrun{
#' countSE <- MOCHA::extractRegion(
#'   SampleTileObj = SampleTileMatrices,
#'   cellPopulations = "ALL",
#'   region = "chr1:18137866-38139912",
#'   numCores = 30,
#'   sampleSpecific = FALSE
#' )
#' }
#' 
#' @export
#'

extractRegion <- function(SampleTileObj,
                          region,
                          cellPopulations = "ALL",
                          groupColumn = NULL,
                          subGroups = NULL,
                          sampleSpecific = FALSE,
                          approxLimit = 100000,
                          binSize = 250,
                          numCores = 1,
                          verbose = FALSE) {
  . <- idx <- score <- NULL

  cellNames <- names(SummarizedExperiment::assays(SampleTileObj))
  metaFile <- SummarizedExperiment::colData(SampleTileObj)
  outDir <- SampleTileObj@metadata$Directory

  if (is.na(outDir)) {
    stop("Error: Missing coverage file directory.")
  }

  if (is.character(region)) {

    # Convert to GRanges
    regionGRanges <- MOCHA::StringsToGRanges(region)
  } else if (class(region)[1] == "GRanges") {
    regionGRanges <- region
  } else {
    stop("Wrong region input type. Input must either be a string, or a GRanges location.")
  }


  if (all(toupper(cellPopulations) == "ALL")) {
    cellPopulations <- cellNames
  }
  if (!all(cellPopulations %in% cellNames)) {
    stop("Some or all cell populations provided are not found.")
  }


  # Pull out a list of samples by group.

  if (!is.null(subGroups) & !is.null(groupColumn)) {
    # If the user defined a list of subgroup(s) within the groupColumn from the metadata, then it subsets to just those samples
    subSamples <- lapply(subGroups, function(x) metaFile[metaFile[, groupColumn] %in% x, "Sample"])
    names(subSamples) <- subGroups
  } else if (!is.null(groupColumn)) {

    # If no subGroup defined, then it'll form a list of samples across all labels within the groupColumn
    subGroups <- unique(metaFile[, groupColumn])

    subSamples <- lapply(subGroups, function(x) metaFile[metaFile[, groupColumn] %in% x, "Sample"])
  } else {

    # If neither groupColumn nor subGroup is defined, then it forms one list of all sample names
    subGroups <- "All"
    subSamples <- list("All" = metaFile[, "Sample"])
  }

  # Determine if binning is needed to simplify things
  if (GenomicRanges::end(regionGRanges) - GenomicRanges::start(regionGRanges) > approxLimit) {
    if (verbose) {
      message(stringr::str_interp("Size of region exceeds ${approxLimit}bp. Binning data over ${binSize}bp windows."))
    }
    binnedData <- regionGRanges %>%
      plyranges::tile_ranges(., binSize) %>%
      dplyr::mutate(idx = c(1:length(.)))
  }


  # Pull up the cell types of interest, and filter for samples and subset down to region of interest
  cellPopulation_Files <- lapply(cellPopulations, function(x) {
    originalCovGRanges <- readRDS(paste(outDir, "/", x, "_CoverageFiles.RDS", sep = ""))
    cellPopSubsampleCov <- parallel::mclapply(subSamples, function(y) {
      sampleNames <- names(originalCovGRanges)
      keepSamples <- grepl(paste(y, collapse = "|"), sampleNames)
      subSampleList <- originalCovGRanges[keepSamples]


      # If we don't want sample-specific, then average and get a dataframe out.
      # else, get a sample-specific count dataframe out.
      if (!sampleSpecific) {
        sampleCount <- sum(keepSamples)
        filterCounts <- lapply(subSampleList, function(z) {
          tmpGR <- plyranges::join_overlap_intersect(z, regionGRanges)

          if (GenomicRanges::end(regionGRanges) - GenomicRanges::start(regionGRanges) > approxLimit) {
            tmpGR <- plyranges::join_overlap_intersect(tmpGR, binnedData) %>%
              plyranges::group_by(idx) %>%
              plyranges::reduce_ranges(score = mean(score)) %>%
              dplyr::ungroup()
          }
          tmpGR
        })

        mergedCounts <- IRanges::stack(methods::as(filterCounts, "GRangesList")) %>%
          plyranges::join_overlap_intersect(regionGRanges) %>%
          plyranges::compute_coverage(weight = .$score / sampleCount) %>%
          plyranges::join_overlap_intersect(regionGRanges)
        if (GenomicRanges::end(regionGRanges) - GenomicRanges::start(regionGRanges) > approxLimit) {
          mergedCounts %>% plyranges::join_overlap_intersect(binnedData)
        } else {
          mergedCounts
        }
      } else {
        tmpCounts <- lapply(subSampleList, function(z) {
          tmpGR <- plyranges::join_overlap_intersect(z, regionGRanges)
          if (GenomicRanges::end(regionGRanges) - GenomicRanges::start(regionGRanges) > approxLimit) {
            tmpGR <- plyranges::join_overlap_intersect(tmpGR, binnedData) %>%
              plyranges::group_by(idx) %>%
              plyranges::reduce_ranges(score = mean(score)) %>%
              dplyr::ungroup()
          }
          tmpGR
        })

        unlist(tmpCounts)
      }
    }, mc.cores = numCores)
    names(cellPopSubsampleCov) <- subGroups
    cellPopSubsampleCov
  })
  names(cellPopulation_Files) <- cellPopulations
  allGroups <- unlist(cellPopulation_Files)

  ## Generate a data.frame for export.

  if (GenomicRanges::end(regionGRanges) - GenomicRanges::start(regionGRanges) > approxLimit) {
    allGroupsDF <- parallel::mclapply(seq_along(allGroups), function(x) {
      subGroupdf <- as.data.frame(allGroups[[x]])
      subGroupdf$Groups <- rep(names(allGroups)[x], length(allGroups[[x]]))

      covdf <- subGroupdf[, c("seqnames", "start", "score", "Groups")]
      colnames(covdf) <- c("chr", "Locus", "Counts", "Groups")

      covdf
    }, mc.cores = numCores)
  } else {
    allGroupsDF <- parallel::mclapply(seq_along(allGroups), function(x) {
      tmp <- allGroups[[x]] %>% plyranges::tile_ranges(width = 1)
      tmp$score <- allGroups[[x]]$score[tmp$partition]
      tmp$Groups <- rep(names(allGroups)[x], length(tmp))

      covdf <- as.data.frame(tmp)[, c("seqnames", "start", "score", "Groups")]
      colnames(covdf) <- c("chr", "Locus", "Counts", "Groups")

      covdf
    }, mc.cores = numCores)
  }

  names(allGroupsDF) <- names(allGroups)

  countSE <- SummarizedExperiment::SummarizedExperiment(allGroupsDF,
    metadata = SampleTileObj@metadata
  )

  return(countSE)
}
