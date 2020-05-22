#' Draw a kneeplot
#'
#' Creates the classic kneeplot.
#'
#' @param empty_drops.out The object returned by the `DropletUtils::emptyDrops()` function
#'   that has been run on a raw (unfiltered) counts matrix.
#' @param fdr.threshold An FDR threshold for which to colour putative cell-containing
#'   droplets versus putative empty droplets.
#' 
#' @return A ggplot2 object with plenty of customization.
#'
#' @import ggplot2
#' @importFrom DropletUtils emptyDrops
#'
#' @seealso [calc_preqc_metrics()] for key calculated preliminary metrics used in this plot.
#'
#' @examples
#' suppressMessages(library(SingleCellExperiment))
#' sce.path <- system.file('extdata', 'sce', 'raw.rds', package = 'AssembleSC')
#' sce <- readRDS(sce.path)
#' 
#' empty_drops.out <- DropletUtils::emptyDrops(counts(sce), lower = 500)
#'
#' draw_kneeplot(empty_drops.out, fdr.threshold = 0.001)
#' 
#' @export
draw_kneeplot <- function(empty_drops.out, fdr.threshold = 0.001) {
    ## Add extra metadata to emptydrops output
    e.df <- as.data.frame(empty_drops.out)
    e.df$putative.status <- ifelse(e.df$FDR > fdr.threshold, 'Empty (estimated)', 'Cell')
    e.df$putative.status <- ifelse(is.na(e.df$FDR), 'Empty (below set lower limit)', e.df$putative.status)
    e.df$rank <- rank(-e.df$Total, ties.method = 'first')

    ## Calc and grab key metrics
    metrics <- calc_preqc_metrics(e.df, fdr.threshold)
    knee <- metrics$knee
    infl <- metrics$inflection
    n.gt_knee <- metrics$n.gt_knee
    n.putative_cell <- metrics$n.putative_cell
    
    ## Plotter
    e.df %>%
        ggplot(aes(x = rank,
                   y = Total,
                   colour = putative.status)) +
        geom_point(size = 0.5) +
        geom_hline(yintercept = knee, linetype = 'dashed', colour = 'red') +
        geom_hline(yintercept = infl, linetype = 'dashed', colour = 'green') +
        geom_text(aes(x = 5, y = knee, label = paste('Knee:', knee)), nudge_y = 0.05, colour = 'red') +
        geom_text(aes(x = 5, y = infl, label = paste('Inflection:', infl)), nudge_y = 0.05, colour = 'green') +
        theme(legend.position = 'bottom',
              aspect.ratio = 1) +
        scale_y_log10() + scale_x_log10() +
        labs(caption = paste0('Inflection (dashed green line); Knee (dashed red line)\n',
                              'Droplets with # UMIs greater than knee: ', n.gt_knee, '\n',
                              'Droplets that are putative cells: ', n.putative_cell),
             x = 'Rank', y = 'Total UMIs')
}



#' Calculate preliminary quality control metrics
#'
#' Calculates some metrics that are nice to have to assess quality of an
#' scRNA-seq sample, based on the output of the `DropletUtils::emptyDrops()`
#' function.
#'
#' @inheritParams draw_kneeplot
#'
#' @return A list with the following metrics:
#'   * knee : the calculated knee point by `DropletUtils::barcodeRanks()`
#'   * inflection : the calculated inflection point by `DropletUtils::barcodeRanks()`
#'   * n.gt_knee : number of droplets with number of UMIs greater than the knee point
#'   * n.putative_cell : number of droplets that are putative cells based on
#'       specified FDR threshold
#'
#' @importFrom DropletUtils emptyDrops
#' 
#' @examples
#' suppressMessages(library(SingleCellExperiment))
#' sce.path <- system.file('extdata', 'sce', 'raw.rds', package = 'AssembleSC')
#' sce <- readRDS(sce.path)
#' 
#' empty_drops.out <- DropletUtils::emptyDrops(counts(sce), lower = 500)
#'
#' calc_preqc_metrics(empty_drops.out, fdr.threshold = 0.001)
#'
#' @export
calc_preqc_metrics <- function(empty_drops.out, fdr.threshold = 0.001) {
    e.df <- empty_drops.out
    
    ## Estimate the likely lower limit that had been set for barcodeRanks
    likely.lower <- max(e.df$Total[is.na(e.df$FDR)])    

    ## Ripped out from DropletUtils::barcodeRanks()
    ## Calculate the knee and inflection point
    .barcodeRanks_ripped <- function(totals, lower, df = 20, fit.bounds = NULL) {
        ##totals <- unname(colSums(m))
        o <- order(totals, decreasing = TRUE)
        stuff <- rle(totals[o])
        run.rank <- cumsum(stuff$lengths) - (stuff$lengths - 1)/2
        run.totals <- stuff$values
        keep <- run.totals > lower
        if (sum(keep) < 3) {
            stop("insufficient unique points for computing knee/inflection points")
        }
        y <- log10(run.totals[keep])
        x <- log10(run.rank[keep])
        d1n <- diff(y)/diff(x)
        right.edge <- which.min(d1n)
        left.edge <- which.max(d1n[seq_len(right.edge)])
        if (is.null(fit.bounds)) {
            new.keep <- left.edge:right.edge
        }
        else {
            new.keep <- y > log10(fit.bounds[1]) & y < log10(fit.bounds[2])
        }
        fit <- smooth.spline(x[new.keep], y[new.keep], df = df)
        d1 <- predict(fit, deriv = 1)$y
        d2 <- predict(fit, deriv = 2)$y
        curvature <- d2/(1 + d1^2)^1.5
        knee <- 10^(y[which.min(curvature)])
        inflection <- 10^(y[right.edge])
        return(list(knee = knee, inflection = inflection))
    }

    bro <- .barcodeRanks_ripped(e.df$Total, lower = likely.lower)
    knee <- bro$knee
    infl <- bro$inflection

    ## Calculate some basic metrics
    n.gt_knee <- sum(e.df$Total > knee, na.rm = TRUE)
    n.putative_cell <- sum(e.df$FDR < fdr.threshold, na.rm = TRUE)

    list(knee = knee,
         inflection = infl,
         n.gt_knee = n.gt_knee,
         n.putative_cell = n.putative_cell)
}
