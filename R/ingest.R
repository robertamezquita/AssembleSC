#' Ingest of scRNA-seq Alignment Results
#'
#' Reads in the resulting count matrices following alignment of scRNA-seq data.
#' Supported aligners are detailed below, as each has their own way
#' of organizing results.
#'
#' @details
#'
#' Below is a list of aligners that are supported or with planned support,
#' and details on how they are run and their respective organization of results:
#'
#' * kb-python (kallisto|bustools) (standard):
#'   * Run via kb-python workflow via `kb count` with `--filter bustools` enabled.
#'   * Alignment results are organized as follows:
#'     * `{sample}/{counts_filtered | counts_unfiltered}/{cells_x_genes.barcodes.txt, cells_x_genes.genes.txt, cells_x_genes.mtx}`
#'   * Alignment statistics are stored in: `{sample}/run_info.json`
#' * kb-python (kallisto|bustools) (--lamanno) [not supported yet]:
#'   * Run via kb-python workflow via `kb count` with `--lamanno` and `--filter bustools` flags enabled.
#'   * Alignment results are organized as follows:
#'     * `{sample}/{counts_filtered | counts_unfiltered}/{cells_x_genes.barcodes.txt, cells_x_genes.genes.txt, cells_x_genes.mtx}`
#'   * Alignment statistics are stored in: `{sample}/run_info.json`
#' * STARsolo:
#'   * Run via STAR in STARsolo mode with with the `--soloFeatures Gene GeneFull SJ Velocyto` and `--soloCellFilter` enabled.
#'   * Alignment results are organized as follows:
#'     * `{sample}/Solo.out/{mode}/{filtered | raw}/{barcodes.tsv, features.tsv, matrix.mtx}`
#'   * Alignment statistics are stored in: `{sample}/Log.final.out`
#'   * Note that running STARsolo with Velocyto feature output only provides a raw, unfiltered output. Not yet supported.
#' * kallisto:
#'   * Run via `kallisto quant`.
#'   * Alignment results are organized as follows:
#'     * `{sample}/{abundance.tsv | abundance.h5}`
#'   * Alignment statistics are stored in: `{sample}/run_info.json`
#'   * See `[ingest_kallisto()]` for a custom function designed for kallisto runs, since the input is
#'     quite different compared to STARsolo/kb-python.
#'
#' Another consideration is the sequencing capture strategy. The reason this is important is
#' for downstream normalization of the counts, specifically whether to adjust counts for transcript-length.
#' In many cases with scRNA-seq, the strategy is 3' tagged, which does not require any gene/transcript
#' length correction as the counts do not have a length bias. However, as is often the case for
#' plate-based protocols such as SMARTseq2 with full-length sequencing, then indeed the counts
#' should be adjusted. The parameter `capture` designates whether the strategy was '3prime_tagged' (default),
#' which does not perform length-correction, versus 'full_length', which does.
#'
#'
#' @param dir A string with the path to the sample directory where a scRNA-seq experiment results 
#'   from an aligner are stored. Required. For the `kallisto` aligner, must be the directory containing all
#'   cells (samples) from the scRNA-seq run.
#' @param aligner A string specifying a supported aligner from one of the following:
#'   'STARsolo', 'kb-python'. For 'kallisto', see `[ingest_kallisto()]` for more details. Required.
#' @param mode A string specifying the mode with which the aligner was run.
#'   'standard' refers to kallisto/kb-python default (and is the same as 'Gene' for
#'   `aligner = STARsolo`), 'lamanno' refers to a kb-python specific mode, and 'Gene',
#'   'GeneFull', 'and 'Velocyto' refer to STARsolo specific modes. Default is
#'   'standard'. 
#' @param type A string specifying whether to use the raw or aligner-filtered count
#'   matrices. Default is 'raw'. Not applicable for when aligner is 'kallisto'.
#' @param strategy A string specifying the sequencing protocol's quantification strategy. Default is '3prime_tagged'.
#' @param .id A string that serves to uniquify samples by prepending to the barcode.
#'   Default is NULL, which means no id is prepended to barcodes. If .id = 'auto',
#'   the .id string is deduced from the basename of the `dir` argument.
#' @param .sep The separator between `.id` and the droplet/cell barcode, default is a full stop.
#' @param tx2gene A two column data.frame where the first column is the transcript id, and the
#'   second column is the corresponding gene id. Required for when the aligner is set to kallisto.
#' 
#' @return An annotated `SingleCellExperiment` object.
#'
#' @importFrom rlang abort inform is_null
#' @importFrom fs dir_exists path path_file
#' @importFrom readr read_tsv
#' @importFrom Matrix readMM t
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#'
#' @seealso ingest_kallisto
#' 
#' @family ingest
#' 
#' @examples
#' dat.path <- system.file('extdata/aligners/STARsolo/Sample_A1/', package = 'AssembleSC')
#' 
#' ingest(dat.path,
#'        aligner = 'STARsolo', mode = 'Gene',
#'        type = 'filtered',
#'        .id = 'A1', .sep = '.')
#'
#' @export
ingest <- function(dir = NULL,
                   aligner = NULL,
                   mode = 'standard',
                   type = 'raw',
                   strategy = '3prime_tagged',
                   .id = NULL, .sep = '.',
                   tx2gene = NULL) {
    ## Suported params
    .aligners <- c('STARsolo', 'kb-python', 'kallisto')
    .modes <- c('standard', 'Gene', 'GeneFull', 'Velocyto', 'lamanno')
    .modes_STAR <- c('standard', 'Gene', 'GeneFull', 'Velocyto')
    .modes_kbp <- c('standard', 'lamanno')    
    .types <- c('raw', 'filtered')
    .strategies <- c('3prime_tagged', 'full_length')
    
    ## Sanity checks -----------------------------------------------------------
    ## Check if arg is supported
    ## Existence of files and if args are supported
    if (is_null(dir) || !dir_exists(dir)) {
        abort('Directory does not exist or not specified.')
    }
    if (is_null(aligner)) {
        abort('Please specify aligner.')
    }
    if (!(aligner %in% .aligners)) {
        abort('Specified aligner not supported.')
    }
    if (!is_null(mode) & !(mode %in% .modes)) {
        abort('Specified mode not supported.')
    }
    if (!is_null(type) & !(type %in% .types)) {
        abort('Type of output not supported.')
    }

    ## Checking if aligner + mode combo supported
    if (aligner == 'STARsolo' & !(mode %in% .modes_STAR)) {
        abort('Specified mode not supported with STARsolo aligner.')
    }
    if (aligner == 'kb-python' & !(mode %in% .modes_kbp)) {
        abort('Specified mode not supported with kb-python aligner.')
    }
    if (aligner == 'kallisto' & mode != 'standard') {
        abort('Specified mode not supported with kallisto aligner.')
    }

    ## Check strategy
    ## TODO: Check if full-length strategy affects kb-python/STARsolo..
    if (!strategy %in% .strategies) {
        abort('Specified strategy not supported.')
    }

    ## Check that tx2gene is provided for kallisto
    if (aligner == 'kallisto' & is_null(tx2gene)) {
        abort('A data.frame mapping transcripts-to-gene must be provided in `tx2gene` parameter.')
    }
    
    ## Construct path to outs per aligner --------------------------------------
    ## NOTE: files always arranged in order: barcodes, genes, mtx

    ## kb-python -------
    ##   + standard ---
    if (aligner == 'kb-python' & mode == 'standard') {
        .suffix <- paste0('cells_x_genes', c('.barcodes.txt', '.genes.txt', '.mtx'))

        if (mode == 'standard') {
            if (type == 'raw') {
                outs.path <- path(dir, 'counts_unfiltered', .suffix)
            }
            if (type == 'filtered') {
                outs.path <- path(dir, 'counts_filtered', .suffix)
            }
        }
    }

    ##   + lamanno ---
    ## TODO: add support for this..
    if (aligner == 'kb-python' & mode == 'lamanno') {
        .suffix <- paste0(rep(c('spliced', 'unspliced'), each = 3),
                          c('.barcodes.txt', '.genes.txt', '.mtx'))
        if (type == 'raw') {
            outs.path <- path(dir, 'counts_unfiltered', .suffix)
        }
        if (type == 'filtered') {
            outs.path <- path(dir, 'counts_filtered', .suffix)
        }

        abort('kb-python count run with --lamanno not yet supported.')
    }

    ## STARsolo --------
    ## - logic: switch modes/aligners/types automagically
    ##   - concatenate together at end since they nicely form `outs.path`
    if (aligner == 'STARsolo' & mode == 'standard') {
        inform('STARsolo with mode set to `standard` is equivalent to `Gene`, setting mode to `Gene` ..')
        mode <- 'Gene'
    }

    ## TODO: add support for STARsolo + Velocyto
    if (aligner == 'STARsolo' & mode == 'Velocyto') {
        abort('STARsolo run in Velocyto mode not yet supported.')
    }
    
    if (aligner == 'STARsolo' & mode == 'Velocyto' & type == 'filtered') {
        inform('STARsolo with Velocyto only provides raw output, switching type to `filtered` ..')
        type <- 'raw'
    }

    if (aligner == 'STARsolo') {
        .suffix <- c('barcodes.tsv', 'features.tsv', 'matrix.mtx')
        outs.path <- path(dir, 'Solo.out', mode, type, .suffix)
    }

    ## kallisto --------
    ## TODO: very different input and downstream processing..not supporting for now
    if (aligner == 'kallisto') {
        ##        abort('Standalone kallisto aligner output not yet supported.')
        inform('Running `ingest_kallisto()` under the hood...')        
        obj <- ingest_kallisto(dir.all = dir, .id = .id,
                               strategy = strategy,
                               tx2gene = tx2gene)
        return(obj)
    }


    ## Read matrices -----------------------------------------------------------
    ## Check if .id is set to auto mode
    if (.id == 'auto') {
        .id <- path_file(dir)
    }
    
    .construct_sce <- function(outs.path, .id = .id, .sep = .sep) {
        ## Read in primary data ------------------------------------------------
        ## Construct barcodes +/- .id_.sep_.barcode
        barcodes <- read_tsv(outs.path[1],
                             col_names = FALSE, col_types = 'c',
                             progress = FALSE)
        colnames(barcodes)[1] <- 'barcode'

        if (!is_null(.id)) {
            barcodes$.id <- .id
            barcodes$barcode <- paste0(.id, .sep, barcodes$barcode)
        }

        ## NOTE: overloading col_types in case >1 gene identifiers are present
        ## only the first one is used though as rowname
        features <- suppressWarnings(
            readr::read_tsv(outs.path[2],
                            col_names = FALSE, col_types = 'ccc',
                            progress = FALSE)
        )
        colnames(features)[1] <- '.row_id'

        mtx <- Matrix::readMM(outs.path[3])

        ## Sanity checks on features/barcodes matching matrix dims -------------
        ## Append features/barcodes as rownames/colnames
        ## NOTE: this may error if features is not entirely unique..
        ## - highly depends on being ENSEMBL identifiers or the like
        ## .. first check barcodes..
        if (nrow(mtx) == nrow(barcodes)) {
            inform('Transposing matrix to features (rows) x barcodes (cols).')
            mtx <- t(mtx) # transpose since its likely cells x genes
        } else if (ncol(mtx) != nrow(barcodes)) {
            abort(paste('Number of unique barcodes does not match either dimension',
                        'of the counts matrix. Check results.'))
        }

        ## .. next check features..
        if (nrow(mtx) != nrow(features)) {
            ## check for uniqueness of first identifier (column) of feature
            abort(paste('Number of unique features does not equal feature space of counts matrix.',
                        'Check feature annotation'))
        }

        ## Annotate matrix object ----------------------------------------------
        mtx <- as.matrix(mtx) # TODO: fix sparse support..?
        rownames(mtx) <- features[, 1, drop = TRUE]
        colnames(mtx) <- barcodes$barcode

        ## Combine into SCE ----------------------------------------------------
        SingleCellExperiment(assays = list(counts = mtx),
                             rowData = features,
                             colData = barcodes)
    }

    .construct_sce(outs.path, .id, .sep)
}



#' Ingest kallisto output
#'
#' Shortcircuits the regular `ingest()` function to return a `SingleCellExperiment`
#' object derived from regular kallisto output. This is used typically for processing
#' plate-based data, where each well corresponds to a single cell, and each cell
#' is in its own directory, e.g. sample folder.
#'
#' @details
#'
#' This function utilizes under the hood the `tximport::tximport()` function to
#' read in the counts data, performing scaling, gene summarisation, and the like
#' as specified. For more details, see `?tximport::tximport()`.
#'
#' Note that depending on the sequencing protocol, different counts-normalization
#' strategies may need to be used. For example, for full-length 5'->3' sequencing,
#' gene/transcript length normalization is a common step, as opposed to 3' targeted
#' strategy which does not require this step.
#'
#' For 3' tagged RNA-seq, correcting the counts for gene length will induce a bias
#' in your analysis, because the counts do not have a length bias. Thus, it is
#' recommended to use the original counts, as a counts matrix.
#' For full-length sequencing, its recommended to set `countsFromAbundance` to
#' 
#' For more information, see:
#' [Downstream DGE in Bioconductor from tximport](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#downstream_dge_in_bioconductor) and [Swimming downstream](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6178912/).
#'
#' @param dir.all A string with the path pointing to the directory *above* the individual sample directories,
#'   the one containing *all* samples (cells) from the run; this is due to between-sample
#'   normalization that needs to be performed prior to analysis. Required.
#' @param .id A character vector specifying the sample ids. Default is 'auto', which deduces the
#'   sample id based on the folder names within the directory specified by `dir.all`.
#' @param strategy A string specifying the sequencing protocol's strategy, used to instruct
#'   the normalization approach. Default is '3prime_tagged'.
#' @param tx2gene A two-column data.frame linking transcript id (column 1) to gene id (column 2).
#'   Required for gene-level summarisation, when `txOut = FALSE`.
#'
#' @importFrom fs path path_file dir_ls
#' @importFrom rlang abort is_null
#' @importFrom readr read_tsv
#' @import tximport
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' 
#' @family ingest
#'
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', 'aligners', 'kallisto', 'Sample_A1', package = 'AssembleSC')
#' t2g.file <- system.file('extdata', 'metadata', 'tx2gene_human.tsv', package = 'AssembleSC')
#' t2g <- readr::read_tsv(t2g.file, col_names = FALSE)
#'
#' ingest_kallisto(dir)
#' }
#'
#' @export
ingest_kallisto <- function(dir.all, .id = 'auto',
                            strategy = 'full_length',
                            tx2gene = NULL) {

    .strategies <- c('full_length', '3prime_tagged')

    ## Sanity checks -----------------------------------------------------------
    ## TODO: add support for transcripts out
    if (is_null(tx2gene)) { #& txOut == FALSE) {
        abort('For summarizing to the gene-level, `tx2gene` is required.`')
    }

    if (is_null(strategy) | !(strategy %in% .strategies)) {
        abort('For ingesting kallisto results, parameter `strategy` must be specified.')
    }

    ## Construct paths ----------------------------------------------------------
    folders <- dir_ls(dir.all)
    outs.path <- path(dir_ls(dir.all), 'abundance.h5')

    if (is_null(.id) | .id == 'auto') {
        .id <- path_file(folders)
    } else if (!is_null(.id) & (length(.id) != length(folders))) {
        abort('Provided `.id` not the same length as number of sample folders in directory. Check input.')
    }

    ## Settings based on strategy strategy --------------------------------------
    if (strategy == 'full_length') {
        .countsFromAbundance <- 'lengthScaledTPM'
    }

    if (strategy == '3prime_tagged') {
        .countsFromAbundance <- 'no'
    }
    
    
    ## Importation --------------------------------------------------------------
    txi <- tximport(outs.path, type = 'kallisto',
                    txIn = TRUE, txOut = FALSE,
                    tx2gene = tx2gene,
                    countsFromAbundance = .countsFromAbundance)
#suppressMessages(#    )


    ## Adapted from: bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
    cts <- round(txi$counts)
    cts <- cts[tx2gene[, 2, drop = TRUE], ] # reorder to be the same 


    ## Organize assays; add eff.length if full-length to allow downstream offset calc
    .assays <- list(counts = cts)
    if (strategy == 'full_length') {
        eff.length <- txi$length
        eff.length <- eff.length[tx2gene[, 2, drop = TRUE], ] # reorder to be the same         
        .assays <- c(.assays, list(eff.length = eff.length))
    }

    ## Create coldata & rowdata
    coldat <- data.frame(.id = .id)
    rowdat <- tx2gene

    ## Create sce with extra info
    sce <- SingleCellExperiment(assays = .assays,
                                colData = coldat,
                                rowData = rowdat)
    colnames(sce) <- .id
    colnames(rowData(sce))[1:2] <- c('transcript_id', 'gene_id')
    metadata(sce)$countsFromAbundance <- .countsFromAbundance

    return(sce)
}
    

    
    ## dir.all = '/fh/fast/gottardo_r/ramezqui_working/analysis/_HVTN/data-raw/kallisto_new_v2'
    ## outs.path <- outs.path[1:5]
    ## .id <- .id[1:5]

    ## if (strategy == 'full_length') {
    ##     eff.length.mat <- txi$length

    ##     ## Add sample/cell id
    ##     colnames(cts) <- colnames(eff.length.mat) <- .id

    ##     ## Obtaining per-observation scaling factors for length, adjusted to avoid
    ##     ## changing the magnitude of the counts.
    ##     normMat <- eff.length.mat / exp(rowMeans(log(eff.length.mat)))
    ##     normCts <- cts / normMat
        
    ##     ## Computing effective library sizes from scaled counts, to account for
    ##     ## composition biases between samples.
    ##     eff.lib <- edgeR::calcNormFactors(normCts, method = norm.method) * colSums(normCts)
        
    ##     ## Combining effective library sizes with the length factors, and calculating
    ##     ## offsets for a log-link GLM.
    ##     normMat.ll <- sweep(normMat, 2, eff.lib, "*")
    ##     normMat.ll <- log(normMat.ll)
        
    ##     ## Creating a DGEList object for offset 
    ##     y <- edgeR::DGEList(cts)
    ##     y <- edgeR::scaleOffset(y, normMat.ll)

    ##     ## Create SCE
    ##     coldat <- data.frame(.id = .id)
    ##     sce <- SingleCellExperiment(assays = list(counts = y$counts,
    ##                                               offset = y$offset,
    ##                                               eff.length = eff.length.mat),
    ##                                 colData = coldat)
    ##     sce$lib.size <- as.integer(y$samples$lib.size)
    ##     sce$eff.lib.size <- as.integer(eff.lib)
        
    ##     ## Add CPM/logCPM (adapted from csaw::calculateCPM)
    ##     assay(sce, 'cpm') <- counts(sce) / exp(assay(sce, 'offset')) * 1e6
        
    ##     ## Add log-CPM (adapted from csaw::calculateCPM)
    ##     ap <- edgeR::addPriorCount(counts(sce), offset = assay(sce, 'offset'), prior.count = 1)
    ##     log.cpm <- log2(ap$y) - as.matrix(ap$offset) / log(2) + log2(1e6)
    ##     assay(sce, 'logcpm') <- log.cpm
    ## }        

