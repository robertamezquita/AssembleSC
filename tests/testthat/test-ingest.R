## Testing ingest functions ---
test_that("Ingest sanity checks.", {
    ## dir
    expect_error(ingest(dir = NULL))    
    expect_error(ingest(dir = 'asdf'))
    ## aligner
    expect_error(ingest('..', aligner = NULL))
    expect_error(ingest('..', aligner = 'asdf'))
    ## mode
    expect_error(ingest('..', aligner = 'STARsolo', mode = 'asdf'))
    ## type
    expect_error(ingest('..', aligner = 'STARsolo', mode = 'standard', type = 'asdf'))
})
