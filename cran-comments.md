## Test environments
* local Windows 7 install, R 3.2
* ubuntu 12.04 (on travis-ci), R 3.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTES:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Andres Villegas <andresmauriciovillegas@gmail.com>'
  New submission

  This is our first submission.

* checking sizes of PDF files under 'inst/doc' ... NOTE
  'qpdf' made some significant size reductions:
     compacted 'StMoMoVignette.pdf' from 804Kb to 684Kb
  consider running tools::compactPDF() on these files
  
  Vignette PDF is as small as possible.