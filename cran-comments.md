## Resubmission
This is a resubmission.  Following CRAN's feedback in this 
version I have:

* updated the package DESCRIPTION field to capitalises only proper 
  nouns and avoid the use of abbreviations.
  
* reduced the size of StMoMoVigentte.pdf in inst/doc from 804KB to 406KB.

* updated the examples to reduce the computing time. 

* removed all \donttest occurrences in the examples.

## Test environments
* local Windows 7 install, R 3.2.0
* ubuntu 12.04 (on travis-ci), R 3.2.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Andres Villegas <andresmauriciovillegas@gmail.com>'
  New submission

  This is our first submission.
  
In addition there was the following comment raised by CHECK:

* Possibly mis-spelled words in DESCRIPTION:
  Dowd (11:26)
  
  Dowd is a surname and is spelled correctly.
