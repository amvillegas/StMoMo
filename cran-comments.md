## Resubmission
This is a resubmission.  Following CRAN's feedback in this 
version I have:

* updated the package DESCRIPTION omitting the obvious fact that this is an 
  R package.
  
* updated the examples in the vignette to reduce the building time. 
  I have also removed unnecessary TeX packages in the preamble of
  the vignette which seemed to be taking much of the vignette building
  time.

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
