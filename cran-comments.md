## Test environments
* local R installation, R 4.0.0
* ubuntu 16.04 (on travis-ci), R 4.0.0
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Revised submission

* Added single quotes to the names of software in DESCRIPTION

* Added a new dataset on coral reefs from an ecology study

* Adjusted the multivariate response stan programs to designate the modified cholesky parameters as local

* Moved beta and u parameters from generated quantities to transformed parameters in all stan files

* Added additional set of tests to cover all Stan files, reduced time of other tests. 

* Revised the mvcorrplot example per submission comments, removing \dontrun