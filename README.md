# CSLSevap

CSLSevap is an R package containing multiple evaporation functions useful to the Wisconsin Department of Natural Resources Central Sands Lakes Study (CSLS). Included evapo(transpir)ation approaches include:

 * FAO Penman-Monteith Reference Evapotranspiration
 * McJannet lake evaporation
 * Unmodified Hamond lake evaporation

For most uses, it is only necessary to interact with the "evaporation" function. However all functions are made accessible to the user since they may be useful in other contexts (e.g., for saturation vapor pressure, dew point temperature).

## Installation

To install from within the WI DNR (assumes you're using R Studio and are familiar with installing other packages):

 1. Navigate to the CSLSevap package by opening `CSLSevap.proj` from the CSLSevap folder (Z:\WQWT_PROJECTS\WY_MS_Monitoring\2018 Summer Field Crew\Central Sands Lake Study\Analysis_2019\CSLSevap). This automatically sets your working directory to the CSLSevap folder.

    **Note:** If you want to be able to knit an Rmd file that is saved within this package, you must have the shared network mapped to a letter drive (e.g., "Z:/" for "//central/water/") and you must open the .proj workspace (or set your working directory) via this letter drive. If you don't, knitting will fail with a less-than-helpful error message (that doesn't mention anything about the letter drive requirement).
 2. While your working directory is set to the top level of `CSLSevap`, use the R package `devtools` to install locally.
    ```R
	install.packages("devtools")
	devtools::install()
	```
Potential issues you may run into:

  * You are asked if you would like to update packages, you say "yes" (to some or all) and then it throws a weird error (e.g., "package 'CSLSevap' is not available (for R version 3.6.1)")
	  * Restart R, and try again not updating any packages
      * Restart R, and try updating packages independently of devtools (e.g., in the "Packages" tab of RStudio, click the green icon to "Update", "Select All", and "Install Updates"
	  * Restart R between troubleshooting attempts and whenever you run into errors - may take a few restarts, or a few attempts at updating packages.
	
## Package Organization

This package is organized based on [R package conventions](http://r-pkgs.had.co.nz/) and includes the following components:

 * Code (`R/`)
 * Package metadata (`DESCRIPTION`)
 * Code and data documentation (`man/`)
 * Tests (`tests/`)
 * Vignettes (`vignettes/`)
 * Knitted .Rmd documents and other documentation (`doc/`)

From within RStudio:

 * For a full list of available functions, use `help(package = CSLSevap)`
 * For help on a specific function, use `?function_name` as usual
