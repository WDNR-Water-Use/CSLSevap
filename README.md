# CSLSevap

CSLSevap is an R package containing multiple evaporation functions useful to the Wisconsin Department of Natural Resources Central Sands Lakes Study (CSLS). Included evapo(transpir)ation approaches include:

 * FAO Penman-Monteith Reference Evapotranspiration
 * McJannet lake evaporation
 * Unmodified Hamond lake evaporation

For most uses, it is only necessary to interact with the "evaporation" function. However all functions are made accessible to the user since they may be useful in other contexts (e.g., for saturation vapor pressure, dew point temperature).

## Installation

To install or explore this R package:

  1. Use `devtools` to install the package. This will allow you to
  use functions and explore vignettes.  
  ```
  devtools::install_github("WDNR-Water-Use/CSLSdata", build_vignettes=T)
  ```
  
  2. Fork from github, clone or download ZIP, then open the CSLSevap.prj file in
  R. This will allow you to explore the source code more easily.
