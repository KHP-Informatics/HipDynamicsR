# HipDynamics

### Summary

*HipDynamics* is a R based downstream live-imaging analysis for the [Human Induced Pluripotent Stem Cell Initiative (HIPSCI)][1]. The analysis allows quantification of cell population dynamics after preliminary image analysis of live images by [CellProfiler][2] and incorporation of `plate_result.txt` files generated from standardised end-point assays by the [Operetta][3].


#### Installation
*HipDynamics* comes with a preconfigured directory structure and settings. In order to install *HipDynamics* simply pull this git-repository and everything should be set up for you. 


#### Settings

*HipDynamics* holds all configurable variables within the `hipDynamics.R` file. 

1. The database login-in- and end-point details need to be specified. 
2. The `plate_result.txt` files directory path should be assigned to `path_PR`.
3. Optional settings allows specifying different analysis-methods and data exports.


#### Go-time

In order to run *HipDynamics* you need to enter the R-shell type the following command (provided you are within *HipDynamics*' directory folder and all settings were specified):

```R
source("hipDynamics.R")
```

[1]: http://www.hipsci.org
[2]: http://www.cellprofiler.org
[3]: http://www.perkinelmer.co.uk/pages/020/cellularimaging/products/operetta.xhtml
