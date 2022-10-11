# GPcluster
This R package allows the estimation and prediction for a clustered Gaussian process model, which can deal with the nonstationarity issue and computational issue for a standard Gaussian process model. The estimation methods for the unknown parameters are based on a stochastic EM algorithm. More details can be seen in [Sung, Haaland, Hwang, and Lu (2023)](https://www3.stat.sinica.edu.tw/ss_newpaper/SS-2020-0456_na.pdf).

You can install the package using `install_github` function as follows:
``` r
library(devtools)
install_github("ChihLi/GPcluster")
```
