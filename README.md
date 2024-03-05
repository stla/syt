# `syt`: Young Tableaux

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/syt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/syt/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

___

## Features

- Generation of all SYTx of a given shape

- Count of all SYTx of a given shape

- Uniform sampling of SYTx

- Robinson-Schensted correspondence

- Conversion to and from ballot sequences

- Conversion to and from paths on the Young graph

- Plancherel growth process

![](http://stla.github.io/stlapblog/posts/assets/img/young_yng_path.png)

```
> gprocess2syt(list(1, 2, c(2,1), c(2,2), c(2,2,1)))
[[1]]
[1] 1 2

[[2]]
[1] 3 4

[[3]]
[1] 5
```
