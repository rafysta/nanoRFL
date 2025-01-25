
# nanoRFL (nanoRF light)

<!-- badges: start -->
<!-- badges: end -->

nanoRFLは、nanoRF ( https://github.com/EarnshawLab/nanoRF )を簡単に使うためのwrapper packageです。

Montaño-Gutierrez, L. F., Ohta, S., Kustatscher, G., Earnshaw, W. C. & Rappsilber, J. Nano Random Forests to mine protein complexes and their relationships in quantitative proteomics data. Mol Biol Cell 28, 673–680 (2017).


## Installation
パッケージはRstudioのConsoleで以下のコマンドをタイプすることでインストールされます。

``` r
install.packages("remotes")
remotes::install_github("rafysta/nanoRFL")
```

## 解析の例
Rstudio上で、以下のコマンドを実行してください。

### ライブラリーの読み込み
``` r
library(nanoRFL)
```

### データの読み込み
``` r
FILE_Proteomics_data <- "Proteomics_RNAseq_data_for_nanoRF.xlsx"
rf <- loadProteomicsData(rf, file = FILE_Proteomics_data, column.ratio = c(8:12), column.target = c(8, 9, 11, 12))
```

### Training data setのIDの読み込み
``` r
rf <- loadIDfile(rf, file = "load_training_set_ids20250110.R")
```

### 読み込んだIDlistを表示
``` r
checkIDcondition(rf)
```

### Negative training setのIDを変更
``` r
rf <- setNegative_vars(rf, ids = c("cytoskelton.raw.ids", "mitochondria.raw.ids", "membrane.raw.ids"))
```

### Target setのIDを変更
``` r
rf <- setTarget_vars(rf, ids = c("HMGA.raw.ids", "macroH2A.raw.ids"))
```

### 変更後のIDlistを確認
``` r
checkIDcondition(rf)
```

### run nanoRF
``` r
rf <- runRF(rf)
```

### output table as text file
``` r
outputTable(rf, file = "RF_score_result.txt")
```

### output cut off value
``` r
outputCutOff(rf, file = "RF_score_cutoff.txt")
```

### 予測結果のプロット
``` r
plotPrediction(rf, output_dir = "out/2025-01-24_test/img/")
```


