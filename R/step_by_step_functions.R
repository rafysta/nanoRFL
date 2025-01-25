#!/usr/bin/Rscript
# Step by Step run of nanoRF packages

#' @title loading Proteomics data
#' @param rf nanoRF object
#' @param file proteomics file path
#' @param column.ratio index of fold-change in input file
#' @param column.target index of important column. at least two columns from these should have values
#' @return nanoRF object
#' @export
#' @examples
#' # rf <- loadProteomicsData(rf, file = FILE_Proteomics_data, column.ratio = c(8:12), column.target = c(8, 9, 11, 12))
loadProteomicsData <- function(rf, file, column.ratio = NULL, column.target = NULL){
  package_names <- c("randomForest", "dplyr", "data.table", "ggplot2", "cowplot", "readxl")
  for(package_name in package_names){
    if (!require(package_name, character.only = TRUE)) {
      install.packages(package_name)
    }
  }
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressWarnings(suppressMessages(library(data.table)))
  library(randomForest)

  rf <- list()
  file_extension <- tools::file_ext(file)

  # 拡張子に応じて処理を分岐
  if (file_extension == "xlsx") {
    protein.groups.raw <-  readxl::read_xlsx(file)
  } else if (file_extension == "txt") {
    protein.groups.raw <-  fread(file)
  } else {
    cat("Unsupported file type.")
    return()
  }

  protein.groups <- protein.groups.raw[rowSums(!is.na(protein.groups.raw[, column.target])) >= 2, ]

  rf[["column.ratio"]] <- column.ratio
  rf[["protein.groups"]] <- as.data.frame(protein.groups)
  rf
}

#' @title load ID file
#' @param rf nanoRF object
#' @param file path to the ID file
#' @return nanoRF object
#' @export
#' @examples
#' # rf <- loadTrainingSet(rf, file = "H:/work/Project/061_20250124_Ohta_program/out/2025-01-24_check_original/load_training_set_ids20250110_T.R")
loadIDfile <- function(rf, file){
  temp_env <- new.env()
  sys.source(file, envir = temp_env)

  all_vars <- ls(envir = temp_env)
  raw_ids_vars <- grep("\\.raw\\.ids$", all_vars, value = TRUE)
  negative_vars <- c("cytoskelton.raw.ids", "mitochondria.raw.ids", "membrane.raw.ids")
  raw_ids_vars <- setdiff(raw_ids_vars, negative_vars)

  rf[["all_vars"]] <- all_vars
  rf[["negative_vars"]] <- negative_vars
  rf[["target_vars"]] <- raw_ids_vars
  rf[["env"]] <- temp_env
  return(rf)
}

#' @title check ID conditions
#' @param rf nanoRF object
#' @export
#' @examples
#' # checkIDcondition(rf)
checkIDcondition <- function(rf){
  cat("All IDs:\n")
  cat(paste(rf[["all_vars"]], collapse = ", "))
  cat("\n\n")

  cat("Negative vars:\n")
  cat(paste(rf[["negative_vars"]], collapse = ", "))
  cat("\n\n")

  cat("Target vars:\n")
  cat(paste(rf[["target_vars"]], collapse = ", "))
  cat("\n")
}

#' @title change negative ID list
#' @param rf nanoRF object
#' @param ids vector of negative ID list
#' @export
#' @examples
#' # rf <- setNegative_vars(rf, ids = c("cytoskelton.raw.ids", "mitochondria.raw.ids", "membrane.raw.ids"))
setNegative_vars <- function(rf, ids){
  rf[["negative_vars"]] <- ids
  rf
}

#' @title change target ID list
#' @param rf nanoRF object
#' @param ids vector of target IDs
#' @export
#' @examples
#' # rf <- setTarget_vars(rf, ids = c("HMGA.raw.ids", "macroH2A.raw.ids"))
setTarget_vars <- function(rf, ids){
  rf[["target_vars"]] <- ids
  rf
}

#' @title run nanoRF calculation
#' @param rf nanoRF object
#' @export
#' @examples
#' # rf <- runRF(rf)
runRF <- function(rf) {
  protein.groups <- rf[["protein.groups"]]
  column.ratio <- rf[["column.ratio"]]
  temp_env <- rf[["env"]]

  # negative data
  cytosol.ids <- c()
  for(ss in rf[["negative_vars"]]){
    ids <- get(ss, envir = temp_env)
    cytosol.ids <- c(cytosol.ids, ids)
  }
  cytosol.protIndices <- retrieve.proteins.identified(cytosol.ids, protein.groups)

  results_list <- list()
  cutoff_list <- list()
  for (raw_id_var in rf[["target_vars"]]) {
    prefix <- sub("\\.raw\\.ids$", "", raw_id_var)
    cat("==============================================================\n")
    cat(paste0("Processing: ", prefix, "\n"))
    raw_ids <- get(raw_id_var, envir = temp_env)
    protIndices <- retrieve.proteins.identified(raw_ids, protein.groups)
    tf <- as.factor(merge.vectors(positive = protIndices, negative = cytosol.protIndices))
    rf <- rf.workflow(cbind(tf, protein.groups[, column.ratio]))
    results_list[[prefix]] <- rf$PREDICTORS.SCORES[, c("training.factor", "T")]
    cutoff_list[[prefix]] <- rf$CUTOFF
  }

  resultsTable <- cbind(
    protein.groups[, c(1:3)],
    do.call(cbind, results_list)
  )

  colnamevec <- c(
    "GeneIDs", "proteinGroupID", "Gene_names",
    unlist(lapply(names(results_list), function(prefix) {
      c(paste0(prefix, ".training.factor"), paste0(prefix, ".RF.score"))
    }))
  )
  colnames(resultsTable) <- colnamevec
  rf[["resultsTable"]] <- resultsTable
  rf[["cutoff_list"]] <- cutoff_list
  return(rf)
}


#' @title output table
#' @param rf nanoRF object
#' @param file output file
#' @export
#' @examples
#' # outputTable(rf, file = "RF_score_result.txt")
outputTable <- function(rf, file){
  write.table(rf[["resultsTable"]], file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

#' @title output cutoff
#' @param rf nanoRF object
#' @param file output file
#' @export
#' @examples
#' # outputCutoff(rf, file = "RF_cut_off.txt")
outputCutOff <- function(rf, file){
  cutoff_list <- rf[["cutoff_list"]]
  df <- data.frame(sample = names(cutoff_list), value = unlist(cutoff_list))
  write.table(df, file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

#' @title plot prediction graphs
#' @param rf nanoRF object
#' @param output_dir output directory. path should have
#' @export
#' @examples
#' # plotPrediction(rf, output_dir = "output_directory/")
plotPrediction <- function(rf, output_dir){
  suppressWarnings(suppressMessages(library(ggplot2)))
  suppressWarnings(suppressMessages(library(cowplot)))

  if (!grepl("/$", output_dir)) {
    output_dir <- paste0(output_dir, "/")
  }

  data <- rf[["resultsTable"]]
  cutoff_list <- rf[["cutoff_list"]]

  for(target in names(cutoff_list)){
    score_col <- paste0(target, ".RF.score")
    factor_col <- paste0(target, ".training.factor")
    threshold <- cutoff_list[[target]]
    data$color_group <- ifelse(data[[factor_col]] == "T", "T",
                               ifelse(data[[factor_col]] == "F", "F",
                                      ifelse(data[[score_col]] > threshold, "Candidate", "Other")))

    data <- data %>% dplyr::mutate(color_group = factor(color_group, levels = c("T", "F", "Candidate", "Other")))

    p <- ggplot(data, aes_string(x = score_col, y =1, color = "color_group")) +
      geom_jitter(size=1.5, position = position_jitter(height = 0.2, width = 0, seed = 42), data=subset(data, color_group == "Other")) +
      geom_jitter(size=1.5, position = position_jitter(height = 0.2, width = 0, seed = 42), data=subset(data, color_group == "Candidate")) +
      geom_jitter(size=2, position = position_jitter(height = 0.2, width = 0, seed = 42), data=subset(data, color_group %in% c("T", "F"))) +
      scale_color_manual(values = c("T" = "red", "F" = "blue", "Candidate" = "orange",  "Other" = "gray")) +
      scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, by=0.2)) +
      geom_vline(xintercept = threshold, linetype = "dashed") +
      labs(title = target, x = "RF score", y = NULL, color = NULL) +
      coord_cartesian(xlim=c(0,1)) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
      )

    save_plot(paste0(output_dir, "RF2_score_for_", target, ".png"), p, base_width = 6, base_height = 2, dpi= 200, unit = "in")
    save_plot(paste0(output_dir, "RF2_score_for_", target, ".pdf"), p, base_width = 6, base_height = 2, dpi= 200, unit = "in")
  }


}






