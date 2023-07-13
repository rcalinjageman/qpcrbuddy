# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!  This is pcrbuddy for the SlugLab")
}


process_directory <- function(
    my_directory = NULL
  ) {


  if (is.null(my_directory)) {
    print("Choose a directory of qPCR export files in xls format.  Be sure each has a Results sheet.")
    my_directory <- choose.dir()
  }

  print(my_directory)

  files <- list.files(my_directory, pattern = "xls", full.names=TRUE)
  tbl_output <- NULL

  for (my_file in files) {
    tbl_output <- rbind(
      tbl_output,
      do_qpcr(my_file)
    )
  }

  print(tbl_output)

  new_name <- gsub(
    "\\\\",
    "-",
    my_directory,
  )

  new_name <- gsub(
    ":",
    "-",
    new_name,
  )

  new_name <- paste(
    "qprcbuddy_",
    new_name,
    ".csv",
    sep = ""
  )

  write.csv(tbl_output, new_name)

}


process_file <- function(
    my_file = NULL
  ) {


  if (is.null(my_file)) {
    print("Choose a qPCR export file in xls format.  Be sure it has a Results sheet.")
    my_file <- file.choose()
  }

  tbl_output <- do_qpcr(my_file)

  print(tbl_output)

  new_name <- paste(
    tools::file_path_sans_ext(my_file),
    "_processed.csv",
    sep = ""
  )


  write.csv(tbl_output, new_name)

}


do_qpcr <- function(my_file) {

  tbl_output <- NULL


  # Find the start of the data
  print(
    glue::glue(
      "Finding start of data in {my_file}"
    )
  )
  my_row <- 1

  value <- suppressMessages(
    read_xls(
      path = my_file,
      sheet = "Results",
      range = paste("A", my_row, ":A", my_row, sep = ""),
      col_names = FALSE,
      progress = FALSE
    )
  )
  value <- value[[1]]

  if (is.null(value)) stop(
    glue::glue(
      "Couldn't find Results tab in your file, {my_file}."
    )
  )

  while (value != "Well" & value != "Well Position") {
    my_row <- my_row + 1

    value <- suppressMessages(
      read_xls(
        path = my_file,
        sheet = "Results",
        range = paste("A", my_row, ":A", my_row, sep = ""),
        col_names = FALSE, progress = FALSE
      )
    )

    if (nrow(value) > 0) value <- value[[1]] else value = "Empty"

    if (my_row > 1000 | value == "Sample Name") {
      print(
        glue::glue(
          "Couldn't find column 'Well' or 'Well Position' in the Results tab in your file, {my_file}.  Searched to row {my_row}."
        )
      )
      return(NULL)
    }

  }

  start_row <- my_row
  print(
    glue::glue(
      "Found start of results table at {start_row} of {my_file}."
    )
  )

  while (value != "Empty") {
    my_row <- my_row + 1

    value <- suppressMessages(
      read_xls(
        path = my_file,
        sheet = "Results",
        range = paste("A", my_row, ":A", my_row, sep = ""),
        col_names = FALSE, progress = FALSE
      )
    )

    if (nrow(value) > 0) value <- value[[1]] else value = "Empty"

  }

  end_row <- my_row - 1
  print(
    glue::glue(
      "Found end of results table at {end_row} of {my_file}"
    )
  )

  my_data <- suppressMessages(
    read_xls(
      path = my_file,
      sheet = "Results",
      range = paste("A", start_row, ":AL", end_row, sep = ""),
      col_names = TRUE
    )
  )

  print(head(my_data))

  # Clean data
  if ("Cq" %in% colnames(my_data)) {
    my_data$CT <- my_data$Cq
  }

  if ("Undetermined" %in% my_data$CT) {
    my_data[my_data$CT == "Undetermined", ]$CT <- "40"
  }

  my_data$CT <- as.numeric(my_data$CT)

  if ("Target" %in% colnames(my_data)) {
    my_data$target_name <- my_data$Target
  } else {
    my_data$target_name <- my_data$`Target Name`
  }

  if ("Sample" %in% colnames(my_data)) {
    my_data$sample_name <- my_data$Sample
  } else {
    my_data$sample_name <- my_data$`Sample Name`
  }

  my_data <- my_data[!is.na(my_data$target_name), ]
  my_data <- my_data[!is.na(my_data$sample_name), ]

  # Summarize by sample and condition and h4
  file_data <- my_data
  my_data <- file_data %>%
    group_by(sample_name, target_name) %>%
    summarise(mean_ct = mean(CT), count = n(), spread = max(CT) - min(CT))


  # Get z series, condition, rt-, and flags
  my_data$condition <- gsub(
    pattern = "([zZ]\\d+)(.)",
    replacement = "\\2",
    x = my_data$sample_name
  )

  my_data$z_series <- str_replace(
    my_data$sample_name,
    my_data$condition,
    ""
  )

  my_data$condition <- trimws(my_data$condition)

  my_data$is_rt_minus <- FALSE
  my_data[contains("rt-", my_data$condition, ignore.case = TRUE), ]$is_rt_minus <- TRUE


  my_data$is_h4 <- FALSE
  my_data[contains("H4", my_data$target_name, ignore.case = TRUE), ]$is_h4 <- TRUE

  my_data$is_control <- FALSE
  my_data[contains("c", my_data$condition, ignore.case = TRUE), ]$is_control <- TRUE

  my_data$is_trained <- !my_data$is_control


  # Now process each sample

  samples <- unique(my_data$z_series)

  print(
    glue::glue(
      "Found {length(samples)} samples."
    )
  )

  for (sample in samples) {

    # sample <- "z350"

    genes <- unique(
      my_data[my_data$z_series == sample & !my_data$is_h4, ]$target_name
    )


    print(
      glue::glue(
        "For sample {sample} found {length(genes)}: {paste(genes, collapse = ', ')}."
      )
    )


    for (gene in genes) {
      gene_c <- my_data[my_data$z_series == sample & my_data$target_name == gene & !my_data$is_rt_minus & my_data$is_control, ]
      gene_t <- my_data[my_data$z_series == sample & my_data$target_name == gene & !my_data$is_rt_minus & my_data$is_trained, ]
      h4_c <- my_data[my_data$z_series == sample & my_data$is_h4 & !my_data$is_rt_minus & my_data$is_control, ]
      h4_t <- my_data[my_data$z_series == sample & my_data$is_h4 & !my_data$is_rt_minus & my_data$is_trained, ]
      gene_rtm_c <- my_data[my_data$z_series == sample & my_data$target_name == gene & my_data$is_rt_minus & my_data$is_control, ]
      gene_rtm_t <- my_data[my_data$z_series == sample & my_data$target_name == gene & my_data$is_rt_minus & my_data$is_trained, ]
      h4_rtm_c <- my_data[my_data$z_series == sample & my_data$is_h4 & my_data$is_rt_minus & my_data$is_control, ]
      h4_rtm_t <- my_data[my_data$z_series == sample & my_data$is_h4 & my_data$is_rt_minus & my_data$is_trained, ]

      fix_empty <- function(table) {
        if (nrow(table) == 0) {
          table[1, "mean_ct"] <- NA
          table$count <- NA
          table$spread <- NA
          table$relative_quantity <- NA
        }
        return(table)
      }

      gene_c <- fix_empty(gene_c)
      gene_t <- fix_empty(gene_t)
      h4_c <- fix_empty(h4_c)
      h4_t <- fix_empty(h4_t)
      gene_rtm_c <- fix_empty(gene_rtm_c)
      gene_rtm_t <- fix_empty(gene_rtm_t)
      h4_rtm_t <- fix_empty(h4_rtm_t)
      h4_rtm_c <- fix_empty(h4_rtm_c)


      gene_t$relative_quantity <- gene_t$mean_ct - h4_t$mean_ct
      gene_c$relative_quantity <- gene_c$mean_ct - h4_c$mean_ct
      lfc <- gene_c$relative_quantity[1] - gene_t$relative_quantity[1]

      gene_result <- data.frame(
        gene,
        sample,
        lfc,
        gene_t$relative_quantity,
        gene_c$relative_quantity,
        gene_t$mean_ct,
        h4_t$mean_ct,
        gene_c$mean_ct,
        h4_c$mean_ct,
        gene_t$count,
        gene_t$spread,
        gene_rtm_t$mean_ct,
        h4_t$count,
        h4_t$spread,
        h4_rtm_t$mean_ct,
        gene_c$count,
        gene_c$spread,
        gene_rtm_c$mean_ct,
        h4_c$count,
        h4_c$spread,
        h4_rtm_c$mean_ct,
        file_time = file.info(my_file)$ctime,
        file_name = basename(my_file),
        directory = dirname(my_file)
      )

      tbl_output <- rbind(
        tbl_output,
        gene_result
      )

    }


  }

  return(tbl_output)


}
