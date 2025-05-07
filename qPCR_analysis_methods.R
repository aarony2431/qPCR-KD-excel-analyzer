# Import R packages needed for qPCR analysis
library(dplyr)
library(tidyverse)
library(readxl)
library(gsubfn)
library(purrr)
library(openxlsx)

# Global Variables
plate_height <- 16
plate_width <- 24
plate_headers <- c(
  'Organ',
  'Group',
  'Control Group',
  'Normalize to Group',
  'Animal',
  'Tissue Plate Number',
  'Replicate',
  'FAM Probe',
  'VIC Probe'
)
plate_offsets <- c(0, 1, 2, 2, 1, 2, 1, 1, 1)
benchling.column_titles <- c(
  'Entity',
  'Animal.Sample',
  'Animal.Group',
  'Parent.Animal',
  'Tissue.or.Fluid.Type',
  'Plate.Number'
)
output.plates <- c('CT_VIC', 'CT_FAM', 'dCT', 'ddCT', '% mRNA expression')
output.benchling <- c(
  'Entity',
  'FAM assay',
  'VIC assay',
  'CT_FAM',
  'Ct Threshold_FAM',
  'CT_VIC',
  'Ct Threshold_VIC',
  '% mRNA expression',
  'Replicate'
)




# Function for processing data.frame from raw data into usable data
# Remove excess lines and transform to tibble dataframe
# Consistent headers in line with QS6 format
# Ablility to detect which format the data are in (QS6 vs. QS7)
process_raw_data <- function(df_in) {
  # From raw data
  # QS6 starts at row 47 in "Results"
  # 'CT', 'Ct Threshold', 'Target Name', 'Reporter'
  # QS7 starts at row 22 in "Results"
  # 'Cq', 'Threshold', 'Target', 'Reporter'
  rawdata.column_titles_line <- list('QS6' = 46, 'QS7' = 21)
  rawdata.column_titles <- list(
    'QS6' = c(
      well = 'Well',
      wellpos = 'Well Position',
      ct = 'CT',
      threshold = 'Ct Threshold',
      target = 'Target Name',
      reporter = 'Reporter'
    ),
    'QS7' = c(
      well = 'Well',
      wellpos = 'Well Position',
      ct = 'Cq',
      threshold = 'Threshold',
      target = 'Target',
      reporter = 'Reporter'
    )
  )
  
  # Identify QS system from raw data "Results"
  QS.column_titles_line <- apply(df_in$Results, 1, function(r)
    max(
      sum(r %in% rawdata.column_titles$QS6),
      sum(r %in% rawdata.column_titles$QS7)
    )) %>% which.max()
  QS.system_name <- which(rawdata.column_titles_line == QS.column_titles_line) %>% names()
  QS.column_titles <- rawdata.column_titles[[QS.system_name]]
  
  # Create cleaned raw data "Results" so it has proper headers and only useful data
  # Collapse to only important columns and split reporters
  # Consistent formatting of column names to QS6
  data.rawdata <- df_in$Results[-c(1:rawdata.column_titles_line[[QS.system_name]]), ]
  names(data.rawdata) <- df_in$Results[rawdata.column_titles_line[['QS6']], ]
  QS.system_name <- 'QS6'
  QS.column_titles <- rawdata.column_titles[[QS.system_name]]
  data.rawdata <- pivot_wider(
    data.rawdata,
    id_cols = c(QS.column_titles[['well']], QS.column_titles[['wellpos']]),
    names_from = c(QS.column_titles[['reporter']]),
    values_from = c(QS.column_titles[['ct']], QS.column_titles[['threshold']])
  )
  
  return(data.rawdata)
}

# Extract the data from the list from the platemaps file based on consistent expectations
# Merge all data to a single data table
extract_platemaps_and_study_info <- function(data_in,
                                             plate_height = 16,
                                             plate_width = 24,
                                             plate_headers = c(
                                               'Organ',
                                               'Group',
                                               'Control Group',
                                               'Normalize to Group',
                                               'Animal',
                                               'Tissue Plate Number',
                                               'Replicate',
                                               'FAM Probe',
                                               'VIC Probe'
                                             ),
                                             plate_offsets = c(0, 1, 2, 2, 1, 2, 1, 1, 1)) {
  # Metadata needed for uniqueness:
  # Group
  # Animal
  # Tissue
  # Plate Number
  # Harvest Day
  
  data.study_design <- data_in$`Study Design`
  data.benchling.probes <- data_in$Probes %>% t() %>% as_tibble()
  names(data.benchling.probes) <- data.benchling.probes[1, ]
  data.benchling.probes <- data.benchling.probes[-1, ]
  platemaps <- data_in$Platemap
  
  # Extract the plate based on standard information
  plate_absolute_offsets <- sapply(1:length(plate_offsets), function(n)
    (n - 1) * (plate_height + 1) + sum(plate_offsets[1:n]))
  get_plate <- function(df,
                        offset,
                        name = 'Row',
                        h = 16,
                        w = 24) {
    out <- df[(1 + offset):(1 + offset + h), 1:(w + 1)]
    out[1, 1] <- name
    names(out) <- out[1, ]
    out <- out[-1, ] %>% mutate(across(everything(), as.character))
    return(out)
  }
  # Convert plate to table
  plate_to_table <- function(df, value = 'value', w = 24) {
    out <- pivot_longer(df,
                        as.character(1:w),
                        names_to = 'Column',
                        values_to = value)
    well_positions <- apply(out, 1, function (r)
      paste0(r[1], r[2]))
    out <- bind_cols('Well Position' = well_positions, out)
    return(out)
  }
  # Convert platemaps tab to merged data table
  data.plate_metadata <- lapply(1:length(plate_headers), function(n)
    plate_to_table(get_plate(platemaps, plate_absolute_offsets[n]), value = plate_headers[n])) %>%
    reduce(full_join, by = c('Well Position', 'Row', 'Column'))
  
  # Add leading '00' to group and animals
  data.plate_metadata$Group <- sapply(data.plate_metadata$Group, function(value)
    str_sub(paste0('00', value), start = -3))
  data.plate_metadata$Animal <- sapply(data.plate_metadata$Animal, function(value)
    str_sub(paste0('00', value), start = -3))
  
  # Create local unique IDs
  data.plate_metadata <- bind_cols(data.plate_metadata,
                                   Local_Unique_ID = apply(data.plate_metadata, 1, function(r)
                                     paste(r['Group'], r['Animal'], r['Organ'], r['Tissue Plate Number'], sep = '-')))
  
  return(list(plate_metadata = data.plate_metadata, probes = data.benchling.probes))
}

# Extract Benchling information from data.frame from file
# Create local IDs to match with the local IDs from the platemaps
extract_and_convert_benchling_info <- function(df_in,
                                               benchling.column_titles = c(
                                                 'Entity',
                                                 'Animal.Sample',
                                                 'Animal.Group',
                                                 'Parent.Animal',
                                                 'Tissue.or.Fluid.Type',
                                                 'Plate.Number'
                                               )) {
  df_out <- tibble(Local_Unique_ID = apply(df_in, 1, function(r)
    paste(r['Parent.Animal'], r['Tissue.or.Fluid.Type'], r['Plate.Number'], sep = '-')),
    Entity = df_in$Entity)
  df_out$Local_Unique_ID <- sapply(df_out$Local_Unique_ID, function(r)
    str_sub(r, start = 1 + str_locate(r, '-')[1]))
  
  return(df_out)
}

# Merge Benchling dictionary with raw data and platemaps
merge_data <- function(df_plate_metadata,
                       df_benchling_samples,
                       df_benchling_probes,
                       df_rawdata) {
  # Merge data with benchling ---------------------------------------------------
  df_out <- reduce(list(df_plate_metadata, df_benchling_samples),
                   left_join,
                   by = 'Local_Unique_ID')
  df_out <- reduce(list(df_out, df_rawdata), full_join, by = 'Well Position')
  df_out$`FAM assay` <- apply(df_out, 1, function(r) {
    probe_name <- as.character(r['FAM Probe'])
    if (is.na(probe_name)) {
      return(NA)
    } else {
      return(as.character(df_benchling_probes[probe_name]))
    }
  })
  df_out$`VIC assay` <- apply(df_out, 1, function(r) {
    probe_name <- as.character(r['VIC Probe'])
    if (is.na(probe_name)) {
      return(NA)
    } else {
      return(as.character(df_benchling_probes[probe_name]))
    }
  })
  
  return(df_out)
}

# Data Analysis
analyze_data <- function(df_in,
                         maxCT = 30,
                         minCT = 10) {
  numeric_conversions <- c(
    'Column',
    'Group',
    'Animal',
    'Replicate',
    'Well',
    'Tissue Plate Number',
    'Normalize to Group',
    'CT_FAM',
    'CT_VIC',
    'Ct Threshold_FAM',
    'Ct Threshold_VIC'
  )
  df_out <- df_in %>% mutate(across(numeric_conversions, as.numeric))
  
  df_out$minCT <- rep(minCT, length(df_out$Group))
  df_out$maxCT <- rep(maxCT, length(df_out$Group))
  df_out$CT_flag <- apply(df_out, 1, function(r) {
    fam <- as.numeric(r['CT_FAM'])
    vic <- as.numeric(r['CT_VIC'])
    if (is.na(fam) |
        is.na(vic) |
        !between(fam, minCT, maxCT) |
        !between(vic, minCT, maxCT)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  
  df_out$dCT <- apply(df_out, 1, function(r) {
    if (!as.logical(r['CT_flag'])) {
      as.numeric(r['CT_FAM']) - as.numeric(r['CT_VIC'])
    } else {
      NA
    }
  })
  data.normalization_groups <- df_out %>%
    group_by(`Normalize to Group`) %>%
    filter(`Control Group` != 'NA', !`CT_flag`) %>%
    summarise_at(vars(dCT), list('Control dCT' = mean))
  df_out <- reduce(list(df_out, data.normalization_groups), left_join, by = 'Normalize to Group')
  df_out$ddCT <- apply(df_out, 1, function(r) {
    if (!as.logical(r['CT_flag'])) {
      as.numeric(r['dCT']) - as.numeric(r['Control dCT'])
    } else {
      NA
    }
  })
  df_out$`% mRNA expression` <- apply(df_out, 1, function(r) {
    if (!as.logical(r['CT_flag'])) {
      100 * 2^(-1 * as.numeric(r['ddCT']))
    } else {
      NA
    }
  })
  
  return(df_out)
}

# Create visualizations for easy viewing of data
create_data_pivots <- function(df_in,
                               output.plates = c('CT_VIC', 'CT_FAM', 'dCT', 'ddCT', '% mRNA expression')) {
  output.pivot.plates <- df_in[which(names(df_in) %in% c('Well Position', output.plates))]
  output.pivot.plates$Row <- str_sub_all(output.pivot.plates$`Well Position`,
                                         start = 1,
                                         end = 1)
  output.pivot.plates$Column <- str_sub_all(
    output.pivot.plates$`Well Position`,
    start = 2,
    end = length(output.pivot.plates$`Well Position`)
  ) %>% as.numeric()
  df_out <- lapply(output.plates, function(plate) {
    cols <- c('Row', 'Column', plate)
    out <- output.pivot.plates[which(names(output.pivot.plates) %in% cols)] %>%
      pivot_wider(names_from = 'Column', values_from = plate)
    names(out)[1] <- plate
    return(out)
  })
  
  return(df_out)
}

# Benchling output
benchling_output <- function(df_in,
                             output.benchling = c(
                               'Entity',
                               'FAM assay',
                               'VIC assay',
                               'CT_FAM',
                               'Ct Threshold_FAM',
                               'CT_VIC',
                               'Ct Threshold_VIC',
                               '% mRNA expression',
                               'Replicate'
                             )) {
  output.pivot.benchling <- df_in[-which(is.na(df_in$Entity)), which(names(df_in) %in% output.benchling)] %>% .[, output.benchling]
  output.benchling.pivot_columns <- c('CT_FAM',
                                      'Ct Threshold_FAM',
                                      'CT_VIC',
                                      'Ct Threshold_VIC',
                                      '% mRNA expression')
  output.pivot.benchling <- pivot_wider(
    output.pivot.benchling,
    names_from = 'Replicate',
    values_from = all_of(output.benchling.pivot_columns)
  )
  df_out <- apply(output.pivot.benchling, 1, function(r) {
    out.data <- sapply(output.benchling.pivot_columns, function(col) {
      cols <- paste(rep(col, each = 2), 1:2, sep = '_')
      avg <- data.frame(mean(as.numeric(r[cols])),
                        check.names = FALSE,
                        fix.empty.names = FALSE)
      return(avg)
    }) %>% as_tibble()
    out <- bind_cols(
      Entity = r['Entity'],
      'FAM assay' = r['FAM assay'],
      'VIC assay' = r['VIC assay'],
      out.data
    )
    return(out)
  }) %>% reduce(bind_rows)
  
  return(df_out)
}

# Create Excel file with final data set
create_qPCR_Excel <- function(df_all,
                              df_plates,
                              df_benchling,
                              file = NULL,
                              plate_height = 16,
                              plate_width = 24,
                              output.plates.offset = 2) {
  wb <- createWorkbook()
  
  addWorksheet(wb, 'Data', zoom = 75)
  addWorksheet(wb, 'Plate Data (Visual)', zoom = 75)
  addWorksheet(wb, 'Benchling', zoom = 75)
  
  writeDataTable(wb, 'Data', df_all)
  walk(1:length(df_plates), function(x) {
    startRow <- 1 + (x - 1) * (plate_height + 1 + output.plates.offset)
    writeData(
      wb,
      'Plate Data (Visual)',
      df_plates[[x]],
      startRow = startRow,
      borders = "surrounding",
      borderColour = "black"
    )
    conditionalFormatting(
      wb,
      'Plate Data (Visual)',
      rows = (startRow + 1):(startRow + 1 + plate_height),
      cols = 2:(plate_width + 1),
      type = 'colourScale',
      style = c('#F8696B', '#FFEB84', '#63BE7B')
    )
  })
  writeDataTable(wb, 'Benchling', df_benchling)
  
  if (is.null(file)) {
    file <- paste(format(Sys.time(), '%F %H-%M-%S'),
                      'qPCR_analysis.xlsx',
                      sep = '_')
  }
  saveWorkbook(wb, file = file, overwrite = TRUE)
}


# # Testing
# datafile.path.platemap <- r"(C:\Users\ayu\OneDrive - Avidity Biosciences\Documents\Data\qPCR KD examples\2024-182 KD resources for automated in R\D28 platemaps.xlsx)"
# datafile.path.rawdata <- r"(C:\Users\ayu\OneDrive - Avidity Biosciences\Documents\Data\qPCR KD examples\2024-182 KD resources for automated in R\2024-04-05 AY 049 D28 KD.xlsx)"
# datafile.path.Benchling <- r"(C:\Users\ayu\OneDrive - Avidity Biosciences\Documents\Data\qPCR KD examples\2024-182 KD resources for automated in R\cDNA D28 Benchling.csv)"
# 
# input.platemap <- excel_sheets(datafile.path.platemap) %>% set_names() %>% map(read_excel, path = datafile.path.platemap)
# input.rawdata <- excel_sheets(datafile.path.rawdata) %>% set_names() %>% map(read_excel, path = datafile.path.rawdata)
# input.Benchling <- read.csv(datafile.path.Benchling)
# 
# # Get data from dfs
# df_benchling_samples <- extract_and_convert_benchling_info(input.Benchling)
# df_rawdata <- process_raw_data(input.rawdata)
# df_info <- extract_platemaps_and_study_info(input.platemap)
# df_plate_metadata <- df_info[['plate_metadata']]
# df_benchling_probes <- df_info[['probes']]
# 
# # Do stuff here
# df_out_all <- merge_data(
#   df_plate_metadata = df_plate_metadata,
#   df_benchling_samples = df_benchling_samples,
#   df_benchling_probes = df_benchling_probes,
#   df_rawdata = df_rawdata
# ) %>% analyze_data()
# df_out_plates <- create_data_pivots(df_out_all)
# df_out_benchling <- benchling_output(df_out_all)
# 
# # Excel
# folderpath <- r"(C:\Users\ayu\OneDrive - Avidity Biosciences\Documents\Data\qPCR KD examples\2024-182 KD resources for automated in R)" %>% tools::file_path_as_absolute()
# create_qPCR_Excel(df_out_all, df_out_plates, df_out_benchling, folderpath)
