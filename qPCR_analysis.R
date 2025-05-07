# Import R packages needed for qPCR analysis
library(dplyr)
library(tidyverse)
library(readxl)
library(gsubfn)
library(purrr)
library(openxlsx)

# Test files
datafile.path.platemap <- r"(C:\Users\ayu\OneDrive - Avidity Biosciences\Documents\Data\qPCR KD examples\2024-182 KD resources for automated in R\D28 platemaps.xlsx)"
datafile.path.rawdata <- r"(C:\Users\ayu\OneDrive - Avidity Biosciences\Documents\Data\qPCR KD examples\2024-182 KD resources for automated in R\2024-04-05 AY 049 D28 KD.xlsx)"
datafile.path.Benchling <- r"(C:\Users\ayu\OneDrive - Avidity Biosciences\Documents\Data\qPCR KD examples\2024-182 KD resources for automated in R\cDNA D28 Benchling.csv)"

input.platemap <- excel_sheets(datafile.path.platemap) %>% set_names() %>% map(read_excel, path = datafile.path.platemap)
input.rawdata <- excel_sheets(datafile.path.rawdata) %>% set_names() %>% map(read_excel, path = datafile.path.rawdata)
input.Benchling <- read.csv(datafile.path.Benchling)

# Raw Data Processing ---------------------------------------------------------

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
QS.column_titles_line <- apply(input.rawdata$Results, 1, function(r)
  max(
    sum(r %in% rawdata.column_titles$QS6),
    sum(r %in% rawdata.column_titles$QS7)
  )) %>% which.max()
QS.system_name <- which(rawdata.column_titles_line == QS.column_titles_line) %>% names()
QS.column_titles <- rawdata.column_titles[[QS.system_name]]

# Create cleaned raw data "Results" so it has proper headers and only useful data
# Collapse to only important columns and split reporters
# Consistent formatting of column names to QS6
data.rawdata <- input.rawdata$Results[-c(1:rawdata.column_titles_line[[QS.system_name]]), ]
names(data.rawdata) <- input.rawdata$Results[rawdata.column_titles_line[['QS6']], ]
QS.system_name <- 'QS6'
QS.column_titles <- rawdata.column_titles[[QS.system_name]]
data.rawdata <- pivot_wider(
  data.rawdata,
  id_cols = c(QS.column_titles[['well']], QS.column_titles[['wellpos']]),
  names_from = c(QS.column_titles[['reporter']]),
  values_from = c(QS.column_titles[['ct']], QS.column_titles[['threshold']])
)


# Plate Maps and Study Design Processing --------------------------------------

# Metadata needed for uniqueness:
# Group
# Animal
# Tissue
# Plate Number
# Harvest Day

data.study_design <- input.platemap$`Study Design`
data.benchling.probes <- input.platemap$Probes %>% t() %>% as_tibble()
names(data.benchling.probes) <- data.benchling.probes[1, ]
data.benchling.probes <- data.benchling.probes[-1, ]
platemaps <- input.platemap$Platemap

# Extract the plate based on standard information
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



# Benchling cDNA table --------------------------------------------------------

benchling.column_titles <- c(
  'Entity',
  'Animal.Sample',
  'Animal.Group',
  'Parent.Animal',
  'Tissue.or.Fluid.Type',
  'Plate.Number'
)
data.benchling.dictionary <- tibble(
  Local_Unique_ID = apply(input.Benchling, 1, function(r)
    paste(r['Parent.Animal'], r['Tissue.or.Fluid.Type'], r['Plate.Number'], sep = '-')),
  Entity = input.Benchling$Entity
)
data.benchling.dictionary$Local_Unique_ID <- sapply(data.benchling.dictionary$Local_Unique_ID, function(r)
  str_sub(r, start = 1 + str_locate(r, '-')[1]))

# Merge data with benchling ---------------------------------------------------
data.output <- reduce(list(data.plate_metadata, data.benchling.dictionary),
                      left_join,
                      by = 'Local_Unique_ID')
data.output <- reduce(list(data.output, data.rawdata), full_join, by = 'Well Position')
data.output$`FAM assay` <- apply(data.output, 1, function(r) {
  probe_name <- as.character(r['FAM Probe'])
  if (is.na(probe_name)) {
    return(NA)
  } else {
    return(as.character(data.benchling.probes[probe_name]))
  }
})
data.output$`VIC assay` <- apply(data.output, 1, function(r) {
  probe_name <- as.character(r['VIC Probe'])
  if (is.na(probe_name)) {
    return(NA)
  } else {
    return(as.character(data.benchling.probes[probe_name]))
  }
})

# Data Analysis ---------------------------------------------------------------
numeric_conversions <- c('Column', 'Group', 'Animal', 'Replicate', 'Well', 'Tissue Plate Number', 'Normalize to Group', 'CT_FAM', 'CT_VIC', 'Ct Threshold_FAM', 'Ct Threshold_VIC')
data.output <- data.output %>% mutate(across(numeric_conversions, as.numeric))
data.output$dCT <- apply(data.output, 1, function(r) {
  as.numeric(r['CT_FAM']) - as.numeric(r['CT_VIC'])
})
data.normalization_groups <- data.output %>%
  group_by(`Normalize to Group`) %>%
  filter(c(`Control Group` != 'NA')) %>%
  summarise_at(vars(dCT), list('Control dCT' = mean))
data.output <- reduce(list(data.output, data.normalization_groups), left_join, by = 'Normalize to Group')
data.output$ddCT <- apply(data.output, 1, function(r)
  as.numeric(r['dCT']) - as.numeric(r['Control dCT']))
data.output$`% mRNA expression` <- apply(data.output, 1, function(r)
  100 * 2^(-1 * as.numeric(r['ddCT'])))

# Optional outputs for visualization of data ----------------------------------
output.plates <- c('CT_VIC', 'CT_FAM', 'dCT', 'ddCT', '% mRNA expression')
output.pivot.plates <- data.output[which(names(data.output) %in% c('Well Position', output.plates))]
output.pivot.plates$Row <- str_sub_all(output.pivot.plates$`Well Position`,
                                       start = 1,
                                       end = 1)
output.pivot.plates$Column <- str_sub_all(
  output.pivot.plates$`Well Position`,
  start = 2,
  end = length(output.pivot.plates$`Well Position`)
) %>% as.numeric()
output.plates.data <- lapply(output.plates, function(plate) {
  cols <- c('Row', 'Column', plate)
  out <- output.pivot.plates[which(names(output.pivot.plates) %in% cols)] %>% pivot_wider(names_from = 'Column', values_from = plate)
  names(out)[1] <- plate
  return(out)
})

# Benchling table output ------------------------------------------------------
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
output.pivot.benchling <- data.output[-which(is.na(data.output$Entity)), which(names(data.output) %in% output.benchling)] %>% .[, output.benchling]
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
output.pivot.benchling.pivoted <- apply(output.pivot.benchling, 1, function(r) {
  out.data <- sapply(output.benchling.pivot_columns, function(col) {
    cols <- paste(rep(col, each = 2), 1:2, sep = '_')
    avg <- data.frame(mean(as.numeric(r[cols])),
                      check.names = FALSE,
                      fix.empty.names = FALSE)
    return(avg)
  }) %>% as_tibble()
  out <- bind_cols(Entity = r['Entity'], 'FAM assay' = r['FAM assay'], 'VIC assay' = r['VIC assay'], out.data)
  return(out)
}) %>% reduce(bind_rows)

# Output excel sheet ----------------------------------------------------------
# data.output
# output.plates.data
# output.pivot.benchling.pivoted
# wb <- createWorkbook()
# 
# addWorksheet(wb, 'Data', zoom = 75)
# addWorksheet(wb, 'Plate Data (Visual)', zoom = 75)
# addWorksheet(wb, 'Benchling', zoom = 75)
# 
# writeDataTable(wb, 'Data', data.output)
# output.plates.offset <- 2
# walk(1:length(output.plates.data), function(x) {
#   startRow <- 1 + (x - 1) * (plate_height + 1 + output.plates.offset)
#   writeData(
#     wb,
#     'Plate Data (Visual)',
#     output.plates.data[[x]],
#     startRow = startRow,
#     borders = "surrounding",
#     borderColour = "black"
#   )
#   conditionalFormatting(
#     wb,
#     'Plate Data (Visual)',
#     rows = (startRow + 1):(startRow + 1 + plate_height),
#     cols = 2:(plate_width + 1),
#     type = 'colourScale',
#     style = c('#F8696B', '#FFEB84', '#63BE7B')
#   )
# })
# writeDataTable(wb, 'Benchling', output.pivot.benchling.pivoted)
# 
# output.filepath <- r"(C:\Users\ayu\OneDrive - Avidity Biosciences\Documents\Data\qPCR KD examples\2024-182 KD resources for automated in R\analyzed_output.xlsx)"
# saveWorkbook(wb, output.filepath, overwrite = TRUE)

