library(shiny)
library(dplyr)
library(purrr)
library(readxl)
library(openxlsx)

source('qPCR_analysis_methods.R')

ui <- fluidPage(
  titlePanel('Avidity In Vivo Benchling qPCR KD Analysis'),
  sidebarLayout(
    sidebarPanel(
      tags$p(
        'This tool is intended to standardize the way that qPCR analysis is performed at Avidity.
        Additionally, this tool aims to simplify the qPCR workflow and address the biggest
        pain point--Benchling entries. While you will still need to interact with Benchling,
        here we have provided a means to consistently and easily register you entities in a single
        copy-paste. Thankfully though, you can use this tool ', tags$i('without '), 'Benchling.'
      ),
      hr(),
      tags$strong('Instructions:'),
      tags$p(
        tags$strong('1. '),
        'Use the provided platemap/study design/probe template to create ',
        tags$strong('ALL '),
        'of your platemaps.',
        downloadLink(outputId = 'download_template', label = 'Download platemap template here')
      ),
      tags$p(
        tags$strong('2. '),
        'Use the ',
        tags$strong('# of plates '),
        'field to indicate your number of plates to analyze.'
      ),
      tags$p(
        tags$strong('3. (Optional) '),
        'Upload the ',
        tags$i('.csv '),
        'of the ',
        tags$strong('cDNA entitities '),
        'that you downloaded from the lookup table from Benchling. 
        If you do not have cDNA entitites, skip this step and proceed to the next step.'
      ),
      tags$p(
        tags$strong('4. '),
        'Upload the platemap and raw exported data ',
        tags$strong('for each plate '),
        'to the designated upload spot.'
      ),
      tags$p(
        tags$strong('5. '),
        'Perform a ready check by clicking the ',
        tags$strong('Ready Check '),
        'button.'
      ),
      tags$p(
        tags$strong('6. '),
        'If the ready check is successful, analyze and download your data by clicking the ',
        tags$strong('Analyze and Download '),
        'button.'
      )
    ),
    
    
    mainPanel(
      uiOutput('page_title'),
      hr(),
      uiOutput('plates'),
      uiOutput('pages')
    )
  )
)

server <- function(input, output) {
  default_maxCT <- 30
  default_minCT <- 10
  default_num_plates <- 1
  
  rv <- reactiveValues(
    ready_check = FALSE,
    ready_check_msg = NULL,
    num_plates = default_num_plates,
    maxCT = default_maxCT,
    minCT = default_minCT,
    data.input.platemaps = list(),
    data.input.rawdatas = list(),
    data.input.benchling = NULL,
    filepaths.benchling = NULL,
    filepaths.plates = list(),
    filepaths.rawdatas = list()
  )
  
  rv$current_page_advance_button <- reactive({
    if (!rv$ready_check) {
      actionButton(inputId = 'check_btn',
                   label = 'Ready Check',
                   class = 'btn-primary')
    } else {
      if (rv$num_plates == 1) {
        downloadButton(outputId = 'download_single',
                       label = 'Analyze and Download',
                       class = 'btn-primary')
      } else {
        downloadButton(outputId = 'download_zip',
                       label = 'Analyze and Download',
                       class = 'btn-primary')
      }
    }
  })
  
  output$page_title <- renderUI({
    tagList(fluidRow(
      column(
        12,
        tags$strong('File Uploads'),
        actionButton(
          inputId = 'restart_btn',
          label = 'Restart',
          class = 'btn-danger'
        ),
        actionButton(
          inputId = 'settings_btn',
          label = NULL,
          icon = icon("cog", lib = "glyphicon")
        ),
        rv$current_page_advance_button()
      )
    ), 
    fluidRow(column(
      12,
      tags$p(tags$span(
        style = ifelse(rv$ready_check, 'color:green', 'color:red'),
        rv$ready_check_msg
      ))
    )),
    fluidRow(column(
      12,
      actionButton(
        inputId = 'examples',
        label = 'Click to view upload examples',
        class = 'btn-warning'
      )
    )))
  })
  
  examplesModal <- function() {
    modalDialog(
      titlePanel('Examples of types of file uploads'),
      tabsetPanel(
        tabPanel(
          title = 'Benchling cDNAs',
          hr(),
          tags$iframe(src = 'example_benchling_cdna.png', width = '100%', height = '500px', seamless = TRUE)
        ),
        tabPanel(
          title = 'Platemaps',
          hr(),
          tags$p('The names here that you used for your probes ', tags$strong('MUST '), 
                      'be the same as the shorthand ones listed in your Benchling probes.'),
          br(),
          tags$p('The Organs that you entered ', tags$strong('MUST '), 
                 'be the same as the ones listed for you Benchling cDNAs.'),
          br(),
          tags$iframe(src = 'D28 platemaps.pdf', width = '100%', height = '500px', seamless = TRUE)
        ),
        tabPanel(
          title = 'Study Design',
          hr(),
          tags$iframe(src = 'example_studydesign.png', width = '100%', height = '500px', seamless = TRUE)
        ),
        tabPanel(
          title = 'Benchling Probes',
          hr(),
          tags$p('The shorthand names here that you used for your probes ', tags$strong('MUST '), 
                 'be the same as those listed in your platemaps.'),
          br(),
          tags$iframe(src = 'example_probes.png', width = '100%', height = '500px', seamless = TRUE)
        ),
        tabPanel(
          title = 'Raw Data (QS6 format)',
          hr(),
          tags$iframe(src = 'example_rawdata_qs6.png', width = '100%', height = '500px', seamless = TRUE)
        ),
        tabPanel(
          title = 'Raw Data (QS7 format)',
          hr(),
          tags$iframe(src = 'example_rawdata_qs7.png', width = '100%', height = '500px', seamless = TRUE)
        )
      ),
      footer = modalButton(label = 'OK',),
      easyClose = TRUE,
      size = 'xl'
    )
  }
  
  observeEvent(input$examples, {
    showModal(examplesModal())
  })
  
  settingsModal <- reactive(function() {
    modalDialog(
      titlePanel('Global Settings'),
      uiOutput('CT_limits'),
      footer = tagList(
        actionButton(
          inputId = 'settings_reset',
          label = 'Reset to Defaults',
          class = 'btn-danger'
        ),
        modalButton(label = 'Cancel'),
        actionButton(
          inputId = 'settings_ok',
          label = 'Confirm',
          class = 'btn-primary'
        )
      )
    )
  })
  
  observeEvent(input$settings_btn, {
    showModal(settingsModal()())
  })
  
  observeEvent(input$settings_ok, {
    old_max <- rv$maxCT
    old_min <- rv$minCT
    if (is.na(input$maxCT) |
        is.null(input$maxCT) | !between(input$maxCT, 2, 40)) {
      rv$maxCT <- default_maxCT
    } else {
      rv$maxCT <- input$maxCT
    }
    if (is.na(input$minCT) |
        is.null(input$minCT) | !between(input$minCT, 1, 39)) {
      rv$minCT <- default_minCT
    } else {
      rv$minCT <- input$minCT
    }
    if (!(old_max == rv$maxCT & old_min == rv$minCT)) {
      rv$ready_check <- FALSE
    }
    removeModal()
  })
  
  observeEvent(input$settings_reset, {
    old_max <- rv$maxCT
    old_min <- rv$minCT
    rv$maxCT <- default_maxCT
    rv$minCT <- default_minCT
    if (!(old_max == rv$maxCT & old_min == rv$minCT)) {
      rv$ready_check <- FALSE
    }
    removeModal()
  })
  
  output$CT_limits <- renderUI({
    fluidRow(column(
      3,
      numericInput(
        inputId = 'maxCT',
        label = 'Max CT',
        value = rv$maxCT,
        min = 2,
        max = 40
      )
    ), column(
      3,
      numericInput(
        inputId = 'minCT',
        label = 'Min CT',
        value = rv$minCT,
        min = 1,
        max = 39
      )
    ))
  })
  
  observeEvent(input$restart_btn, {
    rv$ready_check <- FALSE
    rv$ready_check_msg <- NULL
    rv$num_plates <- default_num_plates
    rv$maxCT <- default_maxCT
    rv$minCT <- default_minCT
    rv$data.input.platemaps <- list()
    rv$data.input.rawdatas <- list()
    rv$data.input.benchling <- NULL
    filepaths.benchling <- NULL
    filepaths.plates <- list()
    filepaths.rawdatas <- list()
  })
  
  observe({
    rv$filepaths.benchling <- input$benchling$datapath
    rv$ready_check <- FALSE
  })
  
  observe({
    rv$filepaths.plates <- map(1:rv$num_plates, function(n) {
      input[[paste0('platemap', n)]]$datapath
    })
    rv$ready_check <- FALSE
  })
  
  observe({
    rv$filepaths.rawdatas <- map(1:rv$num_plates, function(n) {
      input[[paste0('rawdata', n)]]$datapath
    })
    rv$ready_check <- FALSE
  })
  
  ready_check_report <- observeEvent(input$check_btn, {
    withProgress(message = 'Performing ready check...', value = 0, {
      checks <- rv$num_plates + 1
      rv$ready_check <- TRUE
      rv$ready_check_msg <- 'Ready check successful! Please click "Analyze and Download"'
      
      if (is.null(input$benchling)) {
        # rv$ready_check <- FALSE
        # rv$ready_check_msg <- 'Please upload a valid Benchling cDNA registry file.'
        rv$data.input.benchling <- NULL
        rv$ready_check_msg <- 'No Benchling entitites detected. Analysis will be performed but a Benchling table will not be output.'
      } else {
        if (tolower(tools::file_ext(input$benchling$datapath)) %in% c('xlsx', 'xls')) {
          rv$data.input.benchling <- read_excel(input$benchling$datapath)
        } else {
          rv$data.input.benchling <- read.csv(input$benchling$datapath)
        }
      }
      incProgress(1 / checks)
      
      walk(1:rv$num_plates, function(x) {
        platemap_file <- input[[paste0('platemap', x)]]
        rawdata_file <- input[[paste0('rawdata', x)]]
        if (is.null(platemap_file) | is.null(rawdata_file)) {
          rv$ready_check <- FALSE
          rv$ready_check_msg <- 'Please ensure all entered plates have both a platemap and a raw data file.'
          rv$data.input.platemaps[[x]] <- list()
          rv$data.input.rawdatas[[x]] <- list()
        } else {
          rv$data.input.platemaps[[x]] <- excel_sheets(platemap_file$datapath) %>% set_names() %>% map(read_excel, path = platemap_file$datapath)
          rv$data.input.rawdatas[[x]] <- excel_sheets(rawdata_file$datapath) %>% set_names() %>% map(read_excel, path = rawdata_file$datapath)
        }
        incProgress(1 / checks,
                    message = paste0('Performing ready check (', x, '/', rv$num_plates, '...'))
      })
    })
  })
  
  output$plates <- renderUI({
    tagList(fluidRow(column(
      12,
      tags$strong('You must use the provided platemap template. ', style = 'color:red'),
      downloadLink(outputId = 'download_template', label = 'Download platemap template here')
    )), br(), fluidRow(column(
      3,
      numericInput(
        inputId = 'num_plates',
        label = '# of plates',
        value = rv$num_plates,
        min = 1,
        step = 1
      )
    ), column(
      9,
      fileInput(
        inputId = 'benchling',
        label = 'Benchling cDNAs',
        accept = c('.csv', '.xlsx', '.xls'),
        placeholder = 'None',
        width = '100%'
      )
    )), hr())
  })
  
  observeEvent(input$num_plates, {
    rv$num_plates <- input$num_plates
  })
  
  pages_uploads <- reactive({
    map(1:rv$num_plates, function(n) {
      fluidRow(column(
        6,
        fileInput(
          inputId = paste0('platemap', n),
          label = paste0('Platemap #', n),
          ,
          accept = c('.xlsx', '.xls'),
          placeholder = 'None',
          width = '100%'
        )
      ), column(
        6,
        fileInput(
          inputId = paste0('rawdata', n),
          label = paste0('Raw Data #', n),
          ,
          accept = c('.xlsx', '.xls'),
          placeholder = 'None',
          width = '100%'
        )
      ))
    })
    
  })
  
  output$pages <- renderUI({
    tagList(pages_uploads())
  })
  
  template_name <- 'platemaps_template_blank.xlsx'
  template_path <- file.path(getwd(), template_name) %>% tools::file_path_as_absolute()
  
  output$download_template <- downloadHandler(
    filename = template_name,
    content = function(file) {
      file.copy(template_path, file)
    }
  )
  
  downloader <- observeEvent(input$download, {
    
  })
  
  output$download_single <- downloadHandler(
    filename = paste(format(Sys.time(), '%F %H-%M-%S'), 'qPCR_analysis.xlsx', sep = '_'),
    content = function(file) {
      withProgress(message = 'Analyzing...', value = 0, {
        df_benchling_samples <- NULL
        if (!is.null(rv$data.input.benchling)) {
          df_benchling_samples <- extract_and_convert_benchling_info(rv$data.input.benchling)
        }
        df_rawdata <- process_raw_data(rv$data.input.rawdatas[[1]])
        df_info <- extract_platemaps_and_study_info(rv$data.input.platemaps[[1]])
        df_plate_metadata <- df_info[['plate_metadata']]
        df_benchling_probes <- df_info[['probes']]
        incProgress(0.2, message = 'Analyzing (1/1)...')
        
        # Do stuff here
        df_out_all <- merge_data(
          df_plate_metadata = df_plate_metadata,
          df_benchling_probes = df_benchling_probes,
          df_benchling_samples = df_benchling_samples,
          df_rawdata = df_rawdata
        ) %>% analyze_data(maxCT = rv$maxCT, minCT = rv$minCT)
        incProgress(0.2, message = 'Analyzing (1/1)...')
        df_out_plates <- create_data_pivots(df_out_all)
        incProgress(0.2, message = 'Analyzing (1/1)...')
        df_out_benchling <- NULL
        if (is.null(df_benchling_samples)) {
          df_out_benchling <- benchling_output(df_in = NULL)
        } else {
          df_out_benchling <- benchling_output(df_in = df_out_all)
        }
        incProgress(0.2, message = 'Analyzing (1/1)...')
        
        # Actual download portion
        create_qPCR_Excel(
          df_all = df_out_all,
          df_plates = df_out_plates,
          df_benchling = df_out_benchling,
          file = file
        )
        incProgress(0.2, message = 'Analyzing (1/1)...')
      })
    }
  )
  
  zip_file_type <- '.zip'
  
  output$download_zip <- downloadHandler(
    filename = paste0(
      format(Sys.time(), '%F %H-%M-%S'),
      '_qPCR_analysis',
      zip_file_type
    ),
    content = function(file) {
      withProgress(message = 'Analyzing...', value = 0, {
        iters_per_plate <- 5
        iters <- rv$num_plates * iters_per_plate + 1
        #go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        
        #loop through the sheets
        df_benchling_samples <- NULL
        if (!is.null(rv$data.input.benchling)) {
          df_benchling_samples <- extract_and_convert_benchling_info(rv$data.input.benchling)
        }
        incProgress(1 / iters)
        files <- map_vec(1:rv$num_plates, function(n) {
          df_rawdata <- process_raw_data(rv$data.input.rawdatas[[n]])
          df_info <- extract_platemaps_and_study_info(rv$data.input.platemaps[[n]])
          df_plate_metadata <- df_info[['plate_metadata']]
          df_benchling_probes <- df_info[['probes']]
          incProgress(1 / iters,
                      message = paste0('Analyzing (', n, '/', rv$num_plates, ')...'))
          
          # Do stuff here
          df_out_all <- merge_data(
            df_plate_metadata = df_plate_metadata,
            df_benchling_samples = df_benchling_samples,
            df_benchling_probes = df_benchling_probes,
            df_rawdata = df_rawdata
          ) %>% analyze_data(maxCT = rv$maxCT, minCT = rv$minCT)
          incProgress(1 / iters,
                      message = paste0('Analyzing (', n, '/', rv$num_plates, ')...'))
          df_out_plates <- create_data_pivots(df_out_all)
          incProgress(1 / iters,
                      message = paste0('Analyzing (', n, '/', rv$num_plates, ')...'))
          df_out_benchling <- NULL
          if (is.null(df_benchling_samples)) {
            df_out_benchling <- benchling_output(df_in = NULL)
          } else {
            df_out_benchling <- benchling_output(df_in = df_out_all)
          }
          incProgress(1 / iters,
                      message = paste0('Analyzing (', n, '/', rv$num_plates, ')...'))
          
          mini_filename <- paste0('qPCR_analysis_Plate_', n, '.xlsx')
          create_qPCR_Excel(
            df_all = df_out_all,
            df_plates = df_out_plates,
            df_benchling = df_out_benchling,
            file = mini_filename
          )
          incProgress(1 / iters,
                      message = paste0('Analyzing (', n, '/', rv$num_plates, ')...'))
          return(mini_filename)
        })
        
        #create the zip file
        if (zip_file_type == '.zip') {
          zip(file, files)
          # system2("zip", args = (paste(file, files, sep = " ")))
        } else if (zip_file_type == '.tar') {
          tar(file, files)
        }
      })
    }
  )
  
  
  
  
}

shinyApp(ui, server)
