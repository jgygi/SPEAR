#' Get best SPEAR weights per response Y
#'@param SPEARobj SPEAR object (returned from \code{run_cv_spear})
#' @examples
#' SPEAR.run_interactive(SPEARobject)
#'@export
SPEAR.run_interactive <- function(SPEARobj){
  if(!"shiny" %in% rownames(installed.packages())){
    stop("*** NOTE: This function requires the 'shiny' package (can be installed with 'install.packages('shiny')'. You do not have this package installed.")
  }
  if(!"shinydashboardPlus" %in% rownames(installed.packages())){
    stop("*** NOTE: This function requires the 'shinydashboardPlus' package (can be installed with 'install.packages('shinydashboardPlus')'. You do not have this package installed.")
  }
  if(!"shinyWidgets" %in% rownames(installed.packages())){
    stop("*** NOTE: This function requires the 'shinyWidgets' package (can be installed with 'install.packages('shinyWidgets')'. You do not have this package installed.")
  }
  if(!"fresh" %in% rownames(installed.packages())){
    stop("*** NOTE: This function requires the 'fresh' package (can be installed with 'install.packages('fresh')'. You do not have this package installed.")
  }
  if(!"plotly" %in% rownames(installed.packages())){
    stop("*** NOTE: This function requires the 'plotly' package (can be installed with 'install.packages('plotly')'. You do not have this package installed.")
  }
  require(shiny)
  require(shinydashboardPlus)
  require(shinydashboard)
  require(shinyWidgets)
  require(plotly)
  
  runApp(list(
    ui = dashboardPage(skin = "purple",
          dashboardHeader(title = NULL), 
          dashboardSidebar(
            tags$head(tags$style(HTML('.content-wrapper { overflow: auto; }'))),
            sidebarMenu(id = "sidebar_tabs",
                        menuItem("Main", tabName = "main", icon = icon("home", lib = "glyphicon"))
            ), # end sidebarMenu
            # dashboardSidebar options
            collapsed = TRUE
          ), # end dashboardSidebar
          dashboardBody(
            tabItems(
              tabItem(tabName = "main",
                tags$style(HTML(" 
                .box.box-solid.box-primary>.box-header{
                  color:#fff;
                  background:#3F4248
                }
                .box.box-solid.box-primary{
                  border-bottom-color:#3F4248;
                  border-left-color:#3F4248;
                  border-right-color:#3F4248;
                  border-top-color:#3F4248;
                }
                .box.box-solid.box-info>.box-header{
                  color:#000000;
                  background:#CECECE
                }
                .box.box-solid.box-info{
                  border-bottom-color:#CECECE;
                  border-left-color:#CECECE;
                  border-right-color:#CECECE;
                  border-top-color:#CECECE;
                }")),
                box(title = "SPEARobject | Overview", width = 12, status = "primary", solidHeader = TRUE, collapsible = FALSE,
                    fluidRow(
                      column(6,
                             span(htmlOutput("table_spearobj_overview"), style = "font-size: 14px; text-align:left;")
                      ),
                      column(6,
                             tableOutput("table_spearobj_datasets")
                      )
                    )
                ),
                box(title = "SPEARobject | Mean CV Loss", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                  fluidRow(
                    column(4, align = "center",
                      tableOutput("table_spearobj_w_loss")
                    ),
                    column(8,
                      plotlyOutput("plot_spearobj_w_loss", width = "100%")
                    )
                  )
                ),
                column(12, align = "center", uiOutput("ui_SPEARmodel_choices")),
                conditionalPanel(condition = 'input.var_spearmodel_w != "Choose:"',
                  box(title = "SPEARmodel | Factor Scores", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                      fluidRow(
                        column(4, align = "center", 
                               radioGroupButtons(
                                 inputId = "var_spearmodel_factorscores_forecast", label = "forecast", 
                                 choices = c("in.sample", "out.of.sample"), 
                                 justified = TRUE, status = "primary",
                                 selected = "out.of.sample", size = "sm",
                                 individual = FALSE
                               )),
                        column(4, align = "center",
                               radioGroupButtons(
                                 inputId = "var_spearmodel_factorscores_grouptype", label = "Groups?", 
                                 choices = c("None", "Feature", "Response"), 
                                 justified = TRUE, status = "primary",
                                 selected = "None", size = "sm",
                                 individual = FALSE
                               )),
                        column(4, align = "center", 
                               uiOutput("ui_var_spearmodel_factorscores_groups"))
                      ),
                      fluidRow(
                        column(12,
                               plotOutput("plot_spearmodel_factorscores", width = "100%")
                        )
                      )
                  ), 
                  box(title = "Coefficients", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                      fluidRow(
                        column(4,
                               h3("holder")
                        ),
                        column(4,
                               plotOutput("plot_spearmodel_loadings", width = "100%")
                        ),
                        column(4,
                               plotOutput("plot_spearmodel_coefficients", width = "100%")
                        )
                      )
                  ) 
                )
              )
            ) #end tabItems
          ) #end dashboardBody
    ),
    server = function(input, output, session) {
      
      values <- reactiveValues()
      
      bold <- function(string){
        return(paste0("<b>",string,"</b>"))
      }
      color <- function(string, color){
        return(paste0("<font color = \"", color, "\">", string, "</font>"))
      }
      
      SPEARmodel <- reactive({
        if(input$var_spearmodel_w == "Choose:"){
          return(NULL)
        }
        SPEARmodel <- SPEAR::get_SPEAR_model(SPEARobj, w = input$var_spearmodel_w)
        values$SPEARmodel_is_generated <- TRUE
        return(SPEARmodel)
      })
      
      output$ui_SPEARmodel_choices <- renderUI({
        radioGroupButtons(
          inputId = "var_spearmodel_w", label = "SPEARmodel Weight (w)?", 
          choices = c("Choose:", SPEARobj$params$weights), 
          justified = TRUE, status = "primary",
          selected = "Choose:", size = "normal",
          individual = FALSE
        )
      })
      
      output$ui_var_spearmodel_factorscores_groups <- renderUI({
        shiny::validate(shiny::need(input$var_spearmodel_factorscores_grouptype != "None", ""))
        
        choicelist <- list()
        if(input$var_spearmodel_factorscores_grouptype == "Feature"){
          for(dataset in names(SPEARobj$data$xlist)){
            choicelist[[dataset]] <- colnames(SPEARobj$data$xlist[[dataset]])
          }
        }
        if(input$var_spearmodel_factorscores_grouptype == "Response"){
          for(response in colnames(SPEARobj$data$Y)){
            choicelist[[response]] <- response
          }
        }
        
        selectizeInput(
          inputId = "var_spearmodel_factorscores_xaxis", label = "Groups?", 
          choices = choicelist
        )
      })
      
      output$table_spearobj_overview <- renderUI({
        output <- "SPEARobject"
        return(output)
      })
      
      output$table_spearobj_datasets <- renderTable({
        dt <- data.frame(Dataset = character(0),
                         Samples = numeric(0),
                         Features = numeric(0))
        
        for(i in 1:length(SPEARobj$data$xlist)){
          dt <- rbind(dt, c(names(SPEARobj$data$xlist)[i], nrow(SPEARobj$data$xlist[[i]]), ncol(SPEARobj$data$xlist[[i]])))
        }
        colnames(dt) <- c("Dataset", "Samples", "Features")
        dt$Dataset <- sapply(dt$Dataset, function(dataset){
          return(bold(color(dataset, SPEARobj$params$colors$X[[dataset]])))
        })
        return(dt)
      }, sanitize.text.function = function(x) x)
      
      output$plot_spearobj_w_loss <- renderPlotly({
        g <- SPEAR::SPEAR.plot_cv_loss(SPEARobj,
                                       show.w.labels = TRUE, show.sd = TRUE, plot.per.response = TRUE, show.min.w.line = TRUE, show.overall = TRUE
                                       ) + ggplot2::ggtitle(NULL) + ggplot2::xlab("w")
        # GGPlotly:
        ggplotly(
          g
        ) %>% config(displayModeBar = F) %>% layout(xaxis=list(fixedrange=TRUE), yaxis=list(fixedrange=TRUE)) %>% layout(showlegend = FALSE)
      })
      
      output$table_spearobj_w_loss <- renderTable({
        mat <- SPEARobj$cv.eval$cvm
        mat <- cbind(SPEARobj$params$weights, mat)
        colnames(mat) <- c("W", colnames(SPEARobj$data$Y))
        return(mat)
      })
      
      output$plot_spearmodel_loadings <- renderPlot({
        SPEARmodel <- SPEARmodel()
        shiny::validate(shiny::need(!is.null(values$SPEARmodel_is_generated), ""))
        g <- SPEAR::SPEAR.plot_factor_coefficients(SPEARmodel, type = "regression")
        g
      })
      
      output$plot_spearmodel_coefficients <- renderPlot({
        SPEARmodel <- SPEARmodel()
        shiny::validate(shiny::need(!is.null(values$SPEARmodel_is_generated), ""))
        g <- SPEAR::SPEAR.plot_factor_coefficients(SPEARmodel, type = "projection")
        g
      })
      
      output$plot_spearmodel_factorscores <- renderPlot({
        SPEARmodel <- SPEARmodel()
        shiny::validate(shiny::need(!is.null(values$SPEARmodel_is_generated), ""))
        shiny::validate(shiny::need(!is.null(input$var_spearmodel_factorscores_grouptype), ""))
        if(input$var_spearmodel_factorscores_grouptype == "None"){
          custom.groups <- NULL
        } else if(input$var_spearmodel_factorscores_grouptype == "Feature"){
          custom.groups <- SPEARmodel$data$X[,which(colnames(SPEARmodel$data$X) == input$var_spearmodel_factorscores_xaxis)]
          names(custom.groups) <- rownames(SPEARmodel$data$Y)
        } else if(input$var_spearmodel_factorscores_grouptype == "Response"){
          custom.groups <- SPEARmodel$data$Y[,which(colnames(SPEARmodel$data$Y) == input$var_spearmodel_factorscores_xaxis)]
          names(custom.groups) <- rownames(SPEARmodel$data$Y)
        }
        g <- SPEAR::SPEAR.plot_factor_scores(SPEARmodel, groups = custom.groups, forecast = input$var_spearmodel_factorscores_forecast)
        g
      })
      
      
      
      
      
    }
  ))
  
}






