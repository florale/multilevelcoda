#' Interface to \pkg{shinystan}
#'
#' Provide an interface to \pkg{shinystan} for models fitted with \pkg{brms}
#'
#' @aliases launch_shinystan
#'
#' @param object A fitted model object of class \code{brmcoda}.
#' @param ... Optional arguments to pass to
#' \code{\link[shinystan:launch_shinystan]{launch_shinystan}} or 
#' \code{\link[shiny:runApp]{runApp}}.
#'
#' @seealso \code{\link[shinystan:launch_shinystan]{launch_shinystan}}
#'
#' @return An S4 shinystan object
#'
#' @method launch_shinystan brmcoda
#' @importFrom shinystan launch_shinystan
#' @export
launch_shinystan.brmcoda <- function(object, ...) {
  launch_shinystan(object$model, ...)
}

#' multilevelcoda Simulation Study Results
#'
#' Provide the full results for a simulation study testing the performance of \pkg{multilevelcoda}
#'
#' @return An S4 shiny object
#'
#' @importFrom plotly plotlyOutput
#' @importFrom plotly renderPlotly
#' @importFrom shiny br
#' @importFrom shiny column
#' @importFrom shiny conditionalPanel
#' @importFrom shiny fluidPage
#' @importFrom shiny fluidRow
#' @importFrom shiny h5
#' @importFrom shiny HTML
#' @importFrom shiny imageOutput
#' @importFrom shiny mainPanel
#' @importFrom shiny navbarPage
#' @importFrom shiny radioButtons
#' @importFrom shiny renderImage
#' @importFrom shiny selectInput
#' @importFrom shiny shinyApp
#' @importFrom shiny sidebarLayout
#' @importFrom shiny sidebarPanel
#' @importFrom shiny shinyOptions
#' @importFrom shiny tags
#' @importFrom bslib bs_add_variables
#' @importFrom bslib bs_add_variables
#' @importFrom bslib bs_theme
#' @importFrom bslib nav_item
#' @importFrom bslib nav_menu
#' @importFrom bslib nav_panel
#' @importFrom bslib nav_spacer
#' @importFrom DT renderDataTable
#' @importFrom DT datatable
#' @importFrom DT dataTableOutput
#' @importFrom DT DTOutput
#' @importFrom utils data
#' @export
multilevelcoda_sim <- function() {

  # shiny input
  sim <- data(sim, envir = environment())

  brmcoda_tab <- sim[["brmcoda_tab"]]
  sub_tab <- sim[["sub_tab"]]

  brmcoda_d3 <- sim[["brmcoda_plot"]][["brmcoda_d3"]]
  brmcoda_d4 <- sim[["brmcoda_plot"]][["brmcoda_d4"]]
  brmcoda_d5 <- sim[["brmcoda_plot"]][["brmcoda_d5"]]

  sub_d3 <- sim[["sub_plot"]][["sub_d3"]]
  sub_d4 <- sim[["sub_plot"]][["sub_d4"]]
  sub_d5 <- sim[["sub_plot"]][["sub_d5"]]

  ## theme
  shinyOptions(bslib = TRUE)
  theme <- bs_theme(
    # version = 5,
    bootswatch = "minty",
    bg = "#fff",
    fg = "#1C1718",
    primary = "#90766B",
    secondary = "#A1B2C2",
    success = "#708885",
    info = "#b76e79",
    warning = "#FAD899",
    danger = "#944C4C",
    # base_font = "Roboto",
    code_font = "Consolas",
    heading_font = NULL,
    font_scale = NULL,
  )
  floras <-
    bs_add_variables(
      theme,
      .where = "defaults",
      "body-color" = "#58504C",
      "panel-inner-border" = "#CFDAE2",
      "navbar-light-brand-color" = "#FAF7F3",
      "navbar-light-active-color" = "#F1DDCF",
      "navbar-light-hover-color" = "#d18d9a",
      "navbar-light-bg" = "#665C58", #90766B
      "navbar-default-link-hover-bg" = "#708885",
      "navbar-light-panel-bg" = "#A1B2C2",

      "table-bg" = "#EFE3E0",   #E5E8E1 DFD7D6
      "table-accent-bg" = "#fff",
      "table-light-striped-bg" = "#FAF7F3",
      "table-light-striped-bg-active" = "#CFDAE2",
      "table-bg-hover" = "#DFD7D6",

      "pagination-color" = "#3d251e",
      "pagination-bg" = "#BEC7B4",
      "pagination-border-color" = "#BEC7B4",

      "pagination-hover-colour" = "#FAF7F3",
      "pagination-hover-bg" = "#708885",
      "pagination-hover-border-color" = "#708885",

      "pagination-active-bg" = "#8DA290",
      "pagination-active-border-color" = "#8DA290",

      "pagination-disabled-bg" = "#5A6367",
      "pagination-disabled-border-color" = "#5A6367",
      "pagination-disabled-color" = "#708885",

      "navbar-light-table-striped-color" = "#1C1718",
      "navbar-light-table-striped-bg" = "#CFDAE2",
      "navbar-light-table-hover-bg" = "#CFDAE2",

      "dropdown-bg" = "#fff",

      # "table-striped-bg" = "#FAF7F3",
      # "table-hover-bg-factor" = "1.25",
      "btn-default-bg" = "#A1B2C2",
      "btn-success-bg" = "#A1B2C2",
      "input-bg" = "#fff",


      "legend-color" = "#DCD5CE",
      "border-primary" = "#3d251e",
      "bg-primary" = "#3d251e",



      "table-active-bg" = "#F0DBDE  !important"
      # "bs-table-hover-bg" = "#CFDAE2",
      #
      # "bs-table-bg-accent" = "#CFDAE2"

      # "inverse-bg" = "#DCD5CE"

    )

  # shiny server  ----------------------
  server <- function(input, output, session) {

    # Summary Statistics -------------------
    ## brmcoda tab ---------------
    output$simsum_brmcoda_table <- DT::renderDataTable(DT::datatable({
      if (input$par_brmcoda == "All") {
        par_brmcoda <- levels(brmcoda_tab$by)
      } else if (input$par_brmcoda == "b_Intercept") {
        par_brmcoda  <- "  b0"
      } else if (input$par_brmcoda == "b_bilr1") {
        par_brmcoda <- "between ilr1 beta"
      } else if (input$par_brmcoda == "b_bilr2") {
        par_brmcoda <- "between ilr2 beta"
      } else if (input$par_brmcoda == "b_bilr3") {
        par_brmcoda <- "between ilr3 beta"
      } else if (input$par_brmcoda == "b_bilr4") {
        par_brmcoda <- "between ilr4 beta"
      } else if (input$par_brmcoda == "b_wilr1") {
        par_brmcoda <- "within ilr1 beta"
      } else if (input$par_brmcoda == "b_wilr2") {
        par_brmcoda <- "within ilr2 beta"
      } else if (input$par_brmcoda == "b_wilr3") {
        par_brmcoda <- "within ilr3 beta"
      } else if (input$par_brmcoda == "b_wilr4") {
        par_brmcoda <- "within ilr4 beta"
      } else if (input$par_brmcoda == "sd_ID_Intercept") {
        par_brmcoda <- "  u0"
      } else if (input$par_brmcoda == "sigma") {
        par_brmcoda <- "  sigma"
      }
      if (input$rint_sd_brmcoda == "medium" & input$res_sd1_brmcoda == "medium") {
        brmcoda_tab <- brmcoda_tab[condition == "base"]
      } else if (input$rint_sd_brmcoda == "medium" & input$res_sd1_brmcoda == "small") {
        brmcoda_tab <- brmcoda_tab[condition == "REbase_RESsmall"]
      } else if (input$rint_sd_brmcoda == "medium" & input$res_sd1_brmcoda == "large") {
        brmcoda_tab <- brmcoda_tab[condition == "REbase_RESlarge"]
      } else if (input$rint_sd_brmcoda == "small" & input$res_sd2_brmcoda == "large") {
        brmcoda_tab <- brmcoda_tab[condition == "REsmall_RESlarge"]
      } else if (input$rint_sd_brmcoda == "large" & input$res_sd3_brmcoda == "small") {
        brmcoda_tab <- brmcoda_tab[condition == "RElarge_RESsmall"]
      }
      if (input$N_brmcoda != "All") {
        brmcoda_tab <- brmcoda_tab[N == input$N_brmcoda]
      }
      if (input$K_brmcoda != "All") {
        brmcoda_tab <- brmcoda_tab[K == input$K_brmcoda]
      }
      if (input$D_brmcoda != "All") {
        brmcoda_tab <- brmcoda_tab[D == input$D_brmcoda]
      }

      brmcoda_tab <-
        brmcoda_tab[by %in% par_brmcoda &
                      Stat == input$stat_brmcoda, .(Stat,
                                                    Estimates,
                                                    MCSE,
                                                    N,
                                                    K,
                                                    D,
                                                    # sd_ID_Intercept,
                                                    sigma,
                                                    OnTarget)]

      # brmcoda_tab[] <- lapply(brmcoda_tab, function(x) if(is.numeric(x)) round(x, 2) else x)

      brmcoda_tab

    }, options = list(
      info = T,
      searching = T,
      paging = T,
      autoWidth = F,
      scrollX = T
    )))

    ## substitution tab ------------
    output$simsum_sub_table <- DT::renderDataTable(DT::datatable({

      if (input$rint_sd_sub == "medium" & input$res_sd1_sub == "medium") {
        sub_tab <- sub_tab[condition == "base"]
      } else if (input$rint_sd_sub == "medium" & input$res_sd1_sub == "small") {
        sub_tab <- sub_tab[condition == "REbase_RESsmall"]
      } else if (input$rint_sd_sub == "medium" & input$res_sd1_sub == "large") {
        sub_tab <- sub_tab[condition == "REbase_RESlarge"]
      } else if (input$rint_sd_sub == "small" & input$res_sd2_sub == "large") {
        sub_tab <- sub_tab[condition == "REsmall_RESlarge"]
      } else if (input$rint_sd_sub == "large" & input$res_sd3_sub == "small") {
        sub_tab <- sub_tab[condition == "RElarge_RESsmall"]
      }
      if (input$N_sub != "All") {
        sub_tab <- sub_tab[N == input$N_sub]
      }
      if (input$K_sub != "All") {
        sub_tab <- sub_tab[K == input$K_sub]
      }
      if (input$D_sub != "All") {
        sub_tab <- sub_tab[D == input$D_sub]
      }
      if (input$D_sub == 3 & input$delta3_sub != "All") {
        sub_tab <- sub_tab[To == input$delta3_sub]
      } else if (input$D_sub == 4 & input$delta4_sub != "All") {
        sub_tab <- sub_tab[To == input$delta4_sub]
      } else if (input$D_sub == 5 & input$delta5_sub != "All") {
        sub_tab <- sub_tab[To == input$delta5_sub]
      }
      if (input$level_sub == "Between") {
        sub_tab <- sub_tab[Level == "between"]
      } else if (input$level_sub == "Within") {
        sub_tab <- sub_tab[Level == "within"]
      }

      sub_tab <-
        sub_tab[Stat == input$stat_sub, .(Stat,
                                          Estimates,
                                          MCSE,
                                          N,
                                          K,
                                          D,
                                          # sd_ID_Intercept,
                                          sigma,
                                          OnTarget)]

      # sub_tab[] <- lapply(sub_tab, function(x) if(is.numeric(x)) round(x, 2) else x)

      sub_tab

    }, options = list(
      info = T,
      searching = T,
      paging = T,
      autoWidth = F,
      scrollX = T
    )))

    # Summary Plots -------------------
    # ranges_brmcoda <- reactiveValues(x = NULL, y = NULL)

    ## brmcoda plot ---------------
    output$simsum_brmcoda_plot <- renderPlotly({

      if (input$rint_sd_brmcoda_plot == "medium" & input$res_sd1_brmcoda_plot == "medium") {
        if (input$D_brmcoda_plot == 3) {
          .par_plot(brmcoda_d3[Stat == input$stat_brmcoda_plot & condition == "base"], shiny = TRUE)
        } else if (input$D_brmcoda_plot == 4) {
          .par_plot(brmcoda_d4[Stat == input$stat_brmcoda_plot & condition == "base"], shiny = TRUE)
        } else if (input$D_brmcoda_plot == 5) {
          .par_plot(brmcoda_d5[Stat == input$stat_brmcoda_plot & condition == "base"], shiny = TRUE)
        }
      } else if (input$rint_sd_brmcoda_plot == "medium" & input$res_sd1_brmcoda_plot == "small") {
        if (input$D_brmcoda_plot == 3) {
          .par_plot(brmcoda_d3[Stat == input$stat_brmcoda_plot & condition == "REbase_RESsmall"], shiny = TRUE)
        } else if (input$D_brmcoda_plot == 4) {
          .par_plot(brmcoda_d4[Stat == input$stat_brmcoda_plot & condition == "REbase_RESsmall"], shiny = TRUE)
        } else if (input$D_brmcoda_plot == 5) {
          .par_plot(brmcoda_d5[Stat == input$stat_brmcoda_plot & condition == "REbase_RESsmall"], shiny = TRUE)
        }
      } else if (input$rint_sd_brmcoda_plot == "medium" & input$res_sd1_brmcoda_plot == "large") {
        if (input$D_brmcoda_plot == 3) {
          .par_plot(brmcoda_d3[Stat == input$stat_brmcoda_plot & condition == "REbase_RESlarge"], shiny = TRUE)
        } else if (input$D_brmcoda_plot == 4) {
          .par_plot(brmcoda_d4[Stat == input$stat_brmcoda_plot & condition == "REbase_RESlarge"], shiny = TRUE)
        } else if (input$D_brmcoda_plot == 5) {
          .par_plot(brmcoda_d5[Stat == input$stat_brmcoda_plot & condition == "REbase_RESlarge"], shiny = TRUE)
        }
      } else if (input$rint_sd_brmcoda_plot == "small" & input$res_sd2_brmcoda_plot == "large") {
        if (input$D_brmcoda_plot == 3) {
          .par_plot(brmcoda_d3[Stat == input$stat_brmcoda_plot & condition == "REsmall_RESlarge"], shiny = TRUE)
        } else if (input$D_brmcoda_plot == 4) {
          .par_plot(brmcoda_d4[Stat == input$stat_brmcoda_plot & condition == "REsmall_RESlarge"], shiny = TRUE)
        } else if (input$D_brmcoda_plot == 5) {
          .par_plot(brmcoda_d5[Stat == input$stat_brmcoda_plot & condition == "REsmall_RESlarge"], shiny = TRUE)
        }
      } else if (input$rint_sd_brmcoda_plot == "large" & input$res_sd3_brmcoda_plot == "small") {
        if (input$D_brmcoda_plot == 3) {
          .par_plot(brmcoda_d3[Stat == input$stat_brmcoda_plot & condition == "RElarge_RESsmall"], shiny = TRUE)
        } else if (input$D_brmcoda_plot == 4) {
          .par_plot(brmcoda_d4[Stat == input$stat_brmcoda_plot & condition == "RElarge_RESsmall"], shiny = TRUE)
        } else if (input$D_brmcoda_plot == 5) {
          .par_plot(brmcoda_d5[Stat == input$stat_brmcoda_plot & condition == "RElarge_RESsmall"], shiny = TRUE)
        }
      }

    })

    ## substitution plot ---------------
    output$simsum_sub_plot <- renderPlotly({

      if (input$rint_sd_sub_plot == "medium" & input$res_sd1_sub_plot == "medium") {
        if (input$D_sub_plot == 3) {
          .par_plot(sub_d3[Stat == input$stat_sub_plot & condition == "base"], shiny = TRUE)
        } else if (input$D_sub_plot == 4) {
          .par_plot(sub_d4[Stat == input$stat_sub_plot & condition == "base"], shiny = TRUE)
        } else if (input$D_sub_plot == 5) {
          .par_plot(sub_d5[Stat == input$stat_sub_plot & condition == "base"], shiny = TRUE)
        }
      } else if (input$rint_sd_sub_plot == "medium" & input$res_sd1_sub_plot == "small") {
        if (input$D_sub_plot == 3) {
          .par_plot(sub_d3[Stat == input$stat_sub_plot & condition == "REbase_RESsmall"], shiny = TRUE)
        } else if (input$D_sub_plot == 4) {
          .par_plot(sub_d4[Stat == input$stat_sub_plot & condition == "REbase_RESsmall"], shiny = TRUE)
        } else if (input$D_sub_plot == 5) {
          .par_plot(sub_d5[Stat == input$stat_sub_plot & condition == "REbase_RESsmall"], shiny = TRUE)
        }
      } else if (input$rint_sd_sub_plot == "medium" & input$res_sd1_sub_plot == "large") {
        if (input$D_sub_plot == 3) {
          .par_plot(sub_d3[Stat == input$stat_sub_plot & condition == "REbase_RESlarge"], shiny = TRUE)
        } else if (input$D_sub_plot == 4) {
          .par_plot(sub_d4[Stat == input$stat_sub_plot & condition == "REbase_RESlarge"], shiny = TRUE)
        } else if (input$D_sub_plot == 5) {
          .par_plot(sub_d5[Stat == input$stat_sub_plot & condition == "REbase_RESlarge"], shiny = TRUE)
        }
      } else if (input$rint_sd_sub_plot == "small" & input$res_sd2_sub_plot == "large") {
        if (input$D_sub_plot == 3) {
          .par_plot(sub_d3[Stat == input$stat_sub_plot & condition == "REsmall_RESlarge"], shiny = TRUE)
        } else if (input$D_sub_plot == 4) {
          .par_plot(sub_d4[Stat == input$stat_sub_plot & condition == "REsmall_RESlarge"], shiny = TRUE)
        } else if (input$D_sub_plot == 5) {
          .par_plot(sub_d5[Stat == input$stat_sub_plot & condition == "REsmall_RESlarge"], shiny = TRUE)
        }
      } else if (input$rint_sd_sub_plot == "large" & input$res_sd3_sub_plot == "small") {
        if (input$D_sub_plot == 3) {
          .par_plot(sub_d3[Stat == input$stat_sub_plot & condition == "RElarge_RESsmall"], shiny = TRUE)
        } else if (input$D_sub_plot == 4) {
          .par_plot(sub_d4[Stat == input$stat_sub_plot & condition == "RElarge_RESsmall"], shiny = TRUE)
        } else if (input$D_sub_plot == 5) {
          .par_plot(sub_d5[Stat == input$stat_sub_plot & condition == "RElarge_RESsmall"], shiny = TRUE)
        }
      }

    })

    output$coda <- renderImage({
      filename <- normalizePath(file.path('./inst',
                                          paste('coda.png')))

      # Return a list containing the filename and alt text
      list(src = filename,
           heigh = "40px", width = "50px")

    }, deleteFile = FALSE)
    # observe(session$setCurrentTheme(
    #   if (isTRUE(input$dark_mode)) dark else light
    # ))
    # observe({
    #   brush <- input$simsum_brmcoda_plot_brush
    #   if (!is.null(brush)) {
    #     ranges_brmcoda$x <- c(brush$xmin, brush$xmax)
    #     ranges_brmcoda$y <- c(brush$ymin, brush$ymax)
    #
    #   } else {
    #     ranges_brmcoda$x <- NULL
    #     ranges_brmcoda$y <- NULL
    #   }
    # })

  }

  # shiny ui -----------------
  ui <- fluidPage(
    theme = floras,
    tags$style(
      HTML(
        "
        @media (min-width: 768px) {
            body > div .container-fluid {
                width: 750px;
            }
        }
        @media (min-width: 992px) {
            body > div > .container-fluid {
                width: 970px;
            }
        }
        @media (min-width: 1200px) {
            body > div .container-fluid {
                width: 1170px;
            }
        }
        body > div > .container-fluid:nth-of-type(1) {
            margin: 0 auto;
            padding-top: 55px;
        }
        body > div > nav .nav.navbar-nav {
            float: right;
        }
        .navbar-inner {
            height: 80px;
        }
        .table.dataTable tbody td.active, .table.dataTable tbody tr.active td {
            background-color: #d18d9a !important;
        }
        .dataTables_wrapper .dataTables_filter input {
                      width: 30px;
                      background-color: #b68f90;
        }
        .well {
            background-color:#CFDAE2;
        }
        "
      )
    ),
    br(),
    br(),
    br(),

    # headerPanel(title = span(img(src = "multilevelcoda-sims-shiny/www/coda.png", height = 35), "test")),

    # .tabs-above > .nav > li[class=active] > a {
    #   background-color: #6171a9;
    #     color: #FFF;
    # }

    # .dataTables_wrapper .dataTables_length {
    #   float: right;}
    # .dataTables_wrapper .dataTables_filter {
    #   float: right;
    #   text-align: right;}
    # theme = shinytheme("sandstone"),
    # theme = "shinythemes/css/sandstone.min.css",
    # theme = shiny::bootstrapLib(),
    # tags$head(includeCSS("multilevelcoda-sim-shiny/www/florastheme.css")),
    # checkboxInput("dark_mode", "Dark mode", FALSE),

    navbarPage(
      theme = floras,
      position = "fixed-top",
      # title(div(img(src = "coda.png")),
      "multilevelcoda Simulation Study",
      # nav_item(img(src = 'coda.png'), fillable = FALSE, height="30px"),
      # imageOutput("coda.png")),
      ## Simulation Summary -----------------------------
      nav_menu(
        "Summary Statistics",
        # icon = icon("table"),

        ### brmcoda ------------------
        nav_panel(
          "Bayesian Compositional Multilevel",
          fluid = TRUE,
          # titlePanel("Simulation Condition"),

          sidebarLayout(
            sidebarPanel(
              width = 3,
              style = "height:740px",
              h5("Simulation Condition"),
              br(),

              selectInput("N_brmcoda",
                          "Number of individuals (N):",
                          c("All",
                            30, 50, 360, 1200)),
              selectInput("K_brmcoda",
                          "Number of days (K):",
                          c("All",
                            3, 5, 7, 14)),
              radioButtons(
                "D_brmcoda",
                "Number of compositional parts (D):",
                c(3, 4, 5),
                selected = 4,
                inline = TRUE
              ),
              selectInput(
                "rint_sd_brmcoda",
                "Random Intercept variance:",
                c("medium", "small", "large")
              ),
              conditionalPanel(
                condition = "input.rint_sd_brmcoda == 'medium'",
                selectInput(
                  "res_sd1_brmcoda",
                  "Residual variance:",
                  c("medium", "small", "large")
                )
              ),
              conditionalPanel(
                condition = "input.rint_sd_brmcoda == 'small'",
                selectInput("res_sd2_brmcoda",
                            "Residual variance:",
                            c("large"))
              ),
              conditionalPanel(
                condition = "input.rint_sd_brmcoda == 'large'",
                selectInput("res_sd3_brmcoda",
                            "Residual variance:",
                            c("small"))
              )
              # selectInput("res_sd",
              #             "Residual variance:",
              #             c("medium", "small", "large"))
            ),
            mainPanel(
              width = 9,
              fluidRow(column(
                12,
                selectInput(
                  "stat_brmcoda",
                  "Performance Measure:",
                  c(
                    "Bias" = "bias",
                    "Coverage" = "cover",
                    "Bias-Eliminated Coverage" = "becover"
                    # "Empirical Standard Error" = "empse",
                    # "Mean-squared Error" = "mse"
                  ),
                  width = "100%"
                )
              )),
              fluidRow(column(
                12,
                selectInput(
                  "par_brmcoda",
                  "Parameter:",
                  c(
                    "All",
                    "b_Intercept",
                    "b_bilr1",
                    "b_bilr2",
                    "b_bilr3",
                    "b_bilr4",
                    "b_wilr1",
                    "b_wilr2",
                    "b_wilr3",
                    "b_wilr4",
                    "sd_ID_Intercept",
                    "sigma"
                  ),
                  width = "100%"
                )
              )),
              br(),
              DT::dataTableOutput("simsum_brmcoda_table")
            )
          )
        ),

        ### substitution ------------------
        nav_panel(
          "Bayesian Compositional Multilevel Substitution",
          fluid = TRUE,
          # titlePanel("Simulation Condition"),
          # icon = icon("arrow-right-arrow-left", lib = "font-awesome"),
          sidebarLayout(
            sidebarPanel(
              width = 3,
              style = "height:740px",
              h5("Simulation Condition"),
              br(),

              selectInput("N_sub",
                          "Number of individuals (N):",
                          c("All",
                            30, 50, 360, 1200)),
              selectInput("K_sub",
                          "Number of days (K):",
                          c("All",
                            3, 5, 7, 14)),
              radioButtons(
                "D_sub",
                "Number of compositional parts (D):",
                c(3, 4, 5),
                selected = 4,
                inline = TRUE
              ),
              selectInput(
                "rint_sd_sub",
                "Random Intercept variance:",
                c("medium", "small", "large")
              ),
              conditionalPanel(
                condition = "input.rint_sd_sub == 'medium'",
                selectInput(
                  "res_sd1_sub",
                  "Residual variance:",
                  c("medium", "small", "large")
                )
              ),
              conditionalPanel(
                condition = "input.rint_sd_sub == 'small'",
                selectInput("res_sd2_sub",
                            "Residual variance:",
                            c("large"))
              ),
              conditionalPanel(
                condition = "input.rint_sd_sub == 'large'",
                selectInput("res_sd3_sub",
                            "Residual variance:",
                            c("small"))
              )
              # selectInput("res_sd",
              #             "Residual variance:",
              #             c("medium", "small", "large"))
            ),
            mainPanel(
              width = 9,

              fluidRow(column(
                12,
                selectInput(
                  "stat_sub",
                  "Performance Measure:",
                  c(
                    "Bias" = "bias",
                    "Coverage" = "cover",
                    "Bias-Eliminated Coverage" = "becover"
                    # "Empirical Standard Error" = "empse",
                    # "Mean-squared Error" = "mse"
                  ),
                  width = "100%"
                )
              )),

              fluidRow(
                column(
                  6,
                  # selectInput(
                  #   "delta_sub",
                  #   "Substitution of:",
                  #   c("All")
                  # ),
                  conditionalPanel(
                    condition = "input.D_sub == '3'",
                    selectInput(
                      "delta3_sub",
                      "Substitution of:",
                      c(
                        "All",
                        "Sleep",
                        "Physical Activity" = "PA",
                        "Sedentary Behaviour" = "SB"
                      ),
                      width = "100%"
                    )
                  ),
                  conditionalPanel(
                    condition = "input.D_sub == '4'",
                    selectInput(
                      "delta4_sub",
                      "Substitution of:",
                      c(
                        "All",
                        "Sleep",
                        "Moderate-Vigorous Physical Activity" = "MVPA",
                        "Light Physical Activity" = "LPA",
                        "Sedentary Behaviour" = "SB"
                      ),
                      width = "100%"
                    )
                  ),
                  conditionalPanel(
                    condition = "input.D_sub == '5'",
                    selectInput(
                      "delta5_sub",
                      "Substitution of:",
                      c(
                        "All",
                        "Sleep" = "TST",
                        "Awake in Bed" = "WAKE",
                        "Moderate-Vigorous Physical Activity" = "MVPA",
                        "Light Physical Activity" = "LPA",
                        "Sedentary Behaviour" = "SB"
                      ),
                      width = "100%"
                    )
                  )
                ),
                column(
                  6,
                  selectInput(
                    "level_sub",
                    "Level:",
                    c("All", "Between", "Within"),
                    selected = "All",
                    # inline = TRUE,
                    width = "100%"
                  )
                )
              )
              ,
              br(),
              DT::DTOutput("simsum_sub_table")
            )
          )
        )


      ),
      ## Simulation Plots -----------------------------
      nav_menu(
        "Summary Plots",
        # icon = icon("chart-line", lib = "font-awesome"),

        ### brmcoda plots ------------------
        nav_panel(
          "Bayesian Compositional Multilevel Parameters",
          fluid = TRUE,

          sidebarLayout(
            sidebarPanel(
              width = 3,
              style = "height:740px",
              h5("Simulation Condition"),
              br(),

              selectInput(
                "N_brmcoda_plot",
                "Number of individuals (N):",
                c("All"
                  # ,
                  # 30, 50, 360, 1200
                )
              ),
              selectInput("K_brmcoda_plot",
                          "Number of days (K):",
                          c("All"
                            # ,
                            # 3, 5, 7, 14
                          )),
              radioButtons(
                "D_brmcoda_plot",
                "Number of compositional parts (D):",
                c(3, 4, 5),
                selected = 4,
                inline = TRUE,
                width = "100%"
              ),

              selectInput(
                "rint_sd_brmcoda_plot",
                "Random Intercept variance:",
                c("medium", "small", "large"),
                width = "100%"
              ),
              conditionalPanel(
                condition = "input.rint_sd_brmcoda_plot == 'medium'",
                selectInput(
                  "res_sd1_brmcoda_plot",
                  "Residual variance:",
                  c("medium", "small", "large"),
                  width = "100%"
                )
              ),
              conditionalPanel(
                condition = "input.rint_sd_brmcoda_plot == 'small'",
                selectInput(
                  "res_sd2_brmcoda_plot",
                  "Residual variance:",
                  c("large"),
                  width = "100%"
                )
              ),
              conditionalPanel(
                condition = "input.rint_sd_brmcoda_plot == 'large'",
                selectInput(
                  "res_sd3_brmcoda_plot",
                  "Residual variance:",
                  c("small"),
                  width = "100%"
                )
              )
            ),
            mainPanel(
              width = 9,

              fluidRow(column(
                12,
                selectInput(
                  "stat_brmcoda_plot",
                  "Performance Measure:",
                  c(
                    "Bias" = "bias",
                    "Coverage" = "cover",
                    "Bias-Eliminated Coverage" = "becover"
                    # "Empirical Standard Error" = "empse",
                    # "Mean-squared Error" = "mse"
                  ),
                  width = "100%"
                )
              )),

              fluidRow(column(
                12,
                selectInput(
                  "par_brmcoda_plot",
                  "Parameter:",
                  c(
                    "All"
                    # ,
                    # "b_Intercept",
                    # "b_bilr1",
                    # "b_bilr2",
                    # "b_bilr3",
                    # "b_bilr4",
                    # "b_wilr1",
                    # "b_wilr2",
                    # "b_wilr3",
                    # "b_wilr4",
                    # "sd_ID_Intercept",
                    # "sigma"
                  ),
                  width = "100%"
                )
              )),
              # br(),
              plotlyOutput("simsum_brmcoda_plot",
                           height = "1300px")
            )
          )
        ),

        ### substitution plots ------------------
        nav_panel(
          "Bayesian Compositional Multilevel Substitution Estimates",
          fluid = TRUE,
          # icon = icon("arrow-right-arrow-left", lib = "font-awesome"),

          sidebarLayout(
            sidebarPanel(
              width = 3,
              style = "height:740px",
              h5("Simulation Condition"),
              br(),

              selectInput(
                "N_sub_plot",
                "Number of individuals (N):",
                c("All"
                  # ,
                  # 30, 50, 360, 1200
                )),
              selectInput("K_suub_plot",
                          "Number of days (K):",
                          c("All"
                            # ,
                            # 3, 5, 7, 14
                          )),
              radioButtons(
                "D_sub_plot",
                "Number of compositional parts (D):",
                c(3, 4, 5),
                selected = 4,
                inline = TRUE,
                width = "100%"
              ),

              selectInput(
                "rint_sd_sub_plot",
                "Random Intercept variance:",
                c("medium", "small", "large"),
                width = "100%"
              ),
              conditionalPanel(
                condition = "input.rint_sd_sub_plot == 'medium'",
                selectInput(
                  "res_sd1_sub_plot",
                  "Residual variance:",
                  c("medium", "small", "large"),
                  width = "100%"
                )
              ),
              conditionalPanel(
                condition = "input.rint_sd_sub_plot == 'small'",
                selectInput("res_sd2_sub_plot",
                            "Residual variance:",
                            c("large"),
                            width = "100%")
              ),
              conditionalPanel(
                condition = "input.rint_sd_sub_plot == 'large'",
                selectInput("res_sd3_sub_plot",
                            "Residual variance:",
                            c("small"), width = "100%")
              )
            ),
            mainPanel(
              width = 9,

              fluidRow(column(
                12,
                selectInput(
                  "stat_sub_plot",
                  "Performance Measure:",
                  c(
                    "Bias" = "bias",
                    "Coverage" = "cover",
                    "Bias-Eliminated Coverage" = "becover"
                    # "Empirical Standard Error" = "empse",
                    # "Mean-squared Error" = "mse"
                  ),
                  width = "100%"
                )
              )),

              fluidRow(
                column(
                  6,
                  selectInput(
                    "delta_sub_plot",
                    "Substitution of:",
                    c("All"),
                    width = "100%"
                  ),
                  # conditionalPanel(
                  #   condition = "input.D_sub_plot == '3'",
                  #   selectInput(
                  #     "delta3_sub_plot",
                  #     "Substitution of:",
                  #     c(
                  #       "Sleep",
                  #       "Physical Activity" = "PA",
                  #       "Sedentary Behaviour" = "SB"
                  #     ),
                  #     width = "100%"
                  #   )
                  # ),
                  # conditionalPanel(
                  #   condition = "input.D_sub_plot == '4'",
                  #   selectInput(
                  #     "delta4_sub_plot",
                  #     "Substitution of:",
                  #     c(
                  #       "Sleep",
                  #       "Moderate-Vigorous Physical Activity" = "MVPA",
                  #       "Light Physical Activity" = "LPA",
                  #       "Sedentary Behaviour" = "SB"
                  #     ),
                  #     width = "100%"
                  #   )
                  # ),
                  # conditionalPanel(
                  #   condition = "input.D_sub_plot == '5'",
                  #   selectInput(
                  #     "delta5_sub_plot",
                  #     "Substitution of:",
                  #     c(
                  #       "Sleep" = "TST",
                  #       "Awake in Bed" = "WAKE",
                  #       "Moderate-Vigorous Physical Activity" = "MVPA",
                  #       "Light Physical Activity" = "LPA",
                  #       "Sedentary Behaviour" = "SB"
                  #     ),
                  #     width = "100%"
                  #   )
                  # )
                ),
                column(
                  6,
                  selectInput(
                    "level_sub_plot",
                    "Level:",
                    c("All"
                      # ,
                      # "Between", "Within"
                    ),
                    selected = "All",
                    # inline = TRUE,
                    width = "100%"
                  )
                )
              )
              ,
              # br(),
              plotlyOutput("simsum_sub_plot",
                           height = "1600px")
            )
          )
        )
      ),

      nav_spacer(),
      # nav_item(icon = icon("arrow-right-arrow-left", lib = "font-awesome")),
      nav_item(imageOutput("coda", height = "50px"), fillable = FALSE),


    )
  )
  # run app ------------------------
  shinyApp(ui = ui, server = server)


}
