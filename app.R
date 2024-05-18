library(shiny)
library(bslib)
library(rsconnect)

custom_theme <- bs_theme(
  version = 5,
  bg = "#0F0F0F",
  fg = "#FDFEEE",
  primary = "#0199F8",
  secondary = "#FF374B",
  success = "#009E73",
  base_font = font_google("Inter"),
  code_font = font_google("JetBrains Mono")
)

varStationary <- function(timeseries) {
  listLambda <- list()
  sum <- 0
  if (min(timeseries) < 0) {
    sum <- round(floor(min(timeseries)), 0)
    timeseries <- timeseries - sum
  }
  lambda <- forecast::BoxCox.lambda(timeseries, lower = -1, upper = 1)
  listLambda <- list(round(lambda, 2))
  while (round(lambda, 2) != 1) {
    if (round(lambda, 2) < -0.75) {
      timeseries <- 1/timeseries
      listLambda[length(listLambda)] = round(-1, 2)
    } else if (round(lambda, 2) < -0.25) {
      timeseries <- 1/sqrt(timeseries)
      listLambda[length(listLambda)] = round(-0.5, 2)
    } else if (round(lambda, 2) < 0.25) {
      timeseries <- log(timeseries)
      listLambda[length(listLambda)] = round(0, 2)
    } else if (round(lambda, 2) < 0.75) {
      timeseries <- sqrt(timeseries)
      listLambda[length(listLambda)] = round(0.5, 2)
    } else {
      listLambda[length(listLambda)] = round(1, 2)
      break
    }
    lambda <- forecast::BoxCox.lambda(timeseries, lower = -1, upper = 1)
    listLambda[length(listLambda)+1] = round(lambda, 2)
  }
  return(list(data = timeseries, lambda = listLambda, summation = -sum))
}

predBack <- function(timeseries, lambda) {
  if (lambda > 0.75) {
    timeseries <- timeseries
  } else if (lambda > 0.25) {
    timeseries <- timeseries^2
  } else if (lambda > -0.25) {
    timeseries <- exp(timeseries)
  } else if (lambda > -0.75) {
    timeseries <- timeseries^-2
  } else {
    timeseries <- 1/timeseries
  }
  return(timeseries)
}

meanStationary <- function(timeseries) {
  adf <- tseries::adf.test(timeseries)
  d = 0
  while (adf$p.value > 0.05) {
    d = d + 1
    timeseries <- diff(timeseries, lag = 1)
    adf <- tseries::adf.test(timeseries)
  }
  return(list(data = timeseries, diffOrder = d))
}

# See above for the definitions of ui and server
ui <- fluidPage(
  title = "byF",
  theme = custom_theme,
  tags$br(),
  tags$div(tags$h6(tags$b("byF"))),
  tags$h2(tags$b("ARIMA")),
  tags$div("Autoregressive Integrated Moving Average (ARIMA) is used for time series forecasting with only one dependent variable. You could simply analyze and predict your data here."),
  tags$br(),
  tags$em(tags$b("Note: Different type of data may require different treatment.")),
  tags$br(),
  tags$br(),
  fileInput("file1", "Choose CSV file", multiple = F, accept = c(".csv")),
  tags$hr(),
  layout_columns(
    card(
      checkboxInput("header", "Header", T),
      radioButtons("dec", "Decimal separator",
                   c('Comma (,)' = ",", 'Punktum (.)'="."), "komma"),
    ),
    card(
      radioButtons("sep", "Separator",
                   c('Comma (,)' = ",", 'Semicolon (;)' = ";", 'Tab (\t)' = "\t"), "Comma")
    )
  ),
  conditionalPanel(
    condition = "output.fileUploaded",
    navset_card_underline(
      nav_panel(
        title = "Data Summary",
        verbatimTextOutput("summary")
      ),
      nav_panel(
        title = "Stationary Test",
        layout_columns(
          card(
            tags$em("Ensure your data is numeric and has at least 10 rows."),
            selectInput(inputId = "var", "Choose Column Variable", choices = c()),
            conditionalPanel(
              condition = "output.varAvail",
              tags$em('You could adjust the number of train dataset here.'),
              uiOutput("slider"), uiOutput("sliderTest"),
              downloadButton("downloadStatVar", "Download Stationary in Varian Data"),
              downloadButton("downloadStat", "Download Stationary Data (in Varian and Mean)")
            )
          ),
          card(
            conditionalPanel(
              condition = "output.varAvail",
              tags$em(tags$b("Negative Value Check")),
              textOutput("intersum"),
              tags$em(tags$b("Stationary in Varian Check")),
              textOutput("interstvar"),
              tags$em('Note: Box-Cox transformation formula available', tags$a(href = 'https://github.com/anisyasafiraa/ARIMA-GSTARIMA_Java-Island-Inflation/blob/main/Box-Cox%20Transformation.png', 'here.')),
              tags$em(tags$b("Stationary in Mean Check")),
              textOutput("interdiff"),
            )
          )
        )
      )
    )
  ),
  conditionalPanel(
    condition = "output.varAvail",
    tags$h4(tags$b("Time Order Identification")),
    layout_columns(
      card(
        tags$b("Plot ACF for MA order identification."),
        tags$div(tags$em("Choose MA order when lag is significant (cut off in ACF).")),
        plotOutput("acf")
      ),
      card(
        tags$b("Plot PACF for AR order identification."),
        tags$div(tags$em("Choose AR order when lag is significant (cut off in PACF).")),
        plotOutput("pacf")
      )
    ),
    tags$h4(tags$b("ARIMA Model")),
    layout_columns(
      card(
        numericInput("ar", "AR Order", 1, min = 0),
        numericInput("ma", "MA Order", 1, min = 0),
        numericInput("d", "Differencing", 0),
        tags$em(textOutput("diff")),
        uiOutput("const"),
        radioButtons("meth", "Estimation Parameter Methods",
                     c('ML'="ML", 'CSS'="CSS"), selected = "CSS")
      ),
      card(
        navset_card_underline(
          nav_panel(
            title = "Model",
            verbatimTextOutput("arima"),
          ),
          nav_panel(
            title = "Residual",
            tags$h6(tags$b("White Noise")),
            tags$div("Residual is white noise when p-value > 0.05."),
            verbatimTextOutput("wn"),
            tags$h6(tags$b("Normal Distribution")),
            selectInput(inputId = "norm", "Choose Normality Test Methods",
                        choices = c("Shapiro-Wilk" = "wilk",
                                    "Kolmogorov-Smirnov" = "ks",
                                    "Jarque-Bera" = "jb",
                                    "Anderson-Darling" = "ad")),
            tags$div("Residual fulfilled normality test assumption when p-value > 0.05."),
            verbatimTextOutput("normTest")
          ),
          nav_panel(
            title = "Forecast",
            numericInput("forecast", "Forecast", value = 12, min = 6, step = 1),
            dataTableOutput("forcResult"),
            tags$em("Download forecasting result down below."),
            downloadButton("downloadForc", "Download")
          ),
        ),
        conditionalPanel(
          condition = "output.testAvail",
          tags$h5(tags$b("Data Test Identification")),
          verbatimTextOutput("testSummary"),
          plotOutput('plot')
        )
      ),
      col_widths = c(2, 10)
    )
  )
)

server <- function(session, input, output) {
  data <- reactive({
    req(input$file1)
    inFile <- input$file1
    
    if(is.null(inFile))
      return(NULL)
    
    df <- read.csv(inFile$datapath, header = input$header,
                   sep = input$sep, dec = input$dec)
    return(df)
  })
  
  output$fileUploaded <- reactive({
    return(!is.null(data()))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  data1 <- reactive({
    df1 <- data()
    updateSelectInput(session,"var",choices=colnames(df1))
    return(df1)
  })
  
  output$varAvail <- reactive({
    return(is.numeric(data1()[,input$var]))
  })
  outputOptions(output, 'varAvail', suspendWhenHidden=FALSE)
  
  output$slider <- renderUI({
    dt <- data()[,input$var]
    sliderInput("train", "Train set", min = 10, max = length(dt), value = length(dt), step = 1)
  })
  
  output$sliderTest <- renderUI({
    dt <- data()[,input$var]
    if (is.null(input$train == length(dt))) {
      return(NULL)
    } else if (input$train < length(dt)) {
      sliderInput("test", "Test set", min = 0, max = length(dt)-input$train,
                  value = length(dt)-input$train, step = 1)
    }
  })
  
  output$testAvail <- reactive({
    return(!is.null(input$test))
  })
  outputOptions(output, 'testAvail', suspendWhenHidden=FALSE)
  
  stationaryVarian <- reactive({
    var <- data1()[1:input$train,input$var]
    dt <- varStationary(var)
    return(dt)
  })
  
  stationary <- reactive({
    dt <- stationaryVarian()
    stationary <- meanStationary(dt$data)
    return(stationary)
  })

  output$summary <- renderPrint({
    df <- data()
    summary(df)
  })
  
  output$intersum <- renderPrint({
    var <- data1()[1:input$train,input$var]
    if (min(var) < 0) {
      paste("This data has negative values thus it needs to be added by",
            -floor(min(var)), ".")
    } else {
      paste("This data has no negative value.")
    }
  })
  
  output$interstvar <- renderPrint({
    var <- stationaryVarian()
    if (length(var$lambda) >= 2) {
      paste("To be stationary in varian, this data should be transformed. List of lambda occured:", list(var$lambda))
    } else if (length(var$lambda) == 1) {
      paste("This data already stationary in varian.")
    }
  })
  
  output$interdiff <- renderPrint({
    var <- stationary()
    if (var$diffOrder >= 1) {
      paste("To be stationary in mean, this data should be differentiated",
            var$diffOrder, "times.")
    } else if (var$diffOrder == 0) {
      paste("This data already stationary in mean.")
    }
  })
  
  output$diff <- renderPrint({
    var <- stationary()
    paste("Diff order for this data should be", var$diffOrder)
  })
  
  output$downloadStatVar <- downloadHandler(
    filename = function(){"stationary in varian.csv"}, 
    content = function(fname){
      dt <- stationaryVarian()
      write.csv(dt$data, fname)
    }
  )
  
  output$downloadStat <- downloadHandler(
    filename = function(){"stationary.csv"}, 
    content = function(fname){
      dt <- stationary()
      write.csv(dt$data, fname)
    }
  )
  
  output$acf <- renderPlot({
    var <- stationary()
    data <- var$data
    stats::acf(data)
  })
  
  output$pacf <- renderPlot({
    var <- stationary()
    data <- var$data
    stats::pacf(data)
  })
  
  output$const <- renderUI({
    if (input$d == 0) {
      checkboxInput("const", "Using Constanta", F)
    }
  })
  
  arimaModel <- reactive({
    dt <- stationaryVarian()
    if (input$d == 0) {
      p <- input$const
    } else if (input$d > 0) {
      p <- FALSE
    }
    model <- forecast::Arima(dt$data, order = c(input$ar, input$d, input$ma),
                             include.constant = p, method = input$meth)
    return(model)
  })
  
  output$arima <- renderPrint({
    summary(arimaModel())
  })
  
  output$forcNum <- renderUI({
    dt <- stationaryVarian()
    if (input$train == length(dt$data)) {
      numericInput("forecast", "Forecast", value = 12, min = 6, step = 1)
    } else if (input$train < length(dt$data)) {
      numericInput("forecast", "Forecast", value = input$test,
                   min = input$test, step = 1)
    }
  })
  
  forc <- reactive({
    columns <- c("Time", "Forecast")
    forcResult <- data.frame(matrix(nrow = 0, ncol = length(columns)))
    colnames(forcResult) <- columns
    
    dt <- stationaryVarian()
    for (i in 1:length(dt$lambda)) {
      time <- 0
      if (is.null(input$test) == TRUE) {
        time <- input$forecast
      } else if (is.null(input$test) == FALSE) {
        if (input$test < input$forecast) {
          time <- input$forecast
        } else if (input$test > input$forecast) {
          time <- input$test
        }
      }
      pred <- predict(arimaModel(), n.ahead = time)$pred
      pred <- predBack(pred, dt$lambda[[length(dt$lambda)]])
      dt$lambda[[length(dt$lambda)]] <- NULL
    }
    
    pred <- pred - dt$summation
    
    for (i in 1:time) {
      time <- length(dt$data) + i
      forcResult[nrow(forcResult) + 1,] <- c(time, round(pred[i], 2))
    }
    return(forcResult)
  })
  
  output$testSummary <- renderPrint({
    res <- forc()
    forc <- as.numeric(res$Forecast)[!is.na(res$Forecast)]
    forcUsed <- forc[1:input$test]
    dt <- data1()[1:input$test,input$var]
    forecast::accuracy(forc, dt)
  })
  
  output$plot <- renderPlot({
    res <- forc()
    forc <- as.numeric(res$Forecast)[!is.na(res$Forecast)]
    forcUsed <- forc[1:input$test]
    dt <- data1()[1:input$test,input$var]
    matplot(1:length(dt), cbind(forc[1:length(dt)], dt), type = "l", lty = 1,
            col = c("red", "blue"), xlab = "Time", 
            ylab = input$var, main = "Comparison Actual and Forecast Data")
    legend("bottomleft", legend = c("Forecast", "Data Test"), 
           col = c("red", "blue"), lty = 1)
  })
  
  output$wn <- renderPrint({
    res <- resid(arimaModel())
    res <- na.omit(res)
    Box.test(res, type = "Ljung-Box")
  })
  
  output$normTest <- renderPrint({
    res <- resid(arimaModel())
    res <- na.omit(res)
    if (input$norm == "wilk") {
      stats::shapiro.test(res)
    } else if (input$norm == "ks") {
      stats::ks.test(res, "pnorm")
    } else if (input$norm == "jb") {
      DescTools::JarqueBeraTest(res)
    } else if (input$norm == "ad") {
      nortest::ad.test(res)
    }
  })
  
  output$forcResult <- renderDataTable(forc()[1:input$forecast,])
  
  output$downloadForc <- downloadHandler(
    filename = function(){"forecast result.csv"}, 
    content = function(fname){
      write.csv(forc(), fname)
    }
  )
}

shinyApp(ui = ui, server = server)