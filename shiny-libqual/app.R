library(shiny)
library(data.table)
library(DT)
library(shinyFiles)

# List files to rate
imgs = list.files("/mnt/AchTeraD/data/BICRO278/MS102/plots/profileplots/500000", pattern=".png", full.names = TRUE)
imgs = imgs[file.exists(imgs) & file.info(imgs)$size > 0]

# Define UI
ui = fluidPage(
  
  # Set title and buttons
  titlePanel("Cell quality annotation"),
  sidebarLayout(
    sidebarPanel(
      actionButton("good", "Good"),
      actionButton("bad", "Bad"),
      actionButton("intermediate", "Intermediate"),
      actionButton("na", "NA")
    ),
    
    # Set image and data.table output
    mainPanel(
      imageOutput("image"),
      br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
      fluidRow(
        column(1, offset=1, actionButton("previous", "Previous")),
        column(1, offset=1, actionButton("next", "Next")),
        br(),
        DT::dataTableOutput("table", width = "75%")
      ),
      
    )
  )
)

# Server-side calculations
server = function(input, output, session) {

  # Initiate data.table
  dt = data.table(sample = imgs, user_quality = "NA")
  output$table = DT::renderDataTable(dt,
                                     rownames = FALSE,
                                     server = FALSE,
                                     extensions = c("Buttons"),
                                     options = list(dom = 'Bfrtip', buttons = c('csv')),
                                     width = "50%")
  # Set parameters
  index = reactiveVal(1)

  # Increase/decrease index after previous/next buttons
  observeEvent(input[["previous"]], {
    index(max(index()-1, 1))
  })
  observeEvent(input[["next"]], {
    index(min(index()+1, length(imgs)))
  })

  # Quality assessment buttons, index increase and output of updated DT
  observeEvent(input[["good"]], {
    i = index()
    dt[i, user_quality := "good"]
    if(index()+1 > length(imgs)) info(dt)
    index(min(index()+1, length(imgs)))
    output$table = renderDataTable(dt,
                                   rownames = FALSE,
                                   server = FALSE,
                                   extensions = c("Buttons"),
                                   options = list(dom = 'Bfrtip', buttons = c('csv')))
    })

  observeEvent(input[["bad"]], {
    i = index()
    dt[i, user_quality := "bad"]
    if(index()+1 > length(imgs)) info(dt)
    index(min(index()+1, length(imgs)))
    output$table = renderDataTable(dt,
                                   rownames = FALSE,
                                   server = FALSE,
                                   extensions = c("Buttons"),
                                   options = list(dom = 'Bfrtip', buttons = c('csv')))  })

  observeEvent(input[["intermediate"]], {
    i = index()
    dt[i, user_quality := "intermediate"]
    if(index()+1 > length(imgs)) info(dt)
    index(min(index()+1, length(imgs)))
    output$table = renderDataTable(dt,
                                   rownames = FALSE,
                                   server = FALSE,
                                   extensions = c("Buttons"),
                                   options = list(dom = 'Bfrtip', buttons = c('csv')))  })
  observeEvent(input[["na"]], {
    i = index()
    dt[i, user_quality := "NA"]
    if(index()+1 > length(imgs)) info(dt)
    index(min(index()+1, length(imgs)))
    output$table = renderDataTable(dt,
                                   rownames = FALSE,
                                   server = FALSE,
                                   extensions = c("Buttons"),
                                   options = list(dom = 'Bfrtip', buttons = c('csv')))  })

  # Output of new image after increase/decrease of index
  output$image = renderImage({
    x = imgs[index()]
    list(src = x, alt = "alternate text", width = "75%", heigth = "75%")
  }, deleteFile = FALSE)
}

# Run the application 
shinyApp(ui = ui, server = server)