library(shiny)
library(data.table)

imgs = list.files("/mnt/AchTeraD/data/BICRO278/MS102/plots/profileplots/500000", pattern=".png", full.names = TRUE)[1:10]
imgs = imgs[file.exists(imgs) & file.info(imgs)$size > 0]

ui = fluidPage(
  
  titlePanel("Slideshow"),
  sidebarLayout(
    sidebarPanel(
      actionButton("good", "Good"),
      actionButton("bad", "Bad"),
      actionButton("intermediate", "Intermediate"),
      br(),
      actionButton("table", "Refresh table")
    ),
    
    mainPanel(
      imageOutput("image"),
      br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), 
      fluidRow(
        column(1, offset=1, actionButton("previous", "Previous")),
        column(1, offset=1, actionButton("next", "Next")),
        br(),
        dataTableOutput("table")
      ),
      
    )
  )
)

server = function(input, output, session) {
  # Initiate data.table
  dt = data.table(sample = imgs, user_quality = "NA")
  output$table = renderDataTable(dt)
  
  # Set parameters
  index = reactiveVal(1)
  
  observeEvent(input[["previous"]], {
    index(max(index()-1, 1))
  })
  observeEvent(input[["next"]], {
    index(min(index()+1, length(imgs)))
  })
  
  observeEvent(input[["good"]], {
    i = index()
    dt[i, user_quality := "good"]
    if(index()+1 > length(imgs)) info(dt)
    index(min(index()+1, length(imgs)))
  })
  
  observeEvent(input[["bad"]], {
    i = index()
    dt[i, user_quality := "bad"]
    if(index()+1 > length(imgs)) info(dt)
    index(min(index()+1, length(imgs)))
  })
  
  observeEvent(input[["intermediate"]], {
    i = index()
    dt[i, user_quality := "intermediate"]
    if(index()+1 > length(imgs)) info(dt)
    index(min(index()+1, length(imgs)))
  })
  
  observeEvent(input[["table"]], {
    output$table = renderDataTable(dt,
                                   extensions = c("Buttons"), 
                                   options = list(dom = 'Bfrtip',
                                                  buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                                   ))
  })
  output$image = renderImage({
    x = imgs[index()] 
    list(src = x, alt = "alternate text", width = "75%", heigth = "75%")
  }, deleteFile = FALSE)
}

# Run the application 
shinyApp(ui = ui, server = server)