library(shiny)
library(viridis)
library(hrbrthemes)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
# ressources Apps : https://github.com/phillyo/intelligentsia/blob/master/ui.R

# ----------------------------------
# load dataset
dataset1 <- readr::read_tsv("samples.tsv") 
 


ui <- navbarPage(title = "uBiome",
                 theme = "style/style.css",
                 footer = includeHTML("footer.html"),
                 fluid = TRUE, 
                 collapsible = TRUE,
                 
                 # ----------------------------------
                 # tab panel 1 - Home
                 tabPanel("Home",
                          includeHTML("home.html"),
                          tags$script(src = "plugins/scripts.js"),
                          tags$head(
                            tags$link(rel = "stylesheet", 
                                      type = "text/css", 
                                      href = "plugins/font-awesome-4.7.0/css/font-awesome.min.css"),
                            tags$link(rel = "icon", 
                                      type = "image/png", 
                                      href = "images/logo_icon.png")
                          )
                 ),
                 # ----------------------------------
                 # tab panel 2 - My profile
                 tabPanel(title = "My profile",
                          navlistPanel(
                            # ----------------------------------
                            # tab panel 2-A - Single sample analysis
                            tabPanel("Single sample", 
                                       sidebarLayout(
                                         sidebarPanel(
                                           selectInput("sample", "sample selection:",
                                                       c("sample 1" = "sample1",
                                                         "sample 2" = "sample2",
                                                         "sample 3" = "sample3",
                                                         "sample 4" = "sample4")),
                                           radioButtons("pType", "Choose plot type:",
                                                        list("barchart", "Table view"))
                                         ),
                                         mainPanel(
                                           textOutput("title_txt"),
                                           conditionalPanel('input.pType=="barchart"', plotOutput("barchart_SSA")),
                                           conditionalPanel('input.pType=="Table view"', tableOutput("tableview_SSA"))
                                         )
                                       )
                                     ),
                            
                            
                            # ----------------------------------
                            # tab panel 2-B - Temportal analysis
                            tabPanel("Temporal analysis", 
                                     sidebarLayout(
                                         sidebarPanel(
                                           selectInput("sample_sel", "Select your samples:",
                                                       c("sample 1" = "sample1",
                                                         "sample 2" = "sample2",
                                                         "sample 3" = "sample3",
                                                         "sample 4" = "sample4"), 
                                                       selected = "sample1",
                                                       multiple= TRUE),
                                           checkboxInput("alphaShow", "show alpha index", FALSE)
                                          ),
                                          mainPanel(
                                            textOutput("txt"),
                                             plotOutput("barchart_TSA"),
                                             tableOutput("tableview_TSA"),
                                             verbatimTextOutput("value")
                                          )
                                      )
                                  )
                            
                            
            )
          )
)
                 

# ----------------------------------
# ----------------------------------
# ----------------------------------
# SERVER SIDE
# ----------------------------------
# ----------------------------------

server <- function(input, output) {
  
  # ----------------------------------
  # tab panel 2-A - Sample analysis
    output$title_txt <- renderText({
      paste("You chose", input$sample)
    })
    output$barchart_SSA <- renderPlot({  
      validate(need(input$pType=="barchart", message=FALSE))
      display_dataset <- dataset1 %>% select(Phylum, sample_demo, input$sample) %>% 
        gather("sample", "Abundance", -Phylum)
      display_dataset %>% ggplot(aes(fill=Phylum, y=Abundance, x=sample)) + 
        geom_bar(position="fill", stat="identity")+
        scale_fill_viridis(discrete = T) +
        ggtitle("my title") +
        theme_ipsum() +
        xlab("xlab")+
        ylab("ylab")
    })
    table_view <- reactive({
      dataset1 %>% select(Phylum, sample_demo, input$sample)
    })
    output$tableview_SSA <- renderTable({
      validate(need(input$pType=="Table view", message=FALSE))
      table_view()
    })
  
    # ----------------------------------
    # tab panel 2-B - Temporal analysis
    output$txt <- renderText({
      samples_names <- paste(input$sample_sel, collapse = ", ")
      paste("You chose", samples_names)
    })
    output$barchart_TSA <- renderPlot({  
      display_dataset <- dataset1 %>% select(Phylum, input$sample_sel) %>% 
        gather("sample", "Abundance", -Phylum)
      display_dataset %>% ggplot(aes(fill=Phylum, y=Abundance, x=sample)) + 
        geom_bar(position="fill", stat="identity")+
        scale_fill_viridis(discrete = T) +
        ggtitle("my title") +
        theme_ipsum() +
        xlab("xlab")+
        ylab("ylab")
    })
    alpha_view <- reactive({
      dataset1 %>% select(Phylum, input$sample_sel)
    })
    output$tableview_TSA <- renderTable({
      validate(need(input$alphaShow, message=FALSE))
      alpha_view()
    })
}

shinyApp(server = server, ui = ui)