library(shiny)
library(viridis)
library(hrbrthemes)
library(treemapify)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)

# ----------------------------------
# load dataset
dataset1 <- readr::read_csv("format_samples.csv") %>% rename("Shannon index"="Shannon_alpha", "Sorensen index"="Sorensen_beta")
info_dataset<- readr::read_csv("sample_info.csv")
dataset1 <- left_join(dataset1, info_dataset, by="sample")


ui <- navbarPage(title = "UABiome",
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
                                                        list("barchart", "Table view", "PieChart view"))
                                         ),
                                         mainPanel(
                                           textOutput("title_txt"),
                                           conditionalPanel('input.pType=="PieChart view"', plotOutput("bubblechart_SSA")),
                                           conditionalPanel('input.pType=="barchart"', plotOutput("barchart_SSA")),
                                           conditionalPanel('input.pType=="Table view"', tableOutput("tableview_SSA")),
                                           tableOutput("tableview_diversity")
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
                                           checkboxInput("alphaShow", "show alpha index", FALSE),
                                           checkboxInput("TempShow", "show dates", FALSE)
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
      paste("You are vizualizing ", input$sample)
    })
    output$barchart_SSA <- renderPlot({  
      validate(need(input$pType=="barchart", message=FALSE))
      display_dataset <- dataset1 %>% select(Genus, sample, Genus_count) %>% 
        filter(sample %in% c(input$sample, "sample_average")) %>% group_by(sample) %>% 
        arrange(desc(Genus_count)) %>% top_n(10)
      
      display_dataset %>% ggplot(aes(fill=Genus, y=Genus_count, x=sample)) + 
        geom_bar(position="fill", stat="identity")+
        scale_y_continuous(labels = scales::percent) +
        scale_fill_viridis(discrete = T) +
        ggtitle("Comparison to average population") +
        ylab("% of the population")+
        theme_minimal()
    })
    table_view <- reactive({
      dataset1 %>% select(Genus, sample, Genus_count) %>% 
        filter(sample==input$sample)
    })
    output$tableview_SSA <- renderTable({
      validate(need(input$pType=="Table view", message=FALSE))
      table_view()
    })
    output$bubblechart_SSA <- renderPlot({ 
      validate(need(input$pType=="PieChart view", message=FALSE))
      display_dataset <- dataset1 %>% select(Genus, sample, Genus_count) %>% 
        filter(sample==input$sample)
      display_dataset %>% ggplot(aes(x="", y=Genus_count, fill=Genus)) +
        geom_bar(stat="identity", width=1, color="white") +
        scale_fill_viridis(discrete = T) +
        coord_polar("y", start=0) +
        theme_void()
    })
    table_div <- reactive({
      dataset1 %>% select(sample,  "Shannon index", "Sorensen index") %>% 
        filter(sample==input$sample) %>% unique()
    })
    output$tableview_diversity <- renderTable({
      table_div()
    })
  
    # ----------------------------------
    # tab panel 2-B - Temporal analysis
    output$txt <- renderText({
      samples_names <- paste(input$sample_sel, collapse = ", ")
      paste("You chose", samples_names)
    })
    output$barchart_TSA <- renderPlot({  
      display_dataset <- dataset1 %>% select(Genus, sample, Genus_count, date) %>% 
        filter(sample %in% c(input$sample_sel))
      if(!input$TempShow){
        display_dataset %>% ggplot(aes(fill=Genus, y=Genus_count, x=sample)) + 
          geom_bar(position="fill", stat="identity")+
          scale_y_continuous(labels = scales::percent) +
          scale_fill_viridis(discrete = T) +
          ggtitle("Comparing your sample") +
          ylab("% of the microbial population") +
          theme_minimal()
      }else{
        display_dataset %>% ggplot(aes(fill=Genus, y=Genus_count, x=date)) + 
          geom_bar(position="fill", stat="identity")+
          scale_y_continuous(labels = scales::percent) +
          scale_fill_viridis(discrete = T) +
          ggtitle("Comparing your sample") +
          ylab("% of the microbial population") +
          theme_minimal()
      }
    })
    alpha_view <- reactive({
      dataset1 %>% select(sample, "Shannon index") %>% filter(sample %in% input$sample_sel) %>% 
        unique() 
    })
    output$tableview_TSA <- renderTable({
      validate(need(input$alphaShow, message=FALSE))
      alpha_view()
    })
}

shinyApp(server = server, ui = ui)