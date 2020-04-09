library(shiny)
library(viridis)
library(hrbrthemes)
library(treemapify)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# ----------------------------------
# load dataset
dataset1 <- read_csv("format_samples.csv") %>% rename("Shannon index"="Shannon_alpha", "Sorensen index"="Sorensen_beta")
info_dataset<- read_csv("sample_info.csv", col_types=cols(date = col_date("%Y-%m-%d")))
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
                 # tab panel 2 - Analysis
                 tabPanel(title = "Analysis",
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
                                                        list("Barchart", "Table view", "PieChart view"))
                                         ),
                                         mainPanel(
                                           conditionalPanel('input.pType=="PieChart view"', plotOutput("bubblechart_SSA")),
                                           conditionalPanel('input.pType=="Barchart"', plotOutput("barchart_SSA")),
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
                                             plotOutput("barchart_TSA"),
                                             tableOutput("tableview_TSA"),
                                             verbatimTextOutput("value")
                                          )
                                      )
                                  )
                            
                            
            )
          ),
          
                # ----------------------------------
                # tab panel 3 - my Profile
                tabPanel(title = "My profile",
                         sidebarLayout(
                           sidebarPanel(
                             textInput("myalias", "New username", "Username 1"),
                             actionButton("update_name", "Update my alias")
                           ),
                           mainPanel(
                             textOutput("alias"),
                             tableOutput("tableview_samples") 
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
    output$barchart_SSA <- renderPlot({  
      validate(need(input$pType=="Barchart", message=FALSE))
      display_dataset <- dataset1 %>% select(Genus, sample, Genus_count) %>% 
        filter(sample %in% c(input$sample, "sample_average")) %>% group_by(sample) %>% 
        arrange(desc(Genus_count)) %>% top_n(10)
      title_txt <- paste("Comparison of ", input$sample, " to average population")
      
      display_dataset %>% ggplot(aes(fill=Genus, y=Genus_count, x=sample)) + 
        geom_bar(position="fill", stat="identity")+
        scale_y_continuous(labels = scales::percent) +
        scale_fill_viridis(discrete = T) +
        ggtitle(title_txt) +
        ylab("% of the population")+
        xlab("")+
        theme_minimal(base_size = 13)
    })
    table_view <- reactive({
      dataset1 %>%  filter(sample==input$sample) %>%
        select(Genus, Genus_count)
    })
    output$tableview_SSA <- renderTable({
      validate(need(input$pType=="Table view", message=FALSE))
      table_view()
    })
    output$bubblechart_SSA <- renderPlot({ 
      validate(need(input$pType=="PieChart view", message=FALSE))
      display_dataset <- dataset1 %>% select(Genus, sample, Genus_count) %>% 
        filter(sample==input$sample)
      
      title_txt <- paste("Population of ", input$sample)
      
      display_dataset %>% ggplot(aes(x="", y=Genus_count, fill=Genus)) +
        geom_bar(stat="identity", width=1, color="white") +
        scale_fill_viridis(discrete = T) +
        coord_polar("y", start=0) +
        ggtitle(title_txt) +
        theme_void(base_size = 13)
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
    output$barchart_TSA <- renderPlot({  
      display_dataset <- dataset1 %>% select(Genus, sample, Genus_count, date) %>% 
        filter(sample %in% c(input$sample_sel))
      
      samples_names <- paste(input$sample_sel, collapse = ", ")
      title_txt2 <- paste("Composition of ", samples_names)
      
      if(!input$TempShow){
        display_dataset %>% ggplot(aes(fill=Genus, y=Genus_count, x=sample)) + 
          geom_bar(position="fill", stat="identity")+
          scale_y_continuous(labels = scales::percent) +
          scale_fill_viridis(discrete = T) +
          ggtitle(title_txt2) +
          ylab("% of the microbial population") +
          xlab("")+
          theme_minimal(base_size = 13)
      }else{
        display_dataset %>% ggplot(aes(fill=Genus, y=Genus_count, x=date)) + 
          geom_bar(position="fill", stat="identity")+
          scale_y_continuous(labels = scales::percent) +
          scale_fill_viridis(discrete = T) +
          ggtitle(title_txt2) +
          ylab("% of the microbial population") +
          xlab("")+
          theme_minimal(base_size = 13)
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
    
    # ----------------------------------
    # tab panel 3 - User profile
    new_username <- reactiveValues(data = "default_1")
    observeEvent(input$update_name, {
      new_username$data <- input$myalias
    })
    output$alias <- renderText({
        paste("Your username :", new_username$data)
    })
    table_samples <- reactive({
      dataset1 %>%  mutate("submission_date"=as.character(dataset1$date)) %>% select(sample, submission_date) %>% 
        unique() %>% filter(!sample=="sample_average") %>%
        rename("your samples"=sample) 
    })
    output$tableview_samples <- renderTable({
      table_samples()
    })

}

shinyApp(server = server, ui = ui)