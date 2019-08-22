#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("DREAM SCTC challenge"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            radioButtons("subchallenge",
                        "Subchallenge",
                        c("Sub1" = 1,
                        "Sub2" = 2,
                        "Sub3" = 3)
        )),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("heatmap")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$heatmap <- renderPlot({
        # generate bins based on input$bins from ui.R
        n_genes = case_when(input$subchallenge == 1 ~ 60,
                            input$subchallenge == 2 ~ 40,
                            input$subchallenge == 3 ~ 20)
        n_max = as.integer(n_genes / 10)

        list.files("post/final_submissions/", paste0(n_genes, "genes"), full.names = TRUE) %>%
            map(~read_csv(.x, col_names = F, n_max = 6, col_types = cols())  %>%
                    select(-1) %>% gather %>% mutate(key = str_replace_all(.x, c( "post/final_submissions" ="",
                     "genes.csv" = "",
                     n_genes = "",
                          "/" = "") ))) %>% bind_rows() %>% arrange(desc(key)) -> challenge.frequencies

        ggplot(challenge.frequencies) + geom_tile(aes(value, key), fill = "red", color = "black") +
            labs(x= "Gene", y = "Team") +
            theme(axis.text.x = element_text( color = "blue", size = 7, angle = 90, hjust = 0.95))

    })
}

# Run the application
shinyApp(ui = ui, server = server)
