shinyUI(fluidPage(
  
  #useShinyjs(),
  
  tags$head(tags$style(
    "body { word-wrap: break-word; }")),
  
  title = "Rwanda HIV Forecasting",
  
  tags$style(HTML("
                  .tabbable > .nav > li > a                  
                  {background-color: black;  color:white}"
  )),
  
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
  
  
  
  
  fluidRow(
    column(width = 12, titlePanel(div(img(src="hiv.png", width = 1200)))
    #column(width = 12, titlePanel(div(img(src="hiv_ki.png", width = 1200)))
)#,
    
    # column(width = 9, 
    #        br(), 
    #        br(),
    #        titlePanel(title = div(
    #          "HIV Forecasting in Rwanda by Population",
    #          style = "font-family: 'Avenir';
    #                    font-weight: 700; line-height: 1.1; font-size: 120%;
    #                    color: #000033;"))
    # )
  ),
  
  
  tabsetPanel(
    
    tabPanel("Run Model",
             sidebarLayout(
               sidebarPanel(width=3,
                            style = "background-color: lightgrey",
                            
                            selectInput("PrEP_Year",
                                        "Year of PrEP deployment:",
                                        2019:2029,
                                        width = "80%"),
                            
                            uiOutput("final_year_to_forecast"),
                            
                            uiOutput("FSW_p_PrEP_start"),
                            
                            uiOutput("FSW_p_PrEP_end"),
                            
                            uiOutput("MSM_p_PrEP_start"),
                            
                            uiOutput("MSM_p_PrEP_end"),
                            
                            sliderInput("N_sims", 
                                        "Number of simulations:",  
                                        min = 50, max = 1000, 
                                        value = 500, step = 50),
                            
                            br(),
                            tags$head(tags$style(type="text/css", "
                                                 #loadmessage {
                                                 position: fixed;
                                                 top: 0px;
                                                 left: 0px;
                                                 width: 100%;
                                                 padding: 5px 0px 5px 0px;
                                                 text-align: center;
                                                 font-weight: bold;
                                                 font-size: 100%;
                                                 color: #000000;
                                                 background-color: #00FFFF;
                                                 z-index: 105;
                                                 }
                                                 ")),
                            
                            uiOutput("run_button"),
                            
                            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                             tags$div("In Progress...", id="loadmessage"))
                            
                            ),
               
               
               
               mainPanel(
                 fluidRow(
                   column(width=3,
                          uiOutput("stat_select")
                   ),
                   column(width=3,
                          uiOutput("sub_pops_select")
                   ),
                   column(width=3,
                          uiOutput("lump_select")
                   )
                 ),
                 
                 fluidRow(
                   column(width=8,
                          plotOutput("Main_Plot"),
                          uiOutput("Plot_down"),
                          tags$head(tags$style(".butt1{background-color: #337ab7;} 
                                               .butt1{color: white;} 
                                               .butt1{font-weight: bold;}
                                               .butt1{font-family: Courier New}"))
                          ),
                   column(width=4,
                          uiOutput("state_select"))
                          )
                   )
               
               
               )
                            ),
    
    
    
    tabPanel("View Tables",
             sidebarLayout(
               sidebarPanel(width=2,
                            style = "background-color: lightgrey",
                            
                            uiOutput("table_stats_select"),
                            uiOutput("table_sub_pops_select"),
                            uiOutput("table_lump_select"),
                            uiOutput("table_state_select"),
                            uiOutput("Table_down"),
                            br(),
                            uiOutput("all_tables_down")
               ),
               
               
               mainPanel(
                 fluidRow(
                   column(width=12,
                          DT::dataTableOutput("mytable")
                   )
                 )
                 
                 
               )
             )
             
             
    ),
    
    
    
    
    
    
    tabPanel("No PrEP",
             sidebarLayout(
               sidebarPanel(width=2,
                            style = "background-color: lightgrey",
                            
                            uiOutput("no_prep_text"),
                            uiOutput("no_PrEP_run_button"),
                            br(),
                            br(),
                            uiOutput("plot_or_table_select"),
                            br(),
                            uiOutput("plot_or_table_down"),
                            
                            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                             tags$div("In Progress...", id="loadmessage"))
               ),
               
               
               mainPanel(
                 fluidRow(
                   column(width=3,
                          uiOutput("stat_select_no_PrEP")
                   ),
                   column(width=3,
                          uiOutput("sub_pops_select_no_PrEP")
                   ),
                   column(width=3,
                          uiOutput("lump_select_no_PrEP")
                   )
                 ),
                 
                 fluidRow(
                   column(width=8,
                          br(),
                          uiOutput("Plot_or_Table_PrEP_no_PrEP")
                   ),
                   column(width=4,
                          uiOutput("state_select_no_PrEP"))
                 )
                 
                 
               )
             )
             
             
    )
    
    
    
    
    
    
    
               )
               )
  
    )

