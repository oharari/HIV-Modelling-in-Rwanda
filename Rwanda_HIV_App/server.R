options(shiny.maxRequestSize = 5000*1024^2) 





shinyServer(
  function(input, output, session){
    
    
    my_PrEP_year = reactive({
      as.numeric(input$PrEP_Year)
    })
    
    
    
    
    my_time_to_PrEP = reactive({
      my_PrEP_year() - 2017
    })
    
    
    
    
    output$final_year_to_forecast = renderUI({
      first = max(2025, my_PrEP_year() + 6)
      selectInput("Final_Year",
                  "Final year to forecast:",
                  first:(first + 20),
                  width = "80%", selected = my_PrEP_year() + 8)
    })
    
    
    
    
    
    my_final_year = reactive({
      as.numeric(input$Final_Year)
    })
    
    
    
    
    
    my_N_years = reactive({
      my_final_year() - 2017
    })
    
    
    
    
    output$FSW_p_PrEP_start = renderUI({
      str = paste(my_PrEP_year(), 
                  "coverage of PrEP among FSW:")
      
      selectInput("p_PrEP_SW_start",
                  str,
                  paste0(5*(0:20), "%"),
                  width = "80%", selected = "40%")
    })
    
    
    
    
    
    my_p_PrEP_SW_start = reactive({
      as.numeric(substr(input$p_PrEP_SW_start,
                        1, nchar(input$p_PrEP_SW_start)-1))/100
    })
    
    
    
    
    
    
    output$FSW_p_PrEP_end = renderUI({
      str = paste(input$Final_Year, 
                  "coverage of PrEP among FSW:")
      
      temp = as.numeric(substr(input$p_PrEP_SW_start,
                               1, nchar(input$p_PrEP_SW_start)-1)
      )/5
      selectInput("p_PrEP_SW_end",
                  str,
                  paste0(5*(temp:20), "%"),
                  width = "80%", selected = "80%")
    })
    
    
    
    
    
    my_p_PrEP_SW_end = reactive({
      as.numeric(substr(input$p_PrEP_SW_end,
                        1, nchar(input$p_PrEP_SW_end)-1))/100
    })
    
    
    
    
    
    output$MSM_p_PrEP_start = renderUI({
      str = paste(my_PrEP_year(), 
                  "coverage of PrEP among MSM:")
      
      selectInput("p_PrEP_MSM_start",
                  str,
                  paste0(5*(0:20), "%"),
                  width = "80%", selected = "20%")
    })
    
    
    
    
    
    
    my_p_PrEP_MSM_start = reactive({
      as.numeric(substr(input$p_PrEP_MSM_start,
                        1, nchar(input$p_PrEP_MSM_start)-1))/100
    })
    
    
    
    
    
    output$MSM_p_PrEP_end = renderUI({
      str = paste(input$Final_Year, "coverage of PrEP among MSM:")
      
      temp = as.numeric(substr(input$p_PrEP_SW_start,
                               1, nchar(input$p_PrEP_MSM_start)-1)
      )/5
      selectInput("p_PrEP_MSM_end",
                  str,
                  paste0(5*(temp:20), "%"),
                  width = "80%", selected = "50%")
    })
    
    
    
    
    
    my_p_PrEP_MSM_end = reactive({
      as.numeric(substr(input$p_PrEP_MSM_end,
                        1, nchar(input$p_PrEP_MSM_end)-1))/100
    })
    
    
    
    
    
    output$run_button = renderUI({
      actionButton("go", "Run", icon("paper-plane"), 
                   style = "color: #fff; 
                   background-color: #337ab7; 
                   border-color: #2e6da4")
    })
    
    
    
    
    
    
    my_N_sim = reactive({
      input$N_sims
    })
    
    
    
    
    
    
    my_sim = eventReactive(input$go, {
      simulation(my_N_sim(), my_N_years(), 
                 p_condom_GP, p_condom_low, 
                 p_condom_high, condom_effic, 
                 p_vaginal, p_anal, N_partners_GP,
                 acts_per_partner_GP, N_partners_MP,
                 acts_per_partner_MP,
                 p_infected_vaginal_hiv_SW, 
                 p_infected_anal_hiv_SW,
                 p_infected_vaginal_ART_SW, 
                 p_infected_anal_ART_SW,
                 p_infected_vaginal_hiv_MP, 
                 p_infected_anal_hiv_MP,
                 p_infected_vaginal_ART_MP, 
                 p_infected_anal_ART_MP,
                 MP_counts, SW_counts, GP_counts,
                 p_test_SW, p_test_GP,
                 p_death_pre_ART,
                 p_start_ART_low,
                 p_start_ART_high,
                 p_death_ART_suppressed,
                 p_ART_adhere_low,
                 p_ART_adhere_high,
                 p_ART_fail,
                 p_death_ART_nonsuppressed,
                 p_switch_ART,
                 p_alternative_ART_fail,
                 my_p_PrEP_SW_start(), 
                 my_p_PrEP_SW_end(),
                 PrEP_effic, my_time_to_PrEP(),
                 SW_retire_rate, SW_recruit_rate,
                 MP_retire_rate, MP_recruit_rate,
                 growth_rate,
                 MSM_counts, N_partners_MSM,
                 p_infected_anal_hiv_SW,
                 p_infected_anal_ART_MSM,
                 acts_per_partner_MSM,
                 my_p_PrEP_MSM_start(), 
                 my_p_PrEP_MSM_end(),
                 PrEP_effic_MSM, p_condom_MSM_low, 
                 p_condom_MSM_high, p_test_MSM,
                 MSM_retire_rate, MSM_recruit_rate)
    })
    
    
    
    
    
    my_sim_results = reactive({
      if(is.null(my_sim())) return()
      
      simulation_analysis(sim=my_sim(), 
                          alpha=.05, 
                          start_year=2017,
                          N_simulations=my_N_sim())
    })
    
    
    
    
    
    output$stat_select = renderUI({
      if(is.null(my_sim_results())) stop()
      
      selectInput("stats_type",
                  "Statistics to display:",
                  c("Proportions",
                    "Incidences",
                    "HIV counts",
                    "ART counts",
                    "Prevalence",
                    "ART proportions",
                    "Counts",
                    "Overall incidence",
                    "Overall prevalence",
                    "Total HIV case count",
                    "Total ART count")
      )
    })
    
    
    
    
    my_stat_type = reactive({
      switch(input$stats_type,
             "Proportions" = {"props"},
             "Incidences" = {"incidences"},
             "HIV counts"= {"HIV_counts"},
             "ART counts"= {"ART_counts"},
             "Prevalence"= {"HIV_props"},
             "ART proportions"= {"ART_props"},
             "Overall incidence" = {"overall_incidence"},
             "Counts" = {"counts"},
             "Total HIV case count" = {"HIV_counts_all"},
             "Total ART count" = {"ART_counts_all"},
             "Overall prevalence" = {"HIV_props_all"},
             "Overall ART proportion" = {"ART_props_all"}
      )
    })
    
    
    
    
    
    is_overall_stat = reactive({
      my_stat_type() %in% c("overall_incidence",
                            "HIV_counts_all",
                            "HIV_props_all",
                            "ART_counts_all",
                            "ART_props_all")
    })
    
    
    
    
    
    output$state_select = renderUI({
      if(is.null(my_sim_results()) || 
         !(my_stat_type() %in% c("props", "counts")))
        return(NULL)
      
      radioButtons("state",
                   "State:",
                   c("HIV negative",
                     "Undiagnosed HIV",
                     "Diagnosed HIV Pre-ART",
                     "Viral load suppressed HIV on ART",
                     "ART Failed/Discontinued")
      )
    })
    
    
    
    
    
    my_state = reactive({
      switch(input$state,
             "HIV negative" = {1},
             "Undiagnosed HIV" = {2},
             "Diagnosed HIV Pre-ART" = {3},
             "Viral load suppressed HIV on ART" = {4},
             "ART Failed/Discontinued" = {5}
      )
    })
    
    
    
    
    
    output$sub_pops_select = renderUI({
      if(is.null(my_sim_results()) || is_overall_stat()) 
        return(NULL)
      
      checkboxGroupInput("sub_pops", 
                         "Populations to display:", 
                         c("GP","FSW","SC","MSM"),
                         selected = c("GP","FSW","SC","MSM"))
    })
    
    
    
    
    
    
    my_sub_pops = reactive({
      #temp = gsub("SC", "MP", input$sub_pops)
      temp = input$sub_pops
      #gsub("FSW", "SW", temp)
      temp
    })
    
    
    
    
    
    output$lump_select = renderUI({
      if(is.null(my_sim_results()) || is_overall_stat()) 
        return(NULL)
      
      checkboxGroupInput("lump", 
                         "Populations to group together:", 
                         input$sub_pops)
    })
    
    
    
    
    
    my_lump = reactive({
      #temp = gsub("SC", "MP", input$lump)
      temp = input$lump
      #gsub("FSW", "SW", temp)
      temp
    })
    
    
    
    
    
    
    my_plot = reactive({
      if(!is_overall_stat()){
        lumped_groups_plot(my_sub_pops(), 
                           my_lump(), 
                           my_sim_results(), 
                           my_stat_type(), 
                           state=my_state())
      } else{
        plot_outcomes_one(my_sim_results(), 
                          my_stat_type(), 
                          state=3, 
                          sub_pop = "GP")
      }
    })
    
    
    
    
    
    output$Main_Plot = renderPlot({
      my_plot()
    })
    
    
    
    
    
    
    output$down_Plot =downloadHandler(
      filename = function(){
        state_str = switch(my_state(),
                           "1" = {"HIV negative"},
                           "2" = {"Undiagnosed HIV"},
                           "3" = {"Diagnosed HIV Pre-ART"},
                           "4" = {"Viral load suppressed HIV on ART"},
                           "5" = {"ART Failed/Discontinued"}
        )
        
        
        stats_type_str = switch(my_stat_type(),
                                "props" = {"Proportion of "},
                                "incidences" = {"Incidence"},
                                "HIV_counts" = {"Number of people living with HIV"},
                                "HIV_props" = {"Proportion of people living with HIV"},
                                "ART_counts" = {"Number of people on antiretroviral therapy"},
                                "ART_props" = {"Proportion of people with HIV who are on ART"},
                                "counts" = {""},
                                "HIV_props_all" = {"Proportion of people living with HIV"},
                                "ART_counts_all" = {"Number of people on antiretroviral therapy"},
                                "ART_props_all" = {"Proportion of people with HIV who are on ART"},
                                "overall_incidence" = {"Incidence"},
                                "HIV_counts_all" = {"Number of people living with HIV"}
        )
        
        stats_type_str2 = ifelse(my_stat_type() == "counts", 
                                 " by year (no. of cases)", 
                                 " by year")
        
        if(my_stat_type() %in% c("incidences", "HIV_counts", "HIV_counts_all",
                                 "ART_counts", "HIV_props", "ART_props",
                                 "HIV_props_all", "ART_counts_all",
                                 "ART_props_all", "overall_incidence")){
          str = paste0(stats_type_str, stats_type_str2)
        } else{
          str = paste0(stats_type_str, state_str, stats_type_str2)
        }
        paste0(str, ".pdf")
      },
      content = function(file){
        ggsave(file, my_plot(), width=12, height=9)
      }
    )
    
    
    
    
    
    output$Plot_down = renderUI({
      if(is.null(my_plot())){
        return(NULL)
      } else{
        downloadButton("down_Plot", "Save to PDF", 
                       class = "butt1")
      }
    })
    
    
    
    
    
    
    
    
    
    
    
    output$table_stats_select =  renderUI({
      if(is.null(my_sim_results())) stop()
      
      selectInput("table_stats_type",
                  "Statistics to display:",
                  c("Proportions",
                    "Incidences",
                    "HIV counts",
                    "ART counts",
                    "Prevalence",
                    "ART proportions",
                    "Counts",
                    "Overall incidence",
                    "Overall prevalence",
                    "Total HIV case count",
                    "Total ART count")
      )
    })
    
    
    
    
    
    my_table_stat_type = reactive({
      switch(input$table_stats_type,
             "Proportions" = {"props"},
             "Incidences" = {"incidences"},
             "HIV counts"= {"HIV_counts"},
             "ART counts"= {"ART_counts"},
             "Prevalence"= {"HIV_props"},
             "ART proportions"= {"ART_props"},
             "Overall incidence" = {"overall_incidence"},
             "Counts" = {"counts"},
             "Total HIV case count" = {"HIV_counts_all"},
             "Total ART count" = {"ART_counts_all"},
             "Overall prevalence" = {"HIV_props_all"},
             "Overall ART proportion" = {"ART_props_all"}
      )
    })
    
    
    
    
    
    
    is_table_overall_stat = reactive({
      my_table_stat_type() %in% c("overall_incidence",
                                  "HIV_counts_all",
                                  "HIV_props_all",
                                  "ART_counts_all",
                                  "ART_props_all")
    })
    
    
    
    
    
    output$table_state_select = renderUI({
      if(is.null(my_sim_results()) || 
         !(my_table_stat_type() %in% c("props", "counts")))
        return(NULL)
      
      radioButtons("table_state",
                   "State:",
                   c("HIV negative",
                     "Undiagnosed HIV",
                     "Diagnosed HIV Pre-ART",
                     "Viral load suppressed HIV on ART",
                     "ART Failed/Discontinued")
      )
    })
    
    
    
    
    
    my_table_state = reactive({
      switch(input$table_state,
             "HIV negative" = {1},
             "Undiagnosed HIV" = {2},
             "Diagnosed HIV Pre-ART" = {3},
             "Viral load suppressed HIV on ART" = {4},
             "ART Failed/Discontinued" = {5}
      )
    })
    
    
    
    
    
    output$table_sub_pops_select = renderUI({
      if(is.null(my_sim_results()) || is_table_overall_stat()) 
        return(NULL)
      
      checkboxGroupInput("table_sub_pops", 
                         "Populations to display:", 
                         c("GP","FSW","SC","MSM"),
                         selected = c("GP","FSW","SC","MSM"))
    })
    
    
    
    
    
    
    my_table_sub_pops = reactive({
      #temp = gsub("SC", "MP", input$table_sub_pops)
      temp = input$table_sub_pops
      #gsub("FSW", "SW", temp)
      temp
    })
    
    
    
    
    
    output$table_lump_select = renderUI({
      if(is.null(my_sim_results()) || is_table_overall_stat()) 
        return(NULL)
      
      checkboxGroupInput("table_lump", 
                         "Population to group together:", 
                         input$table_sub_pops)
    })
    
    
    
    
    
    my_table_lump = reactive({
      #temp = gsub("SC", "MP", input$table_lump)
      temp = input$table_lump
      #gsub("FSW", "SW", temp)
      temp
    })
    
    
    
    
    my_table_to_display = reactive({
      if(!is_table_overall_stat()){
        tab = group_lumping(my_table_sub_pops(), 
                            my_table_lump(), 
                            my_sim_results(), 
                            my_table_stat_type(), 
                            my_table_state())
      } else{
        tab = group_lumping("GP", 
                            "GP", 
                            my_sim_results(), 
                            my_table_stat_type(), 
                            1)
      }
      names(tab)[1] = "Date"
      
      if(!(my_table_stat_type() %in% c("counts", "HIV_counts", "ART_counts",
                                       "HIV_counts_all", "ART_counts_all"))){
        tab[,2] = sprintf("%.2f", tab[,2])
      } else{
        tab[,2] = round(tab[,2])
      }
      if(!is_table_overall_stat()){
        #tab[,3] = gsub("SW", "FSW", tab[,3])
        #tab[,3] = gsub("MP", "SC", tab[,3])
        names(tab)[3] = c("Population")
      }
      
      if(my_table_stat_type() %in% c("counts", "HIV_counts", "ART_counts",
                                     "HIV_counts_all", "ART_counts_all")){
        names(tab)[2] = "# cases" 
      } else if(my_table_stat_type() %in% c("incidences", "overall_incidence")){
        names(tab)[2] = "# cases/1000 PY"
      } else{
        names(tab)[2] = "Proportion (%)"
      }
      
      tab
    })
    
    
    
    
    
    
    output$mytable = DT::renderDataTable({
      data.frame(my_table_to_display(), check.names=FALSE)
    }, server = FALSE, escape = FALSE,  
    options = list( 
      pageLength = my_N_years()*4 + 1,
      preDrawCallback = JS('function() { 
                           Shiny.unbindAll(this.api().table().node()); }'), 
      drawCallback = JS('function() { 
                        Shiny.bindAll(this.api().table().node()); } '),
      columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ), 
    rownames = FALSE)
    
    
    
    
    
    output$down_Table = downloadHandler(
      filename = function(){
        state_str = switch(my_table_state(),
                           "1" = {"HIV negative"},
                           "2" = {"Undiagnosed HIV"},
                           "3" = {"Diagnosed HIV Pre-ART"},
                           "4" = {"Viral load suppressed HIV on ART"},
                           "5" = {"ART Failed/Discontinued"}
        )
        
        
        stats_type_str = switch(my_table_stat_type(),
                                "props" = {"Proportion of "},
                                "incidences" = {"Incidence"},
                                "HIV_counts" = {"Number of people living with HIV"},
                                "HIV_props" = {"Proportion of people living with HIV"},
                                "ART_counts" = {"Number of people on antiretroviral therapy"},
                                "ART_props" = {"Proportion of people with HIV who are on ART"},
                                "counts" = {""},
                                "HIV_props_all" = {"Proportion of people living with HIV"},
                                "ART_counts_all" = {"Number of people on antiretroviral therapy"},
                                "ART_props_all" = {"Proportion of people with HIV who are on ART"},
                                "overall_incidence" = {"Incidence"},
                                "HIV_counts_all" = {"Number of people living with HIV"}
        )
        
        stats_type_str2 = ifelse(my_table_stat_type() == "counts", 
                                 " by year (no. of cases)", 
                                 " by year")
        
        if(my_table_stat_type() %in% c("incidences", "HIV_counts", "HIV_counts_all",
                                       "ART_counts", "HIV_props", "ART_props",
                                       "HIV_props_all", "ART_counts_all",
                                       "ART_props_all", "overall_incidence")){
          str = paste0(stats_type_str, stats_type_str2)
        } else{
          str = paste0(stats_type_str, state_str, stats_type_str2)
        }
        paste0(str, ".xlsx")
      },
      
      content = function(file){
        tempFile = tempfile(fileext = ".xlsx")
        openxlsx::write.xlsx(my_table_to_display(), tempFile, rowNames=F)
        file.rename(tempFile, file)
      }
    )
    
    
    
    
    
    output$Table_down = renderUI({
      if(is.null(my_table_to_display())){
        return(NULL)
      } else{
        downloadButton("down_Table", "Save to Excel", 
                       class = "butt1")
      }
    })
    
    
    
    
    
    output$down_all_tables = downloadHandler(
      filename = function(){"HIV_Output.xlsx"},
      
      content = function(file){
        helper = function(x){
          x = gsub("SW", "FSW", x)
          x
        }
        
        wb = openxlsx::createWorkbook("HIV_Output")
        for(name in names(my_sim_results())){
          tab = my_sim_results()[[name]]
          openxlsx::addWorksheet(wb, sheetName = name)
          writeData(wb, sheet = name, tab)
        }
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
    
    
    
    
    output$all_tables_down = renderUI({
      if(is.null(my_table_to_display())){
        return(NULL)
      } else{
        downloadButton("down_all_tables", 
                       "Save all tables", 
                       class = "butt1")
      }
    })
    
    
    
    
    
    
    
    
    output$no_prep_text = renderText({
      if(is.null(my_sim_results()))
        return(NULL)
      
      HTML("<font size='2'><b>Simulate without PrEP:</b></font>")
    })
    
    
    
    
    
    output$no_PrEP_run_button = renderUI({
      if(is.null(my_sim_results()))
        return(NULL)
      
      actionButton("go_no_PrEP", "Run", icon("paper-plane"), 
                   style = "color: #fff; 
                   background-color: #337ab7; 
                   border-color: #2e6da4")
    })
    
    
    
    
    
    
    
    my_sim_no_PrEP = eventReactive(input$go_no_PrEP, {
      simulation(my_N_sim(), my_N_years(), 
                 p_condom_GP, p_condom_low, 
                 p_condom_high, condom_effic, 
                 p_vaginal, p_anal, N_partners_GP,
                 acts_per_partner_GP, N_partners_MP,
                 acts_per_partner_MP,
                 p_infected_vaginal_hiv_SW, 
                 p_infected_anal_hiv_SW,
                 p_infected_vaginal_ART_SW, 
                 p_infected_anal_ART_SW,
                 p_infected_vaginal_hiv_MP, 
                 p_infected_anal_hiv_MP,
                 p_infected_vaginal_ART_MP, 
                 p_infected_anal_ART_MP,
                 MP_counts, SW_counts, GP_counts,
                 p_test_SW, p_test_GP,
                 p_death_pre_ART,
                 p_start_ART_low,
                 p_start_ART_high,
                 p_death_ART_suppressed,
                 p_ART_adhere_low,
                 p_ART_adhere_high,
                 p_ART_fail,
                 p_death_ART_nonsuppressed,
                 p_switch_ART,
                 p_alternative_ART_fail,
                 0, 
                 0,
                 PrEP_effic, my_time_to_PrEP(),
                 SW_retire_rate, SW_recruit_rate,
                 MP_retire_rate, MP_recruit_rate,
                 growth_rate,
                 MSM_counts, N_partners_MSM,
                 p_infected_anal_hiv_SW,
                 p_infected_anal_ART_MSM,
                 acts_per_partner_MSM,
                 0, 
                 0,
                 PrEP_effic_MSM, p_condom_MSM_low, 
                 p_condom_MSM_high, p_test_MSM,
                 MSM_retire_rate, MSM_recruit_rate)
    })
    
    
    
    
    
    
    
    output$plot_or_table_select = renderUI({
      if(is.null(my_sim_results_no_PrEP())) return()
      
      radioButtons("plot_or_table",
                   "Display:",
                   c("Plot", "Table"))
    })
    
    
    
    
    
    
    my_sim_results_no_PrEP = reactive({
      if(is.null(my_sim_no_PrEP())) return()
      
      simulation_analysis(sim=my_sim_no_PrEP(), 
                          alpha=.05, 
                          start_year=2017,
                          N_simulations=my_N_sim())
    })
    
    
    
    
    
    
    output$stat_select_no_PrEP = renderUI({
      if(is.null(my_sim_results_no_PrEP())) stop()
      
      selectInput("stats_type_no_PrEP",
                  "Statistics to display:",
                  c("Proportions",
                    "Incidences",
                    "HIV counts",
                    "ART counts",
                    "Prevalence",
                    "ART proportions",
                    "Counts",
                    "Overall incidence",
                    "Overall prevalence",
                    "Total HIV case count",
                    "Total ART count")
      )
    })
    
    
    
    
    my_stat_type_no_PrEP = reactive({
      switch(input$stats_type_no_PrEP,
             "Proportions" = {"props"},
             "Incidences" = {"incidences"},
             "HIV counts"= {"HIV_counts"},
             "ART counts"= {"ART_counts"},
             "Prevalence"= {"HIV_props"},
             "ART proportions"= {"ART_props"},
             "Overall incidence" = {"overall_incidence"},
             "Counts" = {"counts"},
             "Total HIV case count" = {"HIV_counts_all"},
             "Total ART count" = {"ART_counts_all"},
             "Overall prevalence" = {"HIV_props_all"},
             "Overall ART proportion" = {"ART_props_all"}
      )
    })
    
    
    
    
    
    is_overall_stat_no_PrEP = reactive({
      my_stat_type_no_PrEP() %in% c("overall_incidence",
                                    "HIV_counts_all",
                                    "HIV_props_all",
                                    "ART_counts_all",
                                    "ART_props_all")
    })
    
    
    
    
    
    output$state_select_no_PrEP = renderUI({
      if(is.null(my_sim_results_no_PrEP()) || 
         !(my_stat_type_no_PrEP() %in% c("props", "counts")))
        return(NULL)
      
      radioButtons("state_no_PrEP",
                   "State:",
                   c("HIV negative",
                     "Undiagnosed HIV",
                     "Diagnosed HIV Pre-ART",
                     "Viral load suppressed HIV on ART",
                     "ART Failed/Discontinued")
      )
    })
    
    
    
    
    
    my_state_no_PrEP = reactive({
      switch(input$state_no_PrEP,
             "HIV negative" = {1},
             "Undiagnosed HIV" = {2},
             "Diagnosed HIV Pre-ART" = {3},
             "Viral load suppressed HIV on ART" = {4},
             "ART Failed/Discontinued" = {5}
      )
    })
    
    
    
    
    
    output$sub_pops_select_no_PrEP = renderUI({
      if(is.null(my_sim_results_no_PrEP()) || is_overall_stat_no_PrEP()) 
        return(NULL)
      
      checkboxGroupInput("sub_pops_no_PrEP", 
                         "Populations to display:", 
                         c("GP","FSW","SC","MSM"),
                         selected = c("GP","FSW","SC","MSM"))
    })
    
    
    
    
    
    
    my_sub_pops_no_PrEP = reactive({
      #temp = gsub("SC", "MP", input$sub_pops_no_PrEP)
      temp = input$sub_pops_no_PrEP
      #gsub("FSW", "SW", temp)
      temp
    })
    
    
    
    
    
    output$lump_select_no_PrEP = renderUI({
      if(is.null(my_sim_results_no_PrEP()) || is_overall_stat_no_PrEP()) 
        return(NULL)
      
      checkboxGroupInput("lump_no_PrEP", 
                         "Populations to group together:", 
                         input$sub_pops_no_PrEP)
    })
    
    
    
    
    
    my_lump_no_PrEP = reactive({
      #temp = gsub("SC", "MP", input$lump_no_PrEP)
      temp = input$lump_no_PrEP
      #gsub("FSW", "SW", temp)
      temp
    })
    
    
    
    
    
    
    my_no_PrEP_plot = reactive({
      if(!is_overall_stat_no_PrEP()){
        With_Without_PrEP_lumped_plot(my_sub_pops_no_PrEP(), 
                                      my_lump_no_PrEP(), 
                                      my_sim_results(), 
                                      my_sim_results_no_PrEP(), 
                                      my_stat_type_no_PrEP(), 
                                      my_state_no_PrEP())
      } else{
        With_Without_PrEP_plot("GP", 
                               my_sim_results(), 
                               my_sim_results_no_PrEP(),  
                               my_stat_type_no_PrEP(), 
                               1)
      }
    })
    
    
    
    
    
    
    output$Main_Plot_no_PrEP = renderPlot({
      my_no_PrEP_plot()
    })
    
    
    
    
    
    output$down_Plot_no_PrEP =downloadHandler(
      filename = function(){
        state_str = switch(my_state_no_PrEP(),
                           "1" = {"HIV negative"},
                           "2" = {"Undiagnosed HIV"},
                           "3" = {"Diagnosed HIV Pre-ART"},
                           "4" = {"Viral load suppressed HIV on ART"},
                           "5" = {"ART Failed/Discontinued"}
        )
        
        
        stats_type_str = switch(my_stat_type_no_PrEP(),
                                "props" = {"Proportion of "},
                                "incidences" = {"Incidence"},
                                "HIV_counts" = {"Number of people living with HIV"},
                                "HIV_props" = {"Proportion of people living with HIV"},
                                "ART_counts" = {"Number of people on antiretroviral therapy"},
                                "ART_props" = {"Proportion of people with HIV who are on ART"},
                                "counts" = {""},
                                "HIV_props_all" = {"Proportion of people living with HIV"},
                                "ART_counts_all" = {"Number of people on antiretroviral therapy"},
                                "ART_props_all" = {"Proportion of people with HIV who are on ART"},
                                "overall_incidence" = {"Incidence"},
                                "HIV_counts_all" = {"Number of people living with HIV"}
        )
        
        stats_type_str2 = ifelse(my_stat_type_no_PrEP() == "counts", 
                                 " by year (no. of cases)", 
                                 " by year")
        
        if(my_stat_type_no_PrEP() %in% c("incidences", "HIV_counts", "HIV_counts_all",
                                         "ART_counts", "HIV_props", "ART_props",
                                         "HIV_props_all", "ART_counts_all",
                                         "ART_props_all", "overall_incidence")){
          str = paste0(stats_type_str, stats_type_str2)
        } else{
          str = paste0(stats_type_str, state_str, stats_type_str2)
        }
        paste0(str, " PrEP vs No PrEP.pdf")
      },
      content = function(file){
        ggsave(file, my_no_PrEP_plot(), width=12, height=9)
      }
    )
    
    
    
    
    
    output$Plot_down_no_PrEP = renderUI({
      if(is.null(my_no_PrEP_plot())){
        return(NULL)
      } else{
        downloadButton("down_Plot_no_PrEP", "Save to PDF", 
                       class = "butt1")
      }
    })
    
    
    
    
    
    
    
    my_table_to_display_PrEP_no_PrEP = reactive({
      if(!is_overall_stat_no_PrEP()){
        tab_no_PrEP = group_lumping(my_sub_pops_no_PrEP(), 
                            my_lump_no_PrEP(), 
                            my_sim_results_no_PrEP(), 
                            my_stat_type_no_PrEP(), 
                            my_state_no_PrEP())
        
        tab_PrEP = group_lumping(my_sub_pops_no_PrEP(), 
                                 my_lump_no_PrEP(), 
                                 my_sim_results(), 
                                 my_stat_type_no_PrEP(), 
                                 my_state_no_PrEP())
      } else{
        tab_no_PrEP = group_lumping("GP", 
                            "GP", 
                            my_sim_results_no_PrEP(), 
                            my_stat_type_no_PrEP(), 
                            1)
        
        tab_PrEP = group_lumping("GP", 
                                    "GP", 
                                    my_sim_results(), 
                                    my_stat_type_no_PrEP(), 
                                    1)
      }
      
      if(!(my_stat_type_no_PrEP() %in% c("counts", "HIV_counts", "ART_counts",
                                       "HIV_counts_all", "ART_counts_all"))){
        tab_PrEP[,2] = sprintf("%.2f", tab_PrEP[,2])
        tab_no_PrEP[,2] = sprintf("%.2f", tab_no_PrEP[,2])
      } else{
        tab_PrEP[,2] = round(tab_PrEP[,2])
        tab_no_PrEP[,2] = round(tab_no_PrEP[,2])
      }
      if(!is_overall_stat_no_PrEP()){
        #tab_no_PrEP[,3] = gsub("SW", "FSW", tab_no_PrEP[,3])
        #tab_no_PrEP[,3] = gsub("MP", "SC", tab_no_PrEP[,3])
        tab = #as.data.frame(
          cbind.data.frame(tab_no_PrEP[,1], tab_no_PrEP[,2], 
                                  tab_PrEP[,2], tab_no_PrEP[,3])
          #)
        names(tab)[4] = c("Population")
      } else{
        tab = #as.data.frame(
          cbind.data.frame(tab_no_PrEP, tab_PrEP[,2])
          #)
      }
      
      if(my_stat_type_no_PrEP() %in% c("counts", "HIV_counts", "ART_counts",
                                     "HIV_counts_all", "ART_counts_all")){
        names(tab)[2] = "# cases w/o PrEP" 
        names(tab)[3] = "# cases w PrEP" 
      } else if(my_stat_type_no_PrEP() %in% c("incidences", "overall_incidence")){
        names(tab)[2] = "# cases/1000 PY w/o PrEP"
        names(tab)[3] = "# cases/1000 PY w PrEP"
      } else{
        names(tab)[2] = "Proportion (%) w/o PrEP"
        names(tab)[3] = "Proportion (%) w PrEP"
      }
      
      names(tab)[1] = "Date"

      tab
    })
    
    
    
    
    
    output$mytable_PrEP_no_PrEP = DT::renderDataTable({
      data.frame(my_table_to_display_PrEP_no_PrEP(), check.names=FALSE)
    }, server = FALSE, escape = FALSE,  
    options = list( 
      pageLength = my_N_years()*4 + 1,
      preDrawCallback = JS('function() { 
                           Shiny.unbindAll(this.api().table().node()); }'), 
      drawCallback = JS('function() { 
                        Shiny.bindAll(this.api().table().node()); } '),
      columnDefs = list(list(className = 'dt-center', targets = "_all"))
    ), 
    rownames = FALSE)
    
    
    
    
    
    
    output$Plot_or_Table_PrEP_no_PrEP = renderUI({
      if(is.null(my_sim_results_no_PrEP())) return()
      
      if(input$plot_or_table == "Plot"){
        plotOutput("Main_Plot_no_PrEP")
      } else{
        DT::dataTableOutput("mytable_PrEP_no_PrEP")
      }
    })
    
    
    
    
    
    
    output$down_Table_PrEP_no_PrEP = downloadHandler(
      filename = function(){
        state_str = switch(my_state_no_PrEP(),
                           "1" = {"HIV negative"},
                           "2" = {"Undiagnosed HIV"},
                           "3" = {"Diagnosed HIV Pre-ART"},
                           "4" = {"Viral load suppressed HIV on ART"},
                           "5" = {"ART Failed/Discontinued"}
        )
        
        
        stats_type_str = switch(my_stat_type_no_PrEP(),
                                "props" = {"Proportion of "},
                                "incidences" = {"Incidence"},
                                "HIV_counts" = {"Number of people living with HIV"},
                                "HIV_props" = {"Proportion of people living with HIV"},
                                "ART_counts" = {"Number of people on antiretroviral therapy"},
                                "ART_props" = {"Proportion of people with HIV who are on ART"},
                                "counts" = {""},
                                "HIV_props_all" = {"Proportion of people living with HIV"},
                                "ART_counts_all" = {"Number of people on antiretroviral therapy"},
                                "ART_props_all" = {"Proportion of people with HIV who are on ART"},
                                "overall_incidence" = {"Incidence"},
                                "HIV_counts_all" = {"Number of people living with HIV"}
        )
        
        stats_type_str2 = ifelse(my_stat_type_no_PrEP() == "counts", 
                                 " by year (no. of cases)", 
                                 " by year")
        
        if(my_stat_type_no_PrEP() %in% c("incidences", "HIV_counts", "HIV_counts_all",
                                       "ART_counts", "HIV_props", "ART_props",
                                       "HIV_props_all", "ART_counts_all",
                                       "ART_props_all", "overall_incidence")){
          str = paste0(stats_type_str, stats_type_str2)
        } else{
          str = paste0(stats_type_str, state_str, stats_type_str2)
        }
        paste0(str, " PrEP vs No PrEP.xlsx")
      },
      
      content = function(file){
        tempFile = tempfile(fileext = ".xlsx")
        openxlsx::write.xlsx(my_table_to_display_PrEP_no_PrEP(), 
                             tempFile, rowNames=F)
        file.rename(tempFile, file)
      }
    )
    
    
    
    
    
    
    output$PrEP_no_PrEP_Table_down = renderUI({
      if(is.null(my_table_to_display_PrEP_no_PrEP())){
        return(NULL)
      } else{
        downloadButton("down_Table_PrEP_no_PrEP", 
                       "Save to Excel", 
                       class = "butt1")
      }
    })
    
    
    
    
    
    
    output$plot_or_table_down = renderUI({
      if(input$plot_or_table == "Plot"){
        uiOutput("Plot_down_no_PrEP")
      } else{
        uiOutput("PrEP_no_PrEP_Table_down")
      }
    })
    
    
    
    
    }
)

