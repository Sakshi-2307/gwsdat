
uiTimeSeries <- function(csite, img_frmt) {
    
  fluidRow(
    shinydashboard::box(width = 3, status = "warning", title = "Settings",
                        
                        selectInput("sample_loc_select_ts", label = "Select Monitoring Well", choices = csite$ui_attr$sample_loc_names,
                                    selected = csite$ui_attr$sample_loc_select_ts, width = "80%", multiple=TRUE),
                        
                        selectInput("solute_select_ts", label = "Substance", choices = csite$ui_attr$solute_names,
                                    selected = csite$ui_attr$solute_select_ts, width = '80%', multiple=TRUE),
                        
                        radioButtons("solute_conc", label = "Solute Conc. Unit",
                                     choices = csite$ui_attr$conc_unit_list, 
                                     selected = csite$ui_attr$conc_unit_selected),
                        
                        checkboxInput("check_threshold", label = "Display threshold", 
                                      value = csite$ui_attr$show_thresh_ts ),
                        
                        checkboxGroupInput("ts_true_options", label = "Time Series Plot Options", 
                                           choices = names(csite$ui_attr$ts_options),
                                           selected = names(which(csite$ui_attr$ts_options == TRUE)))
                        
    ),
    
    shinydashboard::box(width = 9, status = "primary",
                        plotOutput("time_series"),
                        
                        div(style = "display: inline-block;",
                            selectInput("export_format_ts", label = "Image format", 
                                        choices = img_frmt, #csite$ui_attr$img_formats, 
                                        selected = img_frmt[[1]] #csite$ui_attr$img_formats[[1]]
                            )
                        ),
                        
                        div(style = "display: inline-block; vertical-align:top; margin-top: 25px; margin-right: 10px", 
                            downloadButton("save_timeseries_plot", label = "Export")
                        )
                        
    ),
    
	# This draggable panel contains the time slider for the spatial heatmap plot.s
	#print(csite$ui_attr$timepoints),
    absolutePanel(id = "timecontrol_ts", class = "panel panel-default", 
                  fixed = TRUE, draggable = TRUE, top = "auto", 
                  left = "auto", right = 20, bottom = 20,
                  width = 350, height = 140,  
                  
                  div(style = "margin-left: 15px; margin-top: 5px",
                      h4(textOutput("timepoint_ts_idx_label")),
                      sliderInput("timepoint_ts_idx",
                                  label="",
                                  #label = paste0("Time: ", pasteAggLimit(csite$ui_attr$timepoints[csite$ui_attr$timepoint_sp_idx], csite$GWSDAT_Options$Aggby)),
                                  min = 1,
                                  max = length(csite$ui_attr$timepoints),
                                  #max=10,
                                  step = 1,
                                  value = csite$ui_attr$timepoint_sp_idx,
                                  #value = length(csite$ui_attr$timepoints),
                                  animate = animationOptions(loop = TRUE, interval = 1500)
                      ) # ,
                      
                      # This worked nice for passing a vector of dates to values.
                      # However, update does not work and grid is messed up with too many values.
                      #
                      #sliderValues(
                      #  inputId = "timepoint_sp", label = "Time Point", width = "95%",
                      #  values = csite$ui_attr$timepoints, 
                      #  from = csite$ui_attr$timepoint_sp,
                      #  grid = if (length(csite$ui_attr$timepoints) < 20) {TRUE} else {FALSE},
                      #  animate = animationOptions(interval = 1500, loop = TRUE)
                  )
				)
	)
}
