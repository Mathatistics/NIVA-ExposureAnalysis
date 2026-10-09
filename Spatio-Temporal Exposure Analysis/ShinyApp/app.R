## -- Load packages ------------------------------
pkgs <- c(
  "tidytable", "ggplot2", "purrr",
  "shiny", "shinydashboard", "leaflet"
)
for (pkg in pkgs) require(pkg, character.only = TRUE)

## -- Local functions ------------------------------
get_tissue_conc <- function(cw, logp) {
  cf <- cw * 10 ^(0.76 * logp - 0.23)
  return(cf)
}

## -- Load data set at the start of application ------------------------------
exposure_data <- fread("Exposure_data_AEP.csv") %>% 
  tidytable::as_tidytable() %>%
  tidytable::mutate(FishConc = get_tissue_conc(MEASURED_VALUE, XLogP)) %>% 
  tidytable::mutate(Month = lubridate::month(SAMPLE_DATE, label = TRUE)) %>% 
  tidytable::mutate(Wday = lubridate::wday(SAMPLE_DATE, label = TRUE, abbr = FALSE)) %>% 
  tidytable::select(
    SITE_CODE, SITE_NAME,
    STRESSOR_ID, STRESSOR_NAME,
    WaterConc = MEASURED_VALUE, FishConc, 
    LATITUDE, LONGITUDE,
    Month, Wday
  )

## -- Define UI for application ------------------------------
ui <- shiny::bootstrapPage(
  title = "Exposure Data Analysis",
  theme = shinythemes::shinytheme(theme = "lumen"),
  shiny::pageWithSidebar(
    ## -- Header panel ------------------------------
    headerPanel = shiny::headerPanel(
      title = "Exposure Data Analysis"
    ),
    ## -- Sidebar panel ------------------------------
    sidebarPanel = shiny::sidebarPanel(
      width = 3,
      ## -- Input: site var selector ------------------------------
      selectInput(
        inputId = "site_vars", 
        label = "Select site variable",
        choices = exposure_data[, unique(SITE_NAME)],
        selected = exposure_data[, unique(SITE_NAME)],
        multiple = TRUE,
        selectize = TRUE
      ),
      ## -- Input: stressor var selector ------------------------------
      selectInput(
        inputId = "stressor_vars",
        label = "Select stressor variable",
        choices = exposure_data[, unique(STRESSOR_NAME)],
        selected = exposure_data[, .N, by = STRESSOR_NAME][order(-N), STRESSOR_NAME],
        multiple = TRUE,
        selectize = FALSE,
        size = 5
      ),
      ## -- Input: measured variables ------------------------------
      # selectInput(
      #   inputId = "medium_vars",
      #   label = "Select measure variable",
      #   choices = list(
      #     "Water concentration" = "WaterConc",
      #     "Fish concentration" = "FishCond"
      #   ),
      #   selected = "WaterConc",
      #   multiple = TRUE,
      #   selectize = TRUE
      # ),
      ## -- Summary over site ------------------------------
      gt::gt_output("site_summary")
    ),
    ## -- Main Panel ------------------------------
    mainPanel = shiny::mainPanel(
      ## -- Leaflet map: sample location ------------------------------
      shinydashboard::box(
        title = "Sampling sites",
        width = 12,
        leaflet::leafletOutput("map"),
      ),
      conditionalPanel(
        "input.map_marker_click",
        ## -- Tab-panel: summary table ------------------------------
        shiny::tabsetPanel(
          id = "tables",
          selected = "summary_table",
          shiny::tabPanel(
            "Summary Table",
            value = "summary_table",
            icon = shiny::icon("table"),
            shinydashboard::box(
              title = "Average Concentration of stressor",
              width = 12,
              gt::gt_output("summary_table"),
            )
          ),
        ## -- Tab-panel: Data table ------------------------------
          shiny::tabPanel(
            "Data Table",
            value = "data_table",
            icon = shiny::icon("table"),
            shinydashboard::box(
              title = "Concentration of stressor",
              width = 12,
              gt::gt_output("data_table"),
            )
          ),
        ## -- Tab-panel: Plots ------------------------------
          shiny::tabPanel(
            "By Month",
            value = "by_month",
            icon = shiny::icon("chart-simple"),
            shiny::tabsetPanel(
              shiny::tabPanel(
                icon = shiny::icon("chart-line"),
                title = "Line chart",
                plotOutput("lineplot_by_month"),
              ),
              shiny::tabPanel(
                icon = shiny::icon("chart-area"),
                title = "Scatter plot",
                plotOutput("scatterplot_by_month"),
              )
            )
          ),
          shiny::tabPanel(
            "By Weekday",
            value = "by_wday",
            icon = shiny::icon("chart-simple"),
            shiny::tabsetPanel(
              shiny::tabPanel(
                icon = shiny::icon("chart-line"),
                title = "Line chart",
                plotOutput("lineplot_by_wday"),
              ),
              shiny::tabPanel(
                icon = shiny::icon("chart-area"),
                title = "Scatter plot",
                plotOutput("scatterplot_by_wday"),
              )
            )
          )
        )
      )
    ) # End: Main panel
  ) # End: Page with sidebar
) # End: UI

ui2 <- shiny::fluidPage(
  title = "Exposure Data Analysis",
  theme = shinythemes::shinytheme(theme = "flatly"),
)

## -- Define server logic ----------------------
server <- function(input, output) {
  observeEvent(input$map_marker_click, {
    pts <- input$map_marker_click
    output$summary_table <- gt::render_gt({
      exposure_data %>% 
        filter(SITE_CODE == pts[["id"]]) %>% 
        select(SITE_NAME, STRESSOR_NAME, WaterConc, FishConc) %>% 
        group_by(SITE_NAME, STRESSOR_NAME) %>% 
        summarize(across(.fns = function(x) {
          Hmisc::smean.cl.normal(x) %>% 
            as.list() %>% 
            bind_cols() %>% 
            list()
        })) %>% unnest(names_sep = "_") %>%
        gt::gt(rowname_col = "STRESSOR_NAME", groupname_col = "SITE_NAME") %>% 
        gt::opt_vertical_padding(0.25) %>% 
        gt::data_color(
          columns = gt::ends_with("Mean"),
          direction = "column",
          palette = "Blues"
        ) %>% 
        gt::sub_missing() %>% 
        gt::fmt_number() %>% 
        gt::cols_merge(
          c("WaterConc_Lower", "WaterConc_Upper"),
          pattern = "({1}, {2})"
        ) %>% 
        gt::cols_merge(
          c("FishConc_Lower", "FishConc_Upper"),
          pattern = "({1}, {2})"
        ) %>% 
        gt::cols_merge(c("WaterConc_Mean", "WaterConc_Lower")) %>% 
        gt::cols_merge(c("FishConc_Mean", "FishConc_Lower")) %>% 
        gt::tab_spanner_delim("_") %>% 
        gt::cols_label_with(
          fn = \(x) gsub("Mean", "Mean (95% CI)", x)
        ) %>% 
        gt::tab_style(
          style = gt::cell_text(align = "center"),
          locations = gt::cells_column_labels()
        ) %>% 
        gt::text_transform(
          fn = \(x) gsub("(.*)Conc", "\\1 concentration", x),
          locations = gt::cells_column_spanners()
        ) %>% 
        gt::tab_options()
    })
    output$data_table <- gt::render_gt({
      exposure_data %>% 
        filter(SITE_CODE == pts[["id"]]) %>% 
        select(SITE_NAME, STRESSOR_NAME, WaterConc, FishConc) %>% 
        gt::gt(rowname_col = "STRESSOR_NAME", groupname_col = "SITE_NAME") %>% 
        gt::opt_vertical_padding(0.25) %>% 
        gt::data_color(
          columns = gt::ends_with("Conc"),
          direction = "column",
          palette = "Blues"
        )
    })
    output$scatterplot_by_month <- renderPlot({
      exposure_data %>% 
        filter(SITE_CODE == pts[["id"]]) %>% 
        ggplot(aes(WaterConc, FishConc)) +
        facet_wrap(
          facets = vars(Month),
          scales = "free",
        ) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "lm", formula = "y ~ x") +
        scale_x_log10("log(Water concentration)") +
        scale_y_log10("log(Fish concentration)") +
        ggthemes::theme_few() +
        theme(
          panel.grid = element_line(color = "gray95")
        )
        
    }, res = 100, height = 500)
    output$lineplot_by_month <- renderPlot({
      exposure_data %>% 
        filter(SITE_CODE == pts[["id"]]) %>% 
        pivot_longer(cols = ends_with("Conc")) %>% 
        ggplot(aes(Month, value, group = name)) +
        facet_grid(
          rows = vars(name), 
          scales = "free_y",
          labeller = labeller(
            name = c(
              FishConc = "Fish concentration",
              WaterConc = "Water concentration"
            )
          )
        ) +
        stat_summary(
          fun.data = mean_se,
          geom = "ribbon",
          color = "firebrick",
          fill = NA,
          linetype = "dashed"
        ) +
        stat_summary(fun = mean, geom = "line") +
        stat_summary(
          fun = mean, geom = "point", 
          shape = 21, fill = "whitesmoke"
        ) +
        ggthemes::theme_few() +
        theme(
          panel.grid = element_line(color = "gray95")
        ) +
        labs(y = "Concentraton level of stressor")
        
    }, res = 150, height = 500)
    output$scatterplot_by_wday <- renderPlot({
      exposure_data %>% 
        filter(SITE_CODE == pts[["id"]]) %>% 
        ggplot(aes(WaterConc, FishConc)) +
        facet_wrap(
          facets = vars(Wday),
          scales = "free",
        ) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "lm", formula = "y ~ x") +
        scale_x_log10("log(Water concentration)") +
        scale_y_log10("log(Fish concentration)") +
        ggthemes::theme_few() +
        theme(
          panel.grid = element_line(color = "gray95")
        )
        
    }, res = 100, height = 500)
    output$lineplot_by_wday <- renderPlot({
      exposure_data %>% 
        filter(SITE_CODE == pts[["id"]]) %>% 
        pivot_longer(cols = ends_with("Conc")) %>% 
        ggplot(aes(Wday, value, group = name)) +
        facet_grid(
          rows = vars(name), 
          scales = "free_y",
          labeller = labeller(
            name = c(
              FishConc = "Fish concentration",
              WaterConc = "Water concentration"
            )
          )
        ) +
        stat_summary(
          fun.data = mean_se,
          geom = "ribbon",
          color = "firebrick",
          fill = NA,
          linetype = "dashed"
        ) +
        stat_summary(fun = mean, geom = "line") +
        stat_summary(
          fun = mean, geom = "point", 
          shape = 21, fill = "whitesmoke"
        ) +
        ggthemes::theme_few() +
        theme(
          panel.grid = element_line(color = "gray95")
        ) +
        labs(
          x = "Weekday",
          y = "Concentraton level of stressor"
        )
        
    }, res = 150, height = 500)
    
    
  })
  output$map <- leaflet::renderLeaflet({
    leaflet::leaflet(
      data = exposure_data,
      options = leafletOptions(
        zoomControl = FALSE,
        minZoom = 7, 
        maxZoom = 7
      )
    ) %>% 
      leaflet::addProviderTiles(providers$Esri.WorldStreetMap) %>% 
      addMarkers(
        ~LONGITUDE, ~LATITUDE, 
        label = ~SITE_NAME, 
        layerId = ~SITE_CODE
      )
  })
  output$box_site <- renderPlot({
    exposure_data %>% 
      filter(SITE_NAME %in% input$site_vars) %>% 
      filter(STRESSOR_NAME %in% input$stressor_vars) %>% 
      select(SITE_NAME, STRESSOR_NAME, input$medium_vars, Month) %>% 
      ggplot(aes(Month, get(input$medium_vars))) +
      facet_wrap(facets = vars(SITE_NAME)) +
      geom_boxplot(linewidth = 0.25, outliers = FALSE, color = "gray60") +
      geom_point(
        alpha = 0.5, color = "gray40",
        position = position_jitter(width = 0.15)
      ) +
      scale_y_log10() +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      stat_summary(
        fun.data = mean_se,
        geom = "pointrange",
        color = "firebrick",
        shape = 22
      ) +
      ggthemes::theme_few() +
      theme(
        panel.grid = element_line(color = "#f0f0f0")
      ) +
      labs(
        x = "Month",
        y = glue::glue("log({input$medium_vars})")
      )
  }, res = 150)
  output$site_summary <- gt::render_gt({
    exposure_data %>% 
      pivot_longer(
        cols = c(WaterConc, FishConc),
        names_to = "Medium",
        values_to = "Concentration"
      ) %>% 
      group_by(SITE_NAME, Medium) %>%
      summarize(
        N = length(Concentration),
        Mean = mean(Concentration, na.rm = TRUE),
        SD = sd(Concentration, na.rm = TRUE)
      ) %>% 
      gt::gt(rowname_col = "SITE_NAME", groupname_col = "Medium") %>% 
      gt::fmt_number(-c("N")) %>% 
      gt::opt_vertical_padding(0.25) %>% 
      gt::tab_options(
        table.width = "100%",
        container.width = "100%",
        table.font.size = "10px",
        column_labels.font.weight = "bold",
        row_group.font.weight = "bold"
      )
  })
}

## -- Run the application ------------------------------
shinyApp(ui = ui, server = server)
