Spatio-Temporal Exposure Analysis
================
Raju Rimal
5/12/24

<details>
<summary>Local functions</summary>

``` r
## -- Get summary of a numeric vector ------------------------------------
get_summary <- function(x) {
  list(tidytable(
    N = length(x),
    Mean = mean(x, na.rm = TRUE),
    Max = max(x, na.rm = TRUE),
    Min = min(x, na.rm = TRUE),
    StdDev = sd(x, na.rm = TRUE),
    P5 = quantile(x, 0.05, na.rm = TRUE),
    P95 = quantile(x, 0.95, na.rm = TRUE)
  ))
}
```

</details>

## Explore the dataset

<details>
<summary>Import data from excel file</summary>

``` r
## -- Import different sheets from the data file -----------------------------
codebook <- read_excel("Exposure_data_AEP.xlsx", sheet = 1)
raw_data <- read_excel("Exposure_data_AEP.xlsx", sheet = 2)
site_data <- read_excel("Exposure_data_AEP.xlsx", sheet = 3)
```

</details>
<details>
<summary>Fetch chemical information from API</summary>

``` r
## -- Fetch extra chemical information of the stressor from API
## -- The fetched information is merged with exposure data and saved as csv
## -- We also merge the site (location) information
## -- This block of code only runs if the csv file is missing
if (!file.exists("Exposure_data_AEP.csv")) {
  cid_data <- raw_data %>% 
    select(INCHIKEY) %>% 
    unique() %>% 
    map_df(webchem::get_cid, from = "inchikey", match = "first")
  
  cprop_data <- webchem::pc_prop(cid_data[["cid"]])
  
  exposure_data <- cprop_data %>% 
    mutate(cid = as.character(CID)) %>% 
    select(cid, MolecularFormula, MolecularWeight, 
           InChI, IUPACName, XLogP, ExactMass) %>% 
    inner_join(cid_data, by = "cid") %>% 
    rename(INCHIKEY = query) %>% 
    right_join(raw_data, by = "INCHIKEY") %>% 
    left_join(site_data, by = "SITE_CODE")
  
  fwrite(exposure_data, file = "Exposure_data_AEP.csv")
} else {
  exposure_data <- fread("Exposure_data_AEP.csv")
}
```

</details>

<div class="panel-tabset">

## Boxplot

### Spatial variation of MCPA

<details>
<summary>Box plot of MCPA by site</summary>

``` r
make_box <- function(data, x_var, x_lab) {
  plot_data %>% 
  ggplot(aes(get(x_var), MEASURED_VALUE)) +
  geom_boxplot(outliers = FALSE, width = 0.75) +
  geom_point(
    position = position_jitter(width = 0.25),
    shape = 21, size = rel(1.5), stroke = 1,
    color = "gray40", fill = "grey80"
  ) +
  stat_summary(
    aes(color = "Mean (SE)"),
    fun.data = mean_se, 
    geom = "pointrange", 
    fill = "whitesmoke",
    shape = 22,
    stroke = 1.5
  ) +
  scale_y_continuous(breaks = scales::breaks_extended(10)) +
  scale_color_manual(
    name = NULL,
    values = "firebrick",
    breaks = "Mean (SE)"
  ) +
  ggthemes::theme_few() +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1.1, 1.1)
  ) +
  labs(
    x = x_lab,
    y = NULL
  )
}

plot_data <- exposure_data[STRESSOR_ID == 21]

make_box(plot_data, "SITE_NAME", "Site") +
  labs(y = paste(
    "Water concentration (ug/L) of",
    plot_data[, first(STRESSOR_NAME)],
    sep = "\n"
  ))
```

</details>

<img
src="index_files/figure-commonmark/Box%20plot%20MCPA%20by%20site-1.png"
style="width:100.0%" />

### Temporal variation of MCPA

<details>
<summary>Box plot of MCPA by month and week day</summary>

``` r
plot_data <- exposure_data %>% 
  filter(STRESSOR_ID == 21) %>% 
  mutate(
    Month = lubridate::month(SAMPLE_DATE, label = TRUE, abbr = FALSE),
    Wday = lubridate::wday(SAMPLE_DATE, label = TRUE, abbr = FALSE)
  )

by_month <- make_box(plot_data, "Month", "Month")
by_wday <- make_box(plot_data, "Wday", "Week day")

patchwork::wrap_plots(by_month, by_wday, nrow = 2) +
  labs(tag = paste(
      "Water concentration (ug/L) of",
      plot_data[, first(STRESSOR_NAME)],
      sep = "\n"
    )) +
  theme(
    plot.tag = element_text(size = rel(1), angle = 90, vjust = -0.5),
    plot.tag.position = c(0, 1),
    plot.margin = margin(l = 48)
  )
```

</details>

<img
src="index_files/figure-commonmark/Box%20plot%20MCPA%20over%20time-1.png"
style="width:100.0%" />

## Lineplot

### Spatial variation of MCPA

<details>
<summary>Line plot of MCPA by site</summary>

``` r
make_lines <- function(data, x_var, x_lab) {
  ggplot(data, aes(get(x_var), MEASURED_VALUE)) +
    stat_summary(
      fun.data = "mean_se",
      geom = "ribbon",
      color = "firebrick",
      fill = NA,
      linetype = "dashed",
      group = 1
    ) +
    stat_summary(
      fun = "mean",
      geom = "line",
      group = 1
    ) +
    stat_summary(
      fun = "mean",
      geom = "point",
      shape = 21, size = rel(1.5), stroke = 1,
      color = "gray40", fill = "grey80"
    ) +
    scale_y_continuous(breaks = scales::breaks_extended(10)) +
    ggthemes::theme_few() +
    theme(
      legend.position = c(1, 1),
      legend.justification = c(1.1, 1.1)
    ) +
    labs(
      x = x_lab,
      y = paste(
        "Water concentration (ug/L) of",
        plot_data[, first(STRESSOR_NAME)],
        sep = "\n"
      )
    )
}

plot_data <- exposure_data[STRESSOR_ID == 21]

make_lines(plot_data, "SITE_NAME", "Site")
```

</details>

<img
src="index_files/figure-commonmark/Line%20plot%20MCPA%20by%20site-1.png"
style="width:100.0%" />

### Temporal variation of MCPA

<details>
<summary>Line plot of MCPA by month and week day</summary>

``` r
plot_data <- exposure_data %>% 
  filter(STRESSOR_ID == 21) %>% 
  mutate(
    Month = lubridate::month(SAMPLE_DATE, label = TRUE, abbr = FALSE),
    Wday = lubridate::wday(SAMPLE_DATE, label = TRUE, abbr = FALSE)
  )

by_month <- make_lines(plot_data, "Month", "Month")
by_wday <- make_lines(plot_data, "Wday", "Week day")

patchwork::wrap_plots(by_month, by_wday, nrow = 2) +
  labs(tag = paste(
      "Water concentration (ug/L) of",
      plot_data[, first(STRESSOR_NAME)],
      sep = "\n"
    )) +
  theme(
    plot.tag = element_text(size = rel(1), angle = 90, vjust = -0.5),
    plot.tag.position = c(0, 1),
    plot.margin = margin(l = 48)
  )
```

</details>

<img
src="index_files/figure-commonmark/Line%20plot%20MCPA%20over%20time-1.png"
style="width:100.0%" />

</div>

## Spatio-Temporal variation of MCPA

<details>
<summary>Spatio-Temporal variation of MCPA</summary>

``` r
plot_data <- exposure_data %>% 
  filter(STRESSOR_ID == 21) %>% 
  mutate(Month = lubridate::month(SAMPLE_DATE, label = TRUE)) %>% 
  mutate(Wday = lubridate::wday(SAMPLE_DATE, label = TRUE))

make_fitted_plot <- function(data, x_var, x_lab) {
  data %>% 
    ggplot(aes(as.numeric(get(x_var)), MEASURED_VALUE)) +
      geom_smooth(
        method = "glm", 
        color = "firebrick", 
        formula = "y ~ x", 
        method.args = list(family = gaussian(link = "log"))
      ) +
      geom_point() +
      facet_wrap(
        facets = vars(SITE_NAME), 
        nrow = 1, 
        scales = "free"
      ) +
      scale_x_continuous(
        breaks = scales::breaks_width(1),
        labels = \(x) data[, levels(get(x_var))][x]
      ) +
      ggthemes::theme_few() +
      theme(
        panel.grid = element_line(color = "#f0f0f0")
      ) +
      labs(
        x = x_lab,
        y = NULL
      )
}

by_month <- make_fitted_plot(plot_data, "Month", "Sampled Month")
by_wday <- make_fitted_plot(plot_data, "Wday", "Sampled day")

patchwork::wrap_plots(by_month, by_wday, nrow = 2) +
  labs(tag = paste(
      "Water concentration (ug/L) of",
      plot_data[, first(STRESSOR_NAME)],
      sep = "\n"
    )) +
  theme(
    plot.tag = element_text(size = rel(1), angle = 90, vjust = -0.5),
    plot.tag.position = c(0, 1),
    plot.margin = margin(l = 48)
  )
```

</details>

<img
src="index_files/figure-commonmark/Spatio-Temporal%20variation%20of%20MCPA-1.png"
style="width:100.0%" />

## Summary of spatio-temporal exposure data

<details>
<summary>Water concentration (ug/L) summary of stressor</summary>

``` r
summary_table <- exposure_data %>% 
  group_by(SITE_NAME, STRESSOR_NAME) %>% 
  summarise(across("MEASURED_VALUE", get_summary)) %>% 
  unnest(MEASURED_VALUE) %>% 
  arrange(desc(Mean), )
  
summary_table %>% 
  gt::gt(rowname_col = "STRESSOR_NAME", groupname_col = "SITE_NAME") %>% 
  gt::fmt_number(columns = !gt::starts_with("N"), decimals = 3) %>% 
  gt::data_color(
    columns = -1,
    direction = "column",
    palette = "Blues",
    alpha = 0.75
  ) %>% 
  gt::opt_vertical_padding(0.5) %>% 
  gt::tab_options(
    table.font.size = "10pt",
    table.width = "100%",
    row_group.font.weight = "bold",
    column_labels.font.weight = "bold"
  ) %>% 
  gt::sub_missing() %>% 
  gt::tab_stubhead("Stressor Name")
```

</details>
<div id="zieytawzqd" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#zieytawzqd table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#zieytawzqd thead, #zieytawzqd tbody, #zieytawzqd tfoot, #zieytawzqd tr, #zieytawzqd td, #zieytawzqd th {
  border-style: none;
}

#zieytawzqd p {
  margin: 0;
  padding: 0;
}

#zieytawzqd .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 10pt;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 100%;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#zieytawzqd .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#zieytawzqd .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 2px;
  padding-bottom: 2px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#zieytawzqd .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 1px;
  padding-bottom: 3px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#zieytawzqd .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#zieytawzqd .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#zieytawzqd .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#zieytawzqd .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: bold;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 2.5px;
  padding-bottom: 3.5px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#zieytawzqd .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: bold;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#zieytawzqd .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#zieytawzqd .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#zieytawzqd .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 2.5px;
  padding-bottom: 2.5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#zieytawzqd .gt_spanner_row {
  border-bottom-style: hidden;
}

#zieytawzqd .gt_group_heading {
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: bold;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#zieytawzqd .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: bold;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#zieytawzqd .gt_from_md > :first-child {
  margin-top: 0;
}

#zieytawzqd .gt_from_md > :last-child {
  margin-bottom: 0;
}

#zieytawzqd .gt_row {
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#zieytawzqd .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#zieytawzqd .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#zieytawzqd .gt_row_group_first td {
  border-top-width: 2px;
}

#zieytawzqd .gt_row_group_first th {
  border-top-width: 2px;
}

#zieytawzqd .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#zieytawzqd .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#zieytawzqd .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#zieytawzqd .gt_last_summary_row {
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#zieytawzqd .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#zieytawzqd .gt_first_grand_summary_row {
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#zieytawzqd .gt_last_grand_summary_row_top {
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#zieytawzqd .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#zieytawzqd .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#zieytawzqd .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#zieytawzqd .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 2px;
  padding-bottom: 2px;
  padding-left: 5px;
  padding-right: 5px;
}

#zieytawzqd .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#zieytawzqd .gt_sourcenote {
  font-size: 90%;
  padding-top: 2px;
  padding-bottom: 2px;
  padding-left: 5px;
  padding-right: 5px;
}

#zieytawzqd .gt_left {
  text-align: left;
}

#zieytawzqd .gt_center {
  text-align: center;
}

#zieytawzqd .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#zieytawzqd .gt_font_normal {
  font-weight: normal;
}

#zieytawzqd .gt_font_bold {
  font-weight: bold;
}

#zieytawzqd .gt_font_italic {
  font-style: italic;
}

#zieytawzqd .gt_super {
  font-size: 65%;
}

#zieytawzqd .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#zieytawzqd .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#zieytawzqd .gt_indent_1 {
  text-indent: 5px;
}

#zieytawzqd .gt_indent_2 {
  text-indent: 10px;
}

#zieytawzqd .gt_indent_3 {
  text-indent: 15px;
}

#zieytawzqd .gt_indent_4 {
  text-indent: 20px;
}

#zieytawzqd .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Stressor Name">Stressor Name</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="N">N</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Mean">Mean</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Max">Max</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Min">Min</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="StdDev">StdDev</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="P5">P5</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="P95">P95</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="8" class="gt_group_heading" scope="colgroup" id="Heiabekken">Heiabekken</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_1" scope="row" class="gt_row gt_left gt_stub">imidacloprid</th>
<td headers="Heiabekken stub_1_1 N" class="gt_row gt_right" style="background-color: rgba(8,48,107,0.75); color: #FFFFFF;">12</td>
<td headers="Heiabekken stub_1_1 Mean" class="gt_row gt_right" style="background-color: rgba(8,48,107,0.75); color: #FFFFFF;">1.142</td>
<td headers="Heiabekken stub_1_1 Max" class="gt_row gt_right" style="background-color: rgba(8,48,107,0.75); color: #FFFFFF;">5.300</td>
<td headers="Heiabekken stub_1_1 Min" class="gt_row gt_right" style="background-color: rgba(8,48,107,0.75); color: #FFFFFF;">0.210</td>
<td headers="Heiabekken stub_1_1 StdDev" class="gt_row gt_right" style="background-color: rgba(8,48,107,0.75); color: #FFFFFF;">1.384</td>
<td headers="Heiabekken stub_1_1 P5" class="gt_row gt_right" style="background-color: rgba(8,48,107,0.75); color: #FFFFFF;">0.226</td>
<td headers="Heiabekken stub_1_1 P95" class="gt_row gt_right" style="background-color: rgba(8,48,107,0.75); color: #FFFFFF;">3.100</td></tr>
    <tr><th id="stub_1_2" scope="row" class="gt_row gt_left gt_stub">propamocarb</th>
<td headers="Heiabekken stub_1_2 N" class="gt_row gt_right" style="background-color: rgba(9,72,142,0.75); color: #FFFFFF;">11</td>
<td headers="Heiabekken stub_1_2 Mean" class="gt_row gt_right" style="background-color: rgba(127,184,218,0.75); color: #000000;">0.524</td>
<td headers="Heiabekken stub_1_2 Max" class="gt_row gt_right" style="background-color: rgba(38,117,183,0.75); color: #FFFFFF;">3.900</td>
<td headers="Heiabekken stub_1_2 Min" class="gt_row gt_right" style="background-color: rgba(229,239,249,0.75); color: #000000;">0.028</td>
<td headers="Heiabekken stub_1_2 StdDev" class="gt_row gt_right" style="background-color: rgba(21,94,166,0.75); color: #FFFFFF;">1.141</td>
<td headers="Heiabekken stub_1_2 P5" class="gt_row gt_right" style="background-color: rgba(229,239,249,0.75); color: #000000;">0.030</td>
<td headers="Heiabekken stub_1_2 P95" class="gt_row gt_right" style="background-color: rgba(32,112,180,0.75); color: #FFFFFF;">2.340</td></tr>
    <tr><th id="stub_1_3" scope="row" class="gt_row gt_left gt_stub">2-methyl-4-chlorophenoxyacetic acid (MCPA)</th>
<td headers="Heiabekken stub_1_3 N" class="gt_row gt_right" style="background-color: rgba(64,143,196,0.75); color: #FFFFFF;">8</td>
<td headers="Heiabekken stub_1_3 Mean" class="gt_row gt_right" style="background-color: rgba(135,188,220,0.75); color: #000000;">0.503</td>
<td headers="Heiabekken stub_1_3 Max" class="gt_row gt_right" style="background-color: rgba(194,217,238,0.75); color: #000000;">1.400</td>
<td headers="Heiabekken stub_1_3 Min" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Heiabekken stub_1_3 StdDev" class="gt_row gt_right" style="background-color: rgba(139,191,221,0.75); color: #000000;">0.589</td>
<td headers="Heiabekken stub_1_3 P5" class="gt_row gt_right" style="background-color: rgba(240,247,253,0.75); color: #000000;">0.018</td>
<td headers="Heiabekken stub_1_3 P95" class="gt_row gt_right" style="background-color: rgba(138,190,220,0.75); color: #000000;">1.330</td></tr>
    <tr><th id="stub_1_4" scope="row" class="gt_row gt_left gt_stub">pencycuron</th>
<td headers="Heiabekken stub_1_4 N" class="gt_row gt_right" style="background-color: rgba(8,48,107,0.75); color: #FFFFFF;">12</td>
<td headers="Heiabekken stub_1_4 Mean" class="gt_row gt_right" style="background-color: rgba(217,232,245,0.75); color: #000000;">0.179</td>
<td headers="Heiabekken stub_1_4 Max" class="gt_row gt_right" style="background-color: rgba(213,229,244,0.75); color: #000000;">0.930</td>
<td headers="Heiabekken stub_1_4 Min" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.014</td>
<td headers="Heiabekken stub_1_4 StdDev" class="gt_row gt_right" style="background-color: rgba(210,227,243,0.75); color: #000000;">0.263</td>
<td headers="Heiabekken stub_1_4 P5" class="gt_row gt_right" style="background-color: rgba(239,246,252,0.75); color: #000000;">0.019</td>
<td headers="Heiabekken stub_1_4 P95" class="gt_row gt_right" style="background-color: rgba(209,226,243,0.75); color: #000000;">0.605</td></tr>
    <tr><th id="stub_1_5" scope="row" class="gt_row gt_left gt_stub">metribuzin</th>
<td headers="Heiabekken stub_1_5 N" class="gt_row gt_right" style="background-color: rgba(22,95,167,0.75); color: #FFFFFF;">10</td>
<td headers="Heiabekken stub_1_5 Mean" class="gt_row gt_right" style="background-color: rgba(218,232,246,0.75); color: #000000;">0.175</td>
<td headers="Heiabekken stub_1_5 Max" class="gt_row gt_right" style="background-color: rgba(194,217,238,0.75); color: #000000;">1.400</td>
<td headers="Heiabekken stub_1_5 Min" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Heiabekken stub_1_5 StdDev" class="gt_row gt_right" style="background-color: rgba(179,211,232,0.75); color: #000000;">0.431</td>
<td headers="Heiabekken stub_1_5 P5" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.015</td>
<td headers="Heiabekken stub_1_5 P95" class="gt_row gt_right" style="background-color: rgba(196,218,238,0.75); color: #000000;">0.807</td></tr>
    <tr><th id="stub_1_6" scope="row" class="gt_row gt_left gt_stub">fluroxypyr</th>
<td headers="Heiabekken stub_1_6 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Heiabekken stub_1_6 Mean" class="gt_row gt_right" style="background-color: rgba(221,234,247,0.75); color: #000000;">0.160</td>
<td headers="Heiabekken stub_1_6 Max" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.160</td>
<td headers="Heiabekken stub_1_6 Min" class="gt_row gt_right" style="background-color: rgba(33,113,181,0.75); color: #FFFFFF;">0.160</td>
<td headers="Heiabekken stub_1_6 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Heiabekken stub_1_6 P5" class="gt_row gt_right" style="background-color: rgba(50,128,189,0.75); color: #FFFFFF;">0.160</td>
<td headers="Heiabekken stub_1_6 P95" class="gt_row gt_right" style="background-color: rgba(237,245,252,0.75); color: #000000;">0.160</td></tr>
    <tr><th id="stub_1_7" scope="row" class="gt_row gt_left gt_stub">mandipropamid</th>
<td headers="Heiabekken stub_1_7 N" class="gt_row gt_right" style="background-color: rgba(211,228,243,0.75); color: #000000;">3</td>
<td headers="Heiabekken stub_1_7 Mean" class="gt_row gt_right" style="background-color: rgba(223,236,247,0.75); color: #000000;">0.145</td>
<td headers="Heiabekken stub_1_7 Max" class="gt_row gt_right" style="background-color: rgba(234,243,251,0.75); color: #000000;">0.360</td>
<td headers="Heiabekken stub_1_7 Min" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.014</td>
<td headers="Heiabekken stub_1_7 StdDev" class="gt_row gt_right" style="background-color: rgba(220,234,246,0.75); color: #000000;">0.187</td>
<td headers="Heiabekken stub_1_7 P5" class="gt_row gt_right" style="background-color: rgba(239,246,252,0.75); color: #000000;">0.019</td>
<td headers="Heiabekken stub_1_7 P95" class="gt_row gt_right" style="background-color: rgba(226,238,248,0.75); color: #000000;">0.330</td></tr>
    <tr><th id="stub_1_8" scope="row" class="gt_row gt_left gt_stub">clopyralid</th>
<td headers="Heiabekken stub_1_8 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Heiabekken stub_1_8 Mean" class="gt_row gt_right" style="background-color: rgba(230,240,249,0.75); color: #000000;">0.110</td>
<td headers="Heiabekken stub_1_8 Max" class="gt_row gt_right" style="background-color: rgba(243,249,254,0.75); color: #000000;">0.110</td>
<td headers="Heiabekken stub_1_8 Min" class="gt_row gt_right" style="background-color: rgba(107,174,214,0.75); color: #000000;">0.110</td>
<td headers="Heiabekken stub_1_8 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Heiabekken stub_1_8 P5" class="gt_row gt_right" style="background-color: rgba(124,183,217,0.75); color: #000000;">0.110</td>
<td headers="Heiabekken stub_1_8 P95" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.110</td></tr>
    <tr><th id="stub_1_9" scope="row" class="gt_row gt_left gt_stub">boscalid</th>
<td headers="Heiabekken stub_1_9 N" class="gt_row gt_right" style="background-color: rgba(8,48,107,0.75); color: #FFFFFF;">12</td>
<td headers="Heiabekken stub_1_9 Mean" class="gt_row gt_right" style="background-color: rgba(242,248,253,0.75); color: #000000;">0.041</td>
<td headers="Heiabekken stub_1_9 Max" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.081</td>
<td headers="Heiabekken stub_1_9 Min" class="gt_row gt_right" style="background-color: rgba(235,243,251,0.75); color: #000000;">0.022</td>
<td headers="Heiabekken stub_1_9 StdDev" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.018</td>
<td headers="Heiabekken stub_1_9 P5" class="gt_row gt_right" style="background-color: rgba(234,243,251,0.75); color: #000000;">0.024</td>
<td headers="Heiabekken stub_1_9 P95" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.071</td></tr>
    <tr><th id="stub_1_10" scope="row" class="gt_row gt_left gt_stub">prothioconazole-desthio</th>
<td headers="Heiabekken stub_1_10 N" class="gt_row gt_right" style="background-color: rgba(64,143,196,0.75); color: #FFFFFF;">8</td>
<td headers="Heiabekken stub_1_10 Mean" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.036</td>
<td headers="Heiabekken stub_1_10 Max" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.094</td>
<td headers="Heiabekken stub_1_10 Min" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.010</td>
<td headers="Heiabekken stub_1_10 StdDev" class="gt_row gt_right" style="background-color: rgba(243,249,254,0.75); color: #000000;">0.028</td>
<td headers="Heiabekken stub_1_10 P5" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.010</td>
<td headers="Heiabekken stub_1_10 P95" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.078</td></tr>
    <tr><th id="stub_1_11" scope="row" class="gt_row gt_left gt_stub">metalaxyl</th>
<td headers="Heiabekken stub_1_11 N" class="gt_row gt_right" style="background-color: rgba(40,119,184,0.75); color: #FFFFFF;">9</td>
<td headers="Heiabekken stub_1_11 Mean" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.033</td>
<td headers="Heiabekken stub_1_11 Max" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.059</td>
<td headers="Heiabekken stub_1_11 Min" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Heiabekken stub_1_11 StdDev" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.017</td>
<td headers="Heiabekken stub_1_11 P5" class="gt_row gt_right" style="background-color: rgba(245,249,254,0.75); color: #000000;">0.013</td>
<td headers="Heiabekken stub_1_11 P95" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.059</td></tr>
    <tr><th id="stub_1_12" scope="row" class="gt_row gt_left gt_stub">trifloxystrobin</th>
<td headers="Heiabekken stub_1_12 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Heiabekken stub_1_12 Mean" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.026</td>
<td headers="Heiabekken stub_1_12 Max" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.026</td>
<td headers="Heiabekken stub_1_12 Min" class="gt_row gt_right" style="background-color: rgba(231,241,250,0.75); color: #000000;">0.026</td>
<td headers="Heiabekken stub_1_12 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Heiabekken stub_1_12 P5" class="gt_row gt_right" style="background-color: rgba(233,242,250,0.75); color: #000000;">0.026</td>
<td headers="Heiabekken stub_1_12 P95" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.026</td></tr>
    <tr><th id="stub_1_13" scope="row" class="gt_row gt_left gt_stub">metamitron</th>
<td headers="Heiabekken stub_1_13 N" class="gt_row gt_right" style="background-color: rgba(162,204,226,0.75); color: #000000;">5</td>
<td headers="Heiabekken stub_1_13 Mean" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.020</td>
<td headers="Heiabekken stub_1_13 Max" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.025</td>
<td headers="Heiabekken stub_1_13 Min" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.012</td>
<td headers="Heiabekken stub_1_13 StdDev" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.006</td>
<td headers="Heiabekken stub_1_13 P5" class="gt_row gt_right" style="background-color: rgba(245,249,254,0.75); color: #000000;">0.013</td>
<td headers="Heiabekken stub_1_13 P95" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.025</td></tr>
    <tr><th id="stub_1_14" scope="row" class="gt_row gt_left gt_stub">propoxycarbazone</th>
<td headers="Heiabekken stub_1_14 N" class="gt_row gt_right" style="background-color: rgba(191,216,236,0.75); color: #000000;">4</td>
<td headers="Heiabekken stub_1_14 Mean" class="gt_row gt_right" style="background-color: rgba(245,250,255,0.75); color: #000000;">0.020</td>
<td headers="Heiabekken stub_1_14 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.023</td>
<td headers="Heiabekken stub_1_14 Min" class="gt_row gt_right" style="background-color: rgba(240,246,253,0.75); color: #000000;">0.017</td>
<td headers="Heiabekken stub_1_14 StdDev" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.003</td>
<td headers="Heiabekken stub_1_14 P5" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.017</td>
<td headers="Heiabekken stub_1_14 P95" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.023</td></tr>
    <tr><th id="stub_1_15" scope="row" class="gt_row gt_left gt_stub">diflufenican</th>
<td headers="Heiabekken stub_1_15 N" class="gt_row gt_right" style="background-color: rgba(162,204,226,0.75); color: #000000;">5</td>
<td headers="Heiabekken stub_1_15 Mean" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.018</td>
<td headers="Heiabekken stub_1_15 Max" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.034</td>
<td headers="Heiabekken stub_1_15 Min" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.010</td>
<td headers="Heiabekken stub_1_15 StdDev" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.009</td>
<td headers="Heiabekken stub_1_15 P5" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.010</td>
<td headers="Heiabekken stub_1_15 P95" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.031</td></tr>
    <tr><th id="stub_1_16" scope="row" class="gt_row gt_left gt_stub">bixafen</th>
<td headers="Heiabekken stub_1_16 N" class="gt_row gt_right" style="background-color: rgba(191,216,236,0.75); color: #000000;">4</td>
<td headers="Heiabekken stub_1_16 Mean" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.017</td>
<td headers="Heiabekken stub_1_16 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.023</td>
<td headers="Heiabekken stub_1_16 Min" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.010</td>
<td headers="Heiabekken stub_1_16 StdDev" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.006</td>
<td headers="Heiabekken stub_1_16 P5" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.010</td>
<td headers="Heiabekken stub_1_16 P95" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.023</td></tr>
    <tr><th id="stub_1_17" scope="row" class="gt_row gt_left gt_stub">2,6-dichlorobenzamide (BAM)</th>
<td headers="Heiabekken stub_1_17 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Heiabekken stub_1_17 Mean" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.016</td>
<td headers="Heiabekken stub_1_17 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.016</td>
<td headers="Heiabekken stub_1_17 Min" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.016</td>
<td headers="Heiabekken stub_1_17 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Heiabekken stub_1_17 P5" class="gt_row gt_right" style="background-color: rgba(242,248,253,0.75); color: #000000;">0.016</td>
<td headers="Heiabekken stub_1_17 P95" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.016</td></tr>
    <tr><th id="stub_1_18" scope="row" class="gt_row gt_left gt_stub">pyraclostrobin</th>
<td headers="Heiabekken stub_1_18 N" class="gt_row gt_right" style="background-color: rgba(229,239,249,0.75); color: #000000;">2</td>
<td headers="Heiabekken stub_1_18 Mean" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.015</td>
<td headers="Heiabekken stub_1_18 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.020</td>
<td headers="Heiabekken stub_1_18 Min" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Heiabekken stub_1_18 StdDev" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.006</td>
<td headers="Heiabekken stub_1_18 P5" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Heiabekken stub_1_18 P95" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.020</td></tr>
    <tr><th id="stub_1_19" scope="row" class="gt_row gt_left gt_stub">spinosad</th>
<td headers="Heiabekken stub_1_19 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Heiabekken stub_1_19 Mean" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.013</td>
<td headers="Heiabekken stub_1_19 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.013</td>
<td headers="Heiabekken stub_1_19 Min" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.013</td>
<td headers="Heiabekken stub_1_19 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Heiabekken stub_1_19 P5" class="gt_row gt_right" style="background-color: rgba(245,249,254,0.75); color: #000000;">0.013</td>
<td headers="Heiabekken stub_1_19 P95" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.013</td></tr>
    <tr><th id="stub_1_20" scope="row" class="gt_row gt_left gt_stub">1,1'-(2,2-Dichloro-1,1-ethenediyl)bis(4-chlorobenzene)</th>
<td headers="Heiabekken stub_1_20 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Heiabekken stub_1_20 Mean" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.012</td>
<td headers="Heiabekken stub_1_20 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.012</td>
<td headers="Heiabekken stub_1_20 Min" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.012</td>
<td headers="Heiabekken stub_1_20 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Heiabekken stub_1_20 P5" class="gt_row gt_right" style="background-color: rgba(245,250,255,0.75); color: #000000;">0.012</td>
<td headers="Heiabekken stub_1_20 P95" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.012</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="8" class="gt_group_heading" scope="colgroup" id="Skuterudbekken">Skuterudbekken</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_21" scope="row" class="gt_row gt_left gt_stub">2-methyl-4-chlorophenoxyacetic acid (MCPA)</th>
<td headers="Skuterudbekken stub_1_21 N" class="gt_row gt_right" style="background-color: rgba(93,164,208,0.75); color: #FFFFFF;">7</td>
<td headers="Skuterudbekken stub_1_21 Mean" class="gt_row gt_right" style="background-color: rgba(161,203,226,0.75); color: #000000;">0.425</td>
<td headers="Skuterudbekken stub_1_21 Max" class="gt_row gt_right" style="background-color: rgba(194,217,238,0.75); color: #000000;">1.400</td>
<td headers="Skuterudbekken stub_1_21 Min" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Skuterudbekken stub_1_21 StdDev" class="gt_row gt_right" style="background-color: rgba(134,188,220,0.75); color: #000000;">0.605</td>
<td headers="Skuterudbekken stub_1_21 P5" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.013</td>
<td headers="Skuterudbekken stub_1_21 P95" class="gt_row gt_right" style="background-color: rgba(137,190,220,0.75); color: #000000;">1.340</td></tr>
    <tr><th id="stub_1_22" scope="row" class="gt_row gt_left gt_stub">fluroxypyr</th>
<td headers="Skuterudbekken stub_1_22 N" class="gt_row gt_right" style="background-color: rgba(162,204,226,0.75); color: #000000;">5</td>
<td headers="Skuterudbekken stub_1_22 Mean" class="gt_row gt_right" style="background-color: rgba(204,223,241,0.75); color: #000000;">0.258</td>
<td headers="Skuterudbekken stub_1_22 Max" class="gt_row gt_right" style="background-color: rgba(225,237,248,0.75); color: #000000;">0.600</td>
<td headers="Skuterudbekken stub_1_22 Min" class="gt_row gt_right" style="background-color: rgba(148,196,223,0.75); color: #000000;">0.090</td>
<td headers="Skuterudbekken stub_1_22 StdDev" class="gt_row gt_right" style="background-color: rgba(218,232,246,0.75); color: #000000;">0.206</td>
<td headers="Skuterudbekken stub_1_22 P5" class="gt_row gt_right" style="background-color: rgba(157,201,225,0.75); color: #000000;">0.092</td>
<td headers="Skuterudbekken stub_1_22 P95" class="gt_row gt_right" style="background-color: rgba(214,229,244,0.75); color: #000000;">0.530</td></tr>
    <tr><th id="stub_1_23" scope="row" class="gt_row gt_left gt_stub">clopyralid</th>
<td headers="Skuterudbekken stub_1_23 N" class="gt_row gt_right" style="background-color: rgba(191,216,236,0.75); color: #000000;">4</td>
<td headers="Skuterudbekken stub_1_23 Mean" class="gt_row gt_right" style="background-color: rgba(225,237,248,0.75); color: #000000;">0.134</td>
<td headers="Skuterudbekken stub_1_23 Max" class="gt_row gt_right" style="background-color: rgba(239,246,252,0.75); color: #000000;">0.230</td>
<td headers="Skuterudbekken stub_1_23 Min" class="gt_row gt_right" style="background-color: rgba(207,225,242,0.75); color: #000000;">0.051</td>
<td headers="Skuterudbekken stub_1_23 StdDev" class="gt_row gt_right" style="background-color: rgba(236,244,252,0.75); color: #000000;">0.078</td>
<td headers="Skuterudbekken stub_1_23 P5" class="gt_row gt_right" style="background-color: rgba(204,223,241,0.75); color: #000000;">0.058</td>
<td headers="Skuterudbekken stub_1_23 P95" class="gt_row gt_right" style="background-color: rgba(234,242,251,0.75); color: #000000;">0.219</td></tr>
    <tr><th id="stub_1_24" scope="row" class="gt_row gt_left gt_stub">prothioconazole-desthio</th>
<td headers="Skuterudbekken stub_1_24 N" class="gt_row gt_right" style="background-color: rgba(191,216,236,0.75); color: #000000;">4</td>
<td headers="Skuterudbekken stub_1_24 Mean" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.047</td>
<td headers="Skuterudbekken stub_1_24 Max" class="gt_row gt_right" style="background-color: rgba(243,249,254,0.75); color: #000000;">0.110</td>
<td headers="Skuterudbekken stub_1_24 Min" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.014</td>
<td headers="Skuterudbekken stub_1_24 StdDev" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.045</td>
<td headers="Skuterudbekken stub_1_24 P5" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.014</td>
<td headers="Skuterudbekken stub_1_24 P95" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.101</td></tr>
    <tr><th id="stub_1_25" scope="row" class="gt_row gt_left gt_stub">dichlorprop</th>
<td headers="Skuterudbekken stub_1_25 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Skuterudbekken stub_1_25 Mean" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.035</td>
<td headers="Skuterudbekken stub_1_25 Max" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.035</td>
<td headers="Skuterudbekken stub_1_25 Min" class="gt_row gt_right" style="background-color: rgba(222,235,247,0.75); color: #000000;">0.035</td>
<td headers="Skuterudbekken stub_1_25 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Skuterudbekken stub_1_25 P5" class="gt_row gt_right" style="background-color: rgba(224,236,248,0.75); color: #000000;">0.035</td>
<td headers="Skuterudbekken stub_1_25 P95" class="gt_row gt_right" style="background-color: rgba(245,250,255,0.75); color: #000000;">0.035</td></tr>
    <tr><th id="stub_1_26" scope="row" class="gt_row gt_left gt_stub">imidacloprid</th>
<td headers="Skuterudbekken stub_1_26 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Skuterudbekken stub_1_26 Mean" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.021</td>
<td headers="Skuterudbekken stub_1_26 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.021</td>
<td headers="Skuterudbekken stub_1_26 Min" class="gt_row gt_right" style="background-color: rgba(236,244,251,0.75); color: #000000;">0.021</td>
<td headers="Skuterudbekken stub_1_26 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Skuterudbekken stub_1_26 P5" class="gt_row gt_right" style="background-color: rgba(237,245,252,0.75); color: #000000;">0.021</td>
<td headers="Skuterudbekken stub_1_26 P95" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.021</td></tr>
    <tr><th id="stub_1_27" scope="row" class="gt_row gt_left gt_stub">bixafen</th>
<td headers="Skuterudbekken stub_1_27 N" class="gt_row gt_right" style="background-color: rgba(211,228,243,0.75); color: #000000;">3</td>
<td headers="Skuterudbekken stub_1_27 Mean" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.021</td>
<td headers="Skuterudbekken stub_1_27 Max" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.028</td>
<td headers="Skuterudbekken stub_1_27 Min" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.014</td>
<td headers="Skuterudbekken stub_1_27 StdDev" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.007</td>
<td headers="Skuterudbekken stub_1_27 P5" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.015</td>
<td headers="Skuterudbekken stub_1_27 P95" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.027</td></tr>
    <tr><th id="stub_1_28" scope="row" class="gt_row gt_left gt_stub">propiconazole</th>
<td headers="Skuterudbekken stub_1_28 N" class="gt_row gt_right" style="background-color: rgba(229,239,249,0.75); color: #000000;">2</td>
<td headers="Skuterudbekken stub_1_28 Mean" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.021</td>
<td headers="Skuterudbekken stub_1_28 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.024</td>
<td headers="Skuterudbekken stub_1_28 Min" class="gt_row gt_right" style="background-color: rgba(240,246,253,0.75); color: #000000;">0.017</td>
<td headers="Skuterudbekken stub_1_28 StdDev" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.005</td>
<td headers="Skuterudbekken stub_1_28 P5" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.017</td>
<td headers="Skuterudbekken stub_1_28 P95" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.024</td></tr>
    <tr><th id="stub_1_29" scope="row" class="gt_row gt_left gt_stub">boscalid</th>
<td headers="Skuterudbekken stub_1_29 N" class="gt_row gt_right" style="background-color: rgba(211,228,243,0.75); color: #000000;">3</td>
<td headers="Skuterudbekken stub_1_29 Mean" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.016</td>
<td headers="Skuterudbekken stub_1_29 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.019</td>
<td headers="Skuterudbekken stub_1_29 Min" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.014</td>
<td headers="Skuterudbekken stub_1_29 StdDev" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.003</td>
<td headers="Skuterudbekken stub_1_29 P5" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.014</td>
<td headers="Skuterudbekken stub_1_29 P95" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.019</td></tr>
    <tr><th id="stub_1_30" scope="row" class="gt_row gt_left gt_stub">pyroxsulam</th>
<td headers="Skuterudbekken stub_1_30 N" class="gt_row gt_right" style="background-color: rgba(229,239,249,0.75); color: #000000;">2</td>
<td headers="Skuterudbekken stub_1_30 Mean" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.014</td>
<td headers="Skuterudbekken stub_1_30 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.017</td>
<td headers="Skuterudbekken stub_1_30 Min" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Skuterudbekken stub_1_30 StdDev" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.004</td>
<td headers="Skuterudbekken stub_1_30 P5" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Skuterudbekken stub_1_30 P95" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.017</td></tr>
    <tr><th id="stub_1_31" scope="row" class="gt_row gt_left gt_stub">carbendazim</th>
<td headers="Skuterudbekken stub_1_31 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Skuterudbekken stub_1_31 Mean" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.013</td>
<td headers="Skuterudbekken stub_1_31 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.013</td>
<td headers="Skuterudbekken stub_1_31 Min" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.013</td>
<td headers="Skuterudbekken stub_1_31 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Skuterudbekken stub_1_31 P5" class="gt_row gt_right" style="background-color: rgba(245,249,254,0.75); color: #000000;">0.013</td>
<td headers="Skuterudbekken stub_1_31 P95" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.013</td></tr>
    <tr><th id="stub_1_32" scope="row" class="gt_row gt_left gt_stub">mecoprop</th>
<td headers="Skuterudbekken stub_1_32 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Skuterudbekken stub_1_32 Mean" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.013</td>
<td headers="Skuterudbekken stub_1_32 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.013</td>
<td headers="Skuterudbekken stub_1_32 Min" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.013</td>
<td headers="Skuterudbekken stub_1_32 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Skuterudbekken stub_1_32 P5" class="gt_row gt_right" style="background-color: rgba(245,249,254,0.75); color: #000000;">0.013</td>
<td headers="Skuterudbekken stub_1_32 P95" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.013</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="8" class="gt_group_heading" scope="colgroup" id="Mørdrebekken">Mørdrebekken</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_33" scope="row" class="gt_row gt_left gt_stub">2-methyl-4-chlorophenoxyacetic acid (MCPA)</th>
<td headers="Mørdrebekken stub_1_33 N" class="gt_row gt_right" style="background-color: rgba(162,204,226,0.75); color: #000000;">5</td>
<td headers="Mørdrebekken stub_1_33 Mean" class="gt_row gt_right" style="background-color: rgba(187,214,235,0.75); color: #000000;">0.332</td>
<td headers="Mørdrebekken stub_1_33 Max" class="gt_row gt_right" style="background-color: rgba(199,220,239,0.75); color: #000000;">1.300</td>
<td headers="Mørdrebekken stub_1_33 Min" class="gt_row gt_right" style="background-color: rgba(211,228,243,0.75); color: #000000;">0.046</td>
<td headers="Mørdrebekken stub_1_33 StdDev" class="gt_row gt_right" style="background-color: rgba(151,198,223,0.75); color: #000000;">0.546</td>
<td headers="Mørdrebekken stub_1_33 P5" class="gt_row gt_right" style="background-color: rgba(214,230,244,0.75); color: #000000;">0.047</td>
<td headers="Mørdrebekken stub_1_33 P95" class="gt_row gt_right" style="background-color: rgba(167,206,228,0.75); color: #000000;">1.082</td></tr>
    <tr><th id="stub_1_34" scope="row" class="gt_row gt_left gt_stub">fluroxypyr</th>
<td headers="Mørdrebekken stub_1_34 N" class="gt_row gt_right" style="background-color: rgba(229,239,249,0.75); color: #000000;">2</td>
<td headers="Mørdrebekken stub_1_34 Mean" class="gt_row gt_right" style="background-color: rgba(205,224,241,0.75); color: #000000;">0.250</td>
<td headers="Mørdrebekken stub_1_34 Max" class="gt_row gt_right" style="background-color: rgba(233,242,250,0.75); color: #000000;">0.390</td>
<td headers="Mørdrebekken stub_1_34 Min" class="gt_row gt_right" style="background-color: rgba(107,174,214,0.75); color: #000000;">0.110</td>
<td headers="Mørdrebekken stub_1_34 StdDev" class="gt_row gt_right" style="background-color: rgba(219,233,246,0.75); color: #000000;">0.198</td>
<td headers="Mørdrebekken stub_1_34 P5" class="gt_row gt_right" style="background-color: rgba(99,168,211,0.75); color: #FFFFFF;">0.124</td>
<td headers="Mørdrebekken stub_1_34 P95" class="gt_row gt_right" style="background-color: rgba(223,236,247,0.75); color: #000000;">0.376</td></tr>
    <tr><th id="stub_1_35" scope="row" class="gt_row gt_left gt_stub">bentazone</th>
<td headers="Mørdrebekken stub_1_35 N" class="gt_row gt_right" style="background-color: rgba(191,216,236,0.75); color: #000000;">4</td>
<td headers="Mørdrebekken stub_1_35 Mean" class="gt_row gt_right" style="background-color: rgba(213,229,244,0.75); color: #000000;">0.205</td>
<td headers="Mørdrebekken stub_1_35 Max" class="gt_row gt_right" style="background-color: rgba(231,241,250,0.75); color: #000000;">0.430</td>
<td headers="Mørdrebekken stub_1_35 Min" class="gt_row gt_right" style="background-color: rgba(223,236,247,0.75); color: #000000;">0.034</td>
<td headers="Mørdrebekken stub_1_35 StdDev" class="gt_row gt_right" style="background-color: rgba(218,233,246,0.75); color: #000000;">0.201</td>
<td headers="Mørdrebekken stub_1_35 P5" class="gt_row gt_right" style="background-color: rgba(225,237,248,0.75); color: #000000;">0.035</td>
<td headers="Mørdrebekken stub_1_35 P95" class="gt_row gt_right" style="background-color: rgba(221,234,247,0.75); color: #000000;">0.413</td></tr>
    <tr><th id="stub_1_36" scope="row" class="gt_row gt_left gt_stub">clopyralid</th>
<td headers="Mørdrebekken stub_1_36 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Mørdrebekken stub_1_36 Mean" class="gt_row gt_right" style="background-color: rgba(214,230,244,0.75); color: #000000;">0.200</td>
<td headers="Mørdrebekken stub_1_36 Max" class="gt_row gt_right" style="background-color: rgba(240,246,253,0.75); color: #000000;">0.200</td>
<td headers="Mørdrebekken stub_1_36 Min" class="gt_row gt_right" style="background-color: rgba(9,61,126,0.75); color: #FFFFFF;">0.200</td>
<td headers="Mørdrebekken stub_1_36 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Mørdrebekken stub_1_36 P5" class="gt_row gt_right" style="background-color: rgba(8,80,155,0.75); color: #FFFFFF;">0.200</td>
<td headers="Mørdrebekken stub_1_36 P95" class="gt_row gt_right" style="background-color: rgba(235,243,251,0.75); color: #000000;">0.200</td></tr>
    <tr><th id="stub_1_37" scope="row" class="gt_row gt_left gt_stub">propiconazole</th>
<td headers="Mørdrebekken stub_1_37 N" class="gt_row gt_right" style="background-color: rgba(191,216,236,0.75); color: #000000;">4</td>
<td headers="Mørdrebekken stub_1_37 Mean" class="gt_row gt_right" style="background-color: rgba(220,234,246,0.75); color: #000000;">0.163</td>
<td headers="Mørdrebekken stub_1_37 Max" class="gt_row gt_right" style="background-color: rgba(225,237,248,0.75); color: #000000;">0.580</td>
<td headers="Mørdrebekken stub_1_37 Min" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Mørdrebekken stub_1_37 StdDev" class="gt_row gt_right" style="background-color: rgba(208,225,242,0.75); color: #000000;">0.278</td>
<td headers="Mørdrebekken stub_1_37 P5" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.013</td>
<td headers="Mørdrebekken stub_1_37 P95" class="gt_row gt_right" style="background-color: rgba(216,231,245,0.75); color: #000000;">0.498</td></tr>
    <tr><th id="stub_1_38" scope="row" class="gt_row gt_left gt_stub">pencycuron</th>
<td headers="Mørdrebekken stub_1_38 N" class="gt_row gt_right" style="background-color: rgba(229,239,249,0.75); color: #000000;">2</td>
<td headers="Mørdrebekken stub_1_38 Mean" class="gt_row gt_right" style="background-color: rgba(237,244,252,0.75); color: #000000;">0.070</td>
<td headers="Mørdrebekken stub_1_38 Max" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.100</td>
<td headers="Mørdrebekken stub_1_38 Min" class="gt_row gt_right" style="background-color: rgba(217,232,245,0.75); color: #000000;">0.040</td>
<td headers="Mørdrebekken stub_1_38 StdDev" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.042</td>
<td headers="Mørdrebekken stub_1_38 P5" class="gt_row gt_right" style="background-color: rgba(217,232,245,0.75); color: #000000;">0.043</td>
<td headers="Mørdrebekken stub_1_38 P95" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.097</td></tr>
    <tr><th id="stub_1_39" scope="row" class="gt_row gt_left gt_stub">fenpropimorph</th>
<td headers="Mørdrebekken stub_1_39 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Mørdrebekken stub_1_39 Mean" class="gt_row gt_right" style="background-color: rgba(240,246,253,0.75); color: #000000;">0.052</td>
<td headers="Mørdrebekken stub_1_39 Max" class="gt_row gt_right" style="background-color: rgba(245,250,255,0.75); color: #000000;">0.052</td>
<td headers="Mørdrebekken stub_1_39 Min" class="gt_row gt_right" style="background-color: rgba(206,224,242,0.75); color: #000000;">0.052</td>
<td headers="Mørdrebekken stub_1_39 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Mørdrebekken stub_1_39 P5" class="gt_row gt_right" style="background-color: rgba(209,226,243,0.75); color: #000000;">0.052</td>
<td headers="Mørdrebekken stub_1_39 P95" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.052</td></tr>
    <tr><th id="stub_1_40" scope="row" class="gt_row gt_left gt_stub">thiabendazole</th>
<td headers="Mørdrebekken stub_1_40 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Mørdrebekken stub_1_40 Mean" class="gt_row gt_right" style="background-color: rgba(242,248,253,0.75); color: #000000;">0.041</td>
<td headers="Mørdrebekken stub_1_40 Max" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.041</td>
<td headers="Mørdrebekken stub_1_40 Min" class="gt_row gt_right" style="background-color: rgba(216,231,245,0.75); color: #000000;">0.041</td>
<td headers="Mørdrebekken stub_1_40 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Mørdrebekken stub_1_40 P5" class="gt_row gt_right" style="background-color: rgba(219,233,246,0.75); color: #000000;">0.041</td>
<td headers="Mørdrebekken stub_1_40 P95" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.041</td></tr>
    <tr><th id="stub_1_41" scope="row" class="gt_row gt_left gt_stub">prothioconazole-desthio</th>
<td headers="Mørdrebekken stub_1_41 N" class="gt_row gt_right" style="background-color: rgba(191,216,236,0.75); color: #000000;">4</td>
<td headers="Mørdrebekken stub_1_41 Mean" class="gt_row gt_right" style="background-color: rgba(243,249,254,0.75); color: #000000;">0.033</td>
<td headers="Mørdrebekken stub_1_41 Max" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.091</td>
<td headers="Mørdrebekken stub_1_41 Min" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Mørdrebekken stub_1_41 StdDev" class="gt_row gt_right" style="background-color: rgba(242,248,253,0.75); color: #000000;">0.039</td>
<td headers="Mørdrebekken stub_1_41 P5" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Mørdrebekken stub_1_41 P95" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.080</td></tr>
    <tr><th id="stub_1_42" scope="row" class="gt_row gt_left gt_stub">pinoxaden</th>
<td headers="Mørdrebekken stub_1_42 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Mørdrebekken stub_1_42 Mean" class="gt_row gt_right" style="background-color: rgba(243,249,254,0.75); color: #000000;">0.031</td>
<td headers="Mørdrebekken stub_1_42 Max" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.031</td>
<td headers="Mørdrebekken stub_1_42 Min" class="gt_row gt_right" style="background-color: rgba(226,238,248,0.75); color: #000000;">0.031</td>
<td headers="Mørdrebekken stub_1_42 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Mørdrebekken stub_1_42 P5" class="gt_row gt_right" style="background-color: rgba(228,239,249,0.75); color: #000000;">0.031</td>
<td headers="Mørdrebekken stub_1_42 P95" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.031</td></tr>
    <tr><th id="stub_1_43" scope="row" class="gt_row gt_left gt_stub">azoxystrobin</th>
<td headers="Mørdrebekken stub_1_43 N" class="gt_row gt_right" style="background-color: rgba(191,216,236,0.75); color: #000000;">4</td>
<td headers="Mørdrebekken stub_1_43 Mean" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.021</td>
<td headers="Mørdrebekken stub_1_43 Max" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.027</td>
<td headers="Mørdrebekken stub_1_43 Min" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.012</td>
<td headers="Mørdrebekken stub_1_43 StdDev" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.007</td>
<td headers="Mørdrebekken stub_1_43 P5" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.013</td>
<td headers="Mørdrebekken stub_1_43 P95" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.027</td></tr>
    <tr><th id="stub_1_44" scope="row" class="gt_row gt_left gt_stub">bixafen</th>
<td headers="Mørdrebekken stub_1_44 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Mørdrebekken stub_1_44 Mean" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.015</td>
<td headers="Mørdrebekken stub_1_44 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.015</td>
<td headers="Mørdrebekken stub_1_44 Min" class="gt_row gt_right" style="background-color: rgba(242,248,253,0.75); color: #000000;">0.015</td>
<td headers="Mørdrebekken stub_1_44 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Mørdrebekken stub_1_44 P5" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.015</td>
<td headers="Mørdrebekken stub_1_44 P95" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.015</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="8" class="gt_group_heading" scope="colgroup" id="Timebekken">Timebekken</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_45" scope="row" class="gt_row gt_left gt_stub">metribuzin</th>
<td headers="Timebekken stub_1_45 N" class="gt_row gt_right" style="background-color: rgba(229,239,249,0.75); color: #000000;">2</td>
<td headers="Timebekken stub_1_45 Mean" class="gt_row gt_right" style="background-color: rgba(196,218,238,0.75); color: #000000;">0.302</td>
<td headers="Timebekken stub_1_45 Max" class="gt_row gt_right" style="background-color: rgba(225,237,248,0.75); color: #000000;">0.590</td>
<td headers="Timebekken stub_1_45 Min" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.014</td>
<td headers="Timebekken stub_1_45 StdDev" class="gt_row gt_right" style="background-color: rgba(185,213,234,0.75); color: #000000;">0.407</td>
<td headers="Timebekken stub_1_45 P5" class="gt_row gt_right" style="background-color: rgba(217,232,245,0.75); color: #000000;">0.043</td>
<td headers="Timebekken stub_1_45 P95" class="gt_row gt_right" style="background-color: rgba(212,228,244,0.75); color: #000000;">0.561</td></tr>
    <tr><th id="stub_1_46" scope="row" class="gt_row gt_left gt_stub">fluroxypyr</th>
<td headers="Timebekken stub_1_46 N" class="gt_row gt_right" style="background-color: rgba(191,216,236,0.75); color: #000000;">4</td>
<td headers="Timebekken stub_1_46 Mean" class="gt_row gt_right" style="background-color: rgba(201,221,240,0.75); color: #000000;">0.275</td>
<td headers="Timebekken stub_1_46 Max" class="gt_row gt_right" style="background-color: rgba(230,240,250,0.75); color: #000000;">0.460</td>
<td headers="Timebekken stub_1_46 Min" class="gt_row gt_right" style="background-color: rgba(48,126,188,0.75); color: #FFFFFF;">0.150</td>
<td headers="Timebekken stub_1_46 StdDev" class="gt_row gt_right" style="background-color: rgba(226,237,248,0.75); color: #000000;">0.148</td>
<td headers="Timebekken stub_1_46 P5" class="gt_row gt_right" style="background-color: rgba(60,138,194,0.75); color: #FFFFFF;">0.151</td>
<td headers="Timebekken stub_1_46 P95" class="gt_row gt_right" style="background-color: rgba(219,233,246,0.75); color: #000000;">0.440</td></tr>
    <tr><th id="stub_1_47" scope="row" class="gt_row gt_left gt_stub">florasulam</th>
<td headers="Timebekken stub_1_47 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Timebekken stub_1_47 Mean" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.023</td>
<td headers="Timebekken stub_1_47 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.023</td>
<td headers="Timebekken stub_1_47 Min" class="gt_row gt_right" style="background-color: rgba(234,243,251,0.75); color: #000000;">0.023</td>
<td headers="Timebekken stub_1_47 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Timebekken stub_1_47 P5" class="gt_row gt_right" style="background-color: rgba(235,243,251,0.75); color: #000000;">0.023</td>
<td headers="Timebekken stub_1_47 P95" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.023</td></tr>
    <tr><th id="stub_1_48" scope="row" class="gt_row gt_left gt_stub">propiconazole</th>
<td headers="Timebekken stub_1_48 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Timebekken stub_1_48 Mean" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.016</td>
<td headers="Timebekken stub_1_48 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.016</td>
<td headers="Timebekken stub_1_48 Min" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.016</td>
<td headers="Timebekken stub_1_48 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Timebekken stub_1_48 P5" class="gt_row gt_right" style="background-color: rgba(242,248,253,0.75); color: #000000;">0.016</td>
<td headers="Timebekken stub_1_48 P95" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.016</td></tr>
    <tr><th id="stub_1_49" scope="row" class="gt_row gt_left gt_stub">2-methyl-4-chlorophenoxyacetic acid (MCPA)</th>
<td headers="Timebekken stub_1_49 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Timebekken stub_1_49 Mean" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.011</td>
<td headers="Timebekken stub_1_49 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.011</td>
<td headers="Timebekken stub_1_49 Min" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Timebekken stub_1_49 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Timebekken stub_1_49 P5" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.011</td>
<td headers="Timebekken stub_1_49 P95" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.011</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="8" class="gt_group_heading" scope="colgroup" id="Vasshaglona">Vasshaglona</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_50" scope="row" class="gt_row gt_left gt_stub">fenamidone</th>
<td headers="Vasshaglona stub_1_50 N" class="gt_row gt_right" style="background-color: rgba(229,239,249,0.75); color: #000000;">2</td>
<td headers="Vasshaglona stub_1_50 Mean" class="gt_row gt_right" style="background-color: rgba(212,229,244,0.75); color: #000000;">0.208</td>
<td headers="Vasshaglona stub_1_50 Max" class="gt_row gt_right" style="background-color: rgba(234,243,251,0.75); color: #000000;">0.360</td>
<td headers="Vasshaglona stub_1_50 Min" class="gt_row gt_right" style="background-color: rgba(201,221,240,0.75); color: #000000;">0.057</td>
<td headers="Vasshaglona stub_1_50 StdDev" class="gt_row gt_right" style="background-color: rgba(217,231,245,0.75); color: #000000;">0.214</td>
<td headers="Vasshaglona stub_1_50 P5" class="gt_row gt_right" style="background-color: rgba(187,214,235,0.75); color: #000000;">0.072</td>
<td headers="Vasshaglona stub_1_50 P95" class="gt_row gt_right" style="background-color: rgba(225,237,248,0.75); color: #000000;">0.345</td></tr>
    <tr><th id="stub_1_51" scope="row" class="gt_row gt_left gt_stub">propamocarb</th>
<td headers="Vasshaglona stub_1_51 N" class="gt_row gt_right" style="background-color: rgba(162,204,226,0.75); color: #000000;">5</td>
<td headers="Vasshaglona stub_1_51 Mean" class="gt_row gt_right" style="background-color: rgba(231,241,250,0.75); color: #000000;">0.102</td>
<td headers="Vasshaglona stub_1_51 Max" class="gt_row gt_right" style="background-color: rgba(233,242,251,0.75); color: #000000;">0.370</td>
<td headers="Vasshaglona stub_1_51 Min" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.011</td>
<td headers="Vasshaglona stub_1_51 StdDev" class="gt_row gt_right" style="background-color: rgba(225,237,248,0.75); color: #000000;">0.153</td>
<td headers="Vasshaglona stub_1_51 P5" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.012</td>
<td headers="Vasshaglona stub_1_51 P95" class="gt_row gt_right" style="background-color: rgba(227,238,249,0.75); color: #000000;">0.314</td></tr>
    <tr><th id="stub_1_52" scope="row" class="gt_row gt_left gt_stub">pyridafol</th>
<td headers="Vasshaglona stub_1_52 N" class="gt_row gt_right" style="background-color: rgba(162,204,226,0.75); color: #000000;">5</td>
<td headers="Vasshaglona stub_1_52 Mean" class="gt_row gt_right" style="background-color: rgba(234,243,251,0.75); color: #000000;">0.084</td>
<td headers="Vasshaglona stub_1_52 Max" class="gt_row gt_right" style="background-color: rgba(238,245,252,0.75); color: #000000;">0.240</td>
<td headers="Vasshaglona stub_1_52 Min" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.010</td>
<td headers="Vasshaglona stub_1_52 StdDev" class="gt_row gt_right" style="background-color: rgba(233,242,250,0.75); color: #000000;">0.100</td>
<td headers="Vasshaglona stub_1_52 P5" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.012</td>
<td headers="Vasshaglona stub_1_52 P95" class="gt_row gt_right" style="background-color: rgba(234,242,251,0.75); color: #000000;">0.218</td></tr>
    <tr><th id="stub_1_53" scope="row" class="gt_row gt_left gt_stub">metribuzin</th>
<td headers="Vasshaglona stub_1_53 N" class="gt_row gt_right" style="background-color: rgba(64,143,196,0.75); color: #FFFFFF;">8</td>
<td headers="Vasshaglona stub_1_53 Mean" class="gt_row gt_right" style="background-color: rgba(235,243,251,0.75); color: #000000;">0.079</td>
<td headers="Vasshaglona stub_1_53 Max" class="gt_row gt_right" style="background-color: rgba(240,247,253,0.75); color: #000000;">0.190</td>
<td headers="Vasshaglona stub_1_53 Min" class="gt_row gt_right" style="background-color: rgba(234,243,251,0.75); color: #000000;">0.023</td>
<td headers="Vasshaglona stub_1_53 StdDev" class="gt_row gt_right" style="background-color: rgba(238,245,252,0.75); color: #000000;">0.068</td>
<td headers="Vasshaglona stub_1_53 P5" class="gt_row gt_right" style="background-color: rgba(235,243,251,0.75); color: #000000;">0.023</td>
<td headers="Vasshaglona stub_1_53 P95" class="gt_row gt_right" style="background-color: rgba(236,244,252,0.75); color: #000000;">0.179</td></tr>
    <tr><th id="stub_1_54" scope="row" class="gt_row gt_left gt_stub">thiacloprid</th>
<td headers="Vasshaglona stub_1_54 N" class="gt_row gt_right" style="background-color: rgba(229,239,249,0.75); color: #000000;">2</td>
<td headers="Vasshaglona stub_1_54 Mean" class="gt_row gt_right" style="background-color: rgba(237,244,252,0.75); color: #000000;">0.069</td>
<td headers="Vasshaglona stub_1_54 Max" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.120</td>
<td headers="Vasshaglona stub_1_54 Min" class="gt_row gt_right" style="background-color: rgba(238,245,252,0.75); color: #000000;">0.019</td>
<td headers="Vasshaglona stub_1_54 StdDev" class="gt_row gt_right" style="background-color: rgba(237,245,252,0.75); color: #000000;">0.071</td>
<td headers="Vasshaglona stub_1_54 P5" class="gt_row gt_right" style="background-color: rgba(234,243,251,0.75); color: #000000;">0.024</td>
<td headers="Vasshaglona stub_1_54 P95" class="gt_row gt_right" style="background-color: rgba(240,247,253,0.75); color: #000000;">0.115</td></tr>
    <tr><th id="stub_1_55" scope="row" class="gt_row gt_left gt_stub">fluroxypyr</th>
<td headers="Vasshaglona stub_1_55 N" class="gt_row gt_right" style="background-color: rgba(211,228,243,0.75); color: #000000;">3</td>
<td headers="Vasshaglona stub_1_55 Mean" class="gt_row gt_right" style="background-color: rgba(237,245,252,0.75); color: #000000;">0.067</td>
<td headers="Vasshaglona stub_1_55 Max" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.086</td>
<td headers="Vasshaglona stub_1_55 Min" class="gt_row gt_right" style="background-color: rgba(202,222,240,0.75); color: #000000;">0.056</td>
<td headers="Vasshaglona stub_1_55 StdDev" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.017</td>
<td headers="Vasshaglona stub_1_55 P5" class="gt_row gt_right" style="background-color: rgba(205,224,241,0.75); color: #000000;">0.056</td>
<td headers="Vasshaglona stub_1_55 P95" class="gt_row gt_right" style="background-color: rgba(242,248,254,0.75); color: #000000;">0.083</td></tr>
    <tr><th id="stub_1_56" scope="row" class="gt_row gt_left gt_stub">cyazofamid</th>
<td headers="Vasshaglona stub_1_56 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Vasshaglona stub_1_56 Mean" class="gt_row gt_right" style="background-color: rgba(240,246,253,0.75); color: #000000;">0.052</td>
<td headers="Vasshaglona stub_1_56 Max" class="gt_row gt_right" style="background-color: rgba(245,250,255,0.75); color: #000000;">0.052</td>
<td headers="Vasshaglona stub_1_56 Min" class="gt_row gt_right" style="background-color: rgba(206,224,242,0.75); color: #000000;">0.052</td>
<td headers="Vasshaglona stub_1_56 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Vasshaglona stub_1_56 P5" class="gt_row gt_right" style="background-color: rgba(209,226,243,0.75); color: #000000;">0.052</td>
<td headers="Vasshaglona stub_1_56 P95" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.052</td></tr>
    <tr><th id="stub_1_57" scope="row" class="gt_row gt_left gt_stub">2-methyl-4-chlorophenoxyacetic acid (MCPA)</th>
<td headers="Vasshaglona stub_1_57 N" class="gt_row gt_right" style="background-color: rgba(191,216,236,0.75); color: #000000;">4</td>
<td headers="Vasshaglona stub_1_57 Mean" class="gt_row gt_right" style="background-color: rgba(240,247,253,0.75); color: #000000;">0.050</td>
<td headers="Vasshaglona stub_1_57 Max" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.160</td>
<td headers="Vasshaglona stub_1_57 Min" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.012</td>
<td headers="Vasshaglona stub_1_57 StdDev" class="gt_row gt_right" style="background-color: rgba(237,244,252,0.75); color: #000000;">0.073</td>
<td headers="Vasshaglona stub_1_57 P5" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.012</td>
<td headers="Vasshaglona stub_1_57 P95" class="gt_row gt_right" style="background-color: rgba(239,246,252,0.75); color: #000000;">0.138</td></tr>
    <tr><th id="stub_1_58" scope="row" class="gt_row gt_left gt_stub">aclonifen</th>
<td headers="Vasshaglona stub_1_58 N" class="gt_row gt_right" style="background-color: rgba(191,216,236,0.75); color: #000000;">4</td>
<td headers="Vasshaglona stub_1_58 Mean" class="gt_row gt_right" style="background-color: rgba(240,247,253,0.75); color: #000000;">0.049</td>
<td headers="Vasshaglona stub_1_58 Max" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.130</td>
<td headers="Vasshaglona stub_1_58 Min" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.014</td>
<td headers="Vasshaglona stub_1_58 StdDev" class="gt_row gt_right" style="background-color: rgba(240,246,253,0.75); color: #000000;">0.054</td>
<td headers="Vasshaglona stub_1_58 P5" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.015</td>
<td headers="Vasshaglona stub_1_58 P95" class="gt_row gt_right" style="background-color: rgba(240,247,253,0.75); color: #000000;">0.115</td></tr>
    <tr><th id="stub_1_59" scope="row" class="gt_row gt_left gt_stub">boscalid</th>
<td headers="Vasshaglona stub_1_59 N" class="gt_row gt_right" style="background-color: rgba(8,48,107,0.75); color: #FFFFFF;">12</td>
<td headers="Vasshaglona stub_1_59 Mean" class="gt_row gt_right" style="background-color: rgba(240,247,253,0.75); color: #000000;">0.048</td>
<td headers="Vasshaglona stub_1_59 Max" class="gt_row gt_right" style="background-color: rgba(242,248,253,0.75); color: #000000;">0.150</td>
<td headers="Vasshaglona stub_1_59 Min" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.010</td>
<td headers="Vasshaglona stub_1_59 StdDev" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.043</td>
<td headers="Vasshaglona stub_1_59 P5" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.013</td>
<td headers="Vasshaglona stub_1_59 P95" class="gt_row gt_right" style="background-color: rgba(239,246,252,0.75); color: #000000;">0.133</td></tr>
    <tr><th id="stub_1_60" scope="row" class="gt_row gt_left gt_stub">pencycuron</th>
<td headers="Vasshaglona stub_1_60 N" class="gt_row gt_right" style="background-color: rgba(93,164,208,0.75); color: #FFFFFF;">7</td>
<td headers="Vasshaglona stub_1_60 Mean" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.045</td>
<td headers="Vasshaglona stub_1_60 Max" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.170</td>
<td headers="Vasshaglona stub_1_60 Min" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.014</td>
<td headers="Vasshaglona stub_1_60 StdDev" class="gt_row gt_right" style="background-color: rgba(239,246,253,0.75); color: #000000;">0.056</td>
<td headers="Vasshaglona stub_1_60 P5" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.015</td>
<td headers="Vasshaglona stub_1_60 P95" class="gt_row gt_right" style="background-color: rgba(239,246,253,0.75); color: #000000;">0.131</td></tr>
    <tr><th id="stub_1_61" scope="row" class="gt_row gt_left gt_stub">mandipropamid</th>
<td headers="Vasshaglona stub_1_61 N" class="gt_row gt_right" style="background-color: rgba(127,184,218,0.75); color: #000000;">6</td>
<td headers="Vasshaglona stub_1_61 Mean" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.043</td>
<td headers="Vasshaglona stub_1_61 Max" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.130</td>
<td headers="Vasshaglona stub_1_61 Min" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.010</td>
<td headers="Vasshaglona stub_1_61 StdDev" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.045</td>
<td headers="Vasshaglona stub_1_61 P5" class="gt_row gt_right" style="background-color: rgba(245,250,255,0.75); color: #000000;">0.012</td>
<td headers="Vasshaglona stub_1_61 P95" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.110</td></tr>
    <tr><th id="stub_1_62" scope="row" class="gt_row gt_left gt_stub">spinosad</th>
<td headers="Vasshaglona stub_1_62 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Vasshaglona stub_1_62 Mean" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.034</td>
<td headers="Vasshaglona stub_1_62 Max" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.034</td>
<td headers="Vasshaglona stub_1_62 Min" class="gt_row gt_right" style="background-color: rgba(223,236,247,0.75); color: #000000;">0.034</td>
<td headers="Vasshaglona stub_1_62 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Vasshaglona stub_1_62 P5" class="gt_row gt_right" style="background-color: rgba(225,237,248,0.75); color: #000000;">0.034</td>
<td headers="Vasshaglona stub_1_62 P95" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.034</td></tr>
    <tr><th id="stub_1_63" scope="row" class="gt_row gt_left gt_stub">cyprodinil</th>
<td headers="Vasshaglona stub_1_63 N" class="gt_row gt_right" style="background-color: rgba(211,228,243,0.75); color: #000000;">3</td>
<td headers="Vasshaglona stub_1_63 Mean" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.033</td>
<td headers="Vasshaglona stub_1_63 Max" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.071</td>
<td headers="Vasshaglona stub_1_63 Min" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.012</td>
<td headers="Vasshaglona stub_1_63 StdDev" class="gt_row gt_right" style="background-color: rgba(243,248,254,0.75); color: #000000;">0.033</td>
<td headers="Vasshaglona stub_1_63 P5" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.013</td>
<td headers="Vasshaglona stub_1_63 P95" class="gt_row gt_right" style="background-color: rgba(243,249,254,0.75); color: #000000;">0.066</td></tr>
    <tr><th id="stub_1_64" scope="row" class="gt_row gt_left gt_stub">clomazone</th>
<td headers="Vasshaglona stub_1_64 N" class="gt_row gt_right" style="background-color: rgba(211,228,243,0.75); color: #000000;">3</td>
<td headers="Vasshaglona stub_1_64 Mean" class="gt_row gt_right" style="background-color: rgba(243,249,254,0.75); color: #000000;">0.032</td>
<td headers="Vasshaglona stub_1_64 Max" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.063</td>
<td headers="Vasshaglona stub_1_64 Min" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.013</td>
<td headers="Vasshaglona stub_1_64 StdDev" class="gt_row gt_right" style="background-color: rgba(243,249,254,0.75); color: #000000;">0.027</td>
<td headers="Vasshaglona stub_1_64 P5" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.014</td>
<td headers="Vasshaglona stub_1_64 P95" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.059</td></tr>
    <tr><th id="stub_1_65" scope="row" class="gt_row gt_left gt_stub">bentazone</th>
<td headers="Vasshaglona stub_1_65 N" class="gt_row gt_right" style="background-color: rgba(40,119,184,0.75); color: #FFFFFF;">9</td>
<td headers="Vasshaglona stub_1_65 Mean" class="gt_row gt_right" style="background-color: rgba(244,249,254,0.75); color: #000000;">0.026</td>
<td headers="Vasshaglona stub_1_65 Max" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.045</td>
<td headers="Vasshaglona stub_1_65 Min" class="gt_row gt_right" style="background-color: rgba(239,246,252,0.75); color: #000000;">0.018</td>
<td headers="Vasshaglona stub_1_65 StdDev" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.010</td>
<td headers="Vasshaglona stub_1_65 P5" class="gt_row gt_right" style="background-color: rgba(240,246,253,0.75); color: #000000;">0.018</td>
<td headers="Vasshaglona stub_1_65 P95" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.043</td></tr>
    <tr><th id="stub_1_66" scope="row" class="gt_row gt_left gt_stub">fludioxonil</th>
<td headers="Vasshaglona stub_1_66 N" class="gt_row gt_right" style="background-color: rgba(211,228,243,0.75); color: #000000;">3</td>
<td headers="Vasshaglona stub_1_66 Mean" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.021</td>
<td headers="Vasshaglona stub_1_66 Max" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.037</td>
<td headers="Vasshaglona stub_1_66 Min" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.012</td>
<td headers="Vasshaglona stub_1_66 StdDev" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.014</td>
<td headers="Vasshaglona stub_1_66 P5" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.012</td>
<td headers="Vasshaglona stub_1_66 P95" class="gt_row gt_right" style="background-color: rgba(245,250,255,0.75); color: #000000;">0.035</td></tr>
    <tr><th id="stub_1_67" scope="row" class="gt_row gt_left gt_stub">imidacloprid</th>
<td headers="Vasshaglona stub_1_67 N" class="gt_row gt_right" style="background-color: rgba(211,228,243,0.75); color: #000000;">3</td>
<td headers="Vasshaglona stub_1_67 Mean" class="gt_row gt_right" style="background-color: rgba(245,250,254,0.75); color: #000000;">0.020</td>
<td headers="Vasshaglona stub_1_67 Max" class="gt_row gt_right" style="background-color: rgba(246,251,255,0.75); color: #000000;">0.027</td>
<td headers="Vasshaglona stub_1_67 Min" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.016</td>
<td headers="Vasshaglona stub_1_67 StdDev" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.006</td>
<td headers="Vasshaglona stub_1_67 P5" class="gt_row gt_right" style="background-color: rgba(242,248,253,0.75); color: #000000;">0.016</td>
<td headers="Vasshaglona stub_1_67 P95" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.026</td></tr>
    <tr><th id="stub_1_68" scope="row" class="gt_row gt_left gt_stub">chlorfenvinphos</th>
<td headers="Vasshaglona stub_1_68 N" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">1</td>
<td headers="Vasshaglona stub_1_68 Mean" class="gt_row gt_right" style="background-color: rgba(246,250,255,0.75); color: #000000;">0.017</td>
<td headers="Vasshaglona stub_1_68 Max" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.017</td>
<td headers="Vasshaglona stub_1_68 Min" class="gt_row gt_right" style="background-color: rgba(240,246,253,0.75); color: #000000;">0.017</td>
<td headers="Vasshaglona stub_1_68 StdDev" class="gt_row gt_right" style="background-color: rgba(128,128,128,0.75); color: #FFFFFF;">—</td>
<td headers="Vasshaglona stub_1_68 P5" class="gt_row gt_right" style="background-color: rgba(241,247,253,0.75); color: #000000;">0.017</td>
<td headers="Vasshaglona stub_1_68 P95" class="gt_row gt_right" style="background-color: rgba(247,251,255,0.75); color: #000000;">0.017</td></tr>
  </tbody>
  
  
</table>
</div>

## Tissue concentration prediction

Here we will predict the tissue concentration of the chemical in fish
(Cf) on the basis of the mean water concentration (Cw),

<details>
<summary>Calculate fish concentration level</summary>

``` r
get_tissue_conc <- function(cw, logp) {
  cf <- cw * 10 ^(0.76 * logp - 0.23)
  return(cf)
}
exposure_data <- exposure_data %>% 
  mutate(FishConcentration = get_tissue_conc(MEASURED_VALUE, XLogP))
```

</details>
<details>
<summary>Average Fish and water concentration</summary>

``` r
summary_comparison <- exposure_data %>%
  group_by(SITE_NAME, STRESSOR_NAME) %>% 
  summarize(
    across(
      c(MEASURED_VALUE, FishConcentration),
      function(x) {
        Hmisc::smean.cl.normal(x) %>%
          as.list() %>% 
          bind_cols() %>% 
          list()
      }
    )
  ) %>% unnest(names_sep = "_") %>% 
  tidytable::rename_with(
    ~gsub("MEASURED_VALUE", "Water concentration (ug/L)", .x)
  ) %>% 
  tidytable::rename_with(
    ~gsub("FishConcentration", "Fish concentration (ug/Kg)", .x)
  ) %>% 
  arrange(
    desc(`Water concentration (ug/L)_Mean`), 
    desc(`Fish concentration (ug/Kg)_Mean`)
  )

gt::gt(
  data = summary_comparison, 
  rowname_col = "STRESSOR_NAME",
  groupname_col = "SITE_NAME"
) %>% gt::tab_spanner_delim("_", split = "last", limit = 1) %>% 
  gt::fmt_number() %>% 
  gt::sub_missing() %>% 
  gt::cols_merge(
    columns = gt::starts_with("Water") & !gt::ends_with("Mean"),
    pattern = "({1},{2})"
  ) %>% 
  gt::cols_merge(
    columns = gt::starts_with("Fish") & !gt::ends_with("Mean"),
    pattern = "({1},{2})"
  ) %>% 
  gt::cols_merge(
    columns = 3:4,
    pattern = "{1} {2}"
  ) %>% 
  gt::cols_merge(
    columns = 6:7,
    pattern = "{1} {2}"
  ) %>% 
  gt::cols_label_with(
    fn = \(x) gsub("Mean", "Mean (95% CI)", x)
  ) %>% 
  gt::opt_vertical_padding(0.25) %>%
  gt::tab_style(
    style = gt::cell_text(align = "center"),
    locations = gt::cells_column_labels()
  ) %>% 
  gt::tab_options(
    table.width = "100%",
    table.font.size = "10pt",
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold"
  ) %>% 
  gt::data_color(
    columns = -c(1:2), 
    direction = "column",
    fn = scales::col_numeric(
      palette = "Blues",
      domain = NULL,
      alpha = 0.5
    )
  )
```

</details>
<div id="rvzeorrzax" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#rvzeorrzax table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#rvzeorrzax thead, #rvzeorrzax tbody, #rvzeorrzax tfoot, #rvzeorrzax tr, #rvzeorrzax td, #rvzeorrzax th {
  border-style: none;
}

#rvzeorrzax p {
  margin: 0;
  padding: 0;
}

#rvzeorrzax .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 10pt;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 100%;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#rvzeorrzax .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#rvzeorrzax .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 1px;
  padding-bottom: 1px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#rvzeorrzax .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0px;
  padding-bottom: 2px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#rvzeorrzax .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rvzeorrzax .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rvzeorrzax .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rvzeorrzax .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: bold;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 1.25px;
  padding-bottom: 2.25px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#rvzeorrzax .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: bold;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#rvzeorrzax .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#rvzeorrzax .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#rvzeorrzax .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 1.25px;
  padding-bottom: 1.25px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#rvzeorrzax .gt_spanner_row {
  border-bottom-style: hidden;
}

#rvzeorrzax .gt_group_heading {
  padding-top: 2px;
  padding-bottom: 2px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: bold;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#rvzeorrzax .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: bold;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#rvzeorrzax .gt_from_md > :first-child {
  margin-top: 0;
}

#rvzeorrzax .gt_from_md > :last-child {
  margin-bottom: 0;
}

#rvzeorrzax .gt_row {
  padding-top: 2px;
  padding-bottom: 2px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#rvzeorrzax .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#rvzeorrzax .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#rvzeorrzax .gt_row_group_first td {
  border-top-width: 2px;
}

#rvzeorrzax .gt_row_group_first th {
  border-top-width: 2px;
}

#rvzeorrzax .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 2px;
  padding-bottom: 2px;
  padding-left: 5px;
  padding-right: 5px;
}

#rvzeorrzax .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#rvzeorrzax .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#rvzeorrzax .gt_last_summary_row {
  padding-top: 2px;
  padding-bottom: 2px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rvzeorrzax .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 2px;
  padding-bottom: 2px;
  padding-left: 5px;
  padding-right: 5px;
}

#rvzeorrzax .gt_first_grand_summary_row {
  padding-top: 2px;
  padding-bottom: 2px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#rvzeorrzax .gt_last_grand_summary_row_top {
  padding-top: 2px;
  padding-bottom: 2px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#rvzeorrzax .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#rvzeorrzax .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rvzeorrzax .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rvzeorrzax .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 1px;
  padding-bottom: 1px;
  padding-left: 5px;
  padding-right: 5px;
}

#rvzeorrzax .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rvzeorrzax .gt_sourcenote {
  font-size: 90%;
  padding-top: 1px;
  padding-bottom: 1px;
  padding-left: 5px;
  padding-right: 5px;
}

#rvzeorrzax .gt_left {
  text-align: left;
}

#rvzeorrzax .gt_center {
  text-align: center;
}

#rvzeorrzax .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#rvzeorrzax .gt_font_normal {
  font-weight: normal;
}

#rvzeorrzax .gt_font_bold {
  font-weight: bold;
}

#rvzeorrzax .gt_font_italic {
  font-style: italic;
}

#rvzeorrzax .gt_super {
  font-size: 65%;
}

#rvzeorrzax .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#rvzeorrzax .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#rvzeorrzax .gt_indent_1 {
  text-indent: 5px;
}

#rvzeorrzax .gt_indent_2 {
  text-indent: 10px;
}

#rvzeorrzax .gt_indent_3 {
  text-indent: 15px;
}

#rvzeorrzax .gt_indent_4 {
  text-indent: 20px;
}

#rvzeorrzax .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings gt_spanner_row">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1" scope="col" id=""></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="1" scope="col" id="Water concentration (ug/L)">
        <span class="gt_column_spanner">Water concentration (ug/L)</span>
      </th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="1" scope="col" id="Fish concentration (ug/Kg)">
        <span class="gt_column_spanner">Fish concentration (ug/Kg)</span>
      </th>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" style="text-align: center;" scope="col" id="Mean (95% CI)">Mean (95% CI)</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" style="text-align: center;" scope="col" id="Mean (95% CI)">Mean (95% CI)</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Heiabekken">Heiabekken</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_1" scope="row" class="gt_row gt_left gt_stub">imidacloprid</th>
<td headers="Heiabekken stub_1_1 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #08306B; color: #FFFFFF;">1.14 (0.26,2.02)</td>
<td headers="Heiabekken stub_1_1 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">2.73 (0.63,4.83)</td></tr>
    <tr><th id="stub_1_2" scope="row" class="gt_row gt_left gt_stub">propamocarb</th>
<td headers="Heiabekken stub_1_2 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #7FB8DA; color: #000000;">0.52 (−0.24,1.29)</td>
<td headers="Heiabekken stub_1_2 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">2.52 (−1.16,6.20)</td></tr>
    <tr><th id="stub_1_3" scope="row" class="gt_row gt_left gt_stub">2-methyl-4-chlorophenoxyacetic acid (MCPA)</th>
<td headers="Heiabekken stub_1_3 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #87BCDC; color: #000000;">0.50 (0.01,1.00)</td>
<td headers="Heiabekken stub_1_3 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F3F9FE; color: #000000;">28.01 (0.56,55.45)</td></tr>
    <tr><th id="stub_1_4" scope="row" class="gt_row gt_left gt_stub">pencycuron</th>
<td headers="Heiabekken stub_1_4 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #D9E8F5; color: #000000;">0.18 (0.01,0.35)</td>
<td headers="Heiabekken stub_1_4 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #B1D2E7; color: #000000;">468.87 (31.60,906.14)</td></tr>
    <tr><th id="stub_1_5" scope="row" class="gt_row gt_left gt_stub">metribuzin</th>
<td headers="Heiabekken stub_1_5 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #DAE8F6; color: #000000;">0.18 (−0.13,0.48)</td>
<td headers="Heiabekken stub_1_5 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">2.02 (−1.53,5.58)</td></tr>
    <tr><th id="stub_1_6" scope="row" class="gt_row gt_left gt_stub">fluroxypyr</th>
<td headers="Heiabekken stub_1_6 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #DDEAF7; color: #000000;">0.16 (—,—)</td>
<td headers="Heiabekken stub_1_6 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">2.20 (—,—)</td></tr>
    <tr><th id="stub_1_7" scope="row" class="gt_row gt_left gt_stub">mandipropamid</th>
<td headers="Heiabekken stub_1_7 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #DFECF7; color: #000000;">0.15 (−0.32,0.61)</td>
<td headers="Heiabekken stub_1_7 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #ECF4FC; color: #000000;">78.77 (−173.61,331.15)</td></tr>
    <tr><th id="stub_1_8" scope="row" class="gt_row gt_left gt_stub">clopyralid</th>
<td headers="Heiabekken stub_1_8 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #E6F0F9; color: #000000;">0.11 (—,—)</td>
<td headers="Heiabekken stub_1_8 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">3.63 (—,—)</td></tr>
    <tr><th id="stub_1_9" scope="row" class="gt_row gt_left gt_stub">boscalid</th>
<td headers="Heiabekken stub_1_9 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F2F8FD; color: #000000;">0.04 (0.03,0.05)</td>
<td headers="Heiabekken stub_1_9 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #E6F0F9; color: #000000;">128.91 (92.39,165.44)</td></tr>
    <tr><th id="stub_1_10" scope="row" class="gt_row gt_left gt_stub">prothioconazole-desthio</th>
<td headers="Heiabekken stub_1_10 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F3F8FE; color: #000000;">0.04 (0.01,0.06)</td>
<td headers="Heiabekken stub_1_10 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">2.42 (0.89,3.95)</td></tr>
    <tr><th id="stub_1_11" scope="row" class="gt_row gt_left gt_stub">metalaxyl</th>
<td headers="Heiabekken stub_1_11 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F3F8FE; color: #000000;">0.03 (0.02,0.05)</td>
<td headers="Heiabekken stub_1_11 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.32 (0.19,0.45)</td></tr>
    <tr><th id="stub_1_12" scope="row" class="gt_row gt_left gt_stub">trifloxystrobin</th>
<td headers="Heiabekken stub_1_12 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F4F9FE; color: #000000;">0.03 (—,—)</td>
<td headers="Heiabekken stub_1_12 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #ECF4FB; color: #000000;">81.09 (—,—)</td></tr>
    <tr><th id="stub_1_13" scope="row" class="gt_row gt_left gt_stub">metamitron</th>
<td headers="Heiabekken stub_1_13 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F5FAFE; color: #000000;">0.02 (0.01,0.03)</td>
<td headers="Heiabekken stub_1_13 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.04 (0.03,0.05)</td></tr>
    <tr><th id="stub_1_14" scope="row" class="gt_row gt_left gt_stub">propoxycarbazone</th>
<td headers="Heiabekken stub_1_14 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F5FAFF; color: #000000;">0.02 (0.01,0.02)</td>
<td headers="Heiabekken stub_1_14 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.39 (0.29,0.48)</td></tr>
    <tr><th id="stub_1_15" scope="row" class="gt_row gt_left gt_stub">diflufenican</th>
<td headers="Heiabekken stub_1_15 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F6FAFF; color: #000000;">0.02 (0.01,0.03)</td>
<td headers="Heiabekken stub_1_15 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F2F8FE; color: #000000;">33.95 (12.36,55.54)</td></tr>
    <tr><th id="stub_1_16" scope="row" class="gt_row gt_left gt_stub">bixafen</th>
<td headers="Heiabekken stub_1_16 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F6FAFF; color: #000000;">0.02 (0.01,0.03)</td>
<td headers="Heiabekken stub_1_16 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F2F8FD; color: #000000;">37.36 (14.70,60.03)</td></tr>
    <tr><th id="stub_1_17" scope="row" class="gt_row gt_left gt_stub">2,6-dichlorobenzamide (BAM)</th>
<td headers="Heiabekken stub_1_17 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F6FAFF; color: #000000;">0.02 (—,—)</td>
<td headers="Heiabekken stub_1_17 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.04 (—,—)</td></tr>
    <tr><th id="stub_1_18" scope="row" class="gt_row gt_left gt_stub">pyraclostrobin</th>
<td headers="Heiabekken stub_1_18 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F6FAFF; color: #000000;">0.02 (−0.04,0.07)</td>
<td headers="Heiabekken stub_1_18 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F5FAFE; color: #000000;">11.92 (−32.06,55.90)</td></tr>
    <tr><th id="stub_1_19" scope="row" class="gt_row gt_left gt_stub">spinosad</th>
<td headers="Heiabekken stub_1_19 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.01 (—,—)</td>
<td headers="Heiabekken stub_1_19 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #808080; color: #FFFFFF;">— (—,—)</td></tr>
    <tr><th id="stub_1_20" scope="row" class="gt_row gt_left gt_stub">1,1'-(2,2-Dichloro-1,1-ethenediyl)bis(4-chlorobenzene)</th>
<td headers="Heiabekken stub_1_20 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.01 (—,—)</td>
<td headers="Heiabekken stub_1_20 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #08306B; color: #FFFFFF;">1,476.32 (—,—)</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Skuterudbekken">Skuterudbekken</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_21" scope="row" class="gt_row gt_left gt_stub">2-methyl-4-chlorophenoxyacetic acid (MCPA)</th>
<td headers="Skuterudbekken stub_1_21 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #A1CBE2; color: #000000;">0.43 (−0.13,0.99)</td>
<td headers="Skuterudbekken stub_1_21 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F4F9FE; color: #000000;">23.70 (−7.50,54.89)</td></tr>
    <tr><th id="stub_1_22" scope="row" class="gt_row gt_left gt_stub">fluroxypyr</th>
<td headers="Skuterudbekken stub_1_22 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #CCDFF1; color: #000000;">0.26 (0.00,0.51)</td>
<td headers="Skuterudbekken stub_1_22 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">3.55 (0.02,7.07)</td></tr>
    <tr><th id="stub_1_23" scope="row" class="gt_row gt_left gt_stub">clopyralid</th>
<td headers="Skuterudbekken stub_1_23 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #E1EDF8; color: #000000;">0.13 (0.01,0.26)</td>
<td headers="Skuterudbekken stub_1_23 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F6FBFF; color: #000000;">4.42 (0.32,8.51)</td></tr>
    <tr><th id="stub_1_24" scope="row" class="gt_row gt_left gt_stub">prothioconazole-desthio</th>
<td headers="Skuterudbekken stub_1_24 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F1F7FD; color: #000000;">0.05 (−0.03,0.12)</td>
<td headers="Skuterudbekken stub_1_24 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">3.12 (−1.66,7.90)</td></tr>
    <tr><th id="stub_1_25" scope="row" class="gt_row gt_left gt_stub">dichlorprop</th>
<td headers="Skuterudbekken stub_1_25 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F3F8FE; color: #000000;">0.04 (—,—)</td>
<td headers="Skuterudbekken stub_1_25 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F6FAFF; color: #000000;">7.91 (—,—)</td></tr>
    <tr><th id="stub_1_26" scope="row" class="gt_row gt_left gt_stub">imidacloprid</th>
<td headers="Skuterudbekken stub_1_26 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F5FAFE; color: #000000;">0.02 (—,—)</td>
<td headers="Skuterudbekken stub_1_26 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.05 (—,—)</td></tr>
    <tr><th id="stub_1_27" scope="row" class="gt_row gt_left gt_stub">bixafen</th>
<td headers="Skuterudbekken stub_1_27 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F5FAFE; color: #000000;">0.02 (0.00,0.04)</td>
<td headers="Skuterudbekken stub_1_27 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F1F7FD; color: #000000;">45.42 (7.07,83.77)</td></tr>
    <tr><th id="stub_1_28" scope="row" class="gt_row gt_left gt_stub">propiconazole</th>
<td headers="Skuterudbekken stub_1_28 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F5FAFE; color: #000000;">0.02 (−0.02,0.06)</td>
<td headers="Skuterudbekken stub_1_28 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F6FBFF; color: #000000;">5.52 (−6.45,17.49)</td></tr>
    <tr><th id="stub_1_29" scope="row" class="gt_row gt_left gt_stub">boscalid</th>
<td headers="Skuterudbekken stub_1_29 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F6FAFF; color: #000000;">0.02 (0.01,0.02)</td>
<td headers="Skuterudbekken stub_1_29 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F0F7FD; color: #000000;">49.90 (29.40,70.40)</td></tr>
    <tr><th id="stub_1_30" scope="row" class="gt_row gt_left gt_stub">pyroxsulam</th>
<td headers="Skuterudbekken stub_1_30 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F6FBFF; color: #000000;">0.01 (−0.02,0.05)</td>
<td headers="Skuterudbekken stub_1_30 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.16 (−0.28,0.60)</td></tr>
    <tr><th id="stub_1_31" scope="row" class="gt_row gt_left gt_stub">mecoprop</th>
<td headers="Skuterudbekken stub_1_31 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.01 (—,—)</td>
<td headers="Skuterudbekken stub_1_31 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">1.74 (—,—)</td></tr>
    <tr><th id="stub_1_32" scope="row" class="gt_row gt_left gt_stub">carbendazim</th>
<td headers="Skuterudbekken stub_1_32 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.01 (—,—)</td>
<td headers="Skuterudbekken stub_1_32 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.11 (—,—)</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Mørdrebekken">Mørdrebekken</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_33" scope="row" class="gt_row gt_left gt_stub">2-methyl-4-chlorophenoxyacetic acid (MCPA)</th>
<td headers="Mørdrebekken stub_1_33 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #BBD6EB; color: #000000;">0.33 (−0.35,1.01)</td>
<td headers="Mørdrebekken stub_1_33 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F5F9FE; color: #000000;">18.48 (−19.28,56.24)</td></tr>
    <tr><th id="stub_1_34" scope="row" class="gt_row gt_left gt_stub">fluroxypyr</th>
<td headers="Mørdrebekken stub_1_34 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #CDE0F1; color: #000000;">0.25 (−1.53,2.03)</td>
<td headers="Mørdrebekken stub_1_34 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">3.44 (−21.01,27.88)</td></tr>
    <tr><th id="stub_1_35" scope="row" class="gt_row gt_left gt_stub">bentazone</th>
<td headers="Mørdrebekken stub_1_35 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #D5E5F4; color: #000000;">0.21 (−0.11,0.53)</td>
<td headers="Mørdrebekken stub_1_35 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F5FAFE; color: #000000;">16.25 (−9.02,41.51)</td></tr>
    <tr><th id="stub_1_36" scope="row" class="gt_row gt_left gt_stub">clopyralid</th>
<td headers="Mørdrebekken stub_1_36 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #D6E6F4; color: #000000;">0.20 (—,—)</td>
<td headers="Mørdrebekken stub_1_36 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F6FAFF; color: #000000;">6.59 (—,—)</td></tr>
    <tr><th id="stub_1_37" scope="row" class="gt_row gt_left gt_stub">propiconazole</th>
<td headers="Mørdrebekken stub_1_37 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #DCEAF6; color: #000000;">0.16 (−0.28,0.61)</td>
<td headers="Mørdrebekken stub_1_37 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F1F7FD; color: #000000;">43.87 (−75.26,163.01)</td></tr>
    <tr><th id="stub_1_38" scope="row" class="gt_row gt_left gt_stub">pencycuron</th>
<td headers="Mørdrebekken stub_1_38 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #EDF4FC; color: #000000;">0.07 (−0.31,0.45)</td>
<td headers="Mørdrebekken stub_1_38 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #DEEBF7; color: #000000;">183.27 (−814.74,1,181.29)</td></tr>
    <tr><th id="stub_1_39" scope="row" class="gt_row gt_left gt_stub">fenpropimorph</th>
<td headers="Mørdrebekken stub_1_39 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F0F6FD; color: #000000;">0.05 (—,—)</td>
<td headers="Mørdrebekken stub_1_39 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #D2E3F3; color: #000000;">274.16 (—,—)</td></tr>
    <tr><th id="stub_1_40" scope="row" class="gt_row gt_left gt_stub">thiabendazole</th>
<td headers="Mørdrebekken stub_1_40 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F2F8FD; color: #000000;">0.04 (—,—)</td>
<td headers="Mørdrebekken stub_1_40 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">1.92 (—,—)</td></tr>
    <tr><th id="stub_1_41" scope="row" class="gt_row gt_left gt_stub">prothioconazole-desthio</th>
<td headers="Mørdrebekken stub_1_41 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F3F9FE; color: #000000;">0.03 (−0.03,0.09)</td>
<td headers="Mørdrebekken stub_1_41 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">2.17 (−1.93,6.28)</td></tr>
    <tr><th id="stub_1_42" scope="row" class="gt_row gt_left gt_stub">pinoxaden</th>
<td headers="Mørdrebekken stub_1_42 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F3F9FE; color: #000000;">0.03 (—,—)</td>
<td headers="Mørdrebekken stub_1_42 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #EFF6FD; color: #000000;">57.20 (—,—)</td></tr>
    <tr><th id="stub_1_43" scope="row" class="gt_row gt_left gt_stub">azoxystrobin</th>
<td headers="Mørdrebekken stub_1_43 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F5FAFE; color: #000000;">0.02 (0.01,0.03)</td>
<td headers="Mørdrebekken stub_1_43 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F6FAFF; color: #000000;">8.21 (4.05,12.38)</td></tr>
    <tr><th id="stub_1_44" scope="row" class="gt_row gt_left gt_stub">bixafen</th>
<td headers="Mørdrebekken stub_1_44 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F6FBFF; color: #000000;">0.01 (—,—)</td>
<td headers="Mørdrebekken stub_1_44 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F3F8FE; color: #000000;">32.97 (—,—)</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Timebekken">Timebekken</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_45" scope="row" class="gt_row gt_left gt_stub">metribuzin</th>
<td headers="Timebekken stub_1_45 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #C4DAEE; color: #000000;">0.30 (−3.36,3.96)</td>
<td headers="Timebekken stub_1_45 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">3.48 (−38.73,45.69)</td></tr>
    <tr><th id="stub_1_46" scope="row" class="gt_row gt_left gt_stub">fluroxypyr</th>
<td headers="Timebekken stub_1_46 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #C9DDF0; color: #000000;">0.28 (0.04,0.51)</td>
<td headers="Timebekken stub_1_46 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F6FBFF; color: #000000;">3.78 (0.53,7.02)</td></tr>
    <tr><th id="stub_1_47" scope="row" class="gt_row gt_left gt_stub">florasulam</th>
<td headers="Timebekken stub_1_47 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F5FAFE; color: #000000;">0.02 (—,—)</td>
<td headers="Timebekken stub_1_47 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.22 (—,—)</td></tr>
    <tr><th id="stub_1_48" scope="row" class="gt_row gt_left gt_stub">propiconazole</th>
<td headers="Timebekken stub_1_48 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F6FAFF; color: #000000;">0.02 (—,—)</td>
<td headers="Timebekken stub_1_48 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F6FBFF; color: #000000;">4.31 (—,—)</td></tr>
    <tr><th id="stub_1_49" scope="row" class="gt_row gt_left gt_stub">2-methyl-4-chlorophenoxyacetic acid (MCPA)</th>
<td headers="Timebekken stub_1_49 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.01 (—,—)</td>
<td headers="Timebekken stub_1_49 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.61 (—,—)</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="3" class="gt_group_heading" scope="colgroup" id="Vasshaglona">Vasshaglona</th>
    </tr>
    <tr class="gt_row_group_first"><th id="stub_1_50" scope="row" class="gt_row gt_left gt_stub">fenamidone</th>
<td headers="Vasshaglona stub_1_50 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #D4E5F4; color: #000000;">0.21 (−1.72,2.13)</td>
<td headers="Vasshaglona stub_1_50 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #E1EDF8; color: #000000;">160.36 (−1,320.20,1,640.93)</td></tr>
    <tr><th id="stub_1_51" scope="row" class="gt_row gt_left gt_stub">propamocarb</th>
<td headers="Vasshaglona stub_1_51 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #E7F1FA; color: #000000;">0.10 (−0.09,0.29)</td>
<td headers="Vasshaglona stub_1_51 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.49 (−0.43,1.41)</td></tr>
    <tr><th id="stub_1_52" scope="row" class="gt_row gt_left gt_stub">pyridafol</th>
<td headers="Vasshaglona stub_1_52 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #EAF3FB; color: #000000;">0.08 (−0.04,0.21)</td>
<td headers="Vasshaglona stub_1_52 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">3.28 (−1.62,8.18)</td></tr>
    <tr><th id="stub_1_53" scope="row" class="gt_row gt_left gt_stub">metribuzin</th>
<td headers="Vasshaglona stub_1_53 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #EBF3FB; color: #000000;">0.08 (0.02,0.14)</td>
<td headers="Vasshaglona stub_1_53 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.91 (0.25,1.56)</td></tr>
    <tr><th id="stub_1_54" scope="row" class="gt_row gt_left gt_stub">thiacloprid</th>
<td headers="Vasshaglona stub_1_54 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #EDF4FC; color: #000000;">0.07 (−0.57,0.71)</td>
<td headers="Vasshaglona stub_1_54 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">1.92 (−15.83,19.68)</td></tr>
    <tr><th id="stub_1_55" scope="row" class="gt_row gt_left gt_stub">fluroxypyr</th>
<td headers="Vasshaglona stub_1_55 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #EDF5FC; color: #000000;">0.07 (0.03,0.11)</td>
<td headers="Vasshaglona stub_1_55 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.92 (0.34,1.49)</td></tr>
    <tr><th id="stub_1_56" scope="row" class="gt_row gt_left gt_stub">cyazofamid</th>
<td headers="Vasshaglona stub_1_56 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F0F6FD; color: #000000;">0.05 (—,—)</td>
<td headers="Vasshaglona stub_1_56 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">2.90 (—,—)</td></tr>
    <tr><th id="stub_1_57" scope="row" class="gt_row gt_left gt_stub">2-methyl-4-chlorophenoxyacetic acid (MCPA)</th>
<td headers="Vasshaglona stub_1_57 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F0F7FD; color: #000000;">0.05 (−0.07,0.17)</td>
<td headers="Vasshaglona stub_1_57 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">2.79 (−3.72,9.29)</td></tr>
    <tr><th id="stub_1_58" scope="row" class="gt_row gt_left gt_stub">aclonifen</th>
<td headers="Vasshaglona stub_1_58 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F0F7FD; color: #000000;">0.05 (−0.04,0.14)</td>
<td headers="Vasshaglona stub_1_58 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F4F9FE; color: #000000;">22.41 (−16.99,61.80)</td></tr>
    <tr><th id="stub_1_59" scope="row" class="gt_row gt_left gt_stub">boscalid</th>
<td headers="Vasshaglona stub_1_59 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F0F7FD; color: #000000;">0.05 (0.02,0.08)</td>
<td headers="Vasshaglona stub_1_59 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #E3EEF8; color: #000000;">150.75 (64.62,236.87)</td></tr>
    <tr><th id="stub_1_60" scope="row" class="gt_row gt_left gt_stub">pencycuron</th>
<td headers="Vasshaglona stub_1_60 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F1F7FD; color: #000000;">0.05 (−0.01,0.10)</td>
<td headers="Vasshaglona stub_1_60 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #E7F1FA; color: #000000;">118.94 (−15.84,253.72)</td></tr>
    <tr><th id="stub_1_61" scope="row" class="gt_row gt_left gt_stub">mandipropamid</th>
<td headers="Vasshaglona stub_1_61 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F1F7FD; color: #000000;">0.04 (0.00,0.09)</td>
<td headers="Vasshaglona stub_1_61 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F4F9FE; color: #000000;">23.49 (−1.85,48.83)</td></tr>
    <tr><th id="stub_1_62" scope="row" class="gt_row gt_left gt_stub">spinosad</th>
<td headers="Vasshaglona stub_1_62 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F3F8FE; color: #000000;">0.03 (—,—)</td>
<td headers="Vasshaglona stub_1_62 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #808080; color: #FFFFFF;">— (—,—)</td></tr>
    <tr><th id="stub_1_63" scope="row" class="gt_row gt_left gt_stub">cyprodinil</th>
<td headers="Vasshaglona stub_1_63 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F3F8FE; color: #000000;">0.03 (−0.05,0.11)</td>
<td headers="Vasshaglona stub_1_63 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F6FBFF; color: #000000;">4.46 (−6.41,15.32)</td></tr>
    <tr><th id="stub_1_64" scope="row" class="gt_row gt_left gt_stub">clomazone</th>
<td headers="Vasshaglona stub_1_64 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F3F9FE; color: #000000;">0.03 (−0.04,0.10)</td>
<td headers="Vasshaglona stub_1_64 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">1.50 (−1.65,4.64)</td></tr>
    <tr><th id="stub_1_65" scope="row" class="gt_row gt_left gt_stub">bentazone</th>
<td headers="Vasshaglona stub_1_65 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F4F9FE; color: #000000;">0.03 (0.02,0.03)</td>
<td headers="Vasshaglona stub_1_65 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">2.06 (1.48,2.65)</td></tr>
    <tr><th id="stub_1_66" scope="row" class="gt_row gt_left gt_stub">fludioxonil</th>
<td headers="Vasshaglona stub_1_66 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F5FAFE; color: #000000;">0.02 (−0.01,0.06)</td>
<td headers="Vasshaglona stub_1_66 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">1.17 (−0.75,3.09)</td></tr>
    <tr><th id="stub_1_67" scope="row" class="gt_row gt_left gt_stub">imidacloprid</th>
<td headers="Vasshaglona stub_1_67 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F5FAFE; color: #000000;">0.02 (0.00,0.04)</td>
<td headers="Vasshaglona stub_1_67 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">0.05 (0.01,0.08)</td></tr>
    <tr><th id="stub_1_68" scope="row" class="gt_row gt_left gt_stub">chlorfenvinphos</th>
<td headers="Vasshaglona stub_1_68 Water concentration (ug/L)_Mean" class="gt_row gt_right" style="background-color: #F6FAFF; color: #000000;">0.02 (—,—)</td>
<td headers="Vasshaglona stub_1_68 Fish concentration (ug/Kg)_Mean" class="gt_row gt_right" style="background-color: #F7FBFF; color: #000000;">2.27 (—,—)</td></tr>
  </tbody>
  
  
</table>
</div>

## Spatio-temporal plot of water vs tissue concentration

<details>
<summary>Comparison plot water and fish concentration</summary>

``` r
plot_data <- exposure_data %>% 
  mutate(Month = lubridate::month(SAMPLE_DATE, label = TRUE, abbr = FALSE)) %>% 
  mutate(Wday = lubridate::wday(SAMPLE_DATE, label = TRUE, abbr = FALSE)) %>% 
  mutate(FishConcentration = get_tissue_conc(MEASURED_VALUE, XLogP)) %>% 
  select(
    SITE_NAME, STRESSOR_NAME, 
    WaterConc = MEASURED_VALUE,
    FishConc = FishConcentration,
    Month, Wday
  )


make_two_fitted_plot <- function(data, time_var) {
  plot_data %>% 
    ggplot(aes(WaterConc, FishConc)) +
    facet_grid(
      cols = vars(SITE_NAME), 
      rows = vars(get(time_var)),
      scales = "free_x"
    ) +
    geom_point(alpha = 0.5) +
    geom_smooth(
      method = "lm",
      formula = "y ~ x",
      linewidth = 0.5
    ) +
    scale_x_log10(
      breaks = 10^c(seq(-2, 2, 1)),
      labels = scales::label_log()
    ) +
    scale_y_log10(labels = scales::label_log(digits = 2)) +
    ggthemes::theme_few() +
    theme(
      panel.grid = element_line(color = "#f0f0f0")
    ) +
    labs(
      x = "log(Water concentration (ug/L))",
      y = "log(Fish concentration (ug/Kg))"
    )
}
```

</details>

<div class="panel-tabset">

### By month

<details>
<summary>Comparison plot water and fish concentration by month</summary>

``` r
make_two_fitted_plot(plot_data, "Month")
```

</details>

<img
src="index_files/figure-commonmark/Comparison%20plot%20water%20and%20fish%20concentration%20by%20month-1.png"
style="width:100.0%" />

### By weekday

<details>
<summary>Comparison plot water and fish concentration by
weekday</summary>

``` r
make_two_fitted_plot(plot_data, "Wday")
```

</details>

<img
src="index_files/figure-commonmark/Comparison%20plot%20water%20and%20fish%20concentration%20by%20weekday-1.png"
style="width:100.0%" />

</div>

## Presentation of the results in map of Norway

<details>
<summary>Water and fish concentration in map</summary>

``` r
plot_data <- exposure_data %>% 
  group_by(SITE_NAME, LATITUDE, LONGITUDE) %>% 
  rename(WaterConc = MEASURED_VALUE, FishConc = FishConcentration) %>% 
  summarize(across(c(WaterConc, FishConc), function(x) {
    Hmisc::smean.cl.normal(x) %>% 
      as.list() %>% 
      bind_cols() %>% 
      list()
  })) %>% unnest(names_sep = "_") %>% 
  pivot_longer(
    id_cols = c(SITE_NAME, LATITUDE, LONGITUDE), 
    cols = -c(SITE_NAME, LATITUDE, LONGITUDE)
  ) %>% tidytable::separate("name", c("Medium", "Measure"), "_") %>% 
  pivot_wider(names_from = "Measure", values_from = "value")

csmaps::nor_municip_map_b2020_split_dt %>% 
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(fill = NA, color = "darkgrey", linewidth = 0.1) +
  geom_label(
    data = plot_data, 
    aes(x = LONGITUDE, y = LATITUDE, label = SITE_NAME),
    group = 1,
    size = rel(3),
    hjust = -0.2,
    label.r = unit(2.5, "pt")
  ) +
  geom_point(
    data = plot_data, 
    aes(x = LONGITUDE, y = LATITUDE, size = Mean, color = Medium),
    group = 1
  ) +
  scale_color_brewer(
    palette = "Set1", 
    direction = -1,
    labels = c(
      FishConc = "Fish concentration (ug/Kg)",
      WaterConc = "Water concentration (ug/L)"
    )
  ) +
  theme_void() +
  coord_map() +
  labs(
    title = "Monitoring sites of Spatio-temporal exposure data",
    subtitle = "Concentration level of stressor in water and fish tissue",
    caption = "Year: 2019, Country: Norway",
    size = "Mean concentration",
    color = "Medium"
  )
```

</details>

<img
src="index_files/figure-commonmark/Water%20and%20fish%20concentration%20in%20map-1.png"
style="width:100.0%" />
