Spatio-Temporal Exposure Analysis
================
Raju Rimal
5/12/24

## Dataset

In this analysis we will use the exposure data available at [zenodo
((https://doi.org/10.5281/zenodo.11093687))](https://zenodo.org/records/11093687 "Spatio-temporal dataset - Exposure data").
The dataset contains exposure information from 6-location in Norway from
the year 2019.

<details>
<summary>Code</summary>

``` r
csmaps::nor_municip_map_b2020_split_dt %>% 
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(fill = NA, color = "darkgrey", linewidth = 0.1) +
  geom_label(
    data = site_data, 
    aes(x = LONGITUDE, y = LATITUDE, label = SITE_NAME),
    group = 1,
    size = rel(3),
    hjust = -0.1,
    label.r = unit(2.5, "pt")
  ) +
  geom_point(
    data = site_data, 
    aes(x = LONGITUDE, y = LATITUDE),
    group = 1,
    color = "firebrick"
  ) +
  theme_void() +
  coord_map() +
  labs(
    title = "Monitoring sites of Spatio-temporal exposure data",
    subtitle = "Year: 2019, Country: Norway"
  )
```

</details>

<img src="README_files/figure-commonmark/unnamed-chunk-2-1.png"
style="width:100.0%" />

In each of these locations, levels of various stressors were measured in
\[Freshwater\] between dates \[mai 06, 2019 and oktober 28, 2019\]. The
average level of contration of each stressor in these location is as
follows.

## Average level of stressor

<details>
<summary>Code</summary>

``` r
exposure_data %>% 
    group_by(STRESSOR_NAME, SITE_NAME) %>% 
    summarize(across(MEASURED_VALUE, ~list(ggplot2::mean_se(.x)))) %>% 
    unnest() %>% 
    rename(x = y, xmin = ymin, xmax = ymax) %>% 
    ggplot(aes(x, STRESSOR_NAME)) +
    facet_grid(cols = vars(SITE_NAME)) +
    geom_pointrange(
        aes(xmin = xmin, xmax = xmax),
        shape = 21, fill = "whitesmoke",
        stroke = 1
    ) +
    ggthemes::theme_few() +
    theme(
        panel.grid = element_line(color = "#f0f0f0")
    ) +
    labs(
        x = "Contration in water (ug/L)",
        title = "Average contration of stressors in each location",
        y = NULL
    )
```

</details>

<img src="README_files/figure-commonmark/Average%20exposure-1.png"
style="width:100.0%" />

<details>
<summary>Code</summary>

``` r
exposure_data %>% 
  group_by(STRESSOR_NAME, SITE_NAME) %>% 
  summarize(across(MEASURED_VALUE, ~list(ggplot2::mean_se(.x)))) %>% 
  unnest() %>% 
  ggplot(aes(STRESSOR_NAME, SITE_NAME, fill = y)) +
  geom_tile(color = "#f0f0f0", linewidth = 0.1) +
  scale_fill_viridis_c(
    trans = "log", na.value = "white", 
    breaks = scales::breaks_log(10),
    labels = \(x) formatC(x, format = "f", flag = 0, digits = 3),
    guide = guide_colorbar(barwidth = unit(0.8, "npc"))
  ) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  coord_equal() +
  ggthemes::theme_few() +
  theme(
    legend.position = "top"
  ) +
  labs(
    x = NULL,
    y = "Location",
    fill = NULL
  )
```

</details>

<img
src="README_files/figure-commonmark/Average%20exposure%20heatmap-1.png"
style="width:100.0%" />

## Distribution of concentration

Water concentration (ug/L) of these stressor are skewed towards zero
with long right tail. In one of the location (Timebekken), because of
small sample size, the distribution is not clear.

<details>
<summary>Code</summary>

``` r
exposure_data %>% 
  ggplot(aes(MEASURED_VALUE)) +
  facet_grid(cols = vars(SITE_NAME), scales = "free_x") +
  geom_histogram(
    aes(y = after_stat(density)), 
    bins = 15, linewidth = 0.5,
    fill = "whitesmoke", 
    color = "#3d3d3d"
  ) +
  geom_density(color = "forestgreen", linewidth = 1) +
  geom_rug(length = unit(5, "pt"), alpha = 0.5) +
  scale_x_continuous(breaks = scales::breaks_extended(6)) +
  ggthemes::theme_few() +
  theme(panel.grid = element_line(color = "#f0f0f0")) +
  labs(
    x = "Concentration in water (ug/L)",
    y = "Density"
  )
```

</details>

<img
src="README_files/figure-commonmark/Distribution%20of%20measured%20values-1.png"
style="width:100.0%" />

## Chemical properties of stressor

The chemical properties of stressor such as “MolecularWeight” and
partition coefficient (LogP or XLogP) are extracted using
`webchem::pc_prop` command from `webchem` package using compound ID. The
Compound ID were fetched using the command `webchem::get_cid` using the
`INCHIKEY` variable available in the dataset.

## Concentration in fish

## Shiny App for Exploration

A [shiny app](https://therimalaya.shinyapps.io/NIVA-AppExposure/) can
help to explore the data further. To run the shiny app from the local
computer,

<details>
<summary>Code</summary>

``` r
shiny::runGitHub(
  username = "Mathatistics", 
  repo = "NIVA-ExposureAnalysis", 
  subdir = "ShinyApp/",
  ref = "main"
)
```

</details>
