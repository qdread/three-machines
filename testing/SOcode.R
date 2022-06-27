data <- data.frame(
  Trial = c("H6", "H6", "H6", "H6", "H6", "U8", "U8", "U8", "U8", "U8"),
  Year = c(2020L, 2020L, 2020L, 2020L, 2020L, 2020L, 2020L, 2020L, 2020L, 2020L),
  Plot = c(101L, 102L, 103L, 104L, 107L, 311L, 312L, 313L, 314L, 315L),
  G1 = c(55.8, 56.3, 55.7, 57.5, 58.3, 57.4, 57.6, 54.9, 57.1, 54.4),
  G2 = c(55.8, 56.3, 55.7, 57.5, 58.2, 57.4, 57.6, 54.9, 57, 54.4),
  G3 = c(55.8, 56.3, 55.7, 57.5, 58.2, 57.4, 57.6, 54.9, 57, 54.4),
  P1 = c(56.3, 57.2, 58, 57.9, 57.9, 57.3, 57.5, 55.7, 56.3, 55.5),
  P2 = c(56.1, 57.4, 57.8, 58.5, 57.1, 57, 57, 55.3, 56.1, 55.6),
  P3 = c(56.9, 57.2, 57.9, 58.8, 57.2, 57.3, 57.5, 55.6, 56.1, 55.6),
  V1 = c(56.3, 56.3, 57.7, 57.1, 56.2, 56.8, 56.4, 54.9, 55, 54.7),
  V2 = c(56, 56.5, 57.3, 57.6, 56.5, 56.2, 56.6, 54.2, 54.9, 54.5),
  V3 = c(56.2, 56.6, 57, 57.6, 56.6, 56.2, 56.3, 54.8, 54.9, 54.4)
)

library("tidyverse")
#iccmlr::theme_set_amazon()

GiniSd <- function(x, y = NULL, na.rm = TRUE) {
  if (is.null(y)) {
    # Adapted from Hmisc::GiniMd
    if (na.rm) {
      k <- is.na(x)
      if (any(k)) {
        x <- x[!k]
      }
    }
    n <- length(x)
    if (n < 2) {
      return(NA)
    }
    w <- 2 * ((1:n) - (n - 1) / 2)
    sum(w * sort(x - mean(x)))
  } else {
    sum(outer(x, y, function(xi, yi) abs(xi - yi)))
  }
}

Aggregate <- function(variabilities, ...) {
  variabilities %>%
    group_by(...) %>%
    summarise(
      variability = sum(N * variability) / sum(N),
      N = sum(N)
    )
}

Estimate <- function(sample, instrument, measurement) {
  measurements <- tibble(
    sample,
    instrument,
    measurement
  )
  measurements %>%
    group_by(
      sample,
      instrument
    ) %>%
    summarise(
      measurement = median(measurement)
    ) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = sample,
      names_from = instrument,
      values_from = measurement
    )
}

IntraVariability <- function(sample, instrument, measurement) {
  measurements <- tibble(
    sample,
    instrument,
    measurement
  )
  intra_var <- measurements %>%
    group_by(
      sample,
      instrument
    ) %>%
    summarise(
      N = choose(n(), 2),
      variability = GiniSd(measurement) / N
    ) %>%
    ungroup() %>%
    drop_na()
  list(
    intra_var = intra_var,
    by_sample = Aggregate(intra_var, sample),
    by_instrument = Aggregate(intra_var, instrument)
  )
}

InterVariability <- function(sample, instrument, measurement) {
  measurements <- tibble(
    sample,
    instrument,
    measurement
  )
  
  nested <- measurements %>%
    chop(measurement)
  
  inter_var <-
    inner_join(
      nested %>% rename(instrument1 = instrument, x = measurement),
      nested %>% rename(instrument2 = instrument, y = measurement),
      by = "sample"
    ) %>%
    filter(
      instrument1 > instrument2
    ) %>%
    rowwise() %>%
    mutate(
      N = length(x) * length(y),
      variability = GiniSd(x, y) / N
    ) %>%
    select(
      -c(x, y)
    ) %>%
    drop_na()
  list(
    inter_var = inter_var,
    by_sample = Aggregate(inter_var, sample),
    by_instrument_pairs = Aggregate(inter_var, instrument1, instrument2)
  )
}

plot_measurements <- function(measurements, measurement) {
  measurements %>%
    arrange(
      sample, instrument, {{ measurement }}
    ) %>%
    mutate(
      id = cumsum(sample != lag(sample, default = ""))
    ) %>%
    group_by(
      sample, instrument
    ) %>%
    mutate(
      h = cumsum({{ measurement }} == lag({{ measurement }}, default = 0)),
      y = id + 0.1 * h
    ) %>%
    ggplot(
      aes({{ measurement }}, y,
          shape = instrument,
          color = instrument
      )
    ) +
    geom_hline(
      aes(yintercept = y),
      linetype = 3,
      color = "gray",
      size = 0.5,
      data = tibble(y = seq(1, 10, 1))
    ) +
    geom_point(
      stroke = 1,
      size = 2
    ) +
    scale_shape_manual(
      values = c(0, 1, 2)
    ) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
    )
}

measurements <- data %>%
  unite(
    "sample", c(Trial, Year, Plot)
  ) %>%
  pivot_longer(
    -sample,
    names_to = "instrument",
    names_transform = list(instrument = ~ str_sub(., 1, 1))
  )

estimates <- Estimate(
  measurements$sample,
  measurements$instrument,
  measurements$value
)
estimates %>%
  summarise(
    `G - V` = mean(G - V),
    `P - V` = mean(P - V)
  )

intra_var <- IntraVariability(
  measurements$sample,
  measurements$instrument,
  measurements$value
)
intra_var$by_instrument

inter_var <- InterVariability(
  measurements$sample,
  measurements$instrument,
  measurements$value
)
inter_var$by_instrument_pairs

# Make scatterplots

p <- plot_measurements(measurements, value)
p + ggtitle("Raw measurements")

p <- plot_measurements(
  estimates %>% pivot_longer(c(G, P, V), names_to = "instrument", values_to = "estimate"),
  estimate
)
p + ggtitle("Estimated weight")

p <- plot_measurements(intra_var$intra_var, variability)
p + ggtitle("Intra variability")

p <- plot_measurements(
  inter_var$inter_var %>% filter(instrument1 == "V") %>% rename(instrument = instrument2),
  variability
)
p + ggtitle("Inter variability (V - ?)")