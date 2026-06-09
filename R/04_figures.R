# =============================================================================
# 04_figures.R
# Generate and save all figures.
# Requires objects produced by 03_model_fitting.R
# =============================================================================

theme_pub <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# -----------------------------------------------------------------------------
# Figure 1a – Holling model fit: Grazer
# -----------------------------------------------------------------------------
p_grazer <- ggplot(mydata.Grazer, aes(x, y)) +
  geom_point(aes(col = treatment)) +
  geom_line(aes(y = model_pp( params.Grazer[1],                    mydata.Grazer))) +
  geom_line(aes(y = model_pp2(params.Grazer2[1], params.Grazer2[2], mydata.Grazer)), lty = 2) +
  theme_pub +
  labs(
    y     = expression(paste("Ingestion (cells grazer"^{-1}, " h"^{-1}, ")")),
    x     = expression(paste("Prochlorococcus abundance (10"^7, " cells mL"^{-1}, ")")),
    title = "Holling Type I (—) and Type II (- -) fit: Grazer"
  )

ggsave("figures/fig_grazer_fit.png", p_grazer, width = 6, height = 4, dpi = 300)

# -----------------------------------------------------------------------------
# Figure 1b – Holling model fit: Phage
# -----------------------------------------------------------------------------
p_phage <- ggplot(mydata.Phage, aes(x, y)) +
  geom_point(aes(col = treatment)) +
  geom_line(aes(y = model_pp( params.Phage[1],                  mydata.Phage))) +
  geom_line(aes(y = model_pp2(params.Phage2[1], params.Phage2[2], mydata.Phage)), lty = 2) +
  theme_pub +
  labs(
    y     = expression(paste("Lysis (cells phage"^{-1}, " h"^{-1}, ")")),
    x     = expression(paste("Prochlorococcus abundance (10"^7, " cells mL"^{-1}, ")")),
    title = "Holling Type I (—) and Type II (- -) fit: Phage"
  )

ggsave("figures/fig_phage_fit.png", p_phage, width = 6, height = 4, dpi = 300)

# -----------------------------------------------------------------------------
# Figure S3 – Observed vs. theoretical mortality, all treatments
# -----------------------------------------------------------------------------
p_all <- data_w3 %>%
  ggplot(aes(time, loss_Prochlorococcus * abundance_Prochlorococcus, group = series)) +
  geom_line(col = "grey") +
  geom_point(aes(
    y   = theoretical.loss_Prochlorococcus * abundance_Prochlorococcus,
    col = series
  )) +
  facet_wrap(~ treatment) +
  theme_pub +
  labs(
    x     = "Time (h)",
    y     = expression(paste("Prochlorococcus Mortality (10"^7, " cells mL"^{-1}, " h"^{-1}, ")")),
    title = "Observed (grey) vs. Theoretical mortality by treatment"
  )

ggsave("figures/fig_mortality_all.png", p_all, width = 10, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# Figure 2 – Mean ± SD observed vs. modelled mortality (selected treatments)
# -----------------------------------------------------------------------------
selected_treatments <- c(
  "Low grazer - high virus", "High virus", "Low virus",
  "High grazer", "Low grazer"
)

data_mean <- data_w3 %>%
  pivot_longer(
    cols      = -c(time:series),
    names_to  = c(".value", "population"),
    names_sep = "_"
  ) %>%
  filter(population == "Prochlorococcus") %>%
  filter(treatment %in% selected_treatments) %>%
  mutate(abundance = abundance / 1e7) %>%
  group_by(population, time, treatment) %>%
  summarize(
    mean_Observed = mean(loss * abundance,                  na.rm = TRUE),
    mean_Total    = mean(theoretical.loss * abundance,      na.rm = TRUE),
    mean_Lysis    = mean(theoretical.loss.Phage * abundance, na.rm = TRUE),
    mean_Grazing  = mean(theoretical.loss.Grazer * abundance, na.rm = TRUE),
    sd_Observed   = sd(loss * abundance,                    na.rm = TRUE),
    sd_Total      = sd(theoretical.loss * abundance,        na.rm = TRUE),
    sd_Lysis      = sd(theoretical.loss.Phage * abundance,  na.rm = TRUE),
    sd_Grazing    = sd(theoretical.loss.Grazer * abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols      = -c(population, time, treatment),
    names_to  = c(".value", "source"),
    names_sep = "_"
  ) %>%
  mutate(source = factor(source,
    levels = c("Observed", "Grazing", "Lysis", "Total"),
    labels = c("Observed", "Grazing", "Lysis", "Grazing + Lysis")
  ))

p_fig2 <- data_mean %>%
  ggplot(aes(time, mean, group = source)) +
  geom_ribbon(
    data = data_mean %>% filter(source != "Observed"),
    aes(ymin = mean - sd, ymax = mean + sd, fill = source),
    alpha = 0.3, show.legend = FALSE
  ) +
  geom_line(aes(col = source)) +
  geom_point(
    data = data_mean %>% filter(source == "Observed"),
    aes(time, mean)
  ) +
  geom_linerange(
    data = data_mean %>% filter(source == "Observed"),
    aes(time, ymin = mean - sd, ymax = mean + sd),
    lwd = 0.5
  ) +
  scale_color_manual(
    values = c("black", "seagreen3", "deepskyblue3", "red3"),
    name   = "",
    labels = c("Observed", "Grazing", "Lysis", "Grazing + Lysis")
  ) +
  scale_fill_manual(
    values = c("seagreen3", "deepskyblue3", "red3"),
    name   = "",
    labels = c("Grazing", "Lysis", "Grazing + Lysis")
  ) +
  facet_wrap(treatment ~ ., ncol = 1) +
  theme_pub +
  theme(legend.position = c(0.7, 0.05)) +
  labs(
    x = "Time (h)",
    y = expression(paste("Mortality (10"^7, " cells mL"^{-1}, " h"^{-1}, ")"))
  )

ggsave("figures/fig2_mortality_mean.png", p_fig2, width = 6, height = 14, dpi = 300)

message("All figures saved to figures/.")
