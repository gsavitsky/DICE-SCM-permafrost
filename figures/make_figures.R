# Figures from data/outputs (run from this folder: Rscript make_figures.R). Needs readr, dplyr, tidyr, ggplot2.
# Published comparison values are entered below with their sources.
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(ggplot2)})
OUTD <- "../data/outputs"; IND <- "../data/inputs"; FIG <- "."
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"; BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; RED <- "#e34948"; GREY <- "#8a8984"
th <- theme_minimal(base_size = 16) +
  theme(text = element_text(colour = INK), axis.text = element_text(colour = INK2), panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.4), strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(face = "bold", size = 18), plot.subtitle = element_text(colour = INK2, size = 13),
        plot.caption = element_text(colour = INK2, size = 11, hjust = 0), legend.position = "bottom",
        plot.title.position = "plot", plot.caption.position = "plot", plot.margin = margin(10, 18, 8, 10),
        plot.background = element_rect(fill = "white", colour = NA))
save <- function(p, f, w = 10, h = 5.625) ggsave(file.path(FIG, f), p, width = w, height = h, dpi = 300, bg = "white")
rd <- function(f, dir = OUTD) read_csv(file.path(dir, f), show_col_types = FALSE)
MODCOL <- c("RCMIP / MAGICC" = INK, "DFAIR" = BLUE, "SCM" = ORANGE)

# ---- 1. carbon-cycle feedbacks vs CMIP6 ----
A <- rd("arora2020_cmip6_4xco2.csv", IND) %>%
  mutate(beta_total = beta_land + beta_ocean, gamma_total = gamma_land + gamma_ocean, group = "CMIP6 Earth-system model")
C <- rd("c4mip_feedbacks.csv") %>% mutate(group = model, model = ifelse(model == "DFAIR", "DFAIR (DICE 2023)", "SCM"))
pars <- c("beta_land", "beta_ocean", "beta_total", "gamma_land", "gamma_ocean", "gamma_total")
L <- bind_rows(A, C) %>% select(model, group, all_of(pars)) %>% pivot_longer(all_of(pars), names_to = "par") %>%
  mutate(par = factor(par, levels = pars, labels = c("Land β (PgC/ppm)", "Ocean β (PgC/ppm)", "Total β (PgC/ppm)",
                                                     "Land γ (PgC/°C)", "Ocean γ (PgC/°C)", "Total γ (PgC/°C)")),
         model = factor(model, levels = rev(c("SCM", "DFAIR (DICE 2023)", sort(A$model)))))
p <- ggplot(L, aes(value, model, colour = group, shape = group, size = group)) + geom_point() +
  facet_wrap(~par, scales = "free_x", nrow = 2) +
  scale_colour_manual(values = c("CMIP6 Earth-system model" = INK2, "SCM" = ORANGE, "DFAIR" = BLUE), name = NULL) +
  scale_shape_manual(values = c("CMIP6 Earth-system model" = 16, "SCM" = 17, "DFAIR" = 15), name = NULL) +
  scale_size_manual(values = c("CMIP6 Earth-system model" = 2.3, "SCM" = 3.8, "DFAIR" = 3.6), name = NULL) +
  labs(x = NULL, y = NULL, title = "Carbon-cycle feedbacks: SCM and DFAIR vs 11 CMIP6 models", subtitle = "C4MIP 1%/yr CO2 to 4×CO2",
       caption = "CMIP6: Arora et al. (2020), Table A1. DFAIR land = biosphere reservoir; ocean = deep-ocean + mixed-layer reservoirs.") +
  th + theme(axis.text.y = element_text(size = 10), strip.text = element_text(size = 12.5))
save(p, "fig01_feedbacks_vs_cmip6.png", h = 6.4)

# ---- 2-3. SSP CO2 rise vs RCMIP/MAGICC ----
S <- rd("ssp_emissions_driven.csv") %>%
  mutate(scen = recode(scenario, ssp119 = "SSP1-1.9", ssp126 = "SSP1-2.6", ssp245 = "SSP2-4.5", ssp370 = "SSP3-7.0", ssp585 = "SSP5-8.5"),
         who = factor(case_when(grepl("RCMIP", model) ~ "RCMIP / MAGICC", grepl("^DFAIR", model) ~ "DFAIR", TRUE ~ "SCM"), levels = names(MODCOL)),
         perma = grepl("permafrost", model))
ssp_bars <- function(d, title) {
  ggplot(d, aes(scen, d2100, fill = who)) + geom_col(position = position_dodge(0.8), width = 0.75) + geom_hline(yintercept = 0, colour = INK) +
    geom_text(aes(label = round(d2100), vjust = ifelse(d2100 >= 0, -0.4, 1.3)), position = position_dodge(0.8), size = 3.6, colour = INK) +
    scale_fill_manual(values = MODCOL, name = NULL) + scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    labs(x = NULL, y = "CO2 rise 2020–2100 (ppm)", title = title,
         caption = "Same RCMIP v5.1.0 CO2 emissions and non-CO2 forcing in every model; DICE 2023 two-box climate.") +
    th + theme(panel.grid.major.x = element_blank())
}
save(ssp_bars(S %>% filter(who == "RCMIP / MAGICC" | perma), "Emissions-driven SSPs: CO2 rise, carbon cycles with permafrost"), "fig02_ssp_co2_with_permafrost.png")
save(ssp_bars(S %>% filter(who == "RCMIP / MAGICC" | !perma), "Emissions-driven SSPs: CO2 rise, carbon cycles without permafrost"), "fig03_ssp_co2_no_permafrost.png")

# ---- 4. SSP warming vs AR6 ----
# AR6 WG1 SPM Table SPM.1, 2081-2100 vs 1850-1900: best estimate and very likely range
ar6 <- tibble(scen = c("SSP1-1.9", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"), best = c(1.4, 1.8, 2.7, 3.6, 4.4),
              lo = c(1.0, 1.3, 2.1, 2.8, 3.3), hi = c(1.8, 2.4, 3.5, 4.6, 5.7), who = "AR6 assessed (2081–2100)")
W <- bind_rows(ar6, S %>% filter(perma, who != "RCMIP / MAGICC") %>% transmute(scen, best = T2100, who = paste0(as.character(who), " + permafrost (2100)"))) %>%
  mutate(who = factor(who, levels = c("AR6 assessed (2081–2100)", "DFAIR + permafrost (2100)", "SCM + permafrost (2100)")))
p <- ggplot(W, aes(scen, best, colour = who, shape = who)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.9, position = position_dodge(0.5), na.rm = TRUE) +
  geom_point(size = 3.8, position = position_dodge(0.5)) +
  scale_colour_manual(values = c(INK, BLUE, ORANGE), name = NULL) + scale_shape_manual(values = c(15, 16, 17), name = NULL) +
  labs(x = NULL, y = "Warming above pre-industrial (°C)", title = "Warming under emissions-driven SSPs vs AR6",
       caption = "AR6 point = best estimate, whiskers = very likely range (AR6 WG1 SPM Table SPM.1, 2081–2100 vs 1850–1900). Models: single year 2100.") +
  th + theme(panel.grid.major.x = element_blank())
save(p, "fig04_ssp_warming_ar6.png")

# ---- 5-7. permafrost vs published estimates ----
CF <- rd("configurations.csv")
pgc <- function(pat) CF$PgC_per_degC_2100[CF$configuration == pat]
dot <- function(d, xl, title, cap) {
  d <- d %>% mutate(lbl = factor(lbl, levels = rev(lbl)))
  ggplot(d, aes(v, lbl, colour = kind, shape = kind)) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 1, na.rm = TRUE) + geom_point(size = 4) +
    scale_colour_manual(values = c("Published estimate" = INK, "This study" = ORANGE), name = NULL) +
    scale_shape_manual(values = c("Published estimate" = 15, "This study" = 16), name = NULL) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(x = xl, y = NULL, title = title, caption = cap) + th + theme(axis.text.y = element_text(size = 13.5, colour = INK))
}
lt <- function(lbl, v, lo = NA, hi = NA, kind) tibble(lbl = lbl, v = v, lo = lo, hi = hi, kind = kind)
p <- dot(bind_rows(lt("AR6 (Canadell et al. 2021, Box 5.1)", 18, 3.1, 41, "Published estimate"),
                   lt("Schädel et al. 2026: gradual thaw", 28, 10, 36, "Published estimate"),
                   lt("Schädel et al. 2026: thaw + abrupt + fire", 63, 35, 73, "Published estimate"),
                   lt("This study, DFAIR: gradual only", pgc("DICE 2023 + permafrost (gradual only)"), kind = "This study"),
                   lt("This study, DFAIR: gradual + abrupt", pgc("DICE 2023 + permafrost"), kind = "This study"),
                   lt("This study, SCM: gradual only", pgc("DICE-CP (SCM + permafrost, gradual only)"), kind = "This study"),
                   lt("This study, SCM: gradual + abrupt", pgc("DICE-CP (SCM + permafrost)"), kind = "This study")),
         "PgC released per °C of global warming by 2100", "Permafrost sensitivity",
         "AR6: assessed estimate, 5–95% range (CO2 only). Schädel et al. (2026): mean, 25th–75th percentile.\nThis study: Hector-strength permafrost, DICE 2023 optimal path.")
save(p, "fig05_permafrost_sensitivity.png")
rel <- S %>% filter(perma, who != "RCMIP / MAGICC") %>% group_by(scen) %>% summarise(v = mean(release2100), lo = min(release2100), hi = max(release2100))
p <- dot(bind_rows(lt("Schneider von Deimling et al. 2015: RCP2.6", 39, 20, 58, "Published estimate"),
                   lt("Schneider von Deimling et al. 2015: RCP4.5", 60, 28, 92, "Published estimate"),
                   lt("SROCC: RCP8.5 best estimate (± 1 s.d.)", 92, 75, 109, "Published estimate"),
                   lt("SROCC: RCP8.5 maximum", 240, kind = "Published estimate"),
                   lt("Schuur et al. 2022: high, net Arctic CO2 (RCP8.5)", 149, kind = "Published estimate"),
                   lt("Woodard et al. 2021 (Hector): RCP2.6", 100.9, kind = "Published estimate"),
                   lt("Woodard et al. 2021 (Hector): RCP4.5", 120.6, kind = "Published estimate"),
                   lt("Woodard et al. 2021 (Hector): RCP8.5", 142.3, kind = "Published estimate"),
                   rel %>% transmute(lbl = paste("This study:", scen), v, lo, hi, kind = "This study"),
                   lt("This study: DICE optimal path", mean(CF$release_2100[CF$configuration %in% c("DICE 2023 + permafrost", "DICE-CP (SCM + permafrost)")]), kind = "This study")),
         "Cumulative carbon released by 2100 (PgC)", "Permafrost carbon released by 2100",
         "Published: AR6 WG1 Box 5.1 (Schneider von Deimling ranges, point = midpoint); Schuur et al. (2022) Table 1, 2000–2099;\nWoodard et al. (2021) Table 3, Hector net permafrost CO2 by 2100 since pre-industrial.\nThis study: permafrost CO2 released 2020–2100, mean of DFAIR and SCM (whiskers span the two).")
save(p, "fig06_permafrost_release_2100.png")
Z <- rd("zero_emissions_commitment.csv")
zrow <- function(col, yrs, med, lo, hi) bind_rows(
  lt("MacDougall 2021: ZECMIP-style (1000 PgC)", med, lo, hi, "Published estimate"),
  lt("This study: DFAIR", Z[[col]][Z$carbon == "DFAIR"], kind = "This study"),
  lt("This study: SCM", Z[[col]][Z$carbon == "SCM"], kind = "This study")) %>% mutate(when = paste(yrs, "years after emissions stop"))
ZZ <- bind_rows(zrow("permafrost_ZEC50", 50, 0.06, 0.02, 0.14), zrow("permafrost_ZEC100", 100, 0.09, 0.04, 0.21)) %>%
  mutate(when = factor(when, levels = c("50 years after emissions stop", "100 years after emissions stop")), lbl = factor(lbl, levels = rev(unique(lbl))))
p <- ggplot(ZZ, aes(v, lbl, colour = kind, shape = kind)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 1, na.rm = TRUE) + geom_point(size = 4) + facet_wrap(~when) +
  scale_colour_manual(values = c("Published estimate" = INK, "This study" = ORANGE), name = NULL) +
  scale_shape_manual(values = c("Published estimate" = 15, "This study" = 16), name = NULL) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Extra warming from permafrost after emissions stop (°C)", y = NULL, title = "Permafrost contribution to the zero-emissions commitment",
       caption = "MacDougall (2021), Table 2: median [5th–95th]; ZECMIP experiment from pre-industrial, which excludes thaw from earlier warming.\nHe adds 0.04 [0–0.06] °C at 50 yr for that historical \"thaw lag\"; this study starts in 2020 and includes it (1000 PgC total, emissions stop ~2049).") +
  th + theme(axis.text.y = element_text(size = 13.5, colour = INK))
save(p, "fig07_permafrost_zec.png")

# ---- 8-9. SCC by feedback ----
nm <- c(ocean = "Ocean buffer saturation", oceanT = "Ocean warming", fert = "CO2 fertilisation", photoT = "Photosynthesis warming response",
        vegT = "Plant respiration", soilT = "Soil respiration", perma = "Permafrost")
SV <- rd("shapley_values.csv")
scc_bars <- function(metric, title) {
  d <- SV %>% filter(metric == !!metric)
  B <- bind_rows(d %>% transmute(lab = nm[process], v = shapley) %>% arrange(v), tibble(lab = "Net: all feedbacks", v = d$all_on[1] - d$all_off[1])) %>%
    mutate(type = case_when(lab == "Net: all feedbacks" ~ "Net", v > 0 ~ "Reinforcing (raises SCC)", TRUE ~ "Balancing (lowers SCC)"), lab = factor(lab, levels = lab))
  rng <- max(abs(B$v)) * 1.25
  ggplot(B, aes(v, lab, fill = type)) + geom_col(width = 0.68) + geom_vline(xintercept = 0, colour = INK, linewidth = 0.6) +
    geom_text(aes(x = v + sign(v) * rng * 0.02, label = sprintf("%+.1f", v), hjust = ifelse(v >= 0, 0, 1)), colour = INK, size = 4.8) +
    scale_fill_manual(values = c("Reinforcing (raises SCC)" = RED, "Balancing (lowers SCC)" = BLUE, "Net" = GREY),
                      breaks = c("Reinforcing (raises SCC)", "Balancing (lowers SCC)"), name = NULL) +
    scale_x_continuous(limits = c(-rng, rng)) +
    labs(x = "Contribution to the SCC ($ per tCO2)", y = NULL, title = title,
         subtitle = sprintf("DICE-CP. SCC without carbon-cycle feedbacks $%.0f; with all of them $%.0f.", d$all_off[1], d$all_on[1]),
         caption = "Shapley values over all 128 on/off combinations of the seven processes; they sum to the net bar. DICE 2023 optimal controls, welfare to 2200.") +
    th + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 13.5, colour = INK))
}
save(scc_bars("SCC2020", "What each carbon-cycle feedback adds to the SCC: emission in 2020"), "fig08_scc_by_feedback_2020.png")
save(scc_bars("SCC2080", "What each carbon-cycle feedback adds to the SCC: emission in 2080"), "fig09_scc_by_feedback_2080.png")

# ---- 10. pairwise interactions (2020 SCC) ----
lv <- names(nm)
nm2 <- c(ocean = "Ocean buffer", oceanT = "Ocean warming", fert = "CO2 fertilisation", photoT = "Photosynthesis (T)",
         vegT = "Plant respiration", soilT = "Soil respiration", perma = "Permafrost")
I <- rd("shapley_interactions.csv") %>% filter(metric == "SCC2020")
d <- SV %>% filter(metric == "SCC2020")
M <- bind_rows(I %>% select(i, j, value = interaction, sign_consistent), I %>% select(i = j, j = i, value = interaction, sign_consistent),
               d %>% transmute(i = process, j = process, value = shapley, sign_consistent = TRUE)) %>%
  mutate(diag = i == j, i = factor(nm2[i], levels = nm2[lv]), j = factor(nm2[j], levels = rev(nm2[lv])))
p <- ggplot(M, aes(i, j, fill = value)) + geom_tile(colour = "white") + geom_tile(data = M %>% filter(diag), colour = "black", linewidth = 0.6) +
  geom_text(aes(label = ifelse(sign_consistent, sprintf("%+.1f", value), sprintf("%+.1f*", value)),
                fontface = ifelse(sign_consistent, "bold", "plain"), colour = sign_consistent), size = 3.8) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey45"), guide = "none") +
  scale_fill_gradient2(low = BLUE, mid = "white", high = RED, midpoint = 0, name = "$ per tCO2") +
  labs(x = NULL, y = NULL, title = "SCC 2020: own contribution (diagonal) and pairwise interactions",
       caption = "Off-diagonal: Shapley interaction index. * sign changes across the 32 settings of the other five processes (bold: same sign in all 32).") +
  theme_minimal(base_size = 13) + theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(),
                                        plot.caption = element_text(hjust = 0, colour = INK2), plot.title.position = "plot",
                                        plot.caption.position = "plot", plot.background = element_rect(fill = "white", colour = NA))
save(p, "fig10_scc_interactions.png", w = 10, h = 7.6)
cat("figures done\n")
