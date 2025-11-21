#  load the calculated lambdas from cluster
lambda_popsize_bootpara_c1 <- read_csv("data/derived/IPM_bootstrappedLambda_para_cluster_c1.csv") %>% # census 1 cluster run
  dplyr::select(-sp_seed_to_sc1, -1) # delete typo column
lambda_popsize_bootpara_c2 <- read_csv("data/derived/IPM_bootstrappedLambda_para_cluster_c2.csv") %>% # census 2 cluster run
  dplyr::select(-sp_seed_to_sc1, -1)

# separate the name and extract lambda data
lambda_bootpara_c1 <- lambda_popsize_bootpara_c1 %>%
  dplyr::select(-required_sc1, -required_seeds, -gT, -p_seed_to_sc1,) %>%
  mutate(site = str_sub(parameter_name, 8, 9),
         species = str_sub(parameter_name, 1, 6),
         treat_combi = str_sub(parameter_name, 11, 19),
         treat_comp = str_sub(treat_combi, 1, 4),
         treat_warm = str_sub(treat_combi, 6, 9),
         bootstrap = str_sub(parameter_name, 21, nchar(parameter_name)),
         bootstrap = as.integer(str_extract(bootstrap, "\\d+")))
lambda_bootpara_c2 <- lambda_popsize_bootpara_c2 %>%
  dplyr::select(-required_sc1, -required_seeds, gT, p_seed_to_sc1) %>%
  mutate(site = str_sub(parameter_name, 8, 9),
         species = str_sub(parameter_name, 1, 6),
         treat_combi = str_sub(parameter_name, 11, 19),
         treat_comp = str_sub(treat_combi, 1, 4),
         treat_warm = str_sub(treat_combi, 6, 9),
         bootstrap = str_sub(parameter_name, 21, nchar(parameter_name)),
         bootstrap = as.integer(str_extract(bootstrap, "\\d+")))


### calculate CIs --------------------------------------------------------------

### For bias-corrected CIs you will also need the "real" lambda values, so either run the whole code to calculate them or load the separate file.

# calculate uncorrected confidence intervals
lambda.allcombis.bootpara_c1 <- lambda_bootpara_c1 %>%
  group_by(site, species, treat_combi, treat_warm, treat_comp) %>%
  summarize(boot_mean = mean(lambda),
            lower_CI = quantile(lambda, 0.025),
            upper_CI = quantile(lambda, 0.975)) %>%
  ungroup() %>%
  left_join(lambda.allcombis_c1, by = c("site", "species", "treat_combi", "treat_warm", "treat_comp"))
lambda.allcombis.bootpara_c2 <- lambda_bootpara_c2 %>%
  group_by(site, species, treat_combi, treat_warm, treat_comp) %>%
  summarize(boot_mean = mean(lambda),
            lower_CI = quantile(lambda, 0.025),
            upper_CI = quantile(lambda, 0.975)) %>%
  ungroup() %>%
  left_join(lambda.allcombis_c2, by = c("site", "species", "treat_combi", "treat_warm", "treat_comp"))


# calculate bias corrected confidence intervals of "real" lambdas
alpha <- 0.05

lambda.allcombis.bootpara.bcpi_c1 <- lambda_bootpara_c1 %>%
  group_by(species, site, treat_combi, treat_warm, treat_comp) %>%
  summarize(bootstrapped_values = list(lambda), .groups = 'drop') %>%
  ungroup() %>%
  left_join(lambda.allcombis_c1, by = c("species", "site", "treat_combi", "treat_comp", "treat_warm")) %>% ### REAL LAMBDAS NEEDED HERE! ###
  mutate(bcpi_CI = map2(bootstrapped_values, lambda, ~ bcpi(.y, .x, alpha)),  # .y is real_lambda, .x is t (bootstrapped values)
         lower_BCPI = sapply(bcpi_CI, `[[`, 1),  
         upper_BCPI = sapply(bcpi_CI, `[[`, 2),
         mean_bootstrap = sapply(bootstrapped_values, mean)) %>%
  ungroup() %>%
  mutate(treat_warm = str_sub(treat_combi, 6, 9),
         treat_comp = str_sub(treat_combi, 1,4),
         treat_warm_site = paste(site, treat_warm, sep = "_"), 
         site_treat_combi = factor(paste(site, treat_combi, sep = "_"), levels = c("lo_bare.ambi", "lo_vege.ambi", "hi_bare.ambi", "hi_vege.ambi", "hi_bare.warm", "hi_vege.warm")),
         treat_comp_trend = paste(treat_comp, trend, sep = "_")) %>%
  ungroup()
lambda.allcombis.bootpara.bcpi_c2 <- lambda_bootpara_c2 %>%
  group_by(species, site, treat_combi, treat_warm, treat_comp) %>%
  summarize(bootstrapped_values = list(lambda), .groups = 'drop') %>%
  ungroup() %>%
  left_join(lambda.allcombis_c2, by = c("species", "site", "treat_combi", "treat_comp", "treat_warm")) %>%
  mutate(bcpi_CI = map2(bootstrapped_values, lambda, ~ bcpi(.y, .x, alpha)),  # .y is real_lambda, .x is t (bootstrapped values)
         lower_BCPI = sapply(bcpi_CI, `[[`, 1),  
         upper_BCPI = sapply(bcpi_CI, `[[`, 2),
         mean_bootstrap = sapply(bootstrapped_values, mean)) %>%
  ungroup() %>%
  mutate(treat_warm = str_sub(treat_combi, 6, 9),
         treat_comp = str_sub(treat_combi, 1,4),
         treat_warm_site = paste(site, treat_warm, sep = "_"), 
         site_treat_combi = factor(paste(site, treat_combi, sep = "_"), levels = c("lo_bare.ambi", "lo_vege.ambi", "hi_bare.ambi", "hi_vege.ambi", "hi_bare.warm", "hi_vege.warm")),
         treat_comp_trend = paste(treat_comp, trend, sep = "_")) %>%
  ungroup()

# get averages and SDs of PGR for the different treatments
mean(lambda.allcombis_c1[lambda.allcombis_c1$site == "hi" & lambda.allcombis_c1$treat_combi == "bare.ambi",]$lambda)
sd(lambda.allcombis_c1[lambda.allcombis_c1$site == "hi" & lambda.allcombis_c1$treat_combi == "bare.ambi",]$lambda)

mean(lambda.allcombis_c1[lambda.allcombis_c1$site == "hi" & lambda.allcombis_c1$treat_combi == "bare.warm",]$lambda)
sd(lambda.allcombis_c1[lambda.allcombis_c1$site == "hi" & lambda.allcombis_c1$treat_combi == "bare.warm",]$lambda)

mean(lambda.allcombis_c1[lambda.allcombis_c1$site == "lo" & lambda.allcombis_c1$treat_combi == "bare.ambi",]$lambda)
sd(lambda.allcombis_c1[lambda.allcombis_c1$site == "lo" & lambda.allcombis_c1$treat_combi == "bare.ambi",]$lambda)

mean(lambda.allcombis_c1[lambda.allcombis_c1$site == "hi" & lambda.allcombis_c1$treat_combi == "vege.ambi",]$lambda)
sd(lambda.allcombis_c1[lambda.allcombis_c1$site == "hi" & lambda.allcombis_c1$treat_combi == "vege.ambi",]$lambda)

mean(lambda.allcombis_c1[lambda.allcombis_c1$site == "hi" & lambda.allcombis_c1$treat_combi == "vege.warm",]$lambda)
sd(lambda.allcombis_c1[lambda.allcombis_c1$site == "hi" & lambda.allcombis_c1$treat_combi == "vege.warm",]$lambda)



### plotting -------------------------------------------------------------------

# plot lambdas separately for species

# plot lambdas with BCPI and add mean bootstrapped lambda as black dot
ggplot(data = lambda.allcombis.bootpara.bcpi_c1, aes(x = as.factor(site_treat_combi), y = lambda)) + # , col = site
  geom_point(aes(shape = trend), size = 4, col = "red") + 
  geom_point(aes(x = as.factor(site_treat_combi), y = mean_bootstrap), col = "black") +
  geom_errorbar(aes(ymin = lower_BCPI, ymax = upper_BCPI), width = 0.2) +
  facet_wrap(~species, scales = "free") + # , scales = "free"
  #ylim(0, 3) +
  geom_abline(intercept = 1, slope = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "site - treatment combinations", y = "lambda") +
  ggtitle("Parametric bootstrapping (BCPI)") 


# bootstrap mean and "real" lambdas are quite similar

# prepare for plotting

# get decreasing lambdas for hi vege.ambi
a <- lambda.allcombis_c1 %>%
  filter(site == "hi" & treat_combi == "vege.ambi")

# re-order species for plotting
sp_order <- c("daucar", "medlup", "cenjac", "hypper", "plamed", "broere", "salpra", "brapin", "silvul",  "scacol")
lambda.allcombis.bootpara.bcpi_c1$species <- factor(lambda.allcombis.bootpara.bcpi_c1$species, levels = sp_order)
lambda.allcombis.bootpara_c1$species <- factor(lambda.allcombis.bootpara_c1$species, levels = sp_order)

lambda.allcombis.bootpara.bcpi_c2$species <- factor(lambda.allcombis.bootpara.bcpi_c2$species, levels = sp_order)
lambda.allcombis.bootpara_c2$species <- factor(lambda.allcombis.bootpara_c2$species, levels = sp_order)

# add an indication whether lambda is significantly different from 0, make CIs NA for non-significant lambdas with CIs smaller than the points 
# (just so the errorbars don't get printed: they are smaller than the points and show through because of reduced alpha), generate facting 
# variable to get daucar its own facet
lambda.allcombis.bootpara.bcpi_c1 <- lambda.allcombis.bootpara.bcpi_c1 %>%
  mutate(lambda_sig = if_else(lower_BCPI > 1 | upper_BCPI < 1, "sig.", "n-sig."),
         treat_comp_fact = factor(treat_comp, levels = c("bare", "vege")) %>% 
           relevel("vege"))
lambda.allcombis.bootpara_c1 <- lambda.allcombis.bootpara_c1 %>%
  mutate(lambda_sig = if_else(lower_CI > 1 | upper_CI < 1, "sig.", "n-sig."))

lambda.allcombis.bootpara.bcpi_c2 <- lambda.allcombis.bootpara.bcpi_c2 %>%
  mutate(lambda_sig = if_else(lower_BCPI > 1 | upper_BCPI < 1, "sig.", "n-sig."),
         treat_comp_fact = factor(treat_comp, levels = c("bare", "vege")) %>% 
           relevel("vege"))
lambda.allcombis.bootpara_c2 <- lambda.allcombis.bootpara_c2 %>%
  mutate(lambda_sig = if_else(lower_CI > 1 | upper_CI < 1, "sig.", "n-sig."))

# define position_dodge
posd <- position_dodge(width = 0.7) # position_dodge2(width = 0.7, preserve = "single")


# Fig 2: ambient high site, census 1
png("plots/Fig2_PGRhiambi_census1_20251027.png", width = 17, height = 10, units="cm", res=800)

# use upper_BCPI_red here to prevent the error bars showing through for points with reduced opacity --> but only if the errorbars are shorter than the points (manually defined above)
ggplot(data = lambda.allcombis.bootpara.bcpi_c1[lambda.allcombis.bootpara.bcpi_c1$treat_warm == "ambi" & lambda.allcombis.bootpara.bcpi_c1$site == "hi",],
       aes(x = species, y = lambda, col = site_treat_combi, alpha = lambda_sig)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_BCPI, ymax = upper_BCPI, group = treat_comp_fact), width = 0.2, position = posd) +
  geom_point(aes(group = treat_comp_fact), position = posd, size = 4, col = "white", fill = "white", alpha = 1) +
  geom_point(aes(shape = trend, group = treat_comp_fact), position = posd, size = 4, fill = "white", stroke = 1.5) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, face = "italic"),
    #axis.ticks = element_blank(),
    legend.background = element_blank(), 
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    strip.background = element_blank(),
    legend.position = c(0.82, 0.95),  # move legend to top-right corner
    legend.justification = c(1, 1),
    legend.box = "vertical",  # combine legends vertically into a single box
    legend.box.background = element_rect(fill = "white", color = "black"),
    #legend.key.size = unit(1.5, "lines"),
    #legend.key.width = unit(1.5, "lines"),
    legend.spacing = unit(-0.5, "lines")
  ) +
  labs(y = expression(paste("Population growth rate ", lambda))) +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_x_discrete(labels = function(x) species_names[x]) +
  scale_color_manual(
    values = treat_combi_site_col,
    name = "Competition Treatment",  # Legend title
    labels = c("hi_bare.ambi" = "without competition", "hi_vege.ambi" = "with competition"),
    guide = guide_legend(order = 1)) +
  scale_shape_manual(
    values = c(19, 21),
    labels = c("grow" = expression(paste(lambda, " > 1")), "shrink" = expression(paste(lambda, " < 1"))),
    guide = guide_legend(order = 2)) +
  guides(alpha = "none") #+
#facet_wrap(~daucar_facet, scales = "free")

dev.off()

# Fig S1: ambient high site, census 2
png("plots/FigS1_a_PGRhiambi_census2_20251027.png", width = 17, height = 10, units="cm", res=800)

# use upper_BCPI_red here to prevent the error bars showing through for points with reduced opacity --> but only if the errorbars are shorter than the points (manually defined above)
ggplot(data = lambda.allcombis.bootpara.bcpi_c2[lambda.allcombis.bootpara.bcpi_c2$treat_warm == "ambi" & lambda.allcombis.bootpara.bcpi_c2$site == "hi",],
       aes(x = species, y = lambda, col = site_treat_combi, alpha = lambda_sig)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_BCPI, ymax = upper_BCPI, group = treat_comp_fact), width = 0.2, position = posd) +
  geom_point(aes(group = treat_comp_fact), position = posd, size = 4, col = "white", fill = "white", alpha = 1) +
  geom_point(aes(shape = trend, group = treat_comp_fact), position = posd, size = 4, fill = "white", stroke = 1.5) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, face = "italic"),
    #axis.ticks = element_blank(),
    legend.background = element_blank(), 
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    strip.background = element_blank(),
    legend.position = c(0.82, 0.95),  # move legend to top-right corner
    legend.justification = c(1, 1),
    legend.box = "vertical",  # combine legends vertically into a single box
    legend.box.background = element_rect(fill = "white", color = "black"),
    #legend.key.size = unit(1.5, "lines"),
    #legend.key.width = unit(1.5, "lines"),
    legend.spacing = unit(-0.5, "lines")
  ) +
  labs(y = expression(paste("Population growth rate ", lambda))) +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_x_discrete(labels = function(x) species_names[x]) +
  scale_color_manual(
    values = treat_combi_site_col,
    name = "Competition Treatment",  # Legend title
    labels = c("hi_bare.ambi" = "without competition", "hi_vege.ambi" = "with competition"),
    guide = guide_legend(order = 1)) +
  scale_shape_manual(
    values = c(19, 21),
    labels = c("grow" = expression(paste(lambda, " > 1")), "shrink" = expression(paste(lambda, " < 1"))),
    guide = guide_legend(order = 2)) +
  guides(alpha = "none") #+
#facet_wrap(~daucar_facet, scales = "free")

dev.off()


# Fig X: warmed high site, census 1 (combined with lo site ambient later)
hi_warm_c1 <- ggplot(data = lambda.allcombis.bootpara.bcpi_c1[lambda.allcombis.bootpara.bcpi_c1$treat_warm == "warm" & lambda.allcombis.bootpara.bcpi_c1$site == "hi",],
                     aes(x = species, y = lambda, col = site_treat_combi, alpha = lambda_sig)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_BCPI, ymax = upper_BCPI, group = treat_comp_fact), position = posd, width = 0.2, alpha = 1) +
  geom_point(aes(group = treat_comp_fact), position = posd, shape = 21, size = 6,
             fill = "white", colour = "white", alpha = 1, stroke = 0) +
  geom_point(aes(shape = trend, group = treat_comp_fact), position = posd, size = 4, fill = "white", stroke = 1.5) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_blank(), 
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    strip.background = element_blank(),
    legend.position = c(0.82, 0.95),  # move legend to top-right corner
    legend.justification = c(1, 1),
    legend.box = "vertical",  # combine legends vertically into a single box
    legend.box.background = element_rect(fill = "white", color = "black"),
    #legend.key.size = unit(1.5, "lines"),
    #legend.key.width = unit(1.5, "lines"),
    legend.spacing = unit(-0.5, "lines")
  ) +
  labs(y = expression(paste("Population growth rate ", lambda))) +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_x_discrete(labels = function(x) species_names[x]) +
  scale_color_manual(
    values = treat_combi_site_col,
    labels = c("hi_bare.warm" = "without competition", "hi_vege.warm" = "with competition"),
    guide = guide_legend(order = 1)) +
  scale_shape_manual(
    values = c(19, 21),
    labels = c("grow" = expression(paste(lambda, " > 1")), "shrink" = expression(paste(lambda, " < 1"))),
    guide = guide_legend(order = 2)) +
  guides(alpha = "none") # , shape = "none"

#dev.off()

# Fig X: warmed high site, census 2 
hi_warm_c2 <- ggplot(data = lambda.allcombis.bootpara.bcpi_c2[lambda.allcombis.bootpara.bcpi_c2$treat_warm == "warm" & lambda.allcombis.bootpara.bcpi_c2$site == "hi",],
                     aes(x = species, y = lambda, col = site_treat_combi, alpha = lambda_sig, shape = trend)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_BCPI, ymax = upper_BCPI, group = treat_comp_fact), position = posd, width = 0.2, alpha = 1) +
  geom_point(aes(group = treat_comp_fact), position = posd, shape = 21, size = 6,
             fill = "white", colour = "white", alpha = 1, stroke = 0) +
  geom_point(aes(shape = trend, group = treat_comp_fact), position = posd, size = 4, fill = "white", stroke = 1.5) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_blank(), 
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    strip.background = element_blank(),
    legend.position = c(0.82, 0.95),  # move legend to top-right corner
    legend.justification = c(1, 1),
    legend.box = "vertical",  # combine legends vertically into a single box
    legend.box.background = element_rect(fill = "white", color = "black"),
    #legend.key.size = unit(1.5, "lines"),
    #legend.key.width = unit(1.5, "lines"),
    legend.spacing = unit(-0.5, "lines")
  ) +
  labs(y = expression(paste("Population growth rate ", lambda))) +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_x_discrete(labels = function(x) species_names[x]) +
  scale_color_manual(
    values = treat_combi_site_col,
    labels = c("hi_bare.warm" = "without competition", "hi_vege.warm" = "with competition"),
    guide = guide_legend(order = 1)) +
  scale_shape_manual(
    values = c(19, 21),
    labels = c("grow" = expression(paste(lambda, " > 1")), "shrink" = expression(paste(lambda, " < 1"))),
    guide = guide_legend(order = 2)) +
  guides(alpha = "none") # , shape = "none"

# Fig X: ambient low site, census 1 (combined with hi site warm later)
lo_ambi_c1 <- ggplot(data = lambda.allcombis.bootpara.bcpi_c1[lambda.allcombis.bootpara.bcpi_c1$treat_warm == "ambi" & lambda.allcombis.bootpara.bcpi_c1$site == "lo",],
                     aes(x = species, y = lambda, col = site_treat_combi, alpha = lambda_sig, shape = trend)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_BCPI, ymax = upper_BCPI, group = treat_comp_fact), position = posd, width = 0.2, alpha = 1) +
  geom_point(aes(group = treat_comp_fact), position = posd, shape = 21, size = 6,
             fill = "white", colour = "white", alpha = 1, stroke = 0) +
  geom_point(aes(shape = trend, group = treat_comp_fact), position = posd, size = 4, fill = "white", stroke = 1.5) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, face = "italic"),
    #axis.ticks = element_blank(),
    legend.background = element_blank(), 
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    strip.background = element_blank(),
    legend.position = c(0.82, 0.95),  # move legend to top-right corner
    legend.justification = c(1, 1),
    legend.box = "vertical",  # combine legends vertically into a single box
    legend.box.background = element_rect(fill = "white", color = "black"),
    #legend.key.size = unit(1.5, "lines"),
    #legend.key.width = unit(1.5, "lines"),
    legend.spacing = unit(-0.5, "lines")
  ) +
  labs(y = expression(paste("Population growth rate ", lambda))) +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_x_discrete(labels = function(x) species_names[x]) +
  scale_color_manual(
    values = treat_combi_site_col,
    name = "Competition Treatment",  # Legend title
    labels = c("lo_bare.ambi" = "without competition", "lo_vege.ambi" = "with competition"),
    guide = guide_legend(order = 1)) +
  scale_shape_manual(
    values = c(19, 21),
    labels = c("grow" = expression(paste(lambda, " > 1")), "shrink" = expression(paste(lambda, " < 1"))),
    guide = guide_legend(order = 2)) +
  guides(alpha = "none", shape = "none") +
  ylim(0, 30)

# Fig X: ambient lo site, census 2 
lo_ambi_c2 <- ggplot(data = lambda.allcombis.bootpara.bcpi_c2[lambda.allcombis.bootpara.bcpi_c2$treat_warm == "ambi" & lambda.allcombis.bootpara.bcpi_c2$site == "lo",],
                     aes(x = species, y = lambda, col = site_treat_combi, alpha = lambda_sig, shape = trend)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_BCPI, ymax = upper_BCPI, group = treat_comp_fact), position = posd, width = 0.2, alpha = 1) +
  geom_point(aes(group = treat_comp_fact), position = posd, shape = 21, size = 6,
             fill = "white", colour = "white", alpha = 1, stroke = 0) +
  geom_point(aes(shape = trend, group = treat_comp_fact), position = posd, size = 4, fill = "white", stroke = 1.5) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, face = "italic"),
    #axis.ticks = element_blank(),
    legend.background = element_blank(), 
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    strip.background = element_blank(),
    legend.position = c(0.82, 0.95),  # move legend to top-right corner
    legend.justification = c(1, 1),
    legend.box = "vertical",  # combine legends vertically into a single box
    legend.box.background = element_rect(fill = "white", color = "black"),
    #legend.key.size = unit(1.5, "lines"),
    #legend.key.width = unit(1.5, "lines"),
    legend.spacing = unit(-0.5, "lines")
  ) +
  labs(y = expression(paste("Population growth rate ", lambda))) +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_x_discrete(labels = function(x) species_names[x]) +
  scale_color_manual(
    values = treat_combi_site_col,
    name = "Competition Treatment",  # Legend title
    labels = c("lo_bare.ambi" = "without competition", "lo_vege.ambi" = "with competition"),
    guide = guide_legend(order = 1)) +
  scale_shape_manual(
    values = c(19, 21),
    labels = c("grow" = expression(paste(lambda, " > 1")), "shrink" = expression(paste(lambda, " < 1"))),
    guide = guide_legend(order = 2)) +
  guides(alpha = "none", shape = "none") +
  ylim(0, 30)


# plot 4: combined warmed and lo site plots

# combine the plots 

png("plots/FigSX_PGRhiwarmlo_census1_20251027.png", width = 17, height = 17, units="cm", res=800)

hi_warm_c1 + 
  theme(plot.tag = element_text(face = "bold", size = 16)) +
  lo_ambi_c1 + 
  plot_layout(ncol = 1) +   
  plot_annotation(tag_levels = "A") + 
  theme(plot.tag = element_text(face = "bold", size = 16))

dev.off()

png("plots/Fig3_c_PGRlo_census1_20251027.png", , width = 17, height = 10, units="cm", res=800)

lo_ambi_c1 

dev.off()

png("plots/FigS1_bc_PGRhiwarmlo_census2_20251027.png", width = 17, height = 17, units="cm", res=800)

hi_warm_c2 + 
  theme(plot.tag = element_text(face = "bold", size = 16)) +
  lo_ambi_c2 + 
  plot_layout(ncol = 1) +   
  plot_annotation(tag_levels = "A") + 
  theme(plot.tag = element_text(face = "bold", size = 16))

dev.off()

