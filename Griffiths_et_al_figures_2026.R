# Code to reproduce figures for Griffiths et al.
#### PACKAGES ####
library(tidyverse); library(ggthemes); library(ggrepel); library(patchwork)
library(dplyr); library(tidyr); library(vegan); library(ggsci); library(cowplot)
library(rcartocolor); library(ggpubr); library(Hmisc); library(ggbreak)

#### READ & PREP ####
data <- readr::read_csv(
  "Macroplankton_micronekton_South Georgia_WCB_2009_2019.csv"
) %>%
  mutate(Cruise.Event.Net = stringr::str_c(Cruise, Event.Net, sep = "."))

abund_data <- data %>%
  dplyr::group_by(
    Cruise.Event.Net,
    Cruise,
    Event.Net,
    Net.type,
    Start.of.Event,
    Survey.Season,
    Species.code
  ) %>%
  dplyr::summarise(
    Abundm2 = sum(`Abund (m-2)`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols     = c(Cruise.Event.Net, Start.of.Event, Survey.Season),
    names_from  = Species.code,
    values_from = Abundm2,
    values_fill = 0,
    names_sort  = TRUE
  )

site_meta <- readr::read_csv(
  "WCB_meta_Griffiths_et_al_w_sunangle.csv"
)

site_keys <- site_meta %>%
  transmute(Cruise.Event.Net, Area, Site, Code, nautical_twilight) %>%
  distinct(Cruise.Event.Net, .keep_all = TRUE)

abund_data <- abund_data %>%
  left_join(site_keys, by = "Cruise.Event.Net") %>%
  filter(Site != 4) %>%
  relocate(Area:nautical_twilight, .after = Survey.Season)

abund_data <- abund_data %>%
  dplyr::rename(
    Fish_larvae     = Fish.l,
    Limacina_spp.   = Lima,
    Tomopteris_spp. = Tomop,
    T._gaudichaudii = Them.g,
    P._macropa      = Primn
  )

abund_long <- abund_data %>%
  pivot_longer(Amph:Vibi.a, names_to = "Species.code", values_to = "Abun") %>%
  filter(Abun != 0) %>%
  mutate(Code = stringr::str_c(stringr::str_sub(Survey.Season, -2), ".", Code))

# with E. superba
full <- abund_long %>%
  dplyr::count(Code, Area, wt = Abun, name = "Total_abun")
hist(full$Total_abun, breaks = 20)
full_mwu <- wilcox.test(Total_abun ~ Area, data = full, exact = FALSE)

# without E. superba
sub <- abund_long %>%
  dplyr::filter(Species.code != "E._superba") %>%
  dplyr::count(Code, Area, wt = Abun, name = "Total_abun")
hist(sub$Total_abun, breaks = 20)
sub_mwu <- wilcox.test(Total_abun ~ Area, data = sub, exact = FALSE)

area_labels <- c(Offshelf = "Off-shelf", Onshelf = "On-shelf")

#### FIGURE 2 ABUNDANCE ####
full$Area <- factor(full$Area, levels = c("Offshelf", "Onshelf"))
sub$Area  <- factor(sub$Area,  levels = c("Offshelf", "Onshelf"))

common_theme <- ggthemes::theme_few() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(margin = margin(t = 6, b = 6))
  )

ES_vio <- ggplot(full, aes(x = Area, y = Total_abun)) +
  geom_violin(trim = FALSE, aes(color = Area)) +
  geom_boxplot(width = 0.1) +
  scale_x_discrete(breaks = names(area_labels), labels = area_labels) +
  scale_color_npg(breaks = names(area_labels), labels = area_labels, guide = "none") +
  # IMPORTANT: do NOT use ylim() with ggbreak
  scale_y_continuous(limits = c(0, 850)) +
  ggbreak::scale_y_break(c(175, 375), scales = 0.25,
                         ticklabels = c(400, 500, 600, 700, 800)) +
  labs(
    y = expression("Abundance (individuals m"^-2*")"),
    title = expression(paste("A) With ", italic("E. superba")))
  ) +
  annotate("text", x = 2.2, y = 650,
           label = paste("p-value = ", round(full_mwu$p.value, 2))) +
  common_theme +
  theme(
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank()
  )

noES_vio <- ggplot(sub, aes(x = Area, y = Total_abun)) +
  geom_violin(trim = FALSE, aes(color = Area)) +
  geom_boxplot(width = 0.1) +
  scale_x_discrete(breaks = names(area_labels), labels = area_labels) +
  scale_color_npg(breaks = names(area_labels), labels = area_labels, guide = "none") +
  scale_y_continuous(limits = c(0, 150)) +
  labs(
    y = expression("Abundance (individuals m"^-2*")"),
    title = expression(paste("B) Without ", italic("E. superba")))
  ) +
  annotate("text", x = 2.3, y = 150,
           label = paste("p-value = ", round(sub_mwu$p.value, 2))) +
  common_theme

side_by_side <- cowplot::plot_grid(
  ES_vio, noES_vio,
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1)
)

side_by_side

# ggsave("Fig2_ABUND_sidebyside.png", side_by_side, width = 12, height = 6, dpi = 300)

#### FIGURE 3 DIVERSITY ####
# 95% CI helpers (z = 1.96)
upper_conf <- function(mean, sd, n) mean + (1.96 * (sd / sqrt(n)))
lower_conf <- function(mean, sd, n) mean - (1.96 * (sd / sqrt(n)))

# Site-level diversity metrics
site_div <- abund_long %>%
  dplyr::group_by(Site, Area, Survey.Season, Species.code) %>%
  dplyr::summarise(Abun = sum(Abun), .groups = "drop") %>%
  dplyr::group_by(Site, Area, Survey.Season) %>%
  dplyr::summarise(
    shannon  = vegan::diversity(Abun),
    richness = vegan::specnumber(Abun),
    .groups  = "drop"
  )

site_div$Area <- factor(site_div$Area, levels = c("Onshelf","Offshelf"))
all_years <- sort(unique(site_div$Survey.Season))
site_div  <- site_div %>% mutate(Survey.Season = factor(Survey.Season, levels = all_years))

# 95% CIs for Shannon
sh_ci <- site_div %>%
  dplyr::group_by(Survey.Season, Area) %>%
  dplyr::summarise(
    div_mean = mean(shannon, na.rm = TRUE),
    div_sd   = sd(shannon,   na.rm = TRUE),
    n        = sum(!is.na(shannon)),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(
    div_upper_conf = upper_conf(div_mean, div_sd, n),
    div_lower_conf = pmax(0, lower_conf(div_mean, div_sd, n))
  )

# CIs for Richness per Survey.Season x Area
ri_ci <- site_div %>%
  dplyr::group_by(Survey.Season, Area) %>%
  dplyr::summarise(
    rich_mean = mean(richness, na.rm = TRUE),
    rich_sd   = sd(richness,   na.rm = TRUE),
    n         = sum(!is.na(richness)),
    .groups   = "drop"
  ) %>%
  dplyr::mutate(
    rich_upper_conf = upper_conf(rich_mean, rich_sd, n),
    rich_lower_conf = lower_conf(rich_mean, rich_sd, n)
  )

conf_int_full <- site_div %>%
  left_join(sh_ci, by = c("Survey.Season", "Area")) %>%
  left_join(ri_ci, by = c("Survey.Season", "Area"))

conf_int_full$Area <- factor(conf_int_full$Area, levels = c("Onshelf","Offshelf"))
conf_int_full <- conf_int_full %>% mutate(Survey.Season = factor(Survey.Season, levels = all_years))

div_plot <- ggplot() +
  geom_boxplot(subset(site_div, Area == "Offshelf"),
               mapping = aes(x = Survey.Season, y = shannon,
                             color = Area, alpha = 0, group = Survey.Season),
               coef = 0, fatten = NULL, outlier.shape = NA, show.legend = FALSE) +
  geom_boxplot(subset(site_div, Area == "Onshelf"),
               mapping = aes(x = Survey.Season, y = shannon,
                             color = Area, alpha = 0, group = Survey.Season),
               coef = 0, fatten = NULL, outlier.shape = NA, show.legend = FALSE) +
  geom_errorbar(dplyr::distinct(conf_int_full, Survey.Season, Area, div_lower_conf, div_upper_conf),
                mapping = aes(x = Survey.Season, ymin = div_lower_conf, ymax = div_upper_conf, color = Area)) +
  geom_point(dplyr::distinct(conf_int_full, Survey.Season, Area, div_mean),
             mapping = aes(x = Survey.Season, y = div_mean, color = Area)) +
  labs(x = "Year", y = "Diversity") +
  scale_color_npg(breaks = names(area_labels), labels = area_labels) + theme_few() +
  theme(legend.position = "top",
        panel.background = element_rect(fill='transparent'),
        plot.background  = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background     = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'))

rich_plot <- ggplot() +
  geom_boxplot(subset(site_div, Area == "Offshelf"),
               mapping = aes(x = Survey.Season, y = richness,
                             color = Area, alpha = 0, group = Survey.Season),
               coef = 0, fatten = NULL, outlier.shape = NA, show.legend = FALSE) +
  geom_boxplot(subset(site_div, Area == "Onshelf"),
               mapping = aes(x = Survey.Season, y = richness,
                             color = Area, alpha = 0, group = Survey.Season),
               coef = 0, fatten = NULL, outlier.shape = NA, show.legend = FALSE) +
  geom_errorbar(dplyr::distinct(conf_int_full, Survey.Season, Area, rich_lower_conf, rich_upper_conf),
                mapping = aes(x = Survey.Season, ymin = rich_lower_conf, ymax = rich_upper_conf, color = Area)) +
  geom_point(dplyr::distinct(conf_int_full, Survey.Season, Area, rich_mean),
             mapping = aes(x = Survey.Season, y = rich_mean, color = Area)) +
  labs(x = "Year", y = "Richness") +
  scale_color_npg(breaks = names(area_labels), labels = area_labels) + theme_few() +
  theme(legend.position = "top",
        panel.background = element_rect(fill='transparent'),
        plot.background  = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background     = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'))

div_plot + rich_plot

#### Key Species ####
abund_data_myc <- abund_data %>%
  dplyr::mutate(Myctophidae = Elec.a + Gymnos.br + Gymnos.f + Gymnos + Gymnos.bo + Gymnos.n) %>%
  dplyr::relocate(Myctophidae, .before = Mysida) %>%
  dplyr::select(-Elec.a, -Gymnos.br, -Gymnos.f, -Gymnos, -Gymnos.bo, -Gymnos.n)

abund_data_merged <- abund_data_myc %>%
  dplyr::mutate(
    Antarctomysis_spp. = rowSums(dplyr::across(dplyr::any_of(
      c("Antarctomysis_spp.", "Antarctomysis_spp.", "Anta.m")
    )), na.rm = TRUE),
    
    # Fish eggs + larvae -> Fish_larvae
    Fish_larvae = rowSums(dplyr::across(dplyr::any_of(
      c("Fish_larvae", "Fish.l", "Fish.e")
    )), na.rm = TRUE),
    
    # S._thompsoni + Salpa → Salpa_spp.
    Salpa_spp. = rowSums(dplyr::across(dplyr::any_of(
      c("Salpa_spp.", "S._thompsoni", "Salpa")
    )), na.rm = TRUE),
    
    # Noto + Noto.a → Noto_spp.
    Noto_spp. = rowSums(dplyr::across(dplyr::any_of(
      c("Noto_spp.", "Noto", "Noto.a")
    )), na.rm = TRUE),
    
    # Chaetognatha + Sagi + Sagi.m → Chaetognatha
    Chaetognatha = rowSums(dplyr::across(dplyr::any_of(
      c("Chaetognatha", "Sagi", "Sagi.m")
    )), na.rm = TRUE)
  ) %>%
  dplyr::relocate(Salpa_spp., .before = Sibo) %>%
  dplyr::relocate(Noto_spp., .before = Noto.c) %>%
  dplyr::select(-dplyr::any_of(c(
    "Anta.m",
    "S._thompsoni", "Salpa",
    "Noto", "Noto.a",
    "Sagi", "Sagi.m",
    "Fish.l", "Fish.e"
  )))

#### Long format for species columns
abund_longL <- abund_data_merged %>%
  tidyr::pivot_longer(Amph:Vibi.a, names_to = "Species.code", values_to = "Abun") %>%
  dplyr::mutate(
    Species.code = dplyr::recode(Species.code,
                                 "Cnid"   = "Cnidaria_spp.",
                                 "Cnid.m" = "Cnidaria_spp."
    )
  )

#### Year × Area × Species totals with zero-fill
tot_byAYS <- abund_longL %>%
  dplyr::count(Area, Survey.Season, Species.code, wt = Abun, name = "totalAbun") %>%
  tidyr::complete(Area, Survey.Season, Species.code, fill = list(totalAbun = 0))

#### Area-level summaries across years
average <- tot_byAYS %>%
  dplyr::group_by(Area, Species.code) %>%
  dplyr::summarise(
    observations = sum(totalAbun > 0),
    totalAbun    = sum(totalAbun, na.rm = TRUE),
    avgAbun      = mean(totalAbun, na.rm = TRUE),
    .groups      = "drop"
  )

#### Top species choice: top 10 by mean abundance per Area, plus species occurring ≥ 6 years
top10 <- average %>%
  dplyr::arrange(Area, dplyr::desc(avgAbun)) %>%
  dplyr::group_by(Area) %>% dplyr::slice_head(n = 10) %>% dplyr::ungroup()

freq  <- average %>% dplyr::filter(observations >= 8)

key_sp <- union(top10$Species.code, freq$Species.code) %>%
  union("Fish_larvae")

top_full <- tot_byAYS %>%
  dplyr::filter(Species.code %in% key_sp) %>%
  dplyr::mutate(Species.code = forcats::fct_inorder(Species.code))

#### Palette tied to the Species levels
sp_levels <- levels(top_full$Species.code)
base_pal  <- rcartocolor::carto_pal(n = max(7, length(sp_levels)), name = "Safe")
sp_pal    <- setNames(rep(base_pal, length.out = length(sp_levels)), sp_levels)

#### FIGURE 4 NMDS ####

abund_data_merged_nomyc <- abund_data %>%
  dplyr::mutate(
    Antarctomysis_spp. = rowSums(dplyr::across(dplyr::any_of(
      c("Antarctomysis_spp.", "Antarctomysis_spp.", "Anta.m")
    )), na.rm = TRUE),
    
    # Fish eggs + larvae -> Fish_larvae
    Fish_larvae = rowSums(dplyr::across(dplyr::any_of(
      c("Fish_larvae", "Fish.l", "Fish.e")
    )), na.rm = TRUE),
    
    # S._thompsoni + Salpa -> Salpa_spp.
    Salpa_spp. = rowSums(dplyr::across(dplyr::any_of(
      c("Salpa_spp.", "S._thompsoni", "Salpa")
    )), na.rm = TRUE),
    
    # Noto + Noto.a -> Noto_spp.
    Noto_spp. = rowSums(dplyr::across(dplyr::any_of(
      c("Noto_spp.", "Noto", "Noto.a")
    )), na.rm = TRUE),
    
    # Chaetognatha + Sagi + Sagi.m -> Chaetognatha
    Chaetognatha = rowSums(dplyr::across(dplyr::any_of(
      c("Chaetognatha", "Sagi", "Sagi.m")
    )), na.rm = TRUE)
  ) %>%
  dplyr::relocate(Salpa_spp., .before = Sibo) %>%
  dplyr::relocate(Noto_spp.,  .before = Noto.c) %>%
  dplyr::select(-dplyr::any_of(c(
    "Anta.m",
    "S._thompsoni", "Salpa",
    "Noto", "Noto.a",
    "Sagi", "Sagi.m",
    "Fish.l", "Fish.e"
  )))

metadata_cols <- c("Cruise.Event.Net","Start.of.Event","Survey.Season","Area","Site","Code","nautical_twilight")

dat <- abund_data_merged_nomyc %>%
  dplyr::mutate(ID = Cruise.Event.Net)

# IMPORTANT: keep ID in envdat so you can join later
envdat <- dat %>%
  dplyr::select(ID, dplyr::all_of(metadata_cols)) %>%
  dplyr::mutate(Area = factor(Area, levels = c("Offshelf","Onshelf")))

# IMPORTANT: set rownames = ID before metaMDS, and fill NAs with 0
ecomat <- dat %>%
  dplyr::select(-dplyr::all_of(metadata_cols)) %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~replace(., is.na(.), 0)))

ecomat_mat <- as.data.frame(ecomat)
rownames(ecomat_mat) <- dat$ID

set.seed(24324)
mds <- vegan::metaMDS(sqrt(ecomat_mat), distance = "bray", k = 2, trymax = 100, maxit = 500)

NMDS <- as.data.frame(mds$points) %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(envdat %>% dplyr::select(ID, Area), by = "ID") %>%
  dplyr::rename(x = MDS1, y = MDS2)

NMDS.scores <- as.data.frame(vegan::scores(mds, display = "species"))
colnames(NMDS.scores)[1:2] <- c("NMDS1","NMDS2")
NMDS.scores$species <- rownames(NMDS.scores)
NMDS.scores <- dplyr::filter(NMDS.scores, is.finite(NMDS1), is.finite(NMDS2))

extra_key <- c("Myctophidae", "Elec.a", "Gymnos.br", "Gymnos.f", "Gymnos", "Gymnos.bo", "Gymnos.n")
key_sp_chr <- union(as.character(top_full$Species.code), extra_key)
key_sp_chr <- intersect(key_sp_chr, NMDS.scores$species)

iskey  <- dplyr::filter(NMDS.scores, species %in% key_sp_chr)
nonkey <- dplyr::filter(NMDS.scores, !(species %in% key_sp_chr))

grp.on  <- dplyr::slice(dplyr::filter(NMDS, Area == "Onshelf"),  chull(x, y))
grp.off <- dplyr::slice(dplyr::filter(NMDS, Area == "Offshelf"), chull(x, y))
hull.data <- dplyr::bind_rows(grp.on, grp.off)

fig4_nmds <- ggplot() +
  geom_polygon(data = hull.data,
               aes(x = x, y = y, fill = Area, colour = Area),
               alpha = 0.30) +
  geom_point(data = NMDS,
             aes(x = x, y = y, colour = Area),
             shape = 2, size = 3) +
  geom_point(data = nonkey,
             aes(x = NMDS1, y = NMDS2),
             alpha = 0.25, colour = "grey30", inherit.aes = FALSE, show.legend = FALSE) +
  ggrepel::geom_text_repel(
    data = nonkey,
    aes(x = NMDS1, y = NMDS2, label = species),
    inherit.aes   = FALSE, size = 3.5, max.overlaps = Inf,
    box.padding   = 0.35, point.padding = 0.1, segment.colour = "grey50",
    seed = 24324
  ) +
  geom_point(data = iskey,
             aes(x = NMDS1, y = NMDS2),
             inherit.aes = FALSE, colour = "black", size = 2.2) +
  ggrepel::geom_label_repel(
    data = iskey,
    aes(x = NMDS1, y = NMDS2, label = species),
    inherit.aes = FALSE,
    fill = scales::alpha("white", 0.8), colour = "black",
    label.size = 0.2, max.overlaps = Inf,
    box.padding = 0.4, point.padding = 0.2, segment.colour = "grey50",
    seed = 24324
  ) +
  scale_fill_discrete(breaks = names(area_labels), labels = area_labels) +
  scale_colour_discrete(breaks = names(area_labels), labels = area_labels) +
  annotate("text",
           x = max(NMDS$x) - 0.3, y = max(NMDS$y) + 0.5,
           label = paste("Stress =", round(mds$stress, 2))) +
  coord_equal() + theme_bw() +
  theme(
    legend.position   = c(0.9, 0.85),
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    axis.title        = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.background  = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    aspect.ratio      = 1
  )

fig4_nmds

# ggsave("Final Figures/Griffiths_et_al_fig4.png", width = 10, height = 12, dpi = 300)

## ADONIS ##
dist <- vegdist(sqrt(ecomat))
dummy.adonis_area  <- adonis2(formula=dist~envdat$Area)
dummy.adonis_area 

bd_area <- betadisper(dist, as.factor(envdat$Area))
plot(bd_area, main="betadispersion")
boxplot(bd_area, main="betadispersion")
permutest(bd_area)

#### FIGURE 5 KEY SPECIES ABUND ####
desired_levels <- c(
  "E._superba","Thysanoessa_spp.","E._triacantha","E._frigida",
  "T._gaudichaudii","Antarctomysis_spp.","Salpa_spp.",
  "Chaetognatha","Tomopteris_spp.","Euphausia_spp.",
  "Myctophidae","Fish_larvae","P._macropa"
)

label_map <- setNames(
  list(
    expression(italic("E. superba")),
    expression(italic("Thysanoessa spp.")),
    expression(italic("E. triacantha")),
    expression(italic("E. frigida")),
    expression(italic("T. gaudichaudii")),
    expression(italic("Antarctomysis spp.")),
    expression(italic("Salpa spp.")),
    "Chaetognatha",
    expression(italic("Tomopteris spp.")),
    expression(italic("Euphausia spp.")),
    "Myctophidae",
    "Fish larvae",
    expression(italic("Primno macropa"))
  ),
  desired_levels
)

color_map <- setNames(
  c('grey', '#CC6677', '#AA4499', '#661100',
    '#117733', '#332288', '#88CCEE', '#DDCC77',
    '#44AA99', '#882255', '#999933', '#6699CC',
    '#56B4E9', '#F0E442', '#E69F00', '#CC79A7'
  )[seq_along(desired_levels)],
  desired_levels
)

top_full <- top_full %>%
  dplyr::filter(!is.na(Area), !is.na(Survey.Season), !is.na(Species.code)) %>%
  dplyr::mutate(
    Species.code = factor(as.character(Species.code),
                          levels = intersect(desired_levels, unique(as.character(Species.code)))),
    Survey.Season = factor(as.character(Survey.Season),
                           levels = sort(unique(as.character(Survey.Season))))
  ) %>% droplevels()

bm_long <- data %>%
  dplyr::mutate(Cruise.Event.Net = stringr::str_c(Cruise, Event.Net, sep = ".")) %>%
  dplyr::left_join(site_keys, by = "Cruise.Event.Net") %>%
  dplyr::filter(Site != 4) %>%
  dplyr::mutate(
    Species.code = dplyr::recode(
      Species.code,
      "Fish.l" = "Fish_larvae",
      "Lima"   = "Limacina_spp.",
      "Tomop"  = "Tomopteris_spp.",
      "Them.g" = "T._gaudichaudii",
      "Primn"  = "P._macropa"
    ),
    Species.code = dplyr::case_when(
      Species.code %in% c("Antarctomysis_spp.", "Anta.m")                ~ "Antarctomysis_spp.",
      Species.code %in% c("S._thompsoni", "Salpa", "Salpa_spp.")         ~ "Salpa_spp.",
      Species.code %in% c("Noto", "Noto.a", "Noto_spp.")                 ~ "Noto_spp.",
      Species.code %in% c("Sagi", "Sagi.m", "Chaetognatha")              ~ "Chaetognatha",
      Species.code %in% c("Elec.a","Gymnos.br","Gymnos.f","Gymnos","Gymnos.bo","Gymnos.n")
      ~ "Myctophidae",
      TRUE ~ Species.code
    )
  )

present_all <- intersect(
  desired_levels,
  unique(c(as.character(abund_longL$Species.code), as.character(bm_long$Species.code)))
)

colScale_global <- scale_fill_manual(
  name   = "Species.code",
  values = color_map[present_all],
  breaks = present_all,
  labels = label_map[present_all],
  drop   = FALSE
)

# aggregate biomass (g m^-3)
bm_totals <- bm_long %>%
  dplyr::group_by(Area, Survey.Season, Species.code) %>%
  dplyr::summarise(totalBiomass = sum(`Biomass (g m-3)`, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(Species.code %in% present_all) %>%
  dplyr::mutate(
    Species.code = factor(Species.code, levels = present_all),
    Survey.Season = factor(as.character(Survey.Season),
                           levels = sort(unique(as.character(Survey.Season))))
  ) %>% droplevels()

bm_year_sums <- bm_totals %>%
  dplyr::filter(Species.code != "E._superba") %>%
  dplyr::group_by(Area, Survey.Season) %>%
  dplyr::summarise(sumBiomass = sum(totalBiomass, na.rm = TRUE), .groups = "drop")

bm_shared_max <- 100 * ceiling(max(bm_year_sums$sumBiomass, na.rm = TRUE) / 100)
bm_break_by   <- max(100, round(bm_shared_max/5, -2))  # nice even spacing

# Absolute abundance
abs_off <- top_full %>%
  dplyr::filter(Area == "Offshelf") %>%
  ggplot(aes(x = as.character(Survey.Season), y = totalAbun, fill = Species.code)) +
  geom_col() + colScale_global +
  labs(x = "Year", y = bquote(" Abundance (ind. m"^-2*")"), title = "A)") +
  ggthemes::theme_few() + theme(axis.title.x = element_blank())

abs_on <- top_full %>%
  dplyr::filter(Area == "Onshelf") %>%
  ggplot(aes(x = as.character(Survey.Season), y = totalAbun, fill = Species.code)) +
  geom_col() + colScale_global +
  ggbreak::scale_y_break(
    c(135, 450), scales = 0.25, ticklabels = seq(500, 950, by = 150)
  ) +
  labs(x = "Year", y = bquote(" Abundance (ind. m"^-2*")"), title = "      B)") +
  ggthemes::theme_few() + theme(axis.title.x = element_blank(), legend.position = "none")

# Biomass (w/ and w/out krill)
bm_off_full <- bm_totals %>%
  dplyr::filter(Area == "Offshelf") %>%
  ggplot(aes(x = as.character(Survey.Season), y = totalBiomass, fill = Species.code)) +
  geom_col() + colScale_global +
  labs(x = "Year", y = bquote(" Biomass (g m"^-3*")"), title = "Off-shelf Biomass w/ krill") +
  ggthemes::theme_few() + theme(axis.title.x = element_blank())

bm_on_full <- bm_totals %>%
  dplyr::filter(Area == "Onshelf") %>%
  ggplot(aes(x = as.character(Survey.Season), y = totalBiomass, fill = Species.code)) +
  geom_col() + colScale_global +
  labs(x = "Year", y = bquote(" Biomass (g m"^-3*")"), title = "On-shelf Biomass w/ krill") +
  ggthemes::theme_few() + theme(axis.title.x = element_blank())

ggarrange(
  bm_on_full,
  bm_off_full,
  common.legend = TRUE, legend = "bottom"
)

bm_off <- bm_totals %>%
  dplyr::filter(Area == "Offshelf", Species.code != "E._superba") %>%
  ggplot(aes(x = as.character(Survey.Season), y = totalBiomass, fill = Species.code)) +
  geom_col() + colScale_global +
  labs(x = "Year", y = bquote(" Biomass (g m"^-3*")"), title = "C)") +
  ggthemes::theme_few() + theme(axis.title.x = element_blank())

bm_on <- bm_totals %>%
  dplyr::filter(Area == "Onshelf", Species.code != "E._superba") %>%
  ggplot(aes(x = as.character(Survey.Season), y = totalBiomass, fill = Species.code)) +
  geom_col() + colScale_global +
  labs(x = "Year", y = bquote(" Biomass (g m"^-3*")"), title = "D)") +
  ggthemes::theme_few() + theme(axis.title.x = element_blank())

# Relative abundance (excludes E. superba)
sub_rel_off <- top_full %>%
  dplyr::filter(Area == "Offshelf", Species.code != "E._superba") %>%
  ggplot(aes(x = as.character(Survey.Season), y = totalAbun, fill = Species.code)) +
  geom_col(position = "fill") + colScale_global +
  labs(x = "Year", y = "Relative Abundance", title = "E)") +
  ggthemes::theme_few() + theme(axis.title.x = element_blank())

sub_rel_on <- top_full %>%
  dplyr::filter(Area == "Onshelf", Species.code != "E._superba") %>%
  ggplot(aes(x = as.character(Survey.Season), y = totalAbun, fill = Species.code)) +
  geom_col(position = "fill") + colScale_global +
  labs(x = "Year", y = "Relative Abundance", title = "F)") +
  ggthemes::theme_few() + theme(axis.title.x = element_blank())

grid <- ggarrange(
  abs_off + theme(legend.position = "none"),
  print(abs_on  + theme(legend.position = "none")),
  bm_off  + theme(legend.position = "none"),
  bm_on   + theme(legend.position = "none"),
  sub_rel_off + theme(legend.position = "none"),
  sub_rel_on  + theme(legend.position = "none"),
  ncol = 2, nrow = 3,
  common.legend = TRUE, legend = "bottom"
)

# 2) Column titles, left-aligned within each column
headers <- plot_grid(
  ggdraw() + draw_label("     Off-shelf", x = 0.05, hjust = 0, size = 19),
  ggdraw() + draw_label("     On-shelf", x = 0.05, hjust = 0, size = 19),
  ncol = 2
)

# 3) Stack headers over the grid
final <- plot_grid(headers, grid, ncol = 1, rel_heights = c(0.04, 1))

final

# ggsave("Final Figures/Griffiths_et_al_fig5.png", width = 10, height = 14, dpi = 300)
#
#### FIG. 6 PAIRWISE CORRELATION MATRIX #####
off_spec <- union(subset(top10, Area == "Offshelf")$Species.code,
                  subset(freq,  Area == "Offshelf")$Species.code)
on_spec  <- union(subset(top10, Area == "Onshelf")$Species.code,
                  subset(freq,  Area == "Onshelf")$Species.code)
on_spec <- union(on_spec, "Fish_larvae")

top_off <- top_full %>%
  dplyr::filter(Area == "Offshelf", Species.code %in% off_spec) %>%
  dplyr::select(Area, Species.code, Survey.Season, totalAbun) %>%
  tidyr::pivot_wider(names_from = "Survey.Season",
                     values_from  = "totalAbun",
                     values_fill  = 0)

cor_off <- t(as.matrix(top_off[, 3:ncol(top_off)]))
colnames(cor_off) <- top_off$Species.code

top_on <- top_full %>%
  dplyr::filter(Area == "Onshelf", Species.code %in% on_spec) %>%
  dplyr::select(Area, Species.code, Survey.Season, totalAbun) %>%
  tidyr::pivot_wider(names_from = "Survey.Season",
                     values_from  = "totalAbun",
                     values_fill  = 0)

cor_on <- t(as.matrix(top_on[, 3:ncol(top_on)]))
colnames(cor_on) <- top_on$Species.code

cors <- function(df) {
  M <- Hmisc::rcorr(as.matrix(df), type = "pearson")
  purrr::map(M, ~data.frame(.x))
}

formatted_cors <- function(df){
  cors(df) %>%
    purrr::map(~tibble::rownames_to_column(.x, var = "measure1")) %>%
    purrr::map(~tidyr::pivot_longer(.x, -measure1, names_to = "measure2")) %>%
    dplyr::bind_rows(.id = "id") %>%
    tidyr::pivot_wider(names_from = id, values_from = value) %>%
    dplyr::mutate(
      sig_p    = P < .05,
      p_if_sig = ifelse(P < .05, "*", ""),
      r_if_sig = ifelse(P < .05, paste("r=", round(r, 2)),
                        ifelse(r > 0.8 | r < -0.8, round(r, 2), ""))
    ) %>%
    dplyr::arrange(measure1, measure2) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pair = paste(sort(c(measure1, measure2)), collapse = ",")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(pair, .keep_all = TRUE) %>%
    dplyr::mutate(
      p_if_sig = ifelse(is.na(p_if_sig), "", p_if_sig),
      r_if_sig = ifelse(is.na(r_if_sig), "", r_if_sig)
    )
}

off_heat <- formatted_cors(cor_off) %>% tidyr::drop_na() %>%
  dplyr::mutate(Species1 = gsub("\\_", " ", measure1),
                Species2 = gsub("\\_", " ", measure2)) %>%
  ggplot(aes(Species1, Species2, fill = r, label = paste(p_if_sig, "\n", r_if_sig))) +
  geom_raster() + geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title = "A) Off-shelf") +
  ylab(expression(" Abundance (ind. m "^-2*")")) +
  scale_fill_gradient2(mid = "#FBFEF9", low = "#0C6291", high = "#A63446", limits = c(-1, 1)) +
  geom_text() + theme_classic() +
  scale_y_discrete(expand=c(0,0), labels = c(expression(italic("E. frigida")),
                                             expression(italic("E. superba")), 
                                             expression(italic("E. triacantha")), 
                                             "Fish larvae", "Myctophidae", 
                                             expression(italic("P. macropa")), 
                                             expression(italic("Salpa spp.")),
                                             expression(italic("T. gaudichaudii")), 
                                             expression(italic("Thysanoessa spp.")),
                                             expression(italic("Tomopteris spp.")))) + 
  scale_x_discrete(expand=c(0,0), labels = c("Chaetognatha", expression(italic("E. frigida")),
                                             expression(italic("E. superba")),
                                             expression(italic("E. triacantha")),
                                             "Fish larvae", "Myctophidae", 
                                             expression(italic("P. macropa")), 
                                             expression(italic("Salpa spp.")),
                                             expression(italic("T. gaudichaudii")), 
                                             expression(italic("Thysanoessa spp.")))) +
  theme(axis.text.x = element_text(angle = 35, hjust = 0.95), 
        axis.text   = element_text(size = 12.5),
        plot.title  = element_text(size = 17.5))

on_heat <- formatted_cors(cor_on) %>% tidyr::drop_na() %>%
  dplyr::mutate(Species1 = gsub("\\_", " ", measure1),
                Species2 = gsub("\\_", " ", measure2)) %>%
  ggplot(aes(Species1, Species2, fill = r, label = paste(p_if_sig, "\n", r_if_sig))) +
  geom_raster() + geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title = "B) On-shelf") +
  ylab(expression(" Abundance (ind. m "^-2*")")) +
  scale_fill_gradient2(mid = "#FBFEF9", low = "#0C6291", high = "#A63446", limits = c(-1, 1)) +
  geom_text() + theme_classic() +
  scale_y_discrete(expand=c(0,0), labels = c("Chaetognatha", expression(italic("E. frigida")),
                                             expression(italic("E. superba")),
                                             expression(italic("E. triacantha")),
                                             expression(italic("Euphausia spp.")),
                                             "Fish larvae", expression(italic("P. macropa")),
                                             expression(italic("Salpa spp.")),
                                             expression(italic("T. gaudichaudii")),
                                             expression(italic("Thysanoessa spp.")), 
                                             expression(italic("Tomopteris spp."))
  )) + 
  scale_x_discrete(expand=c(0,0), labels = c(expression(italic("Antarctomysis spp.")),
                                             "Chaetognatha", expression(italic("E. frigida")),
                                             expression(italic("E. superba")),
                                             expression(italic("E. triacantha")),
                                             expression(italic("Euphausia spp.")),
                                             "Fish larvae",
                                             expression(italic("P. macropa")), 
                                             expression(italic("Salpa spp.")),
                                             expression(italic("T. gaudichaudii")), 
                                             expression(italic("Thysanoessa spp.")))) +
  theme(axis.text.x = element_text(angle = 35, hjust = 0.95),
        axis.text   = element_text(size = 12.5),
        plot.title  = element_text(size = 17.5))

legend_big <- theme(
  legend.title = element_text(size = 12),
  legend.text  = element_text(size = 12),
  legend.key.width  = unit(0.8, "cm"),  # thickness of colourbar
  legend.key.height = unit(5.5, "cm")   # height of colourbar
)

off_heat <- off_heat +
  guides(fill = guide_colourbar(barheight = unit(5.5, "cm"),
                                barwidth  = unit(0.8, "cm"))) +
  legend_big

on_heat <- on_heat +
  guides(fill = guide_colourbar(barheight = unit(6, "cm"),
                                barwidth  = unit(1, "cm"))) +
  legend_big

ggarrange(off_heat, on_heat, ncol = 1, common.legend = TRUE, legend = "right")

# ggsave("Final Figures/Griffiths_et_al_fig6.png", device = "png", dpi = 300, height = 12, width = 8.5)

# fig. 6 C
top_onoff <- top_full %>%
  dplyr::filter(!Species.code %in% c("Antarctomysis_spp.","Myctophidae","Euphausia_spp.")) %>%
  dplyr::select(Area, Species.code, Survey.Season, totalAbun) %>%
  tidyr::pivot_wider(names_from = "Survey.Season", values_from = "totalAbun", values_fill = 0) %>%
  dplyr::mutate(SppArea = paste0(sub("shelf$", "", Area), "_", Species.code))  # Offshelf→Off_, Onshelf→On_

year_cols   <- grep("^\\d{4}$", names(top_onoff), value = TRUE)
X           <- top_onoff %>% dplyr::select(dplyr::all_of(year_cols)) %>% as.matrix() %>% t()
colnames(X) <- top_onoff$SppArea
M           <- Hmisc::rcorr(X, type = "pearson")

r_long <- M$r %>% as.data.frame() %>%
  tibble::rownames_to_column("measure1") %>%
  tidyr::pivot_longer(-measure1, names_to = "measure2", values_to = "r")
p_long <- M$P %>% as.data.frame() %>%
  tibble::rownames_to_column("measure1") %>%
  tidyr::pivot_longer(-measure1, names_to = "measure2", values_to = "P")

t <- dplyr::left_join(r_long, p_long, by = c("measure1","measure2")) %>%
  dplyr::filter(grepl("^Off_", measure1), grepl("^On_", measure2)) %>%
  dplyr::mutate(measure1 = sub("^Off_", "", measure1),
                measure2 = sub("^On_",  "", measure2)) %>%
  dplyr::filter(measure1 == measure2) %>%
  dplyr::mutate(p_if_sig = ifelse(P < .05, "*", ""),
                r_if_sig = ifelse(P < .05, paste("r=", round(r,2)),
                                  ifelse(abs(r) > 0.8, round(r,2), ""))) %>%
  dplyr::arrange(measure1)

# horizontal Pearson
sp_keys    <- t$measure1
sp_pretty  <- gsub("_", " ", sp_keys)
lab_expr   <- as.expression(sapply(sp_pretty, function(x)
  if (x %in% c("Chaetognatha","Fish larvae")) x else bquote(italic(.(x)))))

off_vsonHOR <- t %>%
  dplyr::transmute(Species = factor(measure1, levels = sp_keys),
                   r, p_if_sig, r_if_sig) %>%
  dplyr::mutate(Species_lbl = gsub("_"," ", as.character(Species))) %>%
  ggplot2::ggplot(ggplot2::aes(x = "", y = Species_lbl, fill = r,
                               label = paste(p_if_sig, "\n", r_if_sig))) +
  ggplot2::geom_tile() +
  ggplot2::labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation",
                title = "C) Off-shelf vs On-shelf") +
  ggplot2::scale_fill_gradient2(mid = "#FBFEF9", low = "#0C6291",
                                high = "#A63446", limits = c(-1, 1)) +
  ggplot2::geom_text() +
  ggplot2::theme_classic() +
  ggplot2::scale_x_discrete(expand = c(0, 0), position = "top") +
  ggplot2::scale_y_discrete(expand = c(0, 0), limits = sp_pretty, labels = lab_expr) +
  ggplot2::coord_flip() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 75, hjust = 0.95),
                 legend.position = "none")

off_vsonHOR

# ggsave("Final Figures/Griffiths_et_al_fig6.png", device = "png", dpi = 300, height = 12, width = 8.5)

#### TABLE 1 – Biomass summary (biomass, g m^-3) ####
BIOMASS_UNIT <- "g m^-3"
clean_taxon  <- function(s) str_replace_all(s, "_", " ")

norm_keys <- function(df) {
  df %>%
    mutate(
      Survey.Season = as.character(Survey.Season),
      Area = recode(as.character(Area), "Offshelf" = "Off-shelf", "Onshelf" = "On-shelf")
    )
}

bm_long <- norm_keys(bm_long)

season_totals <- bm_long %>%
  dplyr::group_by(Survey.Season, Area) %>%
  dplyr::summarise(
    total_abun = sum(`Abund (m-2)`, na.rm = TRUE),
    total_bio  = sum(`Biomass (g m-3)`,  na.rm = TRUE),
    n_tows     = dplyr::n_distinct(Cruise.Event.Net),
    .groups    = "drop"
  )

agg_wt <- bm_long %>%
  dplyr::group_by(Survey.Season, Area, Species.code) %>%
  dplyr::summarise(
    total_abun = sum(`Abund (m-2)`, na.rm = TRUE),
    total_bio  = sum(`Biomass (g m-3)`,  na.rm = TRUE),
    .groups    = "drop"
  )

most_ind <- agg_wt %>%
  slice_max(total_abun, n = 1, with_ties = FALSE, by = c(Survey.Season, Area)) %>%
  transmute(
    Survey.Season, Area,
    `Abundant taxon` = clean_taxon(as.character(Species.code)),
    Abundance        = total_abun
  )

most_bio <- agg_wt %>%
  slice_max(total_bio, n = 1, with_ties = TRUE, by = c(Survey.Season, Area)) %>%
  transmute(
    Survey.Season, Area,
    `Biomass taxon` = clean_taxon(as.character(Species.code)),
    Biomass         = total_bio
  )

by_area <- most_ind %>%
  left_join(most_bio,     by = c("Survey.Season","Area")) %>%
  left_join(season_totals, by = c("Survey.Season","Area")) %>%
  mutate(
    `Abundance share (%)` = 100 * Abundance / pmax(total_abun, 1e-9),
    `Biomass share (%)`   = 100 * Biomass  / pmax(total_bio,  1e-9)
  ) %>%
  select(
    Survey.Season, Area,
    `Abundant taxon`, Abundance, `Abundance share (%)`,
    `Biomass taxon`, Biomass, `Biomass share (%)`,
    n_tows
  )

table1 <- by_area %>%
  mutate(Area = factor(Area, levels = c("Off-shelf","On-shelf"))) %>%
  arrange(Survey.Season, Area) %>%
  pivot_wider(
    id_cols    = Survey.Season,
    names_from = Area,
    values_from = c(`Abundant taxon`, Abundance, `Abundance share (%)`,
                    `Biomass taxon`, Biomass, `Biomass share (%)`, n_tows),
    names_glue = "{.value} - {Area}"
  ) %>%
  rename_with(~ str_replace(.x, "^Abundance - ", "Abundance (ind. m^-2) - ")) %>%
  rename_with(~ str_replace(.x, "^Biomass - ",   paste0("Biomass (", BIOMASS_UNIT, ") - "))) %>%
  mutate(
    across(
      matches("(Abundance \\(ind\\. m\\^-2\\)|Biomass \\(|share %\\))"),
      ~ round(., 1)
    )
  ) %>%
  dplyr::rename(`Survey season` = Survey.Season) %>%
  arrange(`Survey season`)

# write.csv(table1, "Final Figures/Griffiths_et_al_table1.csv")
