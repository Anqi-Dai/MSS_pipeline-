library (tidyverse)
library (ggpubr)
library (RColorBrewer)


do_plot <- function (my_plot, plot_height = 5, plot_width = 5, plot_legend_position = "right", limit_x = TRUE, limit_y = TRUE, save = TRUE, color_offset = 0, plot_palette = conf_palette) {
  
  # Set common aesthetics for all plots generated by this script. This function will also automatically save the plot in the project folder.
  # Arguments:
  # my_plot: a ggplot2 object
  # plot_height and plot_width: desired dimentions for the plot (in inches -- I know. For some reasons ggsave does not support dimentions in pixels...)
  # plot_legend_position: default is right, but may want to change this for some plots or even to hide legend (position = "none")
  
  colors_in_palette <- brewer.pal (brewer.pal.info[plot_palette,1], plot_palette)
  
  color_order <- 1:length(colors_in_palette) + color_offset
  color_order <- color_order + length(colors_in_palette) * (color_order <= 0)
  color_order <- color_order - length(colors_in_palette) * (color_order > length(colors_in_palette))
  
  plot_colors <- colors_in_palette[order (color_order)]
  
  if (!exists ("i")) i <<- 0 else i <<- i + 1
  
  my_plot  <- my_plot +
    theme_bw () +
    theme (plot.title = element_text(size = 16, face = "bold"),
           plot.subtitle = element_text (size = 14, face = "italic"),
           legend.position = plot_legend_position,
           legend.text = element_text (size = 12),
           legend.title = element_text (size = 12, face = "bold"),
           axis.text = element_text (size = 12),
           axis.title = element_text (size = 16)) +
    scale_color_manual (values = plot_colors, aesthetics = c ("colour", "fill"))
  
  if (save) {
    
    ggsave (plot = my_plot, device = png (), height = plot_height, width = plot_width, units = "in",
            filename = paste (conf_plot_prefix, conf_population, " ", fig_n, description1, " ", gsub ("[^[:alnum:] ]", "", paste (description2, Sys.time (), sep = " ")),".png", sep = ""))
    dev.off ()
    
  } else my_plot
  
}


add_atb_info_and_cleanup <- function (samples, patients, antibiotics) {
  
  new_samples <- samples %>%
    filter (sampletype %in% c("stool", "Stool")) %>%
    inner_join (patients, by = "mrn") %>%
    mutate (day = datecollection - hct) %>%
    filter (day %in% conf_day_min:conf_day_max) %>%
    
    left_join (antibiotics %>%
                 filter (institution == "MSKCC") %>%
                 mutate (mrn = as.numeric (patient_id)),
               by = "mrn") %>%
    mutate (s_atb_start = start.y - datecollection,
            s_atb_stop = stop - datecollection,
            p_atb_start = start.y - hct,
            p_atb_stop = stop - hct) %>%
    
    ## filters out samples with potentially inaccurate (too distant) ATB information -- this also drops samples with no ATB info
    filter (p_atb_stop > -150 & p_atb_start < 150) %>%
    mutate (drug_name_clean = if_else (drug_name_clean == "vancomycin",
                                       if_else (route_simple == "oral",
                                                "vancomycin_oral",
                                                "vancomycin_iv"),
                                       if_else (drug_name_clean %in% c("azithromycin",
                                                                       "cefepime",
                                                                       "ciprofloxacin",
                                                                       "imipenem_cilastatin",
                                                                       "ertapenem",
                                                                       "levofloxacin",
                                                                       "linezolid",
                                                                       "meropenem",
                                                                       "metronidazole",
                                                                       "piperacillin_tazobactam",
                                                                       "sulfamethoxazole_trimethoprim"),
                                                drug_name_clean,
                                                "atb_other"))) %>%
    mutate (atb_can_censor = if_else (drug_name_clean %in% conf_atbs_to_censor, TRUE, FALSE),
            
            ## actual censoring by ATB exposure
            censoring_event = if_else (sampleid %in% sampleid[atb_can_censor == TRUE &
                                                                s_atb_start < 0 &
                                                                s_atb_stop > -conf_atb_days_from_stop],
                                       TRUE, FALSE),
            
            # details on antibiotics responsible for censoring sample
            active_atb_cipro      = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                    s_atb_stop > -conf_atb_days_from_stop &
                                                                    drug_name_clean == "ciprofloxacin"],
                                             TRUE, FALSE),
            active_atb_vanco_iv   = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                    s_atb_stop > -conf_atb_days_from_stop &
                                                                    drug_name_clean == "vancomycin_iv"],
                                             TRUE, FALSE),
            active_atb_vanco_po   = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                    s_atb_stop > -conf_atb_days_from_stop &
                                                                    drug_name_clean == "vancomycin_oral"],
                                             TRUE, FALSE),
            active_atb_imipenem   = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                    s_atb_stop > -conf_atb_days_from_stop &
                                                                    drug_name_clean == "imipenem_cilastatin"],
                                             TRUE, FALSE),
            active_atb_meropenem  = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                      s_atb_stop > -conf_atb_days_from_stop &
                                                                      drug_name_clean == "meropenem"],
                                             TRUE, FALSE),
            active_atb_ertapenem  = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                      s_atb_stop > -conf_atb_days_from_stop &
                                                                      drug_name_clean == "ertapenem"],
                                             TRUE, FALSE),
            active_atb_azitro     = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                    s_atb_stop > -conf_atb_days_from_stop &
                                                                    drug_name_clean == "azithromycin"],
                                             TRUE, FALSE),
            active_atb_cefepime   = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                    s_atb_stop > -conf_atb_days_from_stop &
                                                                    drug_name_clean == "cefepime"],
                                             TRUE, FALSE),
            active_atb_levo       = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                    s_atb_stop > -conf_atb_days_from_stop &
                                                                    drug_name_clean == "levofloxacin"],
                                             TRUE, FALSE),
            active_atb_linezolid  = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                    s_atb_stop > -conf_atb_days_from_stop &
                                                                    drug_name_clean == "linezolid"],
                                             TRUE, FALSE),
            active_atb_metro      = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                      s_atb_stop > -conf_atb_days_from_stop &
                                                                      drug_name_clean == "metronidazole"],
                                             TRUE, FALSE),
            active_atb_piptazo    = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                      s_atb_stop > -conf_atb_days_from_stop &
                                                                      drug_name_clean == "piperacillin_tazobactam"],
                                             TRUE, FALSE),
            active_atb_cotrim     = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                    s_atb_stop > -conf_atb_days_from_stop &
                                                                    drug_name_clean == "sulfamethoxazole_trimethoprim"],
                                             TRUE, FALSE),
            active_atb_other      = if_else (sampleid %in% sampleid[s_atb_start < 0 &
                                                                    s_atb_stop > -conf_atb_days_from_stop &
                                                                    drug_name_clean == "atb_other"],
                                             TRUE, FALSE)) %>%
    filter (!duplicated (sampleid)) %>%

    mutate (consistency = if_else (consistency %in% c("semi-formed, jelly", "semi-formed", "Semi-formed"),
                                   "semi-formed",
                                   if_else (consistency %in% c("formed stool", "Formed stool", "stool"),
                                            "formed",
                                            "liquid")) %>% factor (levels = c("formed", "semi-formed", "liquid")),
              intensity = if_else (intensity %in% c("Ablative", "ABLATIVE"),
                                 "MAC",
                                 if_else (intensity %in% c("NONABL", "Nonablative"),
                                          "NMA",
                                          "RIC")) %>% factor (levels = c ("MAC", "RIC", "NMA")),
            source = if_else (source %in% c("Allo BM SBA-E-",
                                            "Allo PBSC CD34+ (CliniMacs/Miltenyi)",
                                            "Allo PBSC CD34+E- Only",
                                            "Allo PBSC CD34+ (CliniMacs/Miltenyi) and BM SBA-E-",
                                            "Allo PBSC alpha/beta T-Cell and PBSC CD34+"),
                              "TCD",
                              if_else (source %in% c("Allo BM RBC Depleted",
                                                     "Allo BM Unmodified"),
                                       "BM",
                                       if_else (source %in% c("Allo PBSC Unmodified Only",
                                                              "ALLO PBSC & BM Unmodified"),
                                                "PBSC",
                                                if_else (source %in% c("Cord Blood",
                                                                       "Allo Double Cord Blood",
                                                                       "Double Cord Blood",
                                                                       "Allo PBSC CD34+/Double Cord Blood",
                                                                       "ALLO PBSC CD34+/Double Cord Blood",
                                                                       "ALLO PBSC CD34+/SINGLE Cord Blood"),
                                                         "Cord",
                                                         source)))) %>% factor (),
            group = if_else (source != "Cord", as.character (intensity), as.character (source)) %>%
              factor (levels = c("MAC", "Cord", "RIC", "NMA"))) %>%
    select (-atb_can_censor) %>%
    arrange (mrn, day)
  
  return (new_samples)
  
}


group_samples <- function (samples, groups, best_days = NULL, paired = TRUE, only_one_per_pt_per_group = TRUE) {
  
  # arguments:
  # samples: data frame with samples info (usually samples_castori_ag) NB: must be cleaned up (specifically need samples$day)
  # groups: list in the form "group name" = vector of days you want to group together; e.g. list ("pre" = -30:-6, "post" = -5:0)
  # best_days: vector of same length as groups, where a best day is specified for each group
  # paired: if TRUE only patients with samples in every group will be kept, and only one sample per patient per group
  # only_one_per_pt_per_group: if TRUE only one sample per patient per group will be kept
  
  # will return samples grouped as desired. Group will be in samples$code
  
  samples$code <- ""
  
  for (i in names(groups)) {
    samples[samples$day %in% groups[[i]],]$code <- names(groups[i])
  }
  
  samples <- samples[samples$code %in% names(groups),]
  samples$code <- factor (samples$code, levels = names(groups))
  
  if (paired == TRUE || only_one_per_pt_per_group == TRUE) {
    
    if (length (groups) != length (best_days)) stop ("best_days must be specified and must have same length as groups")
    
    if (paired == TRUE) for (i in samples$mrn) if (length (unique (samples$code[samples$mrn == i])) != length (groups)) samples <- samples[!samples$mrn == i,]
    
    samples$diff_day <- abs(samples$day - best_days[match (samples$code, names(groups))])
    samples <- samples[order(samples$diff_day),]
    samples <- samples[!duplicated(paste(samples$mrn, samples$code)),]
    
  }
  
  return (samples)
  
}


add_fold_change <- function (samples, group1, group2) {
  
  # arguments:
  # samples: data frame with samples info NB: must be cleaned up and groups must have been defined
  # group1, group2: fold change will be calculated as group2/group1
  
  g1 <- samples[samples$code == group1,]
  g2 <- samples[samples$code == group2,]
  g1$fold_change <- g2$simpson_reciprocal[match (g2$mrn, g1$mrn)]/g1$simpson_reciprocal
  
  return (g1)
  
}