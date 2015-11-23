# Helper Functions to create plots in Figure 10 
# They are pretty long so I moved them to their own script. 
library(rhdf5) 
library(dplyr)

summarise_data_group <- function(df) { 
  df %>% 
    group_by(signame, peakname, type) %>% 
    summarise(avg = mean(enrichment), med = median(enrichment) ) %>%
    ungroup() %>% 
    arrange(dplyr::desc(avg) ) 
}

# Get the signal from the hdf5 files
readh5vals <- function(fn) {
  tmp <- suppressWarnings(h5read(fn, 'coverages', bit64conversion='bit64'))
  vals <- data.frame(t(tmp$block0_values) )
  vals$rowname <- tmp$axis1
  return(vals) 
}


filter_arid_coverage <- function(arid_coverage_file, peak_names) { 
  # will return just the coverage for the correct list of peaks
  vals <- readh5vals(arid_coverage_file) 
  return(vals[vals$rowname %in% peak_names,]) 
}
make_meta <- function(dfr) { 
  # takes the output of the above functions and plots by group (color) and
  # facets by peak type
  plt <- dfr %>%  
    group_by(peak, group,  Pos2) %>%
    summarise(avg  = mean(Value), ci = 1.96*(sd(Value))/sqrt(n() )) %>%
    mutate(ymin = avg - ci, ymax= avg + ci) %>%
    mutate(Position = rep(seq(-2500, 2490, by = 10)) ) %>%
    ggplot(aes(x=Position, y = avg, color = group, fill = group)) + 
    geom_line(size = 1.5) +
    facet_wrap(~peak) + 
    theme_paper() + 
    scale_fill_manual(values = c('grey20', 'red2') )+ 
    scale_color_manual(values  = c('grey20', 'red2'))  + 
    scale_x_continuous(labels = c('-2.5', '0', '2.5'), breaks = c(-2500, 0,  2500) ) + 
    theme(panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = NA, color = 'grey30'), 
          legend.position = 'none') + 
    xlab('Distance from Peak Center (Kb)') + ylab('Relative Enrichment') 
  return(plt)
}

make_signal_table <- function(dfr) { 
  p <- dfr %>% 
    filter(Pos2 > 100, Pos2 < 400) %>% 
    group_by(group, peak ,rowname) %>% 
    summarise(m = mean(Value, na.rm =T)) %>%
    summarise(u = mean(m, na.rm =T), ci = (qnorm(0.975) * sd(m, na.rm =T)) / sqrt(n() ), count = n() ) %>% 
    mutate(ymin = u - ci, 
           ymax = u + ci  ) 
  return(p) 
}


make_bxplot <- function(dfr, retPlot =TRUE, notch=F) { 
  # just like above but output is a boxplot 
  pltdf <- dfr %>%  
    filter(Pos2 > 100, Pos2 < 400 ) %>%
    group_by(group, peak, rowname)  %>%
    summarise(signal = mean(Value) )
  summary_pvals <- pltdf %>% 
    ungroup() %>%
    group_by(peak) %>% 
    do(tidy(wilcox.test(signal ~ group, data = .)) ) %>%
    print()
  if (retPlot == FALSE){ 
    return(summary_pvals)
  }
  else{  
    plt <- pltdf %>% 
      ggplot(aes(x=group, y = signal, color = group) ) + 
      geom_boxplot(notch = notch, width = 0.9) +
      facet_wrap(~peak) + 
      theme_paper() + 
      scale_color_manual(values = c('grey20', 'red2') ) + 
      theme(panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill = NA, color = 'grey30'), 
            legend.position = 'none', 
            axis.text = element_text(color = 'grey30'), 
            axis.ticks = element_blank() ) + 
      xlab('') + ylab('Median Relative Enrichment') 
    
    return(plt)
  }
}

make_lineplot <- function(df) { 
  plt <- df %>%
    ggplot(aes(x = toupper(name), y = u, col = group)) +
    geom_linerange(aes(ymin = ymin, ymax = ymax), position = position_dodge(0.5)) + 
    geom_point(position = position_dodge(0.5)) + 
    facet_wrap(~peak, ncol =1) + theme_paper() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =1, size = 11), 
          axis.ticks = element_blank(), 
          panel.grid.major.y = element_blank(), 
          panel.grid.major.x = element_line(linetype = 'dashed', color = 'grey90'),
          panel.grid.minor.x = element_line(linetype = 'dashed', color = 'grey90'),
          panel.grid.minor.y = element_blank(), 
          panel.background = element_rect(fill = NA, color = 'grey30') ) + 
    scale_color_manual(values = c('grey20', 'red2') )  + 
    ylab('Mean Signal Per Peak')  + 
    xlab('')  
  return(plt)
}

# Function to compare active and repressed arid1b and arid2
df_for_mp_1b_2 <- function(end, pname_list) { 
  # relies on a df of all peak info 
  # pass in the  end of the filenam we want and the lists of peak names to get
  # and a length(4) list of 
  arid1b_act <- pname_list[[1]]
  arid1b_rep <- pname_list[[2]]
  arid2_act  <- pname_list[[3]]
  arid2_rep  <- pname_list[[4]]
  arid1b_fn <- paste('output/encode_coverages/coverages/arid1b_p_', end, sep ='') 
  arid2_fn  <- paste('output/encode_coverages/coverages/arid2_p_', end, sep = '') 
  act_1b_vals <- filter_arid_coverage(arid1b_fn, arid1b_act) %>% 
    mutate(group = 'Active', peak = 'ARID1B Peaks')
  rep_1b_vals <- filter_arid_coverage(arid1b_fn, arid1b_rep) %>%
    mutate(group = 'Repressed', peak = 'ARID1B Peaks') 
  act_2_vals  <- filter_arid_coverage(arid2_fn, arid2_act)  %>%
    mutate(group = 'Active', peak = 'ARID2 Peaks') 
  rep_2_vals <- filter_arid_coverage(arid2_fn, arid2_rep) %>%
    mutate(group = 'Repressed', peak = 'ARID2 Peaks') 
  all <- rbind(act_1b_vals, rep_1b_vals, act_2_vals, rep_2_vals) 
  all <- all %>% left_join(peak_lookup, by = c('rowname' = 'name')  )
  all %<>% gather('Pos', 'Value', -group, -rowname, -peak, -type) %>% tbl_df() %>%
    mutate(Pos2 = as.numeric(gsub(pattern = 'X', replacement = '', x = Pos)  ) )
  return(all)
} 


# Function to compare cooperative and competitive 1a and 2
df_for_mp_1a_2 <- function(end, pname_list) { 
  # relies on a df of all peak info 
  # pass in the  end of the filenam we want and the lists of peak names to get
  # and a length(4) list of 
  arid1a_coop <- pname_list[[1]]
  arid1a_comp <- pname_list[[2]]
  arid2_coop  <- pname_list[[3]]
  arid2_comp  <- pname_list[[4]]
  arid1a_fn <- paste('output/encode_coverages/coverages/arid1a_p_', end, sep ='') 
  arid2_fn  <- paste('output/encode_coverages/coverages/arid2_p_', end, sep = '') 
  coop_1a_vals <- filter_arid_coverage(arid1a_fn, arid1a_coop) %>% 
    mutate(group = 'Cooperative', peak = 'ARID1A Peaks')
  print(nrow(coop_1a_vals) )
  comp_1a_vals <- filter_arid_coverage(arid1a_fn, arid1a_comp) %>%
    mutate(group = 'Competitive', peak = 'ARID1A Peaks') 
  print(nrow(comp_1a_vals) )
  coop_2_vals  <- filter_arid_coverage(arid2_fn, arid2_coop)  %>%
    mutate(group = 'Cooperative', peak = 'ARID2 Peaks') 
  print(nrow(coop_2_vals) )
  comp_2_vals <- filter_arid_coverage(arid2_fn, arid2_comp) %>%
    mutate(group = 'Competitive', peak = 'ARID2 Peaks') 
  print(nrow(comp_2_vals) )
  all <- rbind(coop_1a_vals, comp_1a_vals, coop_2_vals, comp_2_vals)
  print(dim(all) ) 
  all <- all %>% left_join(peak_lookup, by = c('rowname' = 'name') )
  print(dim(all) ) 
  all %<>% gather('Pos', 'Value', -group, -rowname, -peak, -type) %>% tbl_df() %>%
    mutate(Pos2 = as.numeric(gsub(pattern = 'X', replacement = '', x = Pos)  ) )
  print(dim(all) ) 
  return(all)
} 

