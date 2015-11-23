#helpers
theme_figure <- function() {  
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16),
        panel.background = element_rect(fill = NA, color = NA), 
        panel.grid.major = element_line(linetype = 'dashed', color = 'grey80'), 
        legend.position = 'top', 
        legend.direction = 'horizontal', 
        legend.title = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 16), 
        legend.key = element_rect(color = NA, size = 5, fill = NA),
        legend.text = element_text(size = 14), 
        axis.line = element_line(color = 'grey30')
  )
}

# because I sometimes wrote theme_paper instead

theme_paper <-theme_figure <- function() {  
        theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16),
        panel.background = element_rect(fill = NA, color = NA), 
        panel.grid.major = element_line(linetype = 'dashed', color = 'grey80'), 
        legend.position = 'top', 
        legend.direction = 'horizontal', 
        legend.title = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 16), 
        legend.key = element_rect(color = NA, size = 5, fill = NA),
        legend.text = element_text(size = 14),
        axis.line = element_line(color = 'grey30')
  )
}