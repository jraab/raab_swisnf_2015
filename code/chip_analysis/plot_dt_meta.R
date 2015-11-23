library(dplyr)
library(tidyr) 
library(ggplot2) 
library(data.table) 

file <- '../raab_swisnf_2015/output/dt_matfiles/enhancers_arid1a_mat'
arid1a_enh <- read.table(file, sep = '\t', skip = 1, comment.char = '#')  
colnames(arid1a_enh) <- gsub('V', '', colnames(arid1a_enh) ) 
active <- 1:19707
act_1a_enh <- arid1a_enh[active, 7:ncol(arid1a_enh)] 
rep_1a_enh <- arid1a_enh[-active,7:ncol(arid1a_enh)]

act <- gather(act_1a_enh, pos, val) %>% 
  group_by(pos) %>%
  summarise(u = mean(val, na.rm =T) ) %>%
  mutate(type = 'act') 

rep <- gather(rep_1a_enh, pos, val) %>%
  group_by(pos) %>%
  summarise(u = mean(val, na.rm =T))%>%
  mutate(type = 'rep') 

x <- rbind(act, rep) 
x$pos <- as.numeric(as.character(x$pos) ) 
x %>% ggplot(aes(x=pos, y = u, color = type )) + geom_line() 

read_dt <- function(file) { 
 d <- read.table(file, sep = '\t', skip =1, comment.char = '#') 
 d <- d[,7:ncol(d)] 
 colnames(d) <- 1:ncol(d)
 return(d) 
  }

read_enh_data <- function(file){ 
  d <- read_dt(file)
  active <- 1:19707
  d_act <- d[active,]
  d_rep <- d[-active,]
  d_act <- d_act %>% 
    gather(pos, val) %>%
    group_by(pos) %>%
    summarise(u = mean(val, na.rm =T), s = sd(val, na.rm = T), n = n()) %>%
    mutate(type = 'Active', ci_low = u - (qnorm(0.975) * s/sqrt(n)), 
           ci_high = u + (qnorm(0.975) * s/sqrt(n) ), 
           pos = as.numeric(as.character(pos) ) ) 
  d_rep <- d_rep %>% 
    gather(pos, val) %>%
    group_by(pos) %>%
    summarise(u = mean(val, na.rm =T), s = sd(val, na.rm = T), n = n()) %>%
    mutate(type = 'Repressed', ci_low = u - (qnorm(0.975) * s/sqrt(n)), 
           ci_high = u + (qnorm(0.975) * s/sqrt(n) ), 
           pos = as.numeric(as.character(pos) ) ) 
    
  return(rbind(d_act, d_rep) )
}

arid1a <- read_enh_data(file) %>% mutate(arid = 'ARID1A')
arid1b <- read_enh_data('../raab_swisnf_2015/output/dt_matfiles/enhancers_arid1b_mat') %>% mutate(arid = 'ARID1B')
arid2   <- read_enh_data('../raab_swisnf_2015/output/dt_matfiles/enhancers_arid2_mat') %>% mutate(arid = 'ARID2')

all <- rbind(arid1a, arid1b, arid2)

all %>% 
  ggplot(aes(x = pos, y = u, col = type)) + geom_line() + geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = type), alpha = 0.25) + 
  facet_wrap(~arid)

