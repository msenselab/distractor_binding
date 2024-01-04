library(tidyverse)
library(cowplot)
source("functions.R")

data_files <- c('1_Partial repeat_2022_May_05_0904.csv', '2_Partial repeat_2022_May_05_0909.csv',
                '3_Partial repeat_2022_May_05_1022.csv', '4_Partial repeat_2022_May_05_1024.csv',
                '5_Partial repeat_2022_May_05_1137.csv', '6_Partial repeat_2022_May_05_1301.csv',
                '7_Partial repeat_2022_May_05_1535.csv', '8_Partial repeat_2022_May_06_0914.csv',
                '9_Partial repeat_2022_May_06_1027.csv', '10_Partial repeat_2022_May_06_1028.csv',
                '11_Partial repeat_2022_May_06_1033.csv', '12_Partial repeat_2022_May_06_1145.csv',
                '13_Partial repeat_2022_May_06_1258.csv', '14_Partial repeat_2022_May_06_1258.csv',
                '15_Partial repeat_2022_May_09_1413.csv', '16_Partial repeat_2022_May_09_1415.csv',
                '17_Partial repeat_2022_May_09_1643.csv', '18_Partial repeat_2022_May_10_1139.csv',
                '19_Partial repeat_2022_May_10_1155.csv', '20_Partial repeat_2022_May_10_1411.csv',
                '21_Partial repeat_2022_May_11_1023.csv', '22_Partial repeat_2022_May_11_1033.csv',
                '23_Partial repeat_2022_May_11_1501.csv', '24_Partial repeat_2022_May_11_1526.csv')

data <- reduce(lapply(lapply(data_files, readData),data.frame), rbind)

N <- length(unique(data$participant))

dists <- c(0, 1:4, 3:1)

data$col_rep <- factor(data$col_rep, levels=c(0,1), labels=c("Change", "Repeat"))
data$shape_rep <- factor(data$shape_rep, levels=c(0,1), labels=c("Change", "Repeat"))
data$resp_rep <- factor(data$resp_rep, levels=c(0,1), labels=c("Change", "Repeat"))
data$pos_fixed <- factor(data$pos_fixed, levels=c(0,1), 
                         labels=c("Random", "Fixed"))

data <- data %>% mutate(prime_correct = 
                          ifelse(PrimeProbe.thisN == 0, key_resp.corr, 
                                 lag(key_resp.corr)), 
                        probe_correct = 
                          ifelse(PrimeProbe.thisN == 1, key_resp.corr, 
                                 lead(key_resp.corr)), 
                        dist_rep=ifelse(col_rep=="Repeat" & shape_rep == "Repeat",
                                        "Full repeat", 
                                        ifelse(col_rep=="Change" & shape_rep=="Change",
                                               "Full change", "Partial repeat")))

rt_dat_probe <- data %>% group_by(shape_rep, col_rep, resp_rep, pos_fixed, 
                                  participant) %>%
  filter(PrimeProbe.thisN==1, prime_correct==1, probe_correct==1, key_resp.rt < 3, key_resp.rt > 0.2) %>%
  summarize(rt=mean(key_resp.rt)) 

srt_dat_probe <- rt_dat_probe %>% summarize(mRT=mean(rt)*1000, 
                                            se_rt=sd(rt)/sqrt(N-1)*1000)

priming_rt_dat <- rt_dat_probe %>% spread(resp_rep, rt) %>% mutate(rre=1000*(Change - Repeat)) %>%
  select(-c(Change, Repeat))

spriming_rt_dat <- priming_rt_dat %>% summarize(mrre=mean(rre), 
                                                srre=sd(rre)/sqrt(N-1))

pj = position_dodge(width = 0.3)
rt_dat_fig <- srt_dat_probe %>% ggplot(aes(x=resp_rep, y=mRT, color=col_rep, shape=shape_rep,
                             group=interaction(col_rep, shape_rep))) +
  geom_point(position=pj, size = 3) + geom_line(position=pj) + theme_classic() +
  facet_wrap(~pos_fixed) + theme(legend.position = "None") +
  geom_errorbar(aes(ymin=mRT-se_rt, ymax=mRT+se_rt),width=0.3, position=pj) +
  labs(x="Response", y="Average RT (ms)", color="Color", shape="Shape")

pj = position_dodge(width = 0.9)
rt_prime_fig <- spriming_rt_dat %>% ggplot(aes(x=pos_fixed, y=mrre, 
                               fill=interaction(col_rep, shape_rep))) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Position", y="Response priming effect (ms)", fill="Distractor features") + 
  scale_fill_discrete( labels=c("Full change", "Shape change", "Color change", "Full repeat"),
                       type = c("#2A9D8F", "#E9C46A", '#F4A261','#E76F51')) + 
  theme(legend.position = "None")
  

err_dat_probe <- data %>% group_by(shape_rep, col_rep, resp_rep, pos_fixed, 
                                   participant) %>%
  filter(PrimeProbe.thisN==1, prime_correct==1) %>%
  summarize(err=1-mean(key_resp.corr)) 

priming_err_dat <- err_dat_probe %>% spread(resp_rep, err) %>% mutate(rre=(Change - Repeat)) %>%
  select(-c(Change, Repeat))

serr_dat_probe <- err_dat_probe %>% summarize(merr=mean(err), se_err=sd(err)/sqrt(n()-1))

pj = position_dodge(width = 0.3)
err_dat_fig <- serr_dat_probe %>% ggplot(aes(x=resp_rep, y=merr, color=col_rep, shape=shape_rep,
                              group=interaction(col_rep, shape_rep))) + 
  geom_point(position=pj, size=3) + geom_line(position=pj) + theme_classic() + 
  facet_wrap(~pos_fixed) +   theme(legend.position = "bottom") +
  geom_errorbar(aes(ymin = merr - se_err, ymax = merr + se_err), width=0.2, position=pj) + 
  scale_shape_discrete(labels=c("Change", "Rep.")) + 
  scale_color_discrete(labels=c("Change", "Rep.")) +
  labs(x="Response", y="Error rate", color = "Color", shape = "Shape")

spriming_err_dat <- priming_err_dat %>% summarize(mrre=mean(rre), 
                                                  srre=sd(rre)/sqrt(N-1))

pj = position_dodge(width = 0.9)
err_prime_fig <- spriming_err_dat %>% ggplot(aes(x=pos_fixed, y=mrre, 
                                                 fill=interaction(col_rep, shape_rep))) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Position", y="Response priming (error)", fill="Dist. features") + 
  scale_fill_discrete(labels=c("Full C", "SC CR", "CC SR", "Full R"),
                      type = c("#2A9D8F", "#E9C46A", '#F4A261','#E76F51')) + 
  theme(legend.position = "bottom")


comb_fig <- plot_grid(rt_dat_fig, rt_prime_fig, err_dat_fig, err_prime_fig,
                      nrow=2, ncol=2, rel_heights = c(1,0.8), 
                      labels = c("A", "B", "C", "D"))

ggsave("./figures/RT_ER_E2.png", comb_fig, width = 10, height = 7)

ies_dat_probe <- full_join(rt_dat_probe, err_dat_probe) %>% mutate(ies = rt/(1-err))

priming_ies_dat <- ies_dat_probe %>% select(-c(rt, err)) %>% spread(resp_rep, ies) %>%
  mutate(rre=1000*(Change - Repeat)) %>% select(-c(Change, Repeat))

pj = position_dodge(width = 0.9)
ies_prime_fig <- priming_ies_dat %>% filter(col_rep==shape_rep) %>% mutate(feat_rep=shape_rep) %>% group_by(feat_rep, pos_fixed) %>% summarize(mrre=mean(rre), srre=sd(rre)/sqrt(N-1)) %>% ggplot(aes(x=pos_fixed, y=mrre, fill=feat_rep)) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Position", y="Response priming effect (ms)", fill="Distractor features") + 
  scale_fill_discrete(labels=c("Full change", "Full repeat"), type=c("#2A9D8F", "#E76F51")) + 
  theme(legend.position = "bottom") + scale_y_continuous(limits=c(-25,90))

ggsave('figures/ies_priming_e2.png', ies_prime_fig, width=4, height=3.5)

ies_prime_full_fig <-  priming_ies_dat %>% group_by(col_rep, shape_rep, pos_fixed) %>% summarize(mrre=mean(rre), srre=sd(rre)/sqrt(N-1)) %>% ggplot(aes(x=pos_fixed, y=mrre, fill=interaction(col_rep, shape_rep))) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Position", y="Response priming effect (ms)", fill="Distractor features") + 
  scale_fill_discrete(labels=c("Full change", "Shape change", "Color change", "Full repeat"), type=c("#2A9D8F", "#E9C46A", '#F4A261','#E76F51')) + 
  theme(legend.position = "bottom")

ggsave('figures/ies_priming_full_e2.png', ies_prime_full_fig, width=6.5, height=4.5)

pj = position_dodge(width = 0.9)
spriming_err_dat %>% ggplot(aes(x=pos_fixed, y=mrre, 
                                fill=interaction(col_rep, shape_rep))) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Position", y="Response priming effect (error)", fill="Distractor features") + 
  scale_fill_discrete( labels=c("Full change", "Shape change", "Color change", "Full repeat"))

# ---- Feature priming ----

dist_feat_prime_ies <- ies_dat_probe %>%
  group_by(shape_rep, col_rep, pos_fixed, resp_rep, participant) %>% 
  summarize(ies=mean(ies)) %>% filter(shape_rep==col_rep) %>% ungroup() %>%
  select(-shape_rep) %>% spread(col_rep, ies) %>% mutate(diff = (Change-Repeat)*1000) %>%
  select(-c(Change, Repeat))


dist_feat_prime_resp <- dist_feat_prime_ies %>% group_by(resp_rep, participant) %>% summarize(diff=mean(diff)) 


sdist_feat_prime_ies <- dist_feat_prime_ies %>% group_by(pos_fixed, resp_rep) %>%
  summarize(m_feat_prime=mean(diff), se_prime=sd(diff)/sqrt(n()-1)) 


pj = position_dodge(width = 0.9)
fig_dist_feat_prime <- sdist_feat_prime_ies %>%
  ggplot(aes(x=resp_rep, y=m_feat_prime, fill=pos_fixed)) + 
  geom_bar(stat="identity", color="black", position=pj) + theme_classic() +
  geom_errorbar(aes(ymin=m_feat_prime-se_prime, ymax=m_feat_prime+se_prime),
                width=0.3, position=pj) +
  labs(x="Response", y="Distractor feature priming effect (ms)", fill="Position") + 
  scale_fill_discrete(labels=c("Random","Fixed"), type=c("#2A9D8F", '#E76F51')) + 
  theme(legend.position = "bottom") + scale_y_continuous(limits=c(-25,90))

priming_fig <- plot_grid(ies_prime_fig, fig_dist_feat_prime, labels=c("A","B"))

ggsave("figures/all_priming.png", priming_fig, width=7.5, height=4)

# ---- Distance effects ---- 

dist_dat <- data %>% mutate(tar_dist = 
                              dists[abs(prime_tar_pos - probe_tar_pos)+1]) %>%
  group_by(tar_dist, shape_rep, col_rep, resp_rep, participant) %>% 
  filter(PrimeProbe.thisN==1, prime_correct==1, probe_correct==1, 
         pos_fixed=="Random", key_resp.rt < 3) %>%
  summarize(rt=mean(key_resp.rt))

sdist_dat <- dist_dat %>% summarize(mRT=mean(rt)*1000)

dist_priming_dat <- dist_dat %>% spread(resp_rep, rt) %>% mutate(rre=(Change - Repeat)*1000) %>%
  select(-c(Change, Repeat))

sdist_priming_dat <- dist_priming_dat %>% summarize(mrre=mean(rre), 
                                                    srre=sd(rre)/sqrt(N-1))


dist_dat_err <- data %>% mutate(tar_dist = 
                                  dists[abs(prime_tar_pos - probe_tar_pos)+1]) %>%
  group_by(tar_dist, shape_rep, col_rep, resp_rep, participant) %>% 
  filter(PrimeProbe.thisN==1, prime_correct==1,
         pos_fixed=="Random") %>%
  summarize(err=mean(1-key_resp.corr))

sdist_dat_err <- dist_dat_err %>% summarize(merr=mean(err))

dist_priming_err_dat <- dist_dat_err %>% spread(resp_rep, err) %>% mutate(rre=(Change - Repeat)) %>%
  select(-c(Change, Repeat))

sdist_priming_err_dat<- dist_priming_err_dat %>% summarize(mrre=mean(rre), 
                                                           srre=sd(rre)/sqrt(N-1))


ies_dat_dist <- full_join(dist_dat, dist_dat_err) %>% mutate(ies = rt*1000/(1-err))

# remove participant with missing data 

ies_dat_dist <- ies_dat_dist %>% filter(participant != 3)

ies_priming_dist_dat <- ies_dat_dist %>% select(-c(rt, err)) %>% spread(resp_rep, ies) %>% mutate(rre=(Change - Repeat)) %>%
  select(-c(Change, Repeat))

ies_priming_dist_dat$tar_dist <- factor(ies_priming_dist_dat$tar_dist)

ies_priming_dist_fig <- ies_priming_dist_dat  %>% filter(shape_rep==col_rep) %>% mutate(feat_rep=shape_rep) %>% 
  group_by(feat_rep, tar_dist) %>% summarize(mrre=mean(rre), srre=sd(rre)/sqrt(N-1)) %>%
  ggplot(aes(x=tar_dist, y=mrre, fill=feat_rep)) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Inter-trial position change", y="Response priming effect (ms)", fill="Distractor features") +
  scale_fill_discrete(labels=c("Full change", "Full repeat"), type=c("#2A9D8F", "#E76F51")) + 
  theme(legend.position = "bottom")

ggsave('figures/ies_priming_dist_e2.png', ies_priming_dist_fig, width=4, height=3)


