library(tidyverse)
library(cowplot)
source("functions.R")

# ---- 1. load raw data ----
data_files <- c('data/01_BY_Partial repeat_2021_Oct_27_1751.csv', 
                'data/02_HS_Partial repeat_2021_Nov_02_1208.csv',
                'data/03__HAB_Partial repeat_2021_Oct_28_1353.csv',
                'data/04_Partial repeat_2021_Oct_28_1717.csv',
                'data/05_EA_Partial repeat_2021_Oct_29_1343.csv',
                'data/06_GA_Partial repeat_2021_Nov_11_1413.csv',
                'data/07_NS_Partial repeat_2021_Oct_29_1659.csv',
                'data/08_AB_Partial repeat_2021_Oct_30_1308.csv',
                'data/09_BC_Partial repeat_2021_Oct_30_1521.csv',
                'data/10_MI_Partial repeat_2021_Nov_01_1029.csv',
                'data/11_YA_Partial repeat_2021_Nov_12_1341.csv',
                'data/12_MZ_Partial repeat_2021_Nov_04_0935.csv',
                'data/13_NF_Partial repeat_2021_Nov_02_1555.csv',
                'data/14_AE_Partial repeat_2021_Nov_02_1743.csv',
                'data/15_SS_Partial repeat_2021_Nov_03_1011.csv',
                'data/16_IB_Partial repeat_2021_Nov_03_1359.csv',
                'data/17_KH_Partial repeat_2021_Nov_03_1614.csv',
                'data/18_KR_Partial repeat_2021_Nov_04_1417.csv',
                'data/19_JM_Partial repeat_2021_Nov_04_1653.csv',
                'data/20_YC_Partial repeat_2021_Nov_05_0950.csv',
                'data/21_SU_Partial repeat_2021_Nov_05_1142.csv',
                'data/22_XHM_Partial repeat_2021_Nov_05_1343.csv',
                'data/23_PP_Partial repeat_2021_Nov_05_1528.csv',
                'data/24_UA_Partial repeat_2021_Nov_05_1732.csv')

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
# ---- 2. RT and error rate ----
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

# error rates
err_dat_probe <- data %>% group_by(shape_rep, col_rep, resp_rep, pos_fixed, 
                                   participant) %>%
  filter(PrimeProbe.thisN==1, prime_correct==1) %>%
  summarize(err=1-mean(key_resp.corr)) 

priming_err_dat <- err_dat_probe %>% spread(resp_rep, err) %>% mutate(rre=(Change - Repeat)) %>%
  select(-c(Change, Repeat))

serr_dat_probe <- err_dat_probe %>% summarize(merr=mean(err), se_err=sd(err)/sqrt(n()-1))

# IES scores
ies_dat_probe <- full_join(rt_dat_probe, err_dat_probe) %>% mutate(ies = rt/(1-err))
sies_dat_probe <- ies_dat_probe %>% summarize(mIES=mean(ies)*1000, 
                                            se_ies=sd(ies)/sqrt(N)*1000)

priming_ies_dat <- ies_dat_probe %>% select(-c(rt, err)) %>% 
  spread(resp_rep, ies) %>% mutate(rre=1000*(Change - Repeat)) %>%
  select(-c(Change, Repeat))

spriming_ies_dat <- priming_ies_dat %>% summarize(mrre=mean(rre), 
                                                srre=sd(rre)/sqrt(N-1))

# ---- 3. ANOVA  ----
ies_dat_probe %>% ezANOVA(dv=ies, wid=participant, within=.(shape_rep, col_rep, pos_fixed, resp_rep))

priming_ies_dat %>% ezANOVA(dv=rre, wid=participant, within=.(shape_rep, col_rep, pos_fixed))

# ---- 4. Figures ----

pj3 = position_dodge(width = 0.3)
ies_dat_fig <- sies_dat_probe %>% ggplot(aes(x=resp_rep, y=mIES, color=interaction(col_rep, shape_rep),
                             group=interaction(col_rep, shape_rep))) +
  geom_point(position=pj3, size=3) + geom_line(position=pj3) + theme_classic() + 
  geom_errorbar(aes(ymin=mIES-se_ies, ymax=mIES+se_ies),width=0.3, position=pj3) +
  scale_color_discrete(labels=c("Full change", "Shape change", "Color change", "Full repeat"),
                      type = c("#2A9D8F", "#E9C46A", '#F4A261','#E76F51')) +
  labs(x="Response", y="Mean IES (ms)", color="Dist. features")  +
  theme(legend.position = "bottom", strip.background= element_blank()) + facet_wrap(~pos_fixed)
ies_dat_fig

pj = position_dodge(width = 0.9)
ies_prime_fig <- spriming_ies_dat %>% ggplot(aes(x=pos_fixed, y=mrre, 
                               fill=interaction(col_rep, shape_rep))) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Position", y="Response priming (ms)", fill="Dist. features") + 
  scale_fill_discrete(labels=c("F.C.", "S.C.", "C.C.", "F.R."),
                      type = c("#2A9D8F", "#E9C46A", '#F4A261','#E76F51')) +
  theme(legend.position = "bottom")

ies_prime_fig

comb_fig <- plot_grid(ies_dat_fig, ies_prime_fig,
                      ncol=2, rel_widths = c(3,2), 
                      labels = c("a", "b"))
comb_fig

ggsave("./figures/exp1.png", comb_fig, width = 10, height = 5)




priming_ies_dat <- ies_dat_probe %>% select(-c(rt, err)) %>% spread(resp_rep, ies) %>%
  mutate(rre=1000*(Change - Repeat)) %>% select(-c(Change, Repeat))

priming_ies_dat %>% filter(col_rep==shape_rep) %>% mutate(feat_rep=shape_rep) %>%
  ezANOVA(dv=rre, wid=participant, within=.(feat_rep,  pos_fixed))

pj = position_dodge(width = 0.9)
ies_prime_fig <- priming_ies_dat %>% filter(col_rep==shape_rep) %>% mutate(feat_rep=shape_rep) %>% group_by(feat_rep, pos_fixed) %>% summarize(mrre=mean(rre), srre=sd(rre)/sqrt(N-1)) %>% ggplot(aes(x=pos_fixed, y=mrre, fill=feat_rep)) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Position", y="Response priming effect (ms)", fill="Distractor features") + 
  scale_fill_discrete(labels=c("Full change", "Full repeat"), type=c("#2A9D8F", "#E76F51")) + 
  theme(legend.position = "bottom") + scale_y_continuous(limits=c(-20,90))

ggsave('figures/ies_priming_e1.png', ies_prime_fig, width=4, height=3.5)

ies_prime_full_fig <-  priming_ies_dat %>% group_by(col_rep, shape_rep, pos_fixed) %>% 
  summarize(mrre=mean(rre), srre=sd(rre)/sqrt(N-1)) %>% 
  ggplot(aes(x=pos_fixed, y=mrre, fill=interaction(col_rep, shape_rep))) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Position", y="Response priming effect (ms)", fill="Distractor features") + 
  scale_fill_discrete(labels=c("Full C", "SC CR", "CC SR", "Full R"),
                      type=c("#2A9D8F", "#E9C46A", '#F4A261','#E76F51')) + 
  theme(legend.position = "bottom")

ggsave('figures/ies_priming_full_e1.png', ies_prime_full_fig, width=5.5, height=3.5)


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
  theme(legend.position = "bottom") + scale_y_continuous(limits=c(-20,90))

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

sdist_dat %>% ggplot(aes(x=tar_dist, y = mRT, color=resp_rep)) +
  geom_point() + geom_line() + theme_classic() + 
  facet_wrap(col_rep ~ shape_rep) + 
  labs(x="Target position change",  y="Average RT (ms)", color="Response")

dist_priming_dat <- dist_dat %>% spread(resp_rep, rt) %>% mutate(rre=(Change - Repeat)*1000) %>%
  select(-c(Change, Repeat))

sdist_priming_dat <- dist_priming_dat %>% summarize(mrre=mean(rre), 
                                                  srre=sd(rre)/sqrt(N-1))

pj = position_dodge(width = 0.9)
sdist_priming_dat %>% ggplot(aes(x=tar_dist, y=mrre, 
                                fill=interaction(col_rep, shape_rep))) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Inter-trial position change", y="Response priming effect (rt)", fill="Distractor features") + 
  scale_fill_discrete( labels=c("Full change", "Shape change", "Color change", "Full repeat"))


dist_dat_err <- data %>% mutate(tar_dist = 
                              dists[abs(prime_tar_pos - probe_tar_pos)+1]) %>%
  group_by(tar_dist, shape_rep, col_rep, resp_rep, participant) %>% 
  filter(PrimeProbe.thisN==1, prime_correct==1,
         pos_fixed=="Random") %>%
  summarize(err=mean(1-key_resp.corr))

sdist_dat_err <- dist_dat_err %>% summarize(merr=mean(err))

sdist_dat_err %>% ggplot(aes(x=tar_dist, y = merr, color=resp_rep)) +
  geom_point() + geom_line() + theme_classic() + 
  facet_wrap(col_rep ~ shape_rep) + 
  labs(x="Target position change",  y="Error rate", color="Response")

dist_priming_err_dat <- dist_dat_err %>% spread(resp_rep, err) %>% mutate(rre=(Change - Repeat)) %>%
  select(-c(Change, Repeat))

sdist_priming_err_dat<- dist_priming_err_dat %>% summarize(mrre=mean(rre), 
                                                    srre=sd(rre)/sqrt(N-1))

pj = position_dodge(width = 0.9)
sdist_priming_err_dat %>% ggplot(aes(x=tar_dist, y=mrre, 
                                 fill=interaction(col_rep, shape_rep))) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Inter-trial position change", y="Response priming effect (error)", fill="Distractor features") + 
  scale_fill_discrete( labels=c("Full change", "Shape change", "Color change", "Full repeat"))

ies_dat_dist <- full_join(dist_dat, dist_dat_err) %>% mutate(ies = rt*1000/(1-err))

ies_priming_dist_dat <- ies_dat_dist %>% select(-c(rt, err)) %>% spread(resp_rep, ies) %>% mutate(rre=(Change - Repeat)) %>%
  select(-c(Change, Repeat))

ies_priming_dist_dat$tar_dist <- factor(ies_priming_dist_dat$tar_dist)

ies_priming_dist_fig <- ies_priming_dist_dat  %>% filter(shape_rep==col_rep) %>% mutate(feat_rep=shape_rep) %>% 
  group_by(feat_rep, tar_dist) %>% summarize(mrre=mean(rre), srre=sd(rre)/sqrt(N-1)) %>%
  ggplot(aes(x=tar_dist, y=mrre, fill=feat_rep)) +
  geom_bar(stat="identity", position=pj, color="black", width=0.8) + theme_classic() + 
  geom_errorbar(aes(ymin=mrre-srre, ymax=mrre+srre), width=0.3, position=pj) +
  labs(x="Inter-trial position change", y="Response priming effect (ms)", 
       fill="Distractor features") + 
  scale_fill_discrete(labels=c("Full change", "Full repeat"), type=c("#2A9D8F", "#E76F51")) + 
  theme(legend.position = "bottom")

ggsave('figures/ies_priming_dist_e1.png', ies_priming_dist_fig, width=4, height=3)


