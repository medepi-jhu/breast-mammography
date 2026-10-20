
################################################################################
#
# Project: MedEpi Mammography 
#
# Purpose: Visualize preliminary findings and assess outliers
#
# Last Update: 07 May 2025
#
################################################################################

library("tidyverse")

thm <- theme_bw() +
  theme(axis.text = element_text(color="black", size=10),
        axis.title = element_text(color="black", size=12),
        strip.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=8),
        legend.position = "bottom")


# Read in results ---------------------------------------------------------

res1 <- read_csv("../results/rates_state_yr.csv")

res2 <- read_csv("../results/rates_hiv_yr.csv") %>%
  mutate(screen_lower = rate_screen - 1.96*se_screen,
         screen_lower = ifelse(screen_lower<0, 0, screen_lower),
         screen_upper = rate_screen + 1.96*se_screen,
         hiv_base = factor(hiv_base, levels=c(99, 1, 0),
                           labels=c("Overall", "HIV", "No HIV")),
         group = factor(group, levels=c("overall", "black", "white",
                                        "Midwest", "Northeast", "South", "West",
                                        "B_MW", "B_NE", "B_S", "B_W",
                                        "W_MW", "W_NE", "W_S", "W_W"),
                        labels=c("Overall", "Black", "White",
                                 "Midwest", "Northeast", "South", "West",
                                 "Black; Midwest", "Black; Northeast", "Black; South", "Black; West",
                                 "White; Midwest", "White; Northeast", "White; South", "White; West"))) %>%
  filter(year<2020)

res3 <- read_csv("../results/surv_hiv_wt_overall_50.csv") %>%
  mutate(group = "Overall")
res4 <- read_csv("../results/surv_hiv_wt_black_50.csv") %>%
  mutate(group = "Black")
res5 <- read_csv("../results/surv_hiv_wt_white_50.csv") %>%
  mutate(group = "White")

res3 <- bind_rows(res3, res4, res5) %>%
  filter(state!="(s0)") %>%
  mutate(hiv_base = factor(hiv_base, 
                           levels=c(1,0),
                           labels=c("HIV", "No HIV")),
         state = factor(state,
                        levels=c(1,2,3),
                        labels=c("Mammography", "Death", "Breast Cancer")),
         group = factor(group, levels=c("Overall", "Black", "White")))


# Find outliers -----------------------------------------------------------

res.summ <- res1 %>% 
  group_by(analysis) %>%
  summarize(mean_total = mean(log(rate_total), na.rm=T),
            sd_total = sd(log(rate_total), na.rm=T),
            mean_screen = mean(log(rate_screen), na.rm=T),
            sd_screen = sd(log(rate_screen), na.rm=T))

res1.2 <- left_join(res1, res.summ, by="analysis") %>%
  mutate(outlier_total = ifelse(log(rate_total)>(mean_total + 3*sd_total) | 
                                      log(rate_total)<(mean_total - 3*sd_total), rate_total, NA),
         outlier_screen = ifelse(log(rate_screen)>(mean_screen + 3*sd_screen) | 
                                      log(rate_screen)<(mean_screen - 3*sd_screen), rate_screen, NA))

total.outliers <- filter(res1.2, !is.na(outlier_total))
screen.outliers <- filter(res1.2, !is.na(outlier_screen))


# Rate of all mammographies -----------------------------------------------

# All ages
jpeg("../figures/total_all.jpeg", height=6, width=8, units="in", res=300)
ggplot(data=filter(res1.2, analysis=="all"), aes(x=year, y=rate_total, color=state)) +
  labs(title="Rate of all mammography (all ages)", 
       x="\nCalendar Year", y="Rate (per person-year)\n", color="") +
  geom_vline(xintercept=2002, color="black", linetype="dashed") +
  geom_vline(xintercept=2009, color="black", linetype="dashed") +
  geom_vline(xintercept=2016, color="black", linetype="dashed") +
  geom_line(linewidth=0.75) +
  geom_point(aes(y=outlier_total), color="red", shape=1, size=4) +
  scale_x_continuous(breaks=seq(2001, 2021, 2)) +
  scale_color_manual(values=c("#E31A1C", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
                              "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
                              "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
                              "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
                              "green1", "yellow4", "yellow3", "darkorange4", "brown",
                              "orchid4", "yellow1")) +
  guides(color = guide_legend(nrow=3, byrow=T)) +
  thm
dev.off()

# Ages 40+
jpeg("../figures/total_40.jpeg", height=6, width=8, units="in", res=300)
ggplot(data=filter(res1.2, analysis==">40"), aes(x=year, y=rate_total, color=state)) +
  labs(title="Rate of all mammography (ages 40+)", 
       x="\nCalendar Year", y="Rate (per person-year)\n", color="") +
  geom_vline(xintercept=2002, color="black", linetype="dashed") +
  geom_vline(xintercept=2009, color="black", linetype="dashed") +
  geom_vline(xintercept=2016, color="black", linetype="dashed") +
  geom_line(linewidth=0.75) +
  geom_point(aes(y=outlier_total), color="red", shape=1, size=4) +
  scale_x_continuous(breaks=seq(2001, 2021, 2)) +
  scale_color_manual(values=c("#E31A1C", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
                              "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
                              "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
                              "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
                              "green1", "yellow4", "yellow3", "darkorange4", "brown",
                              "orchid4", "yellow1")) +
  guides(color = guide_legend(nrow=3, byrow=T)) +
  thm
dev.off()

# Ages 50+
jpeg("../figures/total_50.jpeg", height=6, width=8, units="in", res=300)
ggplot(data=filter(res1.2, analysis==">50"), aes(x=year, y=rate_total, color=state)) +
  labs(title="Rate of all mammography (ages 50+)", 
       x="\nCalendar Year", y="Rate (per person-year)\n", color="") +
  geom_vline(xintercept=2002, color="black", linetype="dashed") +
  geom_vline(xintercept=2009, color="black", linetype="dashed") +
  geom_vline(xintercept=2016, color="black", linetype="dashed") +
  geom_line(linewidth=0.75) +
  geom_point(aes(y=outlier_total), color="red", shape=1, size=4) +
  scale_x_continuous(breaks=seq(2001, 2021, 2)) +
  scale_color_manual(values=c("#E31A1C", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
                              "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
                              "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
                              "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
                              "green1", "yellow4", "yellow3", "darkorange4", "brown",
                              "orchid4", "yellow1")) +
  guides(color = guide_legend(nrow=3, byrow=T)) +
  thm
dev.off()


# Rate of screening mammographies -----------------------------------------

# All ages
jpeg("../figures/screen_all.jpeg", height=6, width=8, units="in", res=300)
ggplot(data=filter(res1.2, analysis=="all"), aes(x=year, y=rate_screen, color=state)) +
  labs(title="Rate of screening mammography (all ages)", 
       x="\nCalendar Year", y="Rate (per person-year)\n", color="") +
  geom_vline(xintercept=2002, color="black", linetype="dashed") +
  geom_vline(xintercept=2009, color="black", linetype="dashed") +
  geom_vline(xintercept=2016, color="black", linetype="dashed") +
  geom_line(linewidth=0.75) +
  geom_point(aes(y=outlier_screen), color="red", shape=1, size=4) +
  scale_x_continuous(breaks=seq(2001, 2021, 2)) +
  scale_color_manual(values=c("#E31A1C", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
                              "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
                              "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
                              "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
                              "green1", "yellow4", "yellow3", "darkorange4", "brown",
                              "orchid4", "yellow1")) +
  guides(color = guide_legend(nrow=3, byrow=T)) +
  thm
dev.off()

# Ages 40+
jpeg("../figures/screen_40.jpeg", height=6, width=8, units="in", res=300)
ggplot(data=filter(res1.2, analysis=="40-50"), aes(x=year, y=rate_screen, color=state)) +
  labs(title="Rate of screening mammography (ages 40+)", 
       x="\nCalendar Year", y="Rate (per person-year)\n", color="") +
  geom_vline(xintercept=2002, color="black", linetype="dashed") +
  geom_vline(xintercept=2009, color="black", linetype="dashed") +
  geom_vline(xintercept=2016, color="black", linetype="dashed") +
  geom_line(linewidth=0.75) +
  geom_point(aes(y=outlier_screen), color="red", shape=1, size=4) +
  scale_x_continuous(breaks=seq(2001, 2021, 2)) +
  scale_color_manual(values=c("#E31A1C", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
                              "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
                              "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
                              "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
                              "green1", "yellow4", "yellow3", "darkorange4", "brown",
                              "orchid4", "yellow1")) +
  guides(color = guide_legend(nrow=3, byrow=T)) +
  thm
dev.off()

# Ages 50+
jpeg("../figures/screen_50.jpeg", height=6, width=8, units="in", res=300)
ggplot(data=filter(res1.2, analysis==">50"), aes(x=year, y=rate_screen*100, color=state)) +
  labs(title="Rate of screening mammography (ages 50+)", 
       x="\nCalendar Year", y="Rate (per 100 person-years)\n", color="") +
  geom_vline(xintercept=2002, color="black", linetype="dashed") +
  geom_vline(xintercept=2009, color="black", linetype="dashed") +
  geom_vline(xintercept=2016, color="black", linetype="dashed") +
  geom_line(linewidth=0.75) +
  geom_point(aes(y=outlier_screen*100), color="red", shape=1, size=4) +
  scale_x_continuous(breaks=seq(2001, 2021, 2)) +
  scale_color_manual(values=c("#E31A1C", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
                              "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
                              "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
                              "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
                              "green1", "yellow4", "yellow3", "darkorange4", "brown",
                              "orchid4", "yellow1")) +
  guides(color = guide_legend(nrow=3, byrow=T)) +
  thm
dev.off()

