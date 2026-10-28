
################################################################################
#
# Project: MedEpi Mammography 
#
# Purpose: Visualize results
#
# Last Update: 04 Sep 2025
#
################################################################################


library("tidyverse")
library("patchwork")

thm <- theme_bw() +
  theme(title = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size=10, angle=45, hjust=1),
        axis.title = element_text(color="black", size=12),
        strip.text = element_text(color="black", size=12),
        panel.spacing = unit(2,"lines"),
        legend.text = element_text(color="black", size=12),
        legend.position = "bottom")


# Read in results ---------------------------------------------------------

crude <- read_csv("../results/crude_rates.csv")
wt <- read_csv("../results/wt_rates.csv")
prob <- read_csv("../results/prob2yr.csv")


# By HIV status and race --------------------------------------------------

fig1 <- filter(crude, analysis==">50") %>% 
  filter(hiv_base!=99) %>% 
  mutate(group = factor(group, 
                        levels=c("overall", "black", "white", "Midwest", "Northeast", "South", "West"), 
                        labels=c("Overall", "non-Hispanic Black", "non-Hispanic white", "Midwest", "Northeast", "South", "West")),
         hiv_base = factor(hiv_base, levels=c(1, 0), labels=c("With HIV", "Without HIV")),
         screen_lower = rate_screen - 1.96*se_screen,
         screen_upper = rate_screen + 1.96*se_screen)

jpeg("../figures/crude_rates.jpeg", height=5, width=12, units="in", res=300)
ggplot(data=fig1, aes(x=year, color=hiv_base, fill=hiv_base)) +
  labs(x="\nCalendar year", y="Rate (per 100 person-years)", 
       color="", fill="") +
  #geom_vline(xintercept=2002, color="black", linetype="dashed") +
  #geom_vline(xintercept=2009, color="black", linetype="dashed") +
  #geom_vline(xintercept=2016, color="black", linetype="dashed") +
  geom_ribbon(aes(ymin=screen_lower*100,
                  ymax=screen_upper*100), alpha=0.5) +
  geom_line(aes(y=rate_screen*100)) +
  scale_x_continuous(expand=c(0, 0), breaks=seq(2001, 2019, 2)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 40)) +
  scale_color_manual(values=c("#68ACE5", "#002D72")) + 
  scale_fill_manual(values=c("#68ACE5", "#002D72")) + 
  facet_wrap(vars(group), ncol=3) +
  thm
dev.off()


# Probabilities -----------------------------------------------------------

# Crude
fig3 <- prob %>% 
  mutate(hiv = factor(hiv, levels=c(1, 0), labels=c("With HIV", "Without HIV")),
         group = factor(group, 
                        levels=c("Overall", "Black", "White"), 
                        labels=c("Overall", "non-Hispanic Black", "non-Hispanic white"))) %>% 
  filter(analysis=="crude")

p1 <- ggplot(data=fig3, aes(x=year, color=hiv, fill=hiv)) +
  labs(tag="A)", x="\nCalendar year", y="Probability",
       color="", fill="") +
  geom_ribbon(aes(ymin=p_screened_ll, ymax=p_screened_ul), alpha=0.5) +
  geom_line(aes(y=p_screened), linewidth=0.75) +
  scale_x_continuous(expand=c(0, 0), breaks=seq(2003, 2019, 2)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 1)) +
  scale_color_manual(values=c("#68ACE5", "#002D72")) + 
  scale_fill_manual(values=c("#68ACE5", "#002D72")) + 
  facet_wrap(vars(group), ncol=3) +
  thm

# Weighted
fig4 <- prob %>% 
  mutate(hiv = factor(hiv, levels=c(1, 0), labels=c("With HIV", "Without HIV")),
         group = factor(group, 
                        levels=c("Overall", "Black", "White"), 
                        labels=c("Overall", "non-Hispanic Black", "non-Hispanic white"))) %>% 
  filter(analysis=="ipw1+ipsw")


p2 <- ggplot(data=fig4, aes(x=year, color=hiv, fill=hiv)) +
  labs(tag="B)",
       x="\nCalendar year", y="Probability",
       color="", fill="") +
  geom_ribbon(aes(ymin=p_screened_ll, ymax=p_screened_ul), alpha=0.5) +
  geom_line(aes(y=p_screened), linewidth=0.75) +
  scale_x_continuous(expand=c(0, 0), breaks=seq(2003, 2019, 2)) +
  scale_y_continuous(expand=c(0, 0),limits=c(0, 1)) +
  scale_color_manual(values=c("#68ACE5", "#002D72")) + 
  scale_fill_manual(values=c("#68ACE5", "#002D72")) + 
  facet_wrap(vars(group), ncol=3) +
  thm

jpeg("../figures/prob2yr.jpeg", height=9, width=12, units="in", res=300)
p1 + p2 + plot_layout(ncol=1, nrow=2, guides="collect") & thm
dev.off()


