
pairwise_relatednes_l$Pop_period_comparisson

pairwise_relatednes_l %>%
  group_by(Pop_period_comparisson) %>% 
  dplyr::summarise(freq = sum(r>=.99),
                   n = n()) %>% group_by(Pop_period_comparisson)%>%
  mutate(prop = binconf(freq,
                        n,
                        alpha = 0.05,
                        method = "exact")[1],
         lower = binconf(freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[2],
         upper = binconf(freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[3])%>%
  filter(Pop_period_comparisson %in% c("Quibdo_First 10 months",
                                       "Quibdo_Last 9 months",
                                       "Buenaventura_First 10 months",
                                       "Buenaventura_Last 9 months",
                                       "Guapi_First 10 months",
                                       "Guapi_Last 9 months")) %>%
  ggplot(aes(x = factor(Pop_period_comparisson, levels = c("Quibdo_First 10 months",
                                                           "Quibdo_Last 9 months",
                                                           "Buenaventura_First 10 months",
                                                           "Buenaventura_Last 9 months",
                                                           "Guapi_First 10 months",
                                                           "Guapi_Last 9 months")),
             y = prop, fill = factor(Pop_period_comparisson, levels = c("Quibdo_First 10 months",
                                                                        "Quibdo_Last 9 months",
                                                                        "Buenaventura_First 10 months",
                                                                        "Buenaventura_Last 9 months",
                                                                        "Guapi_First 10 months",
                                                                        "Guapi_Last 9 months")
                                     ))) + 
  geom_col(alpha = .85)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
  scale_fill_manual(values = c("firebrick4", "firebrick1", "dodgerblue4", "dodgerblue1", "gold3", "gold1"))+
  theme_bw()+
  labs(y = "Fraccion de pares de muestras\naltamente relacionadas (IBD > 0.99)")+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none")
