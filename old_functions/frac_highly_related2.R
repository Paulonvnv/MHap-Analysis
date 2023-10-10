
pairwise_relatednes_l$Pop_period_3months_comparisson

write.csv(pairwise_relatednes_l, "pairwise_relatednes_l.csv")
write.csv(pairwise_relatedness, "pairwise_relatednes_dist_matrix.csv")


poplevels = c(
  "Quibdo_Aug - Oct, 2020",
  "Quibdo_Nov, 2020 - Jan, 2021",
  "Quibdo_Feb - Apr, 2021",
  "Quibdo_May - Jul, 2021",
  "Quibdo_Aug - Oct, 2021",
  "Quibdo_Nov, 2021 - Feb, 2022",
  
"Buenaventura_Aug - Oct, 2020",
"Buenaventura_Nov, 2020 - Jan, 2021",
"Buenaventura_Feb - Apr, 2021",
"Buenaventura_May - Jul, 2021",
"Buenaventura_Aug - Oct, 2021",
"Buenaventura_Nov, 2021 - Feb, 2022",

"Guapi_Aug - Oct, 2020",
"Guapi_Nov, 2020 - Jan, 2021",
"Guapi_Feb - Apr, 2021",
"Guapi_May - Jul, 2021",
"Guapi_Aug - Oct, 2021",
"Guapi_Nov, 2021 - Feb, 2022"
)

dateslevels = c("Aug - Oct, 2020",
               "Nov, 2020 - Jan, 2021",
               "Feb - Apr, 2021",
               "May - Jul, 2021",
               "Aug - Oct, 2021",
               "Nov, 2021 - Feb, 2022")


popcolors = c(brewer.pal(6,"Reds"),
              brewer.pal(5,"Blues"),
              brewer.pal(7,"YlGn"))

pairwise_relatednes_l %>%
  group_by(Pop_period_3months_comparisson) %>% 
  dplyr::summarise(freq = sum(r>=.99),
                   n = n()) %>% group_by(Pop_period_3months_comparisson)%>%
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
  filter(Pop_period_3months_comparisson %in% poplevels)%>%
  mutate(Pop_comparisson = gsub("_.+","",Pop_period_3months_comparisson),
         Dates = gsub(".+_","",Pop_period_3months_comparisson))%>%
  ggplot(aes(x = factor(Dates, levels = dateslevels),
             y = prop, fill = factor(Pop_period_3months_comparisson, levels = poplevels
             ))) + 
  geom_col(alpha = .85)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
  scale_fill_manual(values = popcolors)+
  theme_bw()+
  labs(y = "Fraccion de pares de muestras\naltamente relacionadas (IBD > 0.99)")+
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 315, vjust = 0),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none")+
  facet_wrap(~factor(Pop_comparisson, levels = c("Quibdo",
                                                 "Buenaventura",
                                                 "Guapi")), ncol = 1)


fx_plot_network(pairwise_long_matrix = pairwise_relatednes_l,
                variable = "r",
                cols = c("Yi", "Yj"),
                threshold = .99,
                metadata = ampseq$metadata,
                sample_id = "samples",
                group_by = "Pop_Sampling_period_3months",
                levels = poplevels,
                colors = popcolors)


