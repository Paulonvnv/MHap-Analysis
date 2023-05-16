
colombia_metadata = read.csv('../../AMPLseq_Gates/AMPLseq_Colombia/metadata_feb3_23.csv')


colombia_metadata %>% group_by(X2.8..Grupo.étnico) %>% dplyr::summarise(Samples = n())

colombia_metadata %<>% mutate(
  samples = Código,
  Population = case_when(
    grepl('SP0101254',samples) ~ 'Buenaventura',
    grepl('SP0112286',samples) ~ 'Quibdo',
    grepl('SP0126288',samples) ~ 'Guapi',
    grepl('SP0152835',samples) ~ 'Tumaco',
    grepl('SP019400', samples) ~ 'Puerto Inirida'
  ),
  date_of_symp = X2.9..Fecha.inicio.de.síntomas,
  date_of_collection = X2.3..Fecha.de.reclutamiento.,
  Age = X2.6..Edad,
  Ethnicity = case_when(
    grepl('afro', X2.8..Grupo.étnico, ignore.case = T) ~ 'Afro-descendant',
    grepl('ind', X2.8..Grupo.étnico, ignore.case = T) ~ 'Indigenous',
    grepl('mestizo', X2.8..Grupo.étnico, ignore.case = T) ~ 'Mestizo',
    grepl('otro', X2.8..Grupo.étnico, ignore.case = T) ~ 'Other'
  ),
  Feber = case_when(
    grepl('sin',X2.10..Síntomas, ignore.case = T) ~ 'No',
    grepl('^fiebre',X2.10..Síntomas, ignore.case = T) ~ 'Yes'
  ),
  Nationality = case_when(
    grepl('col|vol', ignore.case = T, X2.2..Nacionalidad) ~ 'Colombian',
    grepl('ven|jar', ignore.case = T, X2.2..Nacionalidad) ~ 'Venezuelan'
  ),
  LastTwoWeeksTravelHistory = case_when(
    grepl('Si', X2.13..Viajes.en.los.últimos.días.o.máximo.2.semanas, ignore.case = T) ~ 'Yes',
    grepl('No', X2.13..Viajes.en.los.últimos.días.o.máximo.2.semanas, ignore.case = T) ~ 'No'
  ),
  MicResult = X3.2.4...Especies, 
  ParasiteCount = as.numeric(X3.2.5.Recuento.parasitario),
  Parasitemia = 6000*as.numeric(colombia_metadata$X3.2.5.Recuento.parasitario)/as.numeric(colombia_metadata$X3.2.7..Leucocitos.contados)
)

colombia_metadata[!grepl('(2020)|(2021)|(2022)|(2023)', colombia_metadata$date_of_collection),][['date_of_collection']] =
  colombia_metadata[!grepl('(2020)|(2021)|(2022)|(2023)', colombia_metadata$date_of_collection),][['date_of_symp']]

colombia_metadata %<>% mutate(month_of_collection = substr(date_of_collection, 1, 7))

colombia_metadata[colombia_metadata$month_of_collection == '2023-10',][['month_of_collection']] = '2023-01'

# Imputing missing data with the previous observation
colombia_metadata[which(is.na(colombia_metadata$month_of_collection)),][['month_of_collection']] = colombia_metadata[which(is.na(colombia_metadata$month_of_collection)) - 1,][['month_of_collection']]
colombia_metadata[which(is.na(colombia_metadata$date_of_collection)),][['date_of_collection']] = colombia_metadata[which(is.na(colombia_metadata$date_of_collection)) - 1,][['date_of_collection']]



# Grouping data by quarters---

colombia_metadata$quarter_of_collection = paste(substr(colombia_metadata$date_of_collection, 1,4), quarters(as.Date(colombia_metadata$date_of_collection)), sep = "-")


colombia_metadata = colombia_metadata[,-1:-35]



dates = sort(unique(colombia_metadata$month_of_collection))


colombia_metadata %>%
  filter(!is.na(Population), MicResult %in% c('P. vivax')) %>%
  group_by(MicResult) %>%
  summarise(nsamples= n())

colombia_metadata %>%
  filter(!is.na(Population), MicResult %in% c('P. vivax')) %>%
  group_by(Population, MicResult) %>%
  summarise(nsamples= n())



colombia_metadata %>%
  filter(!is.na(Population), MicResult %in% c('P. falciparum', 'P. vivax')) %>%
  group_by(Population, month_of_collection, MicResult) %>%
  summarise(nsamples= n())%>%
  ggplot(aes(x = month_of_collection, y = nsamples, fill = factor(Population,
                                                                  levels = c("Quibdo",
                                                                             "Buenaventura",
                                                                             "Guapi",
                                                                             'Tumaco',
                                                                             'Puerto Inirida'))))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("firebrick3", "dodgerblue3", "gold3", "darkseagreen3", "lightsalmon2"))+
  facet_wrap(MicResult~factor(Population,
                              levels = c("Quibdo",
                                         "Buenaventura",
                                         "Guapi",
                                         'Tumaco',
                                         'Puerto Inirida')),
             strip.position = "top", ncol = 5)+
  labs(title = 'B) # Samples collected over time',
       y = "Collected samples by PCD",
       x = "Months")+
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 270, vjust = 0),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position =  "none")+
  scale_alpha_manual(values = c(1,.6))+
  scale_x_discrete(limits = dates)


colombia_metadata %>%
  filter(!is.na(Population), MicResult %in% c('P. falciparum', 'P. vivax')) %>%
  group_by(Population, Ethnicity, MicResult) %>%
  summarise(nsamples= n())%>%
  ggplot(aes(x = Ethnicity, y = nsamples, fill = Ethnicity))+
  geom_col()+
  theme_bw()+
  #scale_fill_manual(values = c("firebrick3", "dodgerblue3", "gold3"))+
  facet_wrap(MicResult~factor(Population,
                              levels = c("Quibdo",
                                         "Buenaventura",
                                         "Guapi",
                                         'Tumaco',
                                         'Puerto Inirida')),
             strip.position = "top", ncol = 5)+
  labs(y = "Collected samples by PCD",
       x = "Ethnic groups")+
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 270, vjust = 0),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position =  "none")+
  scale_alpha_manual(values = c(1,.6))



colombia_metadata %>%
  filter(!is.na(Population), MicResult %in% c('P. falciparum', 'P. vivax')) %>%
  group_by(Population, Nationality, MicResult) %>%
  summarise(nsamples= n())%>%
  ggplot(aes(x = Nationality, y = nsamples, fill = Nationality))+
  geom_col()+
  theme_bw()+
  #scale_fill_manual(values = c("firebrick3", "dodgerblue3", "gold3"))+
  facet_wrap(MicResult~factor(Population,
                              levels = c("Quibdo",
                                         "Buenaventura",
                                         "Guapi",
                                         'Tumaco',
                                         'Puerto Inirida')),
             strip.position = "top", ncol = 5)+
  labs(y = "Collected samples by PCD",
       x = "Ethnic groups")+
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 270, vjust = 0),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position =  "none")+
  scale_alpha_manual(values = c(1,.6))
