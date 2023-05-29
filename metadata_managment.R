
library(dplyr)
library(XLConnect)
library(magrittr)
library(ggplot2)
library(stringi)
library(RColorBrewer)

# Get Country codes for wbank----

wbank_countries = read.csv('../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world-regions-according-to-the-world-bank.csv')

# Colombian metadata ----

## Upload metadata----
colombia_metadata = read.csv('../../AMPLseq_Gates/AMPLseq_Colombia/metadata_feb3_23.csv')

## Normalizing Metadata----

### Imputing spelling mistakes and missing data, and adding geographic and temporal variables at different scales----

#### Spelling mistakes in the sample code

colombia_metadata$Código = gsub('sp', 'SP', colombia_metadata$Código)

colombia_metadata$Código = gsub('SP01112286', 'SP0112286', colombia_metadata$Código)

colombia_metadata$Código = gsub('SP112286', 'SP0112286', colombia_metadata$Código)

colombia_metadata$Código = gsub('SP01226288', 'SP0126288', colombia_metadata$Código)

## (Geopolitical administrative levels, date of collection, Microscopy results, ethnicity, nationality, travel history,)
colombia_metadata %<>% mutate(
  Sample_id = Código,
  Study = 'Gates_Colombia',
  Subnational_level2 = case_when(
    grepl('SP0101254',Sample_id) ~ 'Buenaventura',
    grepl('SP0112286',Sample_id) ~ 'Quibdó',
    grepl('SP0126288',Sample_id) ~ 'Guapi',
    grepl('SP0152835',Sample_id) ~ 'Tumaco',
    grepl('SP019400', Sample_id) ~ 'Puerto Inírida'
  ),
  Date_of_Symp = X2.9..Fecha.inicio.de.síntomas,
  Date_of_Collection = X2.3..Fecha.de.reclutamiento.,
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
  Parasite = X3.2.4...Especies, 
  ParasiteCount = as.numeric(X3.2.5.Recuento.parasitario),
  Parasitemia = 6000*as.numeric(colombia_metadata$X3.2.5.Recuento.parasitario)/as.numeric(colombia_metadata$X3.2.7..Leucocitos.contados)
)

### Imputing spelling mistakes and missing data, and adding geographic and temporal variables at different scales----

#### Imputing missing collection date with date of symptoms----
colombia_metadata[!grepl('(2020)|(2021)|(2022)|(2023)', colombia_metadata$Date_of_Collection),][['Date_of_Collection']] =
  colombia_metadata[!grepl('(2020)|(2021)|(2022)|(2023)', colombia_metadata$Date_of_Collection),][['Date_of_Symp']]

#### Adding month of collection----
colombia_metadata %<>% mutate(Month_of_Collection = substr(Date_of_Collection, 1, 7))
colombia_metadata[colombia_metadata$Month_of_Collection == '2023-10',][['Month_of_Collection']] = '2023-01'

#### Imputing missing date data with the previous observation----
colombia_metadata[which(is.na(colombia_metadata$Month_of_Collection)),][['Month_of_Collection']] = colombia_metadata[which(is.na(colombia_metadata$Month_of_Collection)) - 1,][['Month_of_Collection']]
colombia_metadata[which(is.na(colombia_metadata$Date_of_Collection)),][['Date_of_Collection']] = colombia_metadata[which(is.na(colombia_metadata$Date_of_Collection)) - 1,][['Date_of_Collection']]

#### Grouping data by quarters----
colombia_metadata$Quarter_of_Collection = paste(substr(colombia_metadata$Date_of_Collection, 1,4), quarters(as.Date(colombia_metadata$Date_of_Collection)), sep = "-")

#### Grouping data by years----
colombia_metadata$Year_of_Collection = substr(colombia_metadata$Date_of_Collection, 1,4)

#### Grouping data by Subnational levels 1----
source('../../AMPLseq_Gates/AMPLseq_Colombia/fx_gadm_sp_loadCountries.R')

col2.spldf <- gadm_sp_loadCountries("COL",level = 2, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds2/")
delimitations_colombia  = data.frame(Subnational_level1 = col2.spldf$spdf$NAME_1, Subnational_level2 = col2.spldf$spdf$NAME_2)

colombia_metadata$Subnational_level1 = NA

for(Subnational_level2 in unique(colombia_metadata$Subnational_level2)){
  
  colombia_metadata[colombia_metadata$Subnational_level2 == Subnational_level2, ][['Subnational_level1']] =
    unique(delimitations_colombia[toupper(stri_trans_general(delimitations_colombia$Subnational_level2, id = "Latin-ASCII")) %in% 
                                    toupper(stri_trans_general(Subnational_level2, id = "Latin-ASCII")),][['Subnational_level1']])
  
}

#### Grouping data by Subnational levels 0 (none official sub regions)----

colombia_metadata %<>% mutate(Subnational_level0 = case_when(
  Subnational_level1 %in% c("Valle del Cauca", "Chocó",
                            "Cauca", "Nariño") ~ 'Pacific Coast',
  Subnational_level1 %in% c("Guainía") ~ 'East'
))

colombia_metadata$Country = 'Colombia'
colombia_metadata$Country_code = 'COL'
colombia_metadata$United_Nations_Region = 'Latin Ameria and Caribbean'
colombia_metadata$United_Nations_Region_code = 'LAC'

## Selecting variables of interest----

colombia_metadata %<>% select(Sample_id, Study,
                              Date_of_Collection, Date_of_Symp, Month_of_Collection, Quarter_of_Collection, Year_of_Collection,
                              Country, Country_code, United_Nations_Region, United_Nations_Region_code, Subnational_level0, Subnational_level1, Subnational_level2,
                              Age, Ethnicity, Nationality, Feber, LastTwoWeeksTravelHistory, 
                              Parasite, ParasiteCount, Parasitemia
                              )

dates_col = sort(unique(colombia_metadata$Month_of_Collection))
pops0_col = unique(colombia_metadata$Subnational_level0)
pops1_col = unique(colombia_metadata$Subnational_level1)[c(2,1,3,4,5)]
pops2_col = unique(colombia_metadata$Subnational_level2)[c(2,1,3,4,5)]

colombia_metadata %>%
  filter(!is.na(Subnational_level2), Parasite %in% c('P. vivax')) %>%
  group_by(Parasite) %>%
  summarise(nsamples= n())

colombia_metadata %>%
  filter(!is.na(Subnational_level0), Parasite %in% c('P. vivax')) %>%
  group_by(Subnational_level0, Year_of_Collection) %>%
  summarise(nsamples= n())



colombia_metadata %>%
  filter(!is.na(Subnational_level0), Parasite %in% c('P. falciparum', 'P. vivax')) %>%
  group_by(Subnational_level0, Month_of_Collection, Parasite) %>%
  summarise(nsamples= n())%>%
  ggplot(aes(x = Month_of_Collection, y = nsamples, fill = factor(Subnational_level0,
                                                                  levels = pops0_col)))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("firebrick3", "dodgerblue3"#, "gold3", "darkseagreen3", "lightsalmon2"
                               ))+
  facet_wrap(Parasite~factor(Subnational_level0,
                              levels = pops0_col),
             strip.position = "top", ncol = 2)+
  labs(title = 'B) # Samples collected over time',
       y = "Collected samples by PCD",
       x = "Months")+
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position =  "none")+
  scale_alpha_manual(values = c(1,.6))+
  scale_x_discrete(limits = dates_col)


colombia_metadata %>%
  filter(!is.na(Subnational_level0), !is.na(Ethnicity), Parasite %in% c('P. falciparum', 'P. vivax')) %>%
  group_by(Subnational_level0, Ethnicity, Parasite) %>%
  summarise(nsamples= n())%>%
  ggplot(aes(x = Ethnicity, y = nsamples, fill = Ethnicity))+
  geom_col()+
  theme_bw()+
  #scale_fill_manual(values = c("firebrick3", "dodgerblue3", "gold3"))+
  facet_wrap(Parasite~factor(Subnational_level0,
                              levels = pops0_col),
             strip.position = "top", ncol = 2,
             scales = 'free_y')+
  labs(y = "Malaria cases by PCD",
       x = "Ethnic groups")+
  theme(#strip.text.x = element_text(size = 8),
        strip.text.x = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 25, vjust = 0.6),
        axis.title = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.position =  "none")+
  scale_alpha_manual(values = c(1,.6))



colombia_metadata %>%
  filter(!is.na(Subnational_level0), !is.na(Ethnicity), Parasite %in% c('P. falciparum', 'P. vivax')) %>%
  group_by(Subnational_level0, Ethnicity, Parasite) %>%
  summarise(nsamples= n()) %>% pivot_wider(names_from = 'Parasite',
                                           values_from = 'nsamples')%>%
  mutate(`Pv/Pf` = `P. vivax`/`P. falciparum`
         )%>%
  pivot_longer(cols = c('P. falciparum', 'P. vivax', 'Pv/Pf'),
               names_to = 'Parasite', values_to = 'nsamples')%>%
  ggplot(aes(x = Ethnicity, y = nsamples, fill = Ethnicity))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("firebrick3", "gray30", "gray30", "gray30"))+
  facet_wrap(Parasite~factor(Subnational_level0,
                             levels = pops0_col),
             strip.position = "top", ncol = 2,
             scales = 'free_y')+
  labs(y = "Pv/Pf ratio & Malaria cases by PCD",
       x = "Ethnic groups")+
  theme(#strip.text.x = element_text(size = 8),
    strip.text.x = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 25, vjust = 0.6),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 10),
    legend.position =  "none")+
  scale_alpha_manual(values = c(1,.6))

library(Hmisc)

colombia_metadata %>%
  filter(!is.na(Subnational_level0), !is.na(Ethnicity), Parasite %in% c('P. falciparum', 'P. vivax')) %>%
  group_by(Subnational_level0, Ethnicity, Parasite) %>%
  summarise(nsamples= n()) %>% pivot_wider(names_from = 'Parasite',
                                           values_from = 'nsamples')%>%
  mutate(`Pv/Pf` = binconf(x = `P. vivax`, n = `P. falciparum`, method = 'wilson')[[1]],
         lower = binconf(x = `P. vivax`, n = `P. falciparum`, method = 'wilson')[[2]],
         upper = binconf(x = `P. vivax`, n = `P. falciparum`, method = 'wilson')[[3]])%>%
  filter(Subnational_level0 == 'Pacific Coast')%>%
  ggplot(aes(x = Ethnicity, y = `Pv/Pf`, fill = Ethnicity))+
  geom_col()+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .5)+
  theme_bw()+
  scale_fill_manual(values = c("firebrick3", "gray30", "gray30", "gray30"))+
  # facet_wrap(Parasite~factor(Subnational_level0,
  #                            levels = pops0_col),
  #            strip.position = "top", ncol = 2,
  #            scales = 'free_y')+
  labs(y = "Pv/Pf ratio",
       x = "Ethnic groups")+
  theme(#strip.text.x = element_text(size = 8),
    strip.text.x = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 25, vjust = 0.6),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 10),
    legend.position =  "none")+
  scale_alpha_manual(values = c(1,.6))





colombia_metadata %>%
  filter(!is.na(Subnational_level2), !is.na(Ethnicity), Parasite %in% c('P. falciparum', 'P. vivax')) %>%
  group_by(Subnational_level2, Ethnicity, Parasite) %>%
  summarise(nsamples= n()) %>% pivot_wider(names_from = 'Parasite',
                                           values_from = 'nsamples')%>%
  mutate(`Pv/Pf` = `P. vivax`/`P. falciparum`)%>%
  pivot_longer(cols = c('P. falciparum', 'P. vivax', 'Pv/Pf'),
               names_to = 'Parasite', values_to = 'nsamples')%>%
  ggplot(aes(x = Ethnicity, y = nsamples, fill = Ethnicity))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("firebrick3", "gray30", "gray30", "gray30"))+
  facet_wrap(Parasite~factor(Subnational_level2,
                             levels = pops2_col),
             strip.position = "top", ncol = 5,
             scales = 'free_y')+
  labs(y = "Pv/Pf ratio & Malaria cases by PCD",
       x = "Ethnic groups")+
  theme(#strip.text.x = element_text(size = 8),
    strip.text.x = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 25, vjust = 0.6),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 10),
    legend.position =  "none")+
  scale_alpha_manual(values = c(1,.6))


colombia_metadata %>%
  filter(!is.na(Subnational_level0), !is.na(Nationality), Parasite %in% c('P. falciparum', 'P. vivax')) %>%
  group_by(Subnational_level0, Nationality, Parasite) %>%
  summarise(nsamples= n())%>%
  ggplot(aes(x = Nationality, y = nsamples, fill = Nationality))+
  geom_col()+
  theme_bw()+
  #scale_fill_manual(values = c("firebrick3", "dodgerblue3", "gold3"))+
  facet_wrap(Parasite~factor(Subnational_level0,
                              levels = pops0_col),
             strip.position = "top", ncol = 2, scales = 'free_y')+
  labs(y = "Collected samples by PCD",
       x = "Nationality")+
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 270, vjust = 0),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position =  "none")+
  scale_alpha_manual(values = c(1,.6))


# Pv4 metadata----

Pv4_metadata = read.table('../Pv4_samples.txt', sep = '\t', header = T)

## Impute Missing continental information----
Pv4_metadata[Pv4_metadata$Country == 'Thailand',][['Population']] = 'WSEA'
Pv4_metadata[Pv4_metadata$Country == 'China',][['Population']] = 'EAS'
Pv4_metadata[Pv4_metadata$Country == 'Bangladesh',][['Population']] = 'WSEA'
Pv4_metadata[Pv4_metadata$Country == 'Bhutan',][['Population']] = 'WSEA'
Pv4_metadata[Pv4_metadata$Country == 'Madagascar',][['Population']] = 'AF'
Pv4_metadata[Pv4_metadata$Country == 'Mauritania',][['Population']] = 'AF'
Pv4_metadata[Pv4_metadata$Country == 'Myanmar',][['Population']] = 'WSEA'
Pv4_metadata[Pv4_metadata$Country == 'North Korea',][['Population']] = 'EAS'
Pv4_metadata[Pv4_metadata$Country == 'Sudan',][['Population']] = 'AF'

## Homologazing variables----
Pv4_metadata %<>% mutate(
  Sample_id = Sample,
  Study = Study,
  Date_of_Collection = NA,
  Date_of_Symp = NA,
  Month_of_Collection = NA,
  Quarter_of_Collection = NA,
  Year_of_Collection = Year,
  Country = Country,
  Country_code = NA,
  United_Nations_Region = NA,
  United_Nations_Region_code = Population,
  Subnational_level0 = NA,
  Subnational_level1 = gsub('^.+: ','',First.level.administrative.division),
  Subnational_level2 = Site,
  Age = NA,
  Ethnicity = NA,
  Nationality = NA,
  Feber = NA,
  LastTwoWeeksTravelHistory = Is.returning.traveller, 
  Parasite = 'P. vivax',
  ParasiteCount = NA,
  Parasitemia = NA
)

## Selecting variables of interest----
Pv4_metadata %<>% select(Sample_id, Study,
                              Date_of_Collection, Date_of_Symp, Month_of_Collection, Quarter_of_Collection, Year_of_Collection,
                              Country, Country_code, United_Nations_Region, United_Nations_Region_code, Subnational_level0, Subnational_level1, Subnational_level2,
                              Age, Ethnicity, Nationality, Feber, LastTwoWeeksTravelHistory, 
                              Parasite, ParasiteCount, Parasitemia
)

## Getting country codes----

for(country in unique(Pv4_metadata$Country)){
  Pv4_metadata[Pv4_metadata$Country == country,][['Country_code']] = wbank_countries[wbank_countries$Entity == country,][['Code']]
}

## Selecting samples of interest----
Pv4_LAC_metadata = Pv4_metadata %>% filter(United_Nations_Region_code == 'LAM')

Pv4_LAC_metadata$United_Nations_Region = 'Latin America and Caribbean'
Pv4_LAC_metadata$United_Nations_Region_code = 'LAC'

# Venezuelan metadata----

venezuela_metadata = NULL

probe0 = XLConnect::loadWorkbook('../metadata/Shipment 1_David Forero_BASE DE DATOS HARVARD - 29Sep2022 - HRP 1-394 final.xlsx')
sheet_names = getSheets(probe0)
for(sheet in 1:length(getSheets(probe0))){
  probe1 = readWorksheet(probe0, sheet=sheet)
  if(names(probe1)[1] == 'Col1'){
    names(probe1) = probe1[1,]
    probe1 = probe1[-1,]
  }
  venezuela_metadata[[sheet_names[sheet]]] = probe1
}

rm(list = c('probe0', 'probe1', 'sheet', 'sheet_names'))

ven_merged_metadata = NULL

for(sheet in names(venezuela_metadata)){
  temp_obj = venezuela_metadata[[sheet]]
  
  for(col in colnames(temp_obj)){
    temp_obj[[col]] = as.character(temp_obj[[col]])
  }
  
  if(sum(grepl('^Parasito$', names(temp_obj))) == 1){
    
    if(sum(grepl('Procedencia', names(temp_obj))) == 1){
      temp_obj = temp_obj[,c('Código', 'Parasito', 'Procedencia', 'Dirección.', 'X.Fecha', 'Edad', 'Raza')]
      temp_obj$Country = 'Venezuela'
      names(temp_obj) = c('Sample_id', 'Parasite', 'Site_of_Infection', 'House_address', 'Date_of_Collection', 'Age', 'Ethnicity','Country')
      ven_merged_metadata = rbind(ven_merged_metadata, temp_obj)
    }else if(sum(grepl('Zona.Probable.de.Infección', names(temp_obj))) == 1){
      temp_obj = temp_obj[,c('Código', 'Parasito', 'Zona.Probable.de.Infección', 'Dirección.', 'X.Fecha', 'Edad', 'Raza')]
      temp_obj$Country = 'Venezuela'
      names(temp_obj) = c('Sample_id', 'Parasite', 'Site_of_Infection', 'House_address', 'Date_of_Collection', 'Age', 'Ethnicity', 'Country')
      ven_merged_metadata = rbind(ven_merged_metadata, temp_obj)
    }
    
  }else if(sum(grepl('Parásito.actual', names(temp_obj))) == 1){
    temp_obj = temp_obj[,c('Código', 'Parásito.actual', 'ZPI..estado', 'ZPI..municipio', 'Fecha.de.captación', 'Edad')]
    temp_obj$Ethnicity = NA
    temp_obj$Country = 'Venezuela'
    names(temp_obj) = c('Sample_id', 'Parasite', 'Site_of_Infection', 'House_address', 'Date_of_Collection', 'Age', 'Ethnicity', 'Country')
    ven_merged_metadata = rbind(ven_merged_metadata, temp_obj)
  }
}

ven_merged_metadata %<>% filter(!is.na(Sample_id))

ven0.spldf <- gadm_sp_loadCountries("VEN",level = 0, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds0/")
ven1.spldf <- gadm_sp_loadCountries("VEN",level = 1, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds1/")
ven2.spldf <- gadm_sp_loadCountries("VEN",level = 2, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds2/")
delimitations_venezuela  = data.frame(ven2.spldf$spdf$NAME_1, ven2.spldf$spdf$NAME_2)

## Data Normalization----

### Date of infection----
ven_merged_metadata[grepl('1900', ven_merged_metadata$Date_of_Collection),][['Date_of_Collection']] = '2022-01-19'

ven_merged_metadata$Date_of_Collection = as.Date(ven_merged_metadata$Date_of_Collection)
ven_merged_metadata$Quarter_of_Collection = paste(substr(ven_merged_metadata$Date_of_Collection, 1,4), quarters(ven_merged_metadata$Date_of_Collection), sep = "-")
ven_merged_metadata$Month_of_Collection = substr(ven_merged_metadata$Date_of_Collection, 1,7)
ven_merged_metadata$Year_of_Collection = substr(ven_merged_metadata$Date_of_Collection, 1,4)

### Geographic location ----
# Get states and municipalities from GADM data base
## Translate Latin american characters and convert everything to capital letters
ven_states = toupper(stri_trans_general(unique(delimitations_venezuela$ven2.spldf.spdf.NAME_1), id = "Latin-ASCII"))
ven_municipalities = toupper(stri_trans_general(unique(delimitations_venezuela$ven2.spldf.spdf.NAME_2), id = "Latin-ASCII"))

## Translate Latin american characters and convert everything to capital letters in our venezuelan data base
ven_merged_metadata[['Site_of_Infection']] = toupper(stri_trans_general(ven_merged_metadata[['Site_of_Infection']], id = "Latin-ASCII"))
ven_merged_metadata[['House_address']] = toupper(stri_trans_general(ven_merged_metadata[['House_address']], id = "Latin-ASCII"))

# Get Municipality information from:
# Site_of_Infection contains information of the probable site of infection
## in most of the record the format is LOCALITY, MUNICIPIO [MUNICIPALITY NAME]
## in other cases the name of the sate is in Site_of_Infection and the municipality name is in House_address
## in other cases it only contains the name of the city like CIUDAD BOLIVAR that we know the belongs to the municipality of HERES
## in other cases Site_of_Infection is NA so we have to extract the information from House_address

ven_merged_metadata %<>% mutate(
  Subnational_level2 = case_when(
    grepl('MUNICIPIO', Site_of_Infection) ~ 
      gsub('^.+MUNICIPIO |^MUNICIPIO ', '', Site_of_Infection),
    grepl('CIUDAD BOLIVAR', Site_of_Infection) ~ 'HERES',
    grepl('^BOLIVAR$', Site_of_Infection) ~ House_address,
    is.na(Site_of_Infection) & grepl('CIUDAD BOLIVAR', House_address) ~ 'HERES',
    is.na(Site_of_Infection) & grepl('MUNICIPIO', House_address) ~ 
      gsub('^.+MUNICIPIO |^MUNICIPIO ', '', House_address),
    grepl('ANZOATEGUI', Site_of_Infection) ~ gsub(', ANZOATEGUI$', '', Site_of_Infection),
    grepl('DELTA AMACURO', Site_of_Infection) ~ gsub(', DELTA AMACURO$', '', Site_of_Infection),
    grepl('ZULIA', Site_of_Infection) ~ gsub(', ZULIA$', '', Site_of_Infection)
  )
)

## Correction some spelling mistakes in Subnational_level2 information----
ven_merged_metadata %<>% mutate(Subnational_level2 = case_when(
  Subnational_level2 == "ANGOSTURA DEL ORINOCO" ~ 'ANGOSTURA',
  Subnational_level2 == 'MAPIRE, EDO ANZOATEGUI' ~ 'ANZOATEGUI',
  Subnational_level2 == "SIFONTES " ~ "SIFONTES",
  !(Subnational_level2 %in% c("ANGOSTURA DEL ORINOCO", 'MAPIRE, EDO ANZOATEGUI')) ~ Subnational_level2
))

## Get Subnational_level1 information ----
bolivar_mun = c('BOLIVAR', toupper(stri_trans_general(unique(delimitations_venezuela[delimitations_venezuela$ven2.spldf.spdf.NAME_1 == 'Bolívar',]$ven2.spldf.spdf.NAME_2), id = "Latin-ASCII")))
anzoategui_mun = c('ANZOATEGUI', toupper(stri_trans_general(unique(delimitations_venezuela[delimitations_venezuela$ven2.spldf.spdf.NAME_1 == 'Anzoátegui',]$ven2.spldf.spdf.NAME_2), id = "Latin-ASCII")))
zulia_mun = c('ZULIA', toupper(stri_trans_general(unique(delimitations_venezuela[delimitations_venezuela$ven2.spldf.spdf.NAME_1 == 'Zulia',]$ven2.spldf.spdf.NAME_2), id = "Latin-ASCII")))
delta_mun = c('DELTA AMACURO', toupper(stri_trans_general(unique(delimitations_venezuela[delimitations_venezuela$ven2.spldf.spdf.NAME_1 == 'Delta Amacuro',]$ven2.spldf.spdf.NAME_2), id = "Latin-ASCII")))

ven_merged_metadata %<>% mutate(Subnational_level1 = case_when(
  Subnational_level2 %in% bolivar_mun ~ 'BOLIVAR',
  Subnational_level2 %in% anzoategui_mun ~ 'ANZOATEGUI',
  Subnational_level2 %in% zulia_mun ~ 'ZULIA',
  Subnational_level2 %in% delta_mun ~ 'DELTA AMACURO'
))

## Filtering samples with complete geographic information----
venezuela_metadata = ven_merged_metadata %>%filter(!(is.na(Subnational_level1)), !(is.na(Subnational_level2)))

## Homologazing variables----
venezuela_metadata %<>%mutate(
  Sample_id = Sample_id,
  Study = 'Gates_Venezuela',
  Date_of_Collection = Date_of_Collection,
  Date_of_Symp = NA,
  Month_of_Collection = Month_of_Collection,
  Quarter_of_Collection = Quarter_of_Collection,
  Year_of_Collection = Year_of_Collection,
  Country = Country,
  Country_code = 'VEN',
  United_Nations_Region = 'Latin America and Caribbean',
  United_Nations_Region_code = 'LAC',
  Subnational_level0 = NA,
  Subnational_level1 = Subnational_level1,
  Subnational_level2 = Subnational_level2,
  Age = Age,
  Ethnicity = case_when(
    grepl('mestiza', Ethnicity, ignore.case = T) ~ 'Mestizo',
    grepl('indigena', Ethnicity, ignore.case = T) ~ 'Indigenous',
    grepl('negra', Ethnicity, ignore.case = T) ~ 'Afro-descendant'
  ),
  Nationality = NA,
  Feber = NA,
  LastTwoWeeksTravelHistory = NA, 
  Parasite = Parasite,
  ParasiteCount = NA, 
  Parasitemia = NA
)

## Selecting variables of interest----
venezuela_metadata %<>% select(Sample_id, Study,
                              Date_of_Collection, Date_of_Symp, Month_of_Collection, Quarter_of_Collection, Year_of_Collection,
                              Country, Country_code, United_Nations_Region, United_Nations_Region_code, Subnational_level0, Subnational_level1, Subnational_level2,
                              Age, Ethnicity, Nationality, Feber, LastTwoWeeksTravelHistory, 
                              Parasite, ParasiteCount, Parasitemia
)



col11 = brewer.pal(11, 'Set3')


venezuela_metadata %>%
  filter(Subnational_level1 == 'BOLIVAR'#, municipality %in% c('ANGOSTURA', 'HERES', 'SIFONTES', 'SUCRE')
  ) %>%
  group_by(Subnational_level2, Quarter_of_Collection)%>%
  summarise(sample_size = n())%>%
  ggplot(aes(x = Quarter_of_Collection, y = sample_size,
             fill = factor(Subnational_level2,
                           levels = c(
                             'ANGOSTURA',
                             'SIFONTES',
                             'SUCRE',
                             'HERES',
                             'EL CALLAO',
                             'PIAR',
                             'GRAN SABANA',
                             'CEDENO',
                             'PADRE PEDRO CHIEN',
                             'ROSCIO',
                             'CARONI'
                           )))) + 
  geom_col(position = 'stack') +
  facet_wrap(Subnational_level2~., ncol = 2)+
  scale_fill_manual(values = col11)+
  theme(axis.text.x = element_text(angle = 45, vjust = .5))+
  labs(fill = 'Municipality')

# Honduras metadata----

honduras_wb = XLConnect::loadWorkbook('WGS PLASMODIUM SAMPLES.xlsx')
honduras_metadata = readWorksheet(honduras_wb, sheet=1)

honduras_metadata = honduras_metadata[!grepl('Col',colnames(honduras_metadata))]

colnames(honduras_metadata) = c('Sample_id', 'Date_of_Collection2', 'Gender', 'Age',
                                'Subnational_level2', 'Subnational_level1', 'Malaria',
                                'Parasite')

honduras_metadata %<>% mutate(
  Date_of_Collection = case_when(
    grepl('00:00:00',Date_of_Collection2) ~ gsub(' 00:00:00$', '', Date_of_Collection2),
    grepl('/', Date_of_Collection2) ~ paste(paste0('20',gsub('^.+/','', Date_of_Collection2)),
                                           gsub('(^\\d+/|/\\d+$)','', Date_of_Collection2),
                                           gsub('/.+$','', Date_of_Collection2),
                                           sep = '-'),
    !grepl('(/|-)', Date_of_Collection2) & 
      !is.na(Date_of_Collection2) &
      Date_of_Collection2 != 'N/R' ~ paste0(Date_of_Collection2, '-01-01')
  ))


honduras_metadata %<>% mutate(
  Date_of_Collection = case_when(
    as.numeric(gsub('(^\\d+-|-\\d+$)','', Date_of_Collection)) > 12 ~ 
      paste(gsub('-.+$','', Date_of_Collection),
            gsub('^.+-','', Date_of_Collection),
            gsub('(^\\d+-|-\\d+$)','', Date_of_Collection),
            sep = '-'),
    as.numeric(gsub('(^\\d+-|-\\d+$)','', Date_of_Collection)) <= 12 ~ Date_of_Collection
    
  ))


honduras_metadata %<>% mutate(
  Date_of_Collection = case_when(
    nchar(gsub('(^\\d+-|-\\d+$)','', Date_of_Collection)) == 1 ~ 
      paste(gsub('-.+$','', Date_of_Collection),
            paste0('0',gsub('(^\\d+-|-\\d+$)','', Date_of_Collection)),
            gsub('^.+-','', Date_of_Collection),
            sep = '-'),
    nchar(gsub('(^\\d+-|-\\d+$)','', Date_of_Collection)) == 2 ~ Date_of_Collection
    
  ))


honduras_metadata %<>% mutate(Date_of_Collection = case_when(
  Date_of_Collection == '20211-11-15' ~ '2021-11-15',
  Date_of_Collection2 == '1711/22' ~ '2022-11-17',
  Date_of_Collection != '20211-11-15' & Date_of_Collection2 != '1711/22' ~ Date_of_Collection
  ))


honduras_metadata %<>% mutate(
  Month_of_Collection = substr(Date_of_Collection, 1, 7),
  Year_of_Collection = substr(Date_of_Collection, 1, 4))

honduras_metadata %>% filter(Parasite == 'P. VIVAX') %>% nrow()

honduras_metadata %>% filter(Parasite == 'P. VIVAX') %>%
  group_by(Parasite,Subnational_level1,Subnational_level2,Year_of_Collection) %>%
  summarise(Cases = n())

names(honduras_metadata)

honduras_metadata %<>% mutate(
  Sample_id = Sample_id,
  Study = 'U19_Honduras',
  Date_of_Collection = Date_of_Collection,
  Date_of_Symp = NA,
  Month_of_Collection = Month_of_Collection,
  Quarter_of_Collection = paste(Year_of_Collection,quarters(as.Date(Date_of_Collection)), sep = '-'),
  Year_of_Collection = Year_of_Collection,
  Country = 'Honduras',
  Country_code = 'HND',
  United_Nations_Region = 'Latin America and Caribbean',
  United_Nations_Region_code = 'LAC',
  Subnational_level0 = NA,
  Subnational_level1 = Subnational_level1,
  Subnational_level2 = Subnational_level2,
  Age = Age,
  Ethnicity = NA,
  Nationality = NA, 
  Feber = NA,
  LastTwoWeeksTravelHistory = NA, 
  Parasite = Parasite,
  ParasiteCount = NA,
  Parasitemia = NA
)


honduras_metadata %<>% select(Sample_id, Study,
                               Date_of_Collection, Date_of_Symp, Month_of_Collection, Quarter_of_Collection, Year_of_Collection,
                               Country, Country_code, United_Nations_Region, United_Nations_Region_code, Subnational_level0, Subnational_level1, Subnational_level2,
                               Age, Ethnicity, Nationality, Feber, LastTwoWeeksTravelHistory, 
                               Parasite, ParasiteCount, Parasitemia
)

# Merge metadata----

metadata = rbind(colombia_metadata,
      venezuela_metadata,
      honduras_metadata,
      Pv4_LAC_metadata)

metadata %<>%filter(grepl('vivax',Parasite, ignore.case = T))

summary_metadata = metadata %>%
  group_by(Country, Subnational_level1, Year_of_Collection)%>%
  summarise(Sample_Size = n())

summary_metadata$Year_of_Collection = as.integer(summary_metadata$Year_of_Collection)

summary_metadata = rbind(summary_metadata,
                         data.frame(Country = 'Panama',
                                    Subnational_level1 = 'Panama',
                                    Year_of_Collection = c(2009, 2018, 2022),
                                    Sample_Size = c(35, 24, 100)))

summary_metadata = rbind(summary_metadata,
                         data.frame(Country = 'Peru',
                                    Subnational_level1 = 'Amazonas',
                                    Year_of_Collection = 2022,
                                    Sample_Size = 90))

summary_metadata = rbind(summary_metadata,
                         data.frame(Country = 'Guatemala',
                                    Subnational_level1 = c('Escuintla', 'Alta Verapaz', 'Izabal'),
                                    Year_of_Collection = 2023,
                                    Sample_Size = 40))


summary_metadata %<>% mutate(Period_of_Collection = case_when(
  Year_of_Collection < 2015 ~ '< 2015',
  Year_of_Collection >= 2015 & Year_of_Collection <= 2020 ~ '2015 - 2020',
  Year_of_Collection > 2020 ~ '> 2020'
))

summary_metadata %<>%
  group_by(Country, Subnational_level1, Period_of_Collection)%>%
  summarise(Sample_Size = sum(Sample_Size))

summary_metadata%<>% mutate(Population = case_when(
  Country =='Brazil' ~ 'Brazil',
  Country == 'Peru' & Subnational_level1 != 'Amazonas' ~ 'Loreto',
  Country == 'Peru' & Subnational_level1 == 'Amazonas' ~ 'Amazonas',
  Country == 'Honduras' ~ 'Gracias a Dios',
  Country == 'Panama' ~ 'Panama',
  Country == 'Venezuela' ~ 'Bolivar',
  Country == 'Colombia' & (Subnational_level1 %in% c('Antioquia', 'Bolivar', 'Cordoba')) ~ 'North & Center Northwestern',
  Country == 'Colombia' & (Subnational_level1 %in% c('Choco', 'Chocó', 'Colombia', 'Narino','Nariño',
                                                    'Risaralda', 'Valle del Cauca', 'Cauca')) ~ 'Pacific Coast',
  Country == 'Colombia' & (Subnational_level1 %in% c('Guainía')) ~ 'East',
  Country == 'Mexico' ~ 'Mexico',
  Country == 'Guatemala' ~ Subnational_level1
  ))

summary_metadata %<>%
  group_by(Country, Population, Period_of_Collection)%>%
  summarise(Sample_Size = sum(Sample_Size))

summary_metadata %<>% filter(!is.na(Population), !is.na(Period_of_Collection))

summary_metadata %<>% mutate(
  lat = case_when(
    Population == 'Brazil' ~ -8,
    Population == 'East' ~ 3,
    Population == 'North & Center Northwestern' ~ 8,
    Population == 'Pacific Coast' ~ 4,
    Population == 'Alta Verapaz' ~ 15.5,
    Population == 'Escuintla' ~ 14.25,
    Population == 'Izabal' ~ 15.5,
    Population == 'Gracias a Dios' ~ 15,
    Population == 'Mexico' ~ 17,
    Population == 'Panama' ~ 8.98,
    Population == 'Amazonas' ~ -5,
    Population == 'Loreto' ~ -3,
    Population == 'Bolivar' ~ 5.5
  ),
  lon = case_when(
    Population == 'Brazil' ~ -65,
    Population == 'East' ~ -69,
    Population == 'North & Center Northwestern' ~ -76,
    Population == 'Pacific Coast' ~ -77,
    Population == 'Alta Verapaz' ~ -90,
    Population == 'Escuintla' ~ -91,
    Population == 'Izabal' ~ -88.5,
    Population == 'Gracias a Dios' ~ -84,
    Population == 'Mexico' ~ -93,
    Population == 'Panama' ~ -79.5,
    Population == 'Amazonas' ~ -78.5,
    Population == 'Loreto' ~ -74,
    Population == 'Bolivar' ~ -63.5
  )
)

countries = data.frame(country = c('Guatemala', 'Honduras', 'Panama', 'Colombia', 'Venezuela', 'Peru', 'Guyana'),
                       code = c('GTM', 'HND', 'PAN', 'COL', 'VEN', 'PER', 'GUY'),
                       lon = c(-90.52271,-87.202438,-79.516670, -74.063644, -66.916664, -75, -58.1667),
                       lat = c(14.628434,14.5,8.983333, 4.624335, 7.500000, -10, 6.8000))

# Plot study sites----

lac = wbank_countries %>% filter(World.Region.according.to.the.World.Bank == 'Latin America and Caribbean')

lac %<>% filter((Entity %in% c('Guatemala', 'Belize','Honduras', 'Panama', 'El Salvador', 'Nicaragua',
                               'Costa Rica', 'Ecuador', #'Brazil', 'French Guiana', 'Suriname',
                             'Colombia', 'Venezuela', 'Guyana', 'Peru')))

lac_sites = lac %>% filter((Entity %in% c('Guatemala', 'Honduras', 'Panama',
                                          'Colombia', 'Venezuela', 'Guyana', 'Peru')))


lac0.spldf <- gadm_sp_loadCountries(lac$Code,level = 0, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds0/")
lac_sites1.spldf <- gadm_sp_loadCountries(lac_sites$Code,level = 1, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds1/")
lac_sites2.spldf <- gadm_sp_loadCountries(lac_sites$Code,level = 2, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds2/")

gtm_sites = gadm_sp_loadCountries('GTM',level = 1, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds1/")
hnd_sites = gadm_sp_loadCountries('HND',level = 1, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds1/")
col_sites = gadm_sp_loadCountries('COL',level = 2, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds2/")
per_sites = gadm_sp_loadCountries('PER',level = 1, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds1/")
per_sites2 = gadm_sp_loadCountries('PER',level = 2, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds2/")
ven_sites = gadm_sp_loadCountries('VEN',level = 2, basefile = "../../AMPLseq_Gates/AMPLseq_Colombia/GADMTools/world.rds.files/rds2/")



countries = data.frame(country = c('Guatemala', 'Honduras', 'Panama', 'Colombia', 'Venezuela', 'Peru', 'Guyana'),
                       code = c('GTM', 'HND', 'PAN', 'COL', 'VEN', 'PER', 'GUY'),
                       lon = c(-90.52271,-87.202438,-79.516670, -74.063644, -66.916664, -75, -58.1667),
                       lat = c(14.628434,14.5,8.983333, 4.624335, 7.500000, -10, 6.8000))


library(scatterpie)
library(tidyr)
summary_metadata_long = summary_metadata %>% pivot_wider(names_from = 'Period_of_Collection',
                                                         values_from = 'Sample_Size') %>% mutate(
                                                           Sample_Size = sum(`< 2015`, `2015 - 2020`, `> 2020`, na.rm = T)
                                                         )

summary_metadata_long$Sample_Size = log(summary_metadata_long$Sample_Size,1.2)/max(log(summary_metadata_long$Sample_Size,1.2))

summary_metadata_long[is.na(summary_metadata_long)] = 0

unique(metadata[metadata$Country == 'Colombia',c('Subnational_level1', 'Subnational_level2')])

unique(metadata[metadata$Country == 'Colombia',][['Subnational_level2']])
unique(metadata[metadata$Country == 'Colombia',][['Subnational_level1']])

ggplot() +
  geom_polygon(data=lac0.spldf$spdf, aes(x=long, y=lat, group = group), fill = "floralwhite", color = "gray80") +
  geom_polygon(data=lac0.spldf[["spdf"]][lac0.spldf[["spdf"]][["NAME_0"]] %in% countries$country,], aes(x=long, y=lat, group = group), fill = col11[1], color = "gray48") +
  
  geom_polygon(data=gtm_sites[["spdf"]][gtm_sites[["spdf"]][["NAME_1"]] %in% c('Escuintla', 'Alta Verapaz', 'Izabal'),], aes(x=long, y=lat, group = group), fill = col11[4], color = "gray48") +
  geom_polygon(data=hnd_sites[["spdf"]][hnd_sites[["spdf"]][["NAME_1"]] == 'Gracias a Dios',], aes(x=long, y=lat, group = group), fill = col11[4], color = "gray48") +
  
  geom_polygon(data=col_sites[["spdf"]][col_sites[["spdf"]][["NAME_2"]] %in% c('Buenaventura', 'Quibdó', 'Guapí',
                                                                                 'Tumaco', 'Puerto Inírida','Tierralta',
                                                                                 'Santa Cecilia', 'Pichimá', 'Cali'),], aes(x=long, y=lat, group = group), fill = col11[4], color = "gray48") +
  
  
  # geom_polygon(data=col_sites[["spdf"]][col_sites[["spdf"]][["NAME_1"]] %in% c('Antioquia', 'Bolivar', 'Cordoba',
  #                                                                             'Chocó', 'Nariño','Guainía',
  #                                                                             'Risaralda', 'Valle del Cauca', 'Cauca')&
  #                                         col_sites[["spdf"]][["NAME_2"]] %in% c('Buenaventura', 'Quibdó', 'Guapi',
  #                                                                                'Tumaco', 'Puerto Inírida','Tierralta',
  #                                                                                'Santa Cecilia', 'Pichimá', 'Cali'),], aes(x=long, y=lat, group = group), fill = col11[4], color = "gray48") +

  geom_polygon(data=per_sites[["spdf"]][per_sites[["spdf"]][["NAME_1"]] %in% c('Amazonas'),], aes(x=long, y=lat, group = group), fill = col11[4], color = "gray48") +
  geom_polygon(data=per_sites2[["spdf"]][per_sites2[["spdf"]][["NAME_2"]] %in% c('Maynas'),], aes(x=long, y=lat, group = group), fill = col11[4], color = "gray48") +
  geom_polygon(data=ven_sites[["spdf"]][ven_sites[["spdf"]][["NAME_1"]] %in% c('Bolívar'),], aes(x=long, y=lat, group = group), fill = col11[4], color = "gray48") +  
#  geom_point(data=countries, aes(x=lon, y=lat), size = .1) +
#  geom_text(data=countries, aes(x=lon , y=lat + .2, label = country)) +
  theme_bw()+
  labs(#title = "Study sites",
    x = "Longitude",
    y = "Latitude")+
  geom_scatterpie(mapping = aes(x=lon, y=lat, group=Population, r = Sample_Size), data = summary_metadata_long, 
                  cols = colnames(summary_metadata_long)[5:7], alpha=.8) + 
  geom_scatterpie_legend(summary_metadata_long$Sample_Size, x=-87.5, y=-12.5)+
  coord_quickmap()

map_sites = last_plot()

View(summary_metadata_long %>% select(Country, Population, `< 2015`, `2015 - 2020`, `> 2020`))

ggplot()+
  geom_scatterpie(mapping = aes(x=lon, y=lat, group=Population, r = Sample_Size), data = summary_metadata_long, 
                  cols = colnames(summary_metadata_long)[5:7], alpha=.8) + coord_equal()
  


ggsave('map_sites.jpg', map_sites,device = 'jpeg')

#+
#scale_x_continuous(limits = c(-80, -65))+
#scale_y_continuous(limits = c(-5, 15))+
#scale_color_manual(values = localities$colors)+
# theme(legend.position = "none",
#       text = element_text(size = 10),
#       axis.text = element_text(size = 10),
#       title = element_text(size = 10))+
# theme(axis.text = element_text(size = 12),
#       axis.title = element_text(size = 12))




