# Step 1: Get all LAC countries names and codes
setwd("~/OneDrive/Github/G6PD_variants/")

countries_isoa3 <-as_tibble(read.csv("/ISO-3166-1-alpha-3.csv"))
# line 17 should be used in case download of all countries rds files are required

# countries.wbank <-as_tibble(read.csv("world-regions-according-to-the-world-bank.csv"))
# 
# lac.codes.wbank <- countries.wbank%>%
#   filter(World.Region.according.to.the.World.Bank == "Latin America and Caribbean")%>%
#   dplyr::select(Code)%>%
#   unlist()
# 
# lac.entity.wbank <- countries.wbank%>%
#   filter(World.Region.according.to.the.World.Bank == "Latin America and Caribbean")%>%
#   dplyr::select(Entity)%>%
#   unlist()

# Download and Generate a SpatialPolygonDataframe for all every LAC country

# All
# Download shape files from https://gadm.org/download_world.html,
# whole world shp files can be downloaded from here https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_shp.zip
# and whole world at different levels from here https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_levels_shp.zip
# Other formats like .rds are not currently supported for download https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_TUR_1_sp.rds

# Create a folder and download .rds for each country (level = 0)
system("mkdir rds0")

for (i in countries_isoa3$code){
  download.file(url = paste0("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_",i,"_0_sp.rds",sep=""),
                destfile = paste0("rds0/",i,"_adm0.rds",sep=""))
  rm(i)
  }

# SpatialPolygonsDataFrame at country level (level = 0)
world_spdf0 <- gadm_sp_loadCountries(countries_isoa3$code,level = 0, basefile = "./rds0/")
lac.spldf.0 <- gadm_sp_loadCountries(lac.codes.wbank,level = 0, basefile = "./rds0/")

# Create a folder and download .rds for each level 1 geographic administration by country (level = 1)
system("mkdir rds1")

for (i in countries_isoa3$code[c(-1,-4,-12,-28,-37,-41,-49,-55,
                                 -56,-75,-84,-98,-106,-121,-137,
                                 -139,-142,-144,-150,-151,-163,-166,-175, -196, -213, -237)]){
  download.file(url = paste0("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_",i,"_1_sp.rds",sep=""),
                destfile = paste0("rds1/",i,"_adm1.rds",sep=""))
}

# Create a folder and download .rds for each level 1 geographic administration by country (level = 2)
system("mkdir rds2")

for (i in countries_isoa3$code[c(-1,-4, -5, -7,  -10, -11, -12,-13, -14,-21, -25,-26, -28,-30,-31,-34, -37,-41,-49,-51, -52, -55,
                                 -56, -57, -58, -62, -69, -75,-78,-82, -84,-91,-92,-95, -97, -98,-106,-107,-111,-113, -114, -121,
                                 -122,-124,-128,-129,-130,-132,-137,-139,-140,-142,-144,-145,-149,-150,-151,-154,-155,-156,-159,
                                 -163,-166, -170,-171,-175, -178,-181,-186, -187, -192,-195,-196,-197,-198,-202, -204,-213,-214,
                                 -216,-221,-222,-224, -225, -228,-233,-237, -238, -240)]){
  download.file(url = paste0("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_",i,"_2_sp.rds",sep=""),
                destfile = paste0("rds2/",i,"_adm2.rds",sep=""))
}


# Create a folder and download .rds for each level 1 geographic administration by country (level = 3)
system("mkdir rds3")

for (i in countries_isoa3$code[c(-1, -2, -4, -5, -7, -9, -10, -11, -12,-13, -14, -15, -17, -20, -21, -24, -25,-26, -27, -28, -29,
                                 -30,-31,-34, -35, -36, -37, -38, -39, -41,-47,-48,-49,-50,-51, -52, -53, -54, -55,
                                 -56, -57, -58,-59, -61,-62,-63, -64, -65,-67, -68, -69, -73, -74, -75,-77, -78,-79,
                                 -81, -82, -83,-84,-86,-87,-88,-89,-91,-92,-93,-94,-95, -96,-97, -98,-99,-100,-102,-104,-106,
                                 -107,-108,-109,-110,-111,-113, -114, -115,-116,-117,-119,-121,-122,-123,-124,-125,
                                 -128,-129,-130,-131,-132,-133,-135,-136,-137,-139,-140,-142,-143,-144,-145,-147,-149,-150,-151,-153,
                                 -154,-155,-156,-158,-159,-160,-161,-163,-164,-165,-166,-167,-168, -170,-171,-172,-175, -178,-179,-180,-181,
                                 -182,-184,-185,-186, -187, -188,-189,-192,-195,-196,-197,-198,-199,-201,-202, -203,-204,-205,-207,-208,-209,
                                 -210,-211,-212,-213,-214,-215,-216,-217,-221,-222,-224, -225,-226, -227,-228,-232,-233,-234,-235,-236,-237, 
                                 -238, -239,-240,-241,-243,-244,-245,-246,-248,-249)]){
  download.file(url = paste0("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_",i,"_3_sp.rds",sep=""),
                destfile = paste0("rds3/",i,"_adm3.rds",sep=""))
}

#se quedo en canada

c(-1, -2, -4, -5, -7, -9, -10, -11, -12,-13, -14, -15, -17, -20, -21, -24, -25,-26, -27, -28, -29,
  -30,-31,-34, -35, -36, -37, -38, -39, -41,-47,-48,-49,-50,-51, -52, -53, -54, -55,
  -56, -57, -58,-59, -61,-62,-63, -64, -65,-67, -68, -69, -73, -74, -75,-77, -78,-79,
  -81, -82, -83,-84,-86,-87,-88,-89,-91,-92,-93,-94,-95, -96,-97, -98,-99,-100,-102,-104,-106,
  -107,-108,-109,-110,-111,-113, -114, -115,-116,-117,-119,-121,-122,-123,-124,-125,
  -128,-129,-130,-131,-132,-133,-135,-136,-137,-139,-140,-142,-143,-144,-145,-147,-149,-150,-151,-153,
  -154,-155,-156,-158,-159,-160,-161,-163,-164,-165,-166,-167,-168, -170,-171,-172,-175, -178,-179,-180,-181,
  -182,-184,-185,-186, -187, -188,-189,-192,-195,-196,-197,-198,-199,-201,-202, -203,-204,-205,-207,-208,-209,
  -210,-211,-212,-213,-214,-215,-216,-217,-221,-222,-224, -225,-226, -227,-228,-232,-233,-234,-235,-236,-237, 
  -238, -239,-240,-241,-243,-244,-245,-246,-248,-249)
