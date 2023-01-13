## ---------------------------------------------------------------------------
## Function : gadm_loadtCountries (constructor)
## Description : load a file from local system or from GADM repository 
##               You just have to specify the countries (ISO3 CODE) of the
##               file name, like "ARG" for Argentina.
##               Optionally you can specify which level you want to have and
##               simplify or not theshapefile (a value less or equal 0.01
##               is recommended)
## Return : This function creates a gadm_sp that contains a 
##          SpatialPolygonsDataFrame object that contains all maps you 
##          specify in "fileNames".
## ---------------------------------------------------------------------------
gadm_sp_loadCountries <- function (fileNames,
                                   level = 0,
                                   basefile=GADM_BASE,
                                   baseurl=GADM_URL,
                                   simplify=NULL)
{
  loadNamespace("sp")  
  
  # Load file and change Prefix ---------------------------------------------
  loadChangePrefix <- function (fileName, level = 0) {
    FILENAME = sprintf("%s_adm%d.rds", fileName,level)
    LOCAL_FILE = sprintf("%s%s", basefile, FILENAME)
    if (!file.exists(LOCAL_FILE)) {
      .OS <- toupper(Sys.info()["sysname"])
      REMOTEFILE = sprintf("gadm36_%s_%d_sp.rds", fileName,level)
      REMOTE_LINK <- sprintf("%s%s", baseurl, REMOTEFILE)
      if (.OS == "WINDOWS") {
        download.file(REMOTE_LINK, LOCAL_FILE, method="wininet",mode="wb")
      } else {
        download.file(REMOTE_LINK, LOCAL_FILE, method = 'auto')
      }
    }
    gadm <- readRDS(LOCAL_FILE)
    if (!is.null(gadm)) {
      theFile <- spChFIDs(gadm, paste(fileName, row.names(gadm), sep = "_"))
      theFile
    }
  }
  
  polygon <- sapply(fileNames, loadChangePrefix, level)
  polyMap <- do.call("rbind", polygon)
  # ---- Simplify polygones if requested by user
  if (!is.null(simplify)) {
    S <- gSimplify(polyMap, simplify, topologyPreserve = TRUE)
    polyMap@polygons <- S@polygons
  }
  
  # ---- Create gadm_sp object
  structure(list("basename"=basefile,
                 "spdf"=polyMap,
                 "level"=level,
                 "holes" = c(),
                 "L360" = FALSE,
                 "stripped" = FALSE,
                 "hasBGND"  = FALSE),
            class = "gadm_sp")
}
