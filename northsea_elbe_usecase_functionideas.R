## ONLY THREE VISUALS SELECTED TO INCLUDE
## DECISION: MAKE FUNCTIONS AT THEME/CHAPTER LEVEL

# Use case: North Sea: Elbe 3.3.1 
# ECRINS subcatchments with 
# OBJECTID 20, 119, 121, 122, 172, 424 
# DE for Germany

# This code calculate inputs weights based on 
# approximated use case catchment extent  
# for the use case 3.3.1
# - and provides dasymetric mapping results 
# at LAU (municipality) levels 
# and at river catchment levels for the area in focus
# the data has been pre-cut to extent in other script

# Chapters: 
### 1: PRESETTINGS  
### 2: FUNCTIONS  
### 3: INPUT DATA  
### 4: PRE-PROCESSING BEFORE WEIGHTING    
### 5: CALCULATE THE WEIGHTING  
### 6: DASYMETRIC REFINEMENT ANALYSIS  
### 7: AREA INTERPOLATION WITH LAU  
### 8: AREA INTERPOLATION WITH CATCHMENTS: ECRINS BASINS  
### 9: AREA INTERPOLATION WITH CATCHMENTS: ECRINS SUBBASINS  

#  Visual 1 ------------------------------------------------------> 
#  Print input weights based on all artificial surface 

#  Visual 2 ------------------------------------------------------> 
#  Map SUBSET of masked census grid based on all overlaps

#  Visual 3 ------------------------------------------------------> 
#  Map SUBSET of masked CORINE grid based on all overlaps
#  BBOX for subset is pre-defined as small_extent and bbox_polygon

#  Visual 4 ------------------------------------------------------> 
#  Print input weights from full overlaps of artificial surface 

#  Visual 5 ------------------------------------------------------> 
#  Map SUBSET of masked census grid based on full overlaps

#  Visual 6 ------------------------------------------------------> 
#  Map SUBSET of masked CORINE grid based on full overlaps
#  BBOX for subset is pre-defined as small_extent and bbox_polygon

###  Visual 7 ------------------------------------------------------> 
#  Print input weights: 

#  Visual 8 ------------------------------------------------------> 
#  Map NUTS population density + ECRINS river catchments + subcatchments: 

#  Visual 9 ------------------------------------------------------> 
#  Map selected CORINE classes with nuts coverage: 

#  Visual 10 ------------------------------------------------------> 
#  Map CORINE-specific input weights: 

#  Visual 11 ------------------------------------------------------> 
#  Map for each NUTS: sum of unique CORINE percent weights: 

#  Visual 12 ------------------------------------------------------> 
#  Map CORINE-specific normalised input weights

#  Visual 13 ------------------------------------------------------> 
#  Map area of intersect polygons: 

#  Visual 14 ------------------------------------------------------> 
#  Map area-based sum as denominator + numerator:  

#  Visual 15 ------------------------------------------------------> 
#  Map full Area-based weight:  

#  Visual 16 ------------------------------------------------------> 
#  Map polygon estimated populations   

###  Visual 17 ------------------------------------------------------> 
#  Map errors of dasymetric refinement in percentage at LAU level:   

#  Visual 18 ------------------------------------------------------> 
#  Map estimated population density at LAU level:   

#  Visual 19 ------------------------------------------------------> 
#  Map observed population density at LAU level:   

#  Visual 20 ------------------------------------------------------> 
#  Map estimated population density at catchment level

###  Visual 21 ------------------------------------------------------> 
#  Map estimated population density at subcatchment level


### 1: PRESETTINGS:

# install libraries: 
#install.packages("areal") 
#install.packages("classInt") 
#install.packages("dplyr") 
#install.packages("eurostat") 
#install.packages("foreign") 
#install.packages("giscoR") 
#install.packages("mapview") 
##install.packages("raster")
##install.packages("rlang") 
##install.packages("rstudioapi")  
##install.packages("sf")
##install.packages("terra")
##install.packages("units")

# HAS TO BE ALIGNED WITH DOCKER IMPORT 
#import libraries
library(areal) # used for areal interpolation
library(classInt) # used for classifying data into intervals for mapview
library(dplyr) # used for select and group and similar
library(eurostat) # used to get eurostat table data
library(foreign) # used for function read.dbf
library(giscoR) # used to get eurostat spatial data
library(mapview) # used to map data
library(raster) #raster package needed for mapview (but slow for calculations)  
library(rlang) # used to pass column information in functions
library(rstudioapi) # used to set workspace  
library(sf) # used to handle spatial vector objects
library(terra) # raster package for raster calculations
library(units) # used to handle units for e.g. area calculations

# PROPABLY NOT IMPORTANT
# Set directory to script directory: 
script_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_directory)

# CAN BE DELETED I THINK
# Set proper visualisation of numeric values
options(scipen = 100, digits = 4) # can be reset by options(scipen=0)
# Avoid switching to darkmatter basemap but keep default lightgrey basemap when plotting: 
mapviewOptions(basemaps.color.shuffle = FALSE)

# SHOULD BE PART OF VISUALISATION IF USED 
# Define colorpalette default: 
colors_default <- colorRampPalette(c("#FFFF00", "#990000"))  
colors_blue <- colorRampPalette(c("darkblue"))  
colors_contrast <- c("#f0f0f0", "#B6B6B6", "#9B9B9B","#676464",  "#474747")   
colors_dif <- c("white", "#f7fbff", "#deebf7", "#B2DDFC", "#9ecae1","#6baed6", "#2171b5","#084594", "#8734ba", "#5e1989", "#450d68", "#170324") #"#c6dbef",#"#4292c6", 
colors_negative <- colorRampPalette(c("#990000", "#FFFF00"))
colors_positive <- colorRampPalette(c("#90EE90", "darkblue"))


### 2: FUNCTIONS:

# DEFINE Function to extract raster values
get_raster_values <- function(raster_input) {
  # Extract raster values as a vector (removing NA values)
  raster_values <- values(raster_input, mat = FALSE)
  raster_values <- raster_values[!is.na(raster_values)]
  return(raster_values)
}

# USED AS PART OF VISUALISATION (IS ALL THE CODE USED IN THE SELECTED FUNCTIONS?)
# DEFINE Function to make classification break values: 
## Classification style: "equal", "quantile", "jenks", etc.
make_systematic_interval <- function(dataset_column_input, amount_intervals, style_type, rounded) {
  # Classify data using n intervals with the style method
  classifications <- classIntervals(dataset_column_input, n = amount_intervals, style = style_type)
  print(style_type)
  print(classifications)
  # Extract the breaks
  breaks <- classifications$brks
  print(breaks)
  # Calculate frequencies of data in each class
  frequencies <- tabulate(findInterval(dataset_column_input, breaks, rightmost.closed = TRUE))
  # Find breaks that result in non-empty classes
  non_empty_breaks <- breaks[frequencies > 0]
  # Round the breaks to divisible by 10 if rounded = TRUE otherwise round to three decimals after comma:
  if (rounded == TRUE) {
    non_empty_breaks <- round(non_empty_breaks, digits = -1)
  } else {
    non_empty_breaks <- round(non_empty_breaks, digits = 3)
  } 
  print(non_empty_breaks)
  # Add min(dataset_column_input) if lower than min(non_empty_breaks)
  if (min(dataset_column_input) < min(non_empty_breaks)) {
    non_empty_breaks[non_empty_breaks == min(non_empty_breaks)] <- min(dataset_column_input)
  }
  # Add max(dataset_column_input) if higher than max(non_empty_breaks)
  if (max(dataset_column_input) > max(non_empty_breaks)) {
    non_empty_breaks[non_empty_breaks == max(non_empty_breaks)] <- max(dataset_column_input)
  }
  non_empty_breaks <- sort(unique(non_empty_breaks))
  print(non_empty_breaks)
  return(non_empty_breaks)}

# VISUALISATION (IS IT USED?)
# DEFINE Function to make separate classifications for positive and negative: 
make_intervals_separate_pos_neg <- function(dataset_column_input, amount_intervals, style_type, rounded) {
  below_zero <- dataset_column_input[dataset_column_input < 0]
  classification_valuesP1 <- make_systematic_interval(below_zero, amount_intervals, style_type, rounded)
  above_zero <- dataset_column_input[dataset_column_input > 0]
  classification_valuesP2 <- make_systematic_interval(above_zero, amount_intervals, style_type, rounded)
  classification_values <- c(classification_valuesP1, classification_valuesP2)
  classification_values <- add_zeroclass_and_ensure_not_empty(dataset_column_input, classification_values)
  return(classification_values)}

## DEFINE Function to ensure interval classes next to value zero are not empty: 
add_zeroclass_and_ensure_not_empty <- function(dataset_column_input, classification_values) {
  # Step 1: Count how many values are in each interval
  freqs <- tabulate(findInterval(dataset_column_input, classification_values, rightmost.closed = TRUE))
  # Step 2: Keep only intervals with data
  non_empty_indices <- which(freqs > 0)
  valid_indices <- sort(unique(c(non_empty_indices, non_empty_indices + 1)))
  classification_values_cleaned <- classification_values[valid_indices]
  # Step 3: Add 0 *only* if it's BETWEEN data-containing intervals
  min_val <- min(dataset_column_input)
  max_val <- max(dataset_column_input)
  # Include 0 only if there are both negative and positive values
  if (min_val < 0 & max_val > 0) {
    # Only add 0 if it's not already present
    if (!0 %in% classification_values_cleaned) {
      classification_values_cleaned <- sort(c(classification_values_cleaned, 0))
    }
    # Remove any breaks in intervals around zero if no data exists there
    idx_zero <- which(classification_values_cleaned == 0)

    # Check the interval below 0
    if (idx_zero > 1) {
      low_limit <- classification_values_cleaned[idx_zero - 1]
      high_limit <- 0
      low_interval <- dataset_column_input[dataset_column_input > low_limit & dataset_column_input < high_limit]
      if (length(low_interval) == 0) {
        classification_values_cleaned <- classification_values_cleaned[-(idx_zero - 1)]
        idx_zero <- which(classification_values_cleaned == 0)
      }
    }
    # Check the interval above 0
    if (idx_zero < length(classification_values_cleaned)) {
      low_limit <- 0
      high_limit <- classification_values_cleaned[idx_zero + 1]
      high_interval <- dataset_column_input[dataset_column_input > low_limit & dataset_column_input < high_limit]
      if (length(high_interval) == 0) {
        classification_values_cleaned <- classification_values_cleaned[-(idx_zero + 1)]
      }
    }
  }
  return(classification_values_cleaned)}

# IS DEFINITELY USED FOR THE AREA INTERPOLATION
# DEFINE Function to area interpolate at river catchment level
get_catchment_data <- function(catch, tid_input, source_data, sid_input, areaKm2_column, colors_default) {
  # aw_interpolate to interpolate to LAU level: 
  catch_with_pop <- aw_interpolate(
    .data = catch, #referred to as target elsewhere
    tid = !!rlang::sym(tid_input), #id column in target
    source = source_data,
    sid = !!rlang::sym(sid_input), #id column in source
    weight = "sum", #should be sum for intensive vs. extensive interpolations
    output = "sf",
    extensive = "estPopCor" #attribute in source to interpolate 
  )
  # replace NA with 0 in estimated pop
  catch_with_pop$estPopCor[is.na(catch_with_pop$estPopCor)] <- 0
  # KEY DENSITY COLUMN: estimated LAU population density in km2:  
  catch_with_pop$densKm2 <- as.numeric(catch_with_pop$estPopCor / catch_with_pop[[areaKm2_column]])
  return(catch_with_pop)
}


# IS HANDLED DIFFERENTLY FROM THE DATA LAKE 
### 3: INPUT DATA: 

## 1) NUTS WITHOUT POPULATION 
## Data is automatically available through the eurostat 
## and giscoR R packages
#Get NUTS3 shapes using eurostat
nuts3_all <- get_eurostat_geospatial(
  output_class = "sf",
  resolution = "01",
  nuts_level = "3",
  year = "2016",
  crs = "3035"
)

## 2) Get Eurostat population table
## Data is automatically available through the eurostat 
## and giscoR R packages
#freq: #unique(poptable$freq) "A" 
#sex: #unique(poptable$sex) "F" "M" "T"
#unit #unique(poptable$unit) "NR" :number;unitOfMeasure
#age #unique(poptable$age) "TOTAL"  "UNK"    "Y10-14" "Y15-19"...
#geo #unique(poptable$geo) 
#TIME_PERIOD: 2018-01-01 
#values: population 
poptable <- get_eurostat("demo_r_pjangrp3") %>%  
  filter(TIME_PERIOD == "2018-01-01",
         nchar(geo)==5, 
         sex == "T",
         age == "TOTAL")

## 3) LAU POP 2018  
## CURRENTLY IT HAS JUST BEEN MANUALLY DOWNLOADED AS VECTOR 
## FROM HERE: https://ec.europa.eu/eurostat/web/gisco/geodata/statistical-units/local-administrative-units
## AND HAS BEEN JOINED TOGETHER IN ARCGIS PRO WITH POPULATION 
## STATISTICS DOWNLOADED AS EXCELSHEET FROM HERE: https://ec.europa.eu/eurostat/web/nuts/local-administrative-units
LAUpop2018 <- st_read("elbe/LAUpop2018DE.gpkg") # SF VECTOR

## 4) CENSUS GRID AS VECTOR
## Eurostat data link that requires manual download: 
## https://ec.europa.eu/eurostat/web/gisco/geodata/population-distribution/geostat
## In use here: Newest census grid version of 16th of June 2024
census_grid_country <- st_read("elbe/censusDE_catchSE.gpkg") # SF VECTOR

## 5) CORINE LAND COVER AS RASTER: 
## Data link visible in DDAS metadata that requires authorisation: 
## DDAS link: https://aquainfra.dev.52north.org/result/eea:e640a22b-8dde-423e-b051-991bcb383f2f
## Data link: https://sdi.eea.europa.eu/data/9b2a406b-320b-43f8-81a7-8d3899c56607
#Get terra raster: 
cor2018_country_r <- rast("elbe/cor2018DE_catchSE.tif") # TERRA RASTER

## 6) CORINE LAND COVER AS VECTOR: 
## Data link is visible in DDAS metadata that requires 
## authorisation: DDAS link: https://aquainfra.dev.52north.org/result/eea:2828681f-1d19-4d53-a815-677c6fd03c80
## Data link: https://sdi.eea.europa.eu/data/2828681f-1d19-4d53-a815-677c6fd03c80
#CORINE input: 
cor2018_country <-  st_read("elbe/corDE_nutsSE.gpkg") # SF VECTOR

## 7) CORINE CATEGORY IDs
## Data link visible in DDAS metadata that requires authorisation: DDAS link: https://aquainfra.dev.52north.org/result/eea:e640a22b-8dde-423e-b051-991bcb383f2f
## Data link: https://sdi.eea.europa.eu/data/9b2a406b-320b-43f8-81a7-8d3899c56607
#Read legend to get CODE_18
corlegend <- read.dbf("elbe/cor2018DE_catchSE.tif.vat.dbf")
# Convert CODE_18 factor to character and then to integer
corlegend$CODE_18 <- as.integer(as.character(corlegend$CODE_18))  # Now it's numeric
corlegend$Value <- as.integer(corlegend$Value)  # Now it's numeric

## 8) CATCHMENT DATA: ECRINS basins and sub basins
## Ecrins manual data download: 
## https://www.eea.europa.eu/en/datahub/datahubitem-view/a9844d0c-6dfb-4c0c-a693-7d991cc82e6e?activeAccordion=814%2C1273 
## Downloaded: "Natural sub basins of Europe - version 0 Dec. 2011": 
catch_ecrins_rough <- st_read("elbe/catch_ecrins_northsea_elbeSE.gpkg") # very rough
## Downloaded also: "Functional elementary catchments (fec) - version 1 Jun. 2012":
catch_ecrins_detailed <- st_read("elbe/catchsub_ecrins_northsea_elbeSE.gpkg") # detailed


### 4: PRE-PROCESSING BEFORE WEIGHTING  
# IF PREPROCESSING IS TURNED INTO ONE FUNCTION IT HAS THREE OVERALL TENDENCIES: 
# IT MAKES GEOMETRY VALID
# IT CUTS DATASETS SPATIALLY TO EXTENTS (COUNTRY AND CATCHMENT)
# 

# FUNCTION: DO PRE-PROCESSING TO THE RIVER CATCHMENT DATA 
# Fix catchment geometry if needed: 
valid_geoms <- st_is_valid(catch_ecrins_detailed)
if (any(!valid_geoms)) {
  catch_ecrins_detailed[!valid_geoms, ] <- st_make_valid(catch_ecrins_detailed[!valid_geoms, ])
}
# Remove empty or NULL geometries
catch_ecrins_detailed <- catch_ecrins_detailed[!st_is_empty(catch_ecrins_detailed), ]
# Extract only POLYGON or MULTIPOLYGON types
catch_ecrins_detailed <- catch_ecrins_detailed[st_geometry_type(catch_ecrins_detailed) %in% c("POLYGON", "MULTIPOLYGON"), ]

# FUNCTION: SELECT AND COMBINES EUROSTAT DATA FOR COUNTRY ---#SELECT COUNTRY FROM EU NUTS DATA; FROM EU TABLE DATA; AND COMBINES IT  
## NUTS: to country 
# Filter the dataset where CNTR_CODE == 
nuts3_country <- nuts3_all %>% #Simple feature collection with 11 features and 11 fields
  filter(CNTR_CODE == "DE")
# NUTS3 pop: to country: 
poptable_country <- poptable[grep("^DE", poptable$geo), ]
## Get target data: NUTS3 population spatial data through merge: 
nuts3pop <- nuts3_country %>%
  left_join(poptable_country, by = c("NUTS_ID" = "geo")) 
#NUTS_ID is NUTS ID in the gisco spatial dataset
#geo is NUTS ID in the Eurostat table

# FUNCTION: GET SPATIAL EXTENT OF THE CATCHMENT IN FOCUS 
# spatial extents
analysis_spatial_catch_extent <- st_union(catch_ecrins_rough)

# Select country NUTS polygons that intersect the catchment
nuts3pop <- nuts3pop[st_intersects(nuts3pop, 
                                   analysis_spatial_catch_extent, 
                                   sparse = FALSE), ]

# Turn the NUTS selection into main spatial extent: 
analysis_spatial_extent <- st_union(nuts3pop)

# Find LAU polygons that intersect the selected nuts
LAUpop2018_sub <- LAUpop2018[st_intersects(LAUpop2018, 
                                             analysis_spatial_extent, 
                                       sparse = FALSE), ]
LAUpop2018 <- st_intersection(LAUpop2018_sub, 
                              analysis_spatial_extent)
# Check LAU for valid geometries
valid_geoms <- st_is_valid(LAUpop2018)
# Try to fix invalid geometries
if (any(!valid_geoms)) {
  LAUpop2018[!valid_geoms, ] <- st_make_valid(LAUpop2018[!valid_geoms, ])
}
# Remove empty or NULL geometries
LAUpop2018 <- LAUpop2018[!st_is_empty(LAUpop2018), ]
# Extract only POLYGON or MULTIPOLYGON types
LAUpop2018 <- LAUpop2018[st_geometry_type(LAUpop2018) %in% c("POLYGON", "MULTIPOLYGON"), ]
# force geometry to be multipolygon
LAUpop2018 <- st_cast(LAUpop2018, "MULTIPOLYGON")


### 5: CALCULATE THE WEIGHTING:
# FUNCTION ON THE WEIGHTING CALCULATION 

## Weighting step 1: 
## WEIGHTING PROCEDURE BASED ON ALL CORINE-CENSUS OVERLAPS

#Find raster values where CODE_18 is between 100 and 199 
#(inclusive lower, exclusive upper)
cor_urban_values <- corlegend$Value[corlegend$CODE_18 >= 100 & corlegend$CODE_18 < 200]

#Logical (TRUE/FALSE vector) identifying which raster cells 
#have values found in the specific list
mask_cor_logical <- cor2018_country_r[] %in% cor_urban_values

#create raster based on mask (urban=True): 
mask_cor_raster <- setValues(rast(cor2018_country_r), 
                             mask_cor_logical)

# Apply mask on corine raster to keep artificial surfaces: 
cor2018_artificial <- mask(cor2018_country_r, 
                           mask_cor_raster, 
                           maskvalue = FALSE)

# Convert census grid to terra vector: 
census_grid_country_vector <- vect(census_grid_country)

# Count number of non-NA CORINE pixels in each census grid polygon
# (census grid is 1000x1000 m and CORINE raster is 100x100 m cellsize)
counts <- extract(cor2018_artificial, 
                  census_grid_country_vector, 
                  # Count non-NA raster values in each polygon: 
                  fun = function(x, ...) sum(!is.na(x)), 
                  # raster cells that touch the polygon: 
                  touches = TRUE) 

# Add the result as a new second census grid column (first column is id): 
census_grid_country_vector$count_corine_pixels <- counts[,2]

# Create raster with count of urban CORINE cells within each 
# census grid polygon (they do not have to be unique classes)  
count_raster <- terra::rasterize(census_grid_country_vector, 
                        cor2018_artificial, 
                        field = "count_corine_pixels")

# Create raster with census grid population masked to 
# only be with values for artificial surface CORINE categories
# (use max value if overlapping polygons within cell):
census_raster <- terra::rasterize(census_grid_country_vector, 
                                  cor2018_artificial, 
                                  field = "T", # population value
                                  fun = "max") # max value if overlaps

# Mask census raster where count is zero
# (this removes 31,1769,631 - 28,2099,421 = 29,670,210 
# which is 9.517 % of population size)
census_raster[count_raster == 0] <- NA

# Divide census values by the CORINE urban class count per census grid
census_raster_correctedF1 <- census_raster / count_raster

# Also mask the CORINE artificial surface grid correspondingly: 
cor2018_country_maskedF1 <- mask(cor2018_artificial, 
                                 census_raster_correctedF1)

## Statistics table 1: 
# Average population per 100x100 m per CORINE urban category:
avg_pop_per_corineF1 <- zonal(census_raster_correctedF1, 
                            cor2018_country_maskedF1, 
                            fun = "mean",  
                            na.rm = TRUE) # ignore NA values
# Summing the mean population density across all urban CORINE classes: 
total_avg_sumF1 <- sum(avg_pop_per_corineF1$T)
# Add sum of all mean cases as column in statistics table: 
avg_pop_per_corineF1$percentF1 <- round(avg_pop_per_corineF1$T / total_avg_sumF1 * 100, 2)


## Weighting step 2: 
## WEIGHTING PROCEDURE BASED ON FULL OVERLAPS WHERE 
## ONE CENSUS GRID IS MADE COMPLETELY UP BY A SPECIFIC 
## URBAN CORINE CATEGORY: 

# Reset converting census grid to terra vector: 
census_grid_country_vector <- vect(census_grid_country)

# Extract unique CORINE classes for each census grid cell
count_unique_corine_classes <- extract(cor2018_country_r, 
                                       census_grid_country_vector,
                                       # Count number of unique CORINE classes in census grid: 
                                       fun = function(x, ...) length(unique(x)), 
                                       # raster cells that touch the polygon:
                                       touches = TRUE)

# Add unique CORINE class count as a new column in the census grid vector (first column is ID)
census_grid_country_vector$unique_corine_classes <- count_unique_corine_classes[,2]  

# Set all values that are NOT 1 to NA in the 'unique_corine_classes' column
census_grid_country_vector$unique_corine_classes[census_grid_country_vector$unique_corine_classes != 1] <- NA

# Create raster with unique CORINE class count 
mask_logical_full_raster <- terra::rasterize(census_grid_country_vector, 
                                      cor2018_artificial, 
                                      field = "unique_corine_classes", 
                                      fun = "min")

# Create raster with census grid population masked to 
# only be with values for artificial surface CORINE categories
# (use max value if overlapping polygons within cell):
census_rasterF2 <- terra::rasterize(census_grid_country_vector, 
                                   cor2018_country_r, 
                                   field = "T", # population values
                                   fun = "max") 

# Divide census values by 100 due to this factor difference 
# in cell size between the census grid and the CORINE raster 
# to adjust the population values to the CORINE raster extent
census_raster_correctedF2 <- census_rasterF2 / 100

# Mask the CORINE cells that are NOT a unique class taking the space 
# of the whole census grid cell: ensuring the full overlap this way
# ensures complete mapping of census population value 
# to CORINE category:
census_raster_correctedF2 <- mask(census_raster_correctedF2, 
                                   mask_logical_full_raster, 
                                   maskvalue = NA)

# Also mask the CORINE artificial surface grid correspondingly: 
cor2018_country_maskedF2 <- mask(cor2018_artificial, 
                                  census_raster_correctedF2)

## Statistics table 2: 
# Average population per CORINE cell per CORINE category
avg_pop_per_corineF2 <- zonal(census_raster_correctedF2, 
                              cor2018_country_maskedF2, 
                              fun = "mean", 
                              na.rm = TRUE)

# summing the mean population density across all CORINE classes: 
total_avg_sumF2 <- sum(avg_pop_per_corineF2$T)

# Add sum of mean as column
avg_pop_per_corineF2$percentF2 <- round(avg_pop_per_corineF2$T / total_avg_sumF2 * 100, 2)


## Statistics table combined: 

# Combine the two tables: 
avg_pop_per_corine_combined <- avg_pop_per_corineF1 %>%
  left_join(avg_pop_per_corineF2, by = c("LABEL3" = "LABEL3"))

# Replace 0 values in 'percentF2' with 0.1 to 
# separate these observations of empty population during full 
# overlap of same CORINE category within census grid 
# - separate from cases with no existing full overlaps 
if (any(avg_pop_per_corine_combined$percentF2 == 0)) {
  avg_pop_per_corine_combined$percentF2[avg_pop_per_corine_combined$percentF2 == 0] <- 0.1
}

# Calculate the mean of 'percentF1' and 'percentF2' for categories 
# with non-NA values - create a new column 'percent' which is the 
# mean of 'percentF1' and 'percentF2' for each row 
avg_pop_per_corine_combined$percent <- apply(avg_pop_per_corine_combined[, c("percentF1", 
                                                                             "percentF2")], 
                                             1, # apply over rows 
                                             function(x) {
  # Check if both 'percentF1' and 'percentF2' are non-NA
  if (all(!is.na(x))) {
    return(mean(x, # take mean
                na.rm = TRUE)) # do not consider NA values
  } else {
    # If one value is NA, keep the non-NA value
    return(na.omit(x)[1])
  }
})

# Merge the input weight with the CORINE legend info
weight_table_full <- avg_pop_per_corine_combined %>%
  left_join(corlegend, by = c("LABEL3" = "LABEL3"))
# Select only CODE_18, percent, and LABEL3 columns 
weight_table <- weight_table_full %>%
  dplyr::select(CODE_18, percent, LABEL3)

# Filter rows where percent >= 1.0 to only keep well-contributing 
# CORINE classes as weights
weight_table_final <- weight_table %>%
  filter(percent >= 1.0)

# WE KEEP THIS ONE
#  Visual 7 ------------------------------------------------------> 
#  Print input weights: 
cat("Final weight table is:\n", paste(capture.output(print(weight_table_final)), collapse = "\n"))
############################


### 6: DASYMETRIC REFINEMENT ANALYSIS:
# TURNED INTO A SUBFUNCTION 

# new AREA COLUMN: calculate NUTS area in km2:
nuts3pop$area_nuts <- set_units((st_area(nuts3pop)), km^2)
# new DENSITY COLUMN: calculate overall population density 
# per Ha for each nuts:  
nuts3pop$nutsDens2018 <- nuts3pop$values / (nuts3pop$area_nuts) 

# Only pick population-related subset of CORINE data based on populated classes: 
cor2018_subset <- cor2018_country[cor2018_country$Code_18 %in% weight_table_final$CODE_18,]
# Convert to numeric: 
cor2018_subset$Code_18 <- as.integer(cor2018_subset$Code_18)

# Get the legend colors for the CORINE categories in focus: 
corlegend_subset <- corlegend[corlegend$CODE_18 %in% weight_table_final$CODE_18,]
corlegend_subset$color <- rgb(corlegend_subset$Red, corlegend_subset$Green, corlegend_subset$Blue)

# Merge the input weight with the adjusted CORINE object
cor_detailed <- cor2018_subset %>%
  left_join(weight_table_final, by = c("Code_18" = "CODE_18"))

# Calculate intersection of relevant CORINE classes and nuts population (first intersection): 
cor_detailed_filtered <- cor_detailed[st_intersects(cor_detailed, analysis_spatial_extent, sparse = FALSE), ] # spatial filter
nuts_cor_intersect <- st_intersection(cor_detailed_filtered, nuts3pop)

# new ID column for all unique CORINE-NUTS combinations
nuts_cor_intersect$nuts_cor_id <- interaction(nuts_cor_intersect$NUTS_ID, nuts_cor_intersect$Code_18)

# Get unique CORINE weights per NUTS ID 
##first is used to avoid double-counting CORINE-NUTS combinations 
## if more polygons of the same CORINE category should exist per NUTS 
table_unique_weights <- nuts_cor_intersect %>%
  as_tibble() %>% #to get it as a table and not a spatial output
  dplyr::select(nuts_cor_id, NUTS_ID, percent) %>%
  mutate(percent = as.numeric(percent)) %>%
  group_by(nuts_cor_id) %>%
  summarize(percUnique = first(percent), 
            NUTS_ID = first(NUTS_ID)) 

# Sum unique input weights per NUTS ID 
table_sum_of_weights <- table_unique_weights %>%
  dplyr::select(NUTS_ID, percUnique) %>%
  group_by(NUTS_ID) %>%
  summarize(percSum = sum(percUnique)) 

# Merge the weight sum with the intersect data
nuts_cor_weights <- nuts_cor_intersect %>%
  left_join(table_sum_of_weights, by = "NUTS_ID")

# Remove areas that are not NUTS: 
nuts_cor_weights <- nuts_cor_weights %>%
  filter(!is.na(NUTS_ID) & NUTS_ID != "")

# Estimate final CORINE category input weight:
nuts_cor_weights$newWeight <- (as.numeric(nuts_cor_weights$percent) / #specific corine weight
                             nuts_cor_weights$percSum * # sum of all corine weights present in that nuts
                             100.) #percentage

# new AREA COLUMN: calculate intersection area in km2:
# be aware: multiple polygons will exist for same CORINE category within each NUTS
nuts_cor_weights$area_cor <- set_units((st_area(nuts_cor_weights)), km^2)

#new column with the normalised corine weight multiplied the corine intersect area:  
nuts_cor_weights$area_corNum <- as.numeric(nuts_cor_weights$area_cor)
nuts_cor_weights$areaWeight <- nuts_cor_weights$newWeight * nuts_cor_weights$area_corNum

# Get sum of area-based weights for each NUTS 
table_area_weights <- nuts_cor_weights %>%
  as_tibble() %>%
  mutate(areaWeight = as.numeric(areaWeight)) %>%
  filter(!is.na(NUTS_ID)) %>%
  group_by(NUTS_ID) %>%
  summarize(areaWeightSum = sum(areaWeight, na.rm = TRUE))
# Get sum of area-based weights for each NUTS 
table_area_weights <- nuts_cor_weights %>%
  as_tibble() %>% #to get it as a table and not a spatial output
  dplyr::select(NUTS_ID, areaWeight) %>%
  group_by(NUTS_ID) %>%
  summarize(areaWeightSum = sum(areaWeight)) 

# Merge the calculated area-based weights data with the processed NUTS data
nuts_cor_final <- nuts_cor_weights %>%
  left_join(table_area_weights, by = "NUTS_ID")

# Merge the calculated area-based weights data with the original NUTS data for visualisation purposes
nuts_cor_visual <- nuts3pop %>%
  left_join(table_area_weights, by = "NUTS_ID")

# New column with full area-based weights:  
nuts_cor_final$areaWeightFull <- nuts_cor_final$areaWeight / nuts_cor_final$areaWeightSum 

# New column with estimated population per corine polygon:  
nuts_cor_final$estPopCor <- as.numeric(nuts_cor_final$values) * as.numeric(nuts_cor_final$areaWeightFull)


### 7: AREA INTERPOLATION WITH LAU
# NEW FUNCTION WHICH USES THE PRODUCED ANCILLARY DATA POLYGONS TO PRODUCE AREA-WEIGHTED RESULTS AT LAU LEVEL

# new ID column for all intersects 
nuts_cor_final$nuts_cor_int_id <- interaction(nuts_cor_final$NUTS_ID, nuts_cor_final$ID_corPol)

# aw_interpolate to interpolate to LAU level: 
estimated_lau_pop <- aw_interpolate(
  .data = LAUpop2018, #referred to as target elsewhere
  tid = "LAU_ID", #id column in target
  source = nuts_cor_final,
  sid = "nuts_cor_int_id", #id column in source
  weight = "sum", #should be sum for intensive vs. extensive interpolations
  output = "sf",
  extensive = "estPopCor" #attribute in source to interpolate 
)

# replace NA with 0 in estimated pop
estimated_lau_pop$estPopCor[is.na(estimated_lau_pop$estPopCor)] <- 0

# Calculate percentage overestimation (marked with plus) or underestimation (marked with minus) of estimated population 
estimated_lau_pop$pop2018dif <- ((-(estimated_lau_pop$POP_2018 - as.numeric(estimated_lau_pop$estPopCor)) / estimated_lau_pop$POP_2018) *100.)

#Check geometries
st_geometry_type(estimated_lau_pop)
# Keep only POLYGON and MULTIPOLYGON geometries
estimated_lau_pop <- estimated_lau_pop[st_geometry_type(estimated_lau_pop) %in% c("POLYGON", "MULTIPOLYGON"), ]
estimated_lau_pop <- st_collection_extract(estimated_lau_pop, "POLYGON")
estimated_lau_pop <- st_cast(estimated_lau_pop, "MULTIPOLYGON")

# Create a new column with absolute values before plotting
estimated_lau_pop$pop2018dif_abs <- abs(estimated_lau_pop$pop2018dif)

# WE USE THIS VISUAL
#  Visual 17 ------------------------------------------------------> 
#  Map errors of dasymetric refinement in percentage at LAU level:   
# Specify final classification and colors 
classification_values <- c(0, 5, 10, 15, 20, 30, 50, 60, 70, 80, 90, 100)
# replace 100 if max(estimated_lau_pop$pop2018dif_abs > 100: 
if (max(estimated_lau_pop$pop2018dif_abs, na.rm = TRUE) > 100) {
  classification_values[classification_values == 100] <- max(estimated_lau_pop$pop2018dif_abs, na.rm = TRUE)
}
map_list <- c(mapview(estimated_lau_pop, # dataset
                      zcol = "pop2018dif_abs", # values
                      col.regions = colors_dif, # input colors
                      alpha.regions = 0.9, # almost solid
                      at = classification_values, # intervals
                      layer.name="Errors in %")) # layer name
Reduce(`+`, map_list)
############################

# KEY DENSITY COLUMN: estimated LAU population density in km2:  
estimated_lau_pop$densKm2Est <- as.numeric(estimated_lau_pop$estPopCor) / estimated_lau_pop$AREA_KM2 


### 8: AREA INTERPOLATION WITH CATCHMENTS: ECRINS BASINS --- CAN BE LEFT OUT 
# NEW FUNCTION WHICH USES THE PRODUCED ANCILLARY DATA POLYGONS TO PRODUCE AREA-WEIGHTED RESULTS AT RIVER SUB-CATCHMENT LEVEL

# KEY DENSITY COLUMN: estimated LAU population density in km2:  
# converts area from hectare to km2 by dividing with 1000000
catch_ecrins_rough$AreaKm2 <- as.numeric(catch_ecrins_rough$Shape_Area / 1000000.)

# Interpolate ecrins catchment 
catch_ecrins_with_pop <- get_catchment_data(catch_ecrins_rough, 
                                     tid_input = "SB", #WSO_ID
                                     source_data = nuts_cor_final, 
                                     sid_input = "nuts_cor_int_id", 
                                     areaKm2_column = "AreaKm2", 
                                     colors_default)

### 9: AREA INTERPOLATION WITH CATCHMENTS: ECRINS SUBBASINS

# KEY DENSITY COLUMN: estimated LAU population density in km2:  
# converts area from hectare to km2 by dividing with 1000000
catch_ecrins_detailed$AreaKm2 <- as.numeric(catch_ecrins_detailed$Shape_Area / 1000000.)

# Interpolate ecrins catchment 
subcatch_ecrins_with_pop <- get_catchment_data(catch_ecrins_detailed, 
                                            tid_input = "OBJECTID", #"WX02ID"
                                            source_data = nuts_cor_final, 
                                            sid_input = "nuts_cor_int_id", 
                                            areaKm2_column = "AreaKm2", 
                                            colors_default)

# FINAL VISUALISATIONS 
#  Visual 21 ------------------------------------------------------> 
#  Map estimated population density at subcatchment level
classification_values <- make_systematic_interval(subcatch_ecrins_with_pop$densKm2, 8, "jenks", TRUE)
classification_values <- c(0.01, classification_values)
classification_values <- sort(classification_values)
map_list <- c(mapview(subcatch_ecrins_with_pop, # dataset
                      zcol = "densKm2", # values
                      col.regions = c("#FFFBE2", colors_default(length(classification_values) - 1)), # input colors 
                      alpha.regions = 0.9, # almost solid
                      at = classification_values, # intervals
                      layer.name = paste("Estimated ecrins catchment population density"))) # layer name
Reduce(`+`, map_list)
############################

### END OF SCRIPT ######################################