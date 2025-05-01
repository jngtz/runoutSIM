library(raster) # dealing with gridded data
#library(rgdal)
library(sp) # spatial data point, polygon, line
library(sf)
library(mapview)


# Functions ##############################################################


gMeanThreshold <- function(predict, response, plot_curve = FALSE) {
  # automatically determine threshold value to classify source areas
  #https://stackoverflow.com/questions/16347507/obtaining-threshold-values-from-a-roc-curve
  perf <- ROCR::performance(ROCR::prediction(predict, response), "tpr", "fpr")
  df <- data.frame(cut = perf@alpha.values[[1]], fpr = perf@x.values[[1]], tpr = perf@y.values[[1]])
  
  gmean <- sqrt(df$tpr * (1 - df$fpr))
  opt_threshold <- df[which.max(gmean), "cut"]
  
  if(plot_curve){
    
    plot(df$fpr, df$tpr, type = "l", ylab = "True positive rate", xlab = "False postive rate")
    points(df[which.max(gmean), "fpr"], df[which.max(gmean), "tpr"], pch = 19, col = "red")
    
  }
  
  return(opt_threshold)
  
}



# Read data for modeling #################################################

#setwd("/home/jason/Data/m_buchhart/")
#setwd("C:/Users/ku52jek/Google Drive/UniProjects/supervise/bachelor/max_buchhart/data")
#setwd("C:\\Users\\jgoetz\\OneDrive - Wilfrid Laurier University\\Documents\\GitProjects\\sedconnect\\Data")
setwd("~/Desktop/sda/sedconnect/data")


# Read Polygons and Source Points using RGDAL package
runout_polygons <- st_read("debris_flow_runout_polygons.shp")
source_points <- st_read("debris_flow_source_points.shp")

runout_polygons

river_channel <- st_read("river_channel.shp")

# Convert to sp objects
runout_polygons <- as_Spatial(st_zm(runout_polygons))
source_points <- as_Spatial(st_zm(source_points))
river_channel <- as_Spatial(st_zm(river_channel))

# Read raster using raster package
dem <- raster("elev.tif")
dem

# Note: This DEM is from the (free) publicly available ALOS PALSAR Radiometric Terrain Corrected
#       (RTC) high-resolution (12.5 m) digital elevation model.

#       ASF DAAC: ALOS PALSAR Radiometric Terrain Corrected high
#         resolution digital elevation model. Includes Material ©JAXA/METI
#         (2011). DOI: 10.5067/Z97HFCNKR6VA.

#       https://asf.alaska.edu/data-sets/derived-data-sets/alos-palsar-rtc/alos-palsar-radiometric-terrain-correction/

# Make a hillslope
slope <- terrain(dem, opt='slope')
aspect <- terrain(dem, opt='aspect')
hs <- hillShade(slope, aspect, 40, 270)

# Plot data
plot(hs, col=grey(0:100/100), legend = FALSE)
plot(dem, col=terrain.colors(25, alpha=0.35), add = TRUE)
plot(runout_polygons, add = TRUE)

# Read terrain attribute data
layer_files <- c("elev.tif",
                 "cslope.tif",
                 "slope.tif",
                 "logcarea.tif",
                 "plan_curv.tif",
                 "prof_curv.tif",
                 "swi.tif")

layers <- stack(layer_files)
plot(layers)

# Note: These terrain attributes were processed using SAGA-GIS.
#       First, a sink filled DEM was computed using the "Fill Sinks
#       (Planchon/Darboux, 2001)" tool. Using this sink-filled DEM,
#       the "SAGA Wetness Index" and the "Slope, Aspect, Curvature"
#       tools were used to compute the terrain attributes.

# Create mask of study area
mask_v <- getValues(layers$elev)
mask_v[!is.na(mask_v)] = 1
mask_area <- setValues(layers$elev, mask_v)

# Mask out runout polygons as well
mask_area <- mask(mask_area,runout_polygons, inverse = TRUE)
mask_area <- mask(mask_area, river_channel, inverse = TRUE)

# Create 
runout_area <- setValues(layers$elev, mask_v)
runout_area <- mask(runout_area, runout_polygons)
runout_area <- mask(runout_area, river_channel, inverse = TRUE)

# Create model sample ####################################################

# to provide data for our spatially predictive models we need to create
# sample of the response variable: slide and non-slides

# set the see of R's random number generator to reproduce the random samples
# if the script is run at a later time
set.seed(1234)

# first, we randomly sample cell locations within slides
smp.slides <- as.data.frame(sampleRandom(runout_area, size = 500, xy = TRUE))
smp.slides[,3] = NULL
plot(runout_area)
points(smp.slides)


# randomly sample cell locations of non-slide cells
smp.noslides <- as.data.frame(sampleRandom(mask_area, size = 500, xy = TRUE))
smp.noslides[,3] = NULL
points(smp.noslides)

# extract predictor variable values from grids
# for the slide samples
df.slides <- extract(layers, smp.slides, sp = TRUE)
df.slides <- as.data.frame(df.slides)
# add xy values
df.slides$x <- smp.slides$x
df.slides$y <- smp.slides$y
# add variable indicating that these are slides
df.slides$slide <- "TRUE"

# for the non-slide sample
df.noslides <- as.data.frame(extract(layers, smp.noslides, sp = TRUE))
df.noslides$x <- smp.noslides$x
df.noslides$y <- smp.noslides$y
df.noslides$slide <- "FALSE"

# combine the slides and no slides data into one dataframe
d <- rbind(df.slides, df.noslides)

# encode the slide variable as a factor
d$slide <- as.factor(d$slide)



# Exploratory analysis ###############################################

# To visualize the potential influencing of the predictor factors
# on slide occurence, we can produce spinograms, which show
# a plot of the conditional probabilites of lanslides occuring

summary(d[d$slide == TRUE,])
summary(d[d$slide == FALSE,])

# Transform curvature
d$plan_curv[d$plan_curv > 0.1] <- 0.1
d$plan_curv[d$plan_curv < -0.1] <- -0.1

d$prof_curv[d$plan_curv > 0.1] <- 0.1
d$prof_curv[d$plan_curv < -0.1] <- -0.1


par(mfrow=c(2,4),mex=0.7,mar=c(5,4,2,3))
with(d, {

  spineplot(slide ~ slope, ylevels = c("TRUE", "FALSE"),
            xlab = "Slope angle", ylab = "", yaxlabels = c("slide", "No landsl."))

  spineplot(slide ~ plan_curv, ylevels = c("TRUE", "FALSE"),
            xlab = "Plan curvature", ylab = "", yaxlabels = c("slide", "No landsl."))

  spineplot(slide ~ prof_curv, ylevels = c("TRUE", "FALSE"),
            xlab = "Profile curvature", ylab = "", yaxlabels = c("slide", "No landsl."))

  spineplot(slide ~ logcarea, ylevels = c("TRUE", "FALSE"),
            xlab = "Log. contributing area", ylab = "", yaxlabels = c("Lsl.", "No slide"))

  spineplot(slide ~ elev, ylevels = c("TRUE", "FALSE"),
            xlab = "Elevation", ylab = "", yaxlabels = c("slide", "No landsl."))

  spineplot(slide ~ cslope, ylevels = c("TRUE", "FALSE"),
            xlab = "Catchment slope", ylab = "", yaxlabels = c("slide", "No landsl."))

  spineplot(slide ~ swi, ylevels = c("TRUE", "FALSE"),
            xlab = "SAGA wetness index", ylab = "", yaxlabels = c("slide", "No landsl."))


} )


# Source prediction w GAM ################################################

library(mgcv)

# make formula
fml.gam <- slide ~ s(elev, k = 4) + s(slope, k = 4) + s(cslope, k = 4) +
  s(plan_curv, k = 4) + s(prof_curv, k = 4) + s(logcarea, k = 4) +
  s(swi, k = 4)

# fit a generalized additive model (GAM)
model.gam <- gam(fml.gam, family = binomial, data = d)

# visualize smoothing functions in the GAM
par(mfrow=c(2,4))
plot(model.gam)


# Model validation ######################################################

library(sperrorest)

gam.results <- sperrorest(formula = fml.gam, data = d, coords = c("x","y"),
                          model_fun = gam,
                          model_args = list(family=binomial),
                          pred_args = list(type="response"),
                          smp_fun = partition_kmeans,
                          smp_args = list(repetition = 1:10, nfold = 5))

gam.auroc.test <- unlist(summary(gam.results$error_rep,level=1)[,"test_auroc"])
mean(gam.auroc.test)
sd(gam.auroc.test)


# Apply prediction to raster ############################################

# Transform any layers
layers$plan_curv[layers$plan_curv > 0.1] <- 0.1
layers$plan_curv[layers$plan_curv < -0.1] <- -0.1

layers$prof_curv[layers$plan_curv > 0.1] <- 0.1
layers$prof_curv[layers$plan_curv < -0.1] <- -0.1


pred.gam<- raster::predict(layers, model.gam, type = "response", progress = TRUE)

# Mask out areas in river channel bead
maskpred.gam <- mask(pred.gam, river_channel, inverse = TRUE)
# export to a raster format
writeRaster(maskpred.gam, "src_pred_mask.tiff", format = "GTiff", overwrite = TRUE)


# Plot prediction map ####################################################

library(viridis) # for virdis colour palette

par(mfrow = c(1,2))
plot(hs, col=grey(0:100/100), legend = FALSE, main = "Release susceptibility")
plot(maskpred.gam, col=viridis(10, alpha=0.4, direction = -1), add = TRUE)
plot(runout_polygons, add = TRUE)
plot(source_points, add = TRUE)


# Automatically classifiy source areas ##########################################

src_class <- maskpred.gam(m)

d$pred <- as.numeric(predict(model.gam, type = "response", newdata = d))


# Determine optimital prediction threshold to classify source areas ##############
opt_threshold <- gMeanThreshold(est_src, obs_src, plot_curve = TRUE)

perf <- ROCR::performance(ROCR::prediction(d$pred, d$slide), "tpr", "fpr")
df <- data.frame(cut = perf@alpha.values[[1]], fpr = perf@x.values[[1]], tpr = perf@y.values[[1]])

gmean <- sqrt(df$tpr * (1 - df$fpr))
opt_threshold <- df[which.max(gmean), "cut"]
opt_threshold

# Make figure of results
par(mfrow = c(1,1))
png(filename = "gmean_source_class_roc.png", res = 300, width = 3, height = 3,
    units = "in", pointsize = 7)

plot(df$fpr, df$tpr, type = "l", ylab = "True positive rate", xlab = "False postive rate")
points(df[which.max(gmean), "fpr"], df[which.max(gmean), "tpr"], pch = 19, col = "red")

dev.off()

# Apply classification

src_class[src_class >= opt_threshold] <- 1
src_class[src_class < opt_threshold] <- 0
src_class[src_class == 0] <- NA

plot(source_class)

writeRaster(src_class, filename = "auto_classified_source_areas.tif", format = "GTiff",
            overwrite = TRUE)

plot(hs, col=grey(0:100/100), legend = FALSE, main = "Release classification")
plot(src_class, alpha = 0.6, add = TRUE)
plot(runout_polygons, add = TRUE)


# Make figures ###############################################################

library(ggplot2) # Plotting maps
library(patchwork) # Tile multiple maps together
library(ggnewscale) #
library(ggspatial)
library(sf)



# Read river polyline
river_channel <- readOGR(".", "river_channel")
river <- readOGR(".", "river_rio_olivares")
river_buff <- buffer(river, 30)
object = rasterize(river_buff, dem)

rivers <- st_read("river_rio_olivares.shp")
runout <- st_read("debris_flow_runout_polygons.shp")


source_cells <- src_class

# For visualization make a hillshade model from the DEM
slope <- terrain(dem, opt='slope')
aspect <- terrain(dem, opt='aspect')
hillshade <- hillShade(slope, aspect, angle=40, direction=270)

# Transform raster data into a data frame for use with ggplot
hillshade_df <- as.data.frame(hillshade, xy = TRUE)
hillshade_df <- hillshade_df[!is.na(hillshade_df[,3]),]


gppGGplot <- function(x){
  #parea_cdf <- rasterCdf(x)
  x <- as.data.frame(x, xy = TRUE)
  x <- x[!is.na(x[,3]),]
  return(x)
}

pred_df<- gppGGplot(maskpred.gam)
class_df<- gppGGplot(src_class)


map.pred <- ggplot() +
  geom_sf() +
  geom_raster(data=hillshade_df, aes(x=x, y=y, fill = layer),
              show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +
  
  #geom_sf(data = boundary, alpha = 0.4, fill = "white") +
  
  new_scale("fill") +
  geom_tile(data=pred_df, aes(x=x, y=y, fill = pred_df[,3]) ) +
  scale_fill_viridis_c(name = "Source area\npred.\n(Prob.)",  alpha = 0.6, direction = -1) +
  
  geom_sf(data = rivers, colour = "#85C1E9", fill = alpha("#85C1E9", 0.5)) +
  #geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  #geom_sf(data = boundary, colour = "#7f8c8d", fill = NA) +
  
  geom_sf(data = runout, color=alpha("black", 0.3), fill = alpha("white",0.2) , size = 0.3) +
  
  xlab("") +
  ylab("") +
  theme_bw() +
  annotation_scale(aes(style = "ticks", location = "br"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) +
  theme(text = element_text(size = 7), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6), axis.text.y = element_text(angle = 90),
        legend.key.width=unit(0.5,"cm"), legend.key.height = unit(0.5, "cm"))

map.class<- ggplot() +
  geom_sf() +
  geom_raster(data=hillshade_df, aes(x=x, y=y, fill = layer),
              show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +
  
  #geom_sf(data = boundary, alpha = 0.4, fill = "white") +
  
  new_scale("fill") +
  geom_tile(data=class_df, aes(x=x, y=y, fill = "") ) +
  scale_fill_manual(name = "Opt. source\nareas", values = "#162654" ) +
  
  geom_sf(data = rivers, colour = "#85C1E9", fill = alpha("#85C1E9", 0.5)) +
  #geom_sf(data = water, colour = "#85C1E9", fill = "#85C1E9") +
  #geom_sf(data = boundary, colour = "#7f8c8d", fill = NA) +
  
  geom_sf(data = runout, color=alpha("black", 0.3), fill = alpha("white",0.2) , size = 0.3) +
  
  xlab("") +
  ylab("") +
  theme_bw() +
  annotation_scale(aes(style = "ticks", location = "br"), text_cex = 0.5,
                   bar_cols = "black", line_width = 0.7) +
  theme(text = element_text(size = 7), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6), axis.text.y = element_text(angle = 90),
        legend.key.width=unit(0.5,"cm"), legend.key.height = unit(0.5, "cm"))

#setwd("C:/Users/ku52jek/Nextcloud/project/results/figures/report")

map.pred + ggtitle("Source area susceptibility model") +
  map.class + ggtitle("Optimal source area classification") + 
  plot_annotation(tag_levels = "A")
ggsave("map__olivares_srcpred_srcclass.png", dpi = 300, width = 7.4, height = 5, units = "in")
