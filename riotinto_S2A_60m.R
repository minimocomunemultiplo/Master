# We will work now on 60 m S2A resolution images, so we will have more bandwidths to deal with and we can elaborate a more detailed spectrum.

#### SETUP #####
setwd("C:/Rio Tinto program/60m script")
install.packages("raster")
install.packages("rasterVis")
install.packages("ggplot2")
install.packages("rgdal")
install.packages("rgeos")
install.packages("jpeg")
install.packages("viridis")
library(raster)
library(rasterVis)
library(ggplot2)
library(rgdal)
library(rgeos)
library(jpeg)
library(viridis)

# The SENTINEL-2A instrument acquires measurements at 12 bits. These measurements are converted to reflectances and stored as 16 bit 
# integers in the S2 product. The Level-2A processing includes a Scene Classification and an Atmospheric Correction applied to 
# Top-Of-Atmosphere (TOA) Level-1C orthoimage products. Level-2A main output is an orthoimage atmospherically corrected, Surface 
# Reflectance product.


####### DATASET ######

rlist <- list.files(pattern="T29SQB_20230710T110621_B") # rlist let's us take all files with a certain name.
rlist <- c(rlist[1:7], "T29SQB_20230710T110621_B8A_60m.jp2", rlist[8:10]) # this reorders the list. Band 8a was in the wrong place.
rlist
import <- lapply(rlist,raster) # lapply applies a function (in this case raster) to all elements from rlist.
rt23_60 <- stack(import) # stack function puts files imported as raster in one variable (memory reference) in R.
rt23_60 #let's show our result:

spectrum_cover <- c("B1 443 nm", "B2 490 nm", "B3 560 nm", "B4 665 nm", 
                    "B5 705 nm", "B6 740 nm", "B7 783 nm", "B8a 865 nm", 
                    "B9 940 nm", "B11 1610 nm", "B12 2190 nm")

plotRGB(rt23_60,3,2,1,stretch="lin") # Check if plotRGB works for true colors
plotRGB(rt23_60,8,4,3,stretch="lin") # Check if plotRGB works for false colors
plot(rt23_60[[8]])

ext <- c(710000,715000,4175500,4180500) # the crop is larger than the one I calculated on the 10 m spatial res.
rt23_60 <- crop(rt23_60,ext)
plotRGB(rt23_60,4,3,2,stretch="lin") # the crop is good enough. We have to stay small, otherwise my pc explodes.
ext_rt23_60 <- extent(rt23_60)
ext_rt23_60

###### DATA CLEAN ######

# NDVI CALCULATION = (Band NIR – Band RED) / (Band NIR + Band RED).
NDVI_60 <- (rt23_60[[8]]-rt23_60[[4]])/(rt23_60[[8]]+rt23_60[[4]])
NDVI_60
plot(NDVI_60, col = viridis(100))

threshold_60 <- 0.35
mask_60 <- NDVI_60 > threshold_60
plot(mask_60)
getValues(mask_60)
rt23_60_masked <- rt23_60 
rt23_60_masked[mask_60==TRUE] <- NA 
plotRGB(rt23_60_masked, 4,3,2, stretch = "lin") # dense vegetation should be masked now.

# Now we need reflectance values. Apply the formula given by Sentinel2A online manual 
# L2A_SR = (L2A_DN + BOA_ADD_OFFSET) / QUANTIFICATION_VALUE

BOA <- -1000
QV <- 10000
rt23_60_masked_SR <- (rt23_60_masked + BOA)/QV
rt23_60_masked_SR
getValues(rt23_60_masked_SR)
plot(rt23_60_masked_SR[[8]]) # plot 8a band for quick n easy look.

rt23_60_masked_SR[rt23_60_masked_SR < 0] <- NA # remove incoherent negative numbers
rt23_60_masked_SR

IR_c <- colorRampPalette(c("white", "yellow", "orange","red","black"))(100) # REMEMBER BLACK IS NOW HIGH REFLECTANCE
par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))
for (i in 1:11) {
  plot(rt23_60_masked_SR[[i]], 
       main = paste("Band", i, spectrum_cover[[i]]),
       col = IR_c)
}

dev.off()



###### SPECTRUM EXTRACTION ######

# there is a very interesting area on the west side of the map. It could be anything, but it strongly reflects SWIR radiation. 
# let's extract coordinates to crop out a smaller section identifying that area a bit closer.
plot(rt23_60_masked_SR[[8]])

ext2 <- c(710000,711500,4177000,4179000) # the crop is larger than the one I calculated on the 10 m spatial res.
rt23_60_interest <- crop(rt23_60_masked_SR,ext2)

# I will plot using all SWIR wavelengths, to check what reflects NIR (bluish) and what 
# reflects other deeper SWIR frequencies.
par(mfrow = c(1, 3))
SWIRbands_60 <- c(8,10,11) # this is needed because plot function is stupid and needs a pre-defined array to take the position of elements.
plot(rt23_60_interest[[SWIRbands_60]], 
     col=IR_c)
dev.off() 
plotRGB(rt23_60_interest,11,10,8,stretch="lin") 
dev.off()

# let's now extract the spectrum of the average. Look carefully at this piece of script, because it is a good thing to use again.
# I'm telling him FOR EVERY ELEMENT in the rt23 object, take the values, remove the NA's and calculate a mean. Then store every 
# result for every element of the object in a DEDICATED POSITION inside a vector (means_rt23...).

for(i in 1:11) {
  means_rt23_60_interest[i] <- mean(getValues(rt23_60_interest[[i]]),na.rm=T)
}
means_rt23_60_interest

# NOW PLOT 

plot(1:length(spectrum_cover), means_rt23_60_interest,
     main = "Value of mean reflectance for the interest area", # Create a plot with x representing the layer and y representing the pixel values
     type = "b", pch = 19, 
     xlab = "wavelength", 
     ylab = "Reflectance Value",
     xaxt = "n")
axis(1, at = 1:length(spectrum_cover), labels = spectrum_cover, las = 2) 

dev.off()
