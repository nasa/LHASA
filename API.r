# 7-day ARI calculation v1.1.1
# 2021-8-2
# Thomas Stanley NASA GSFC/GESTAR/USRA
# Calculates a 7-day antecedent rainfall index

# Load R packages
library(raster)

# Set working directory
setwd('C:/LHASA')
# List all daily rainfall files
files <- list.files(path='IMERG', pattern='*liquid.tif', full.names=TRUE)

# Set antecedent rainfall window
ARI.window <- 7
# Set IDW exponent
exponent <- 2
# Calculate weights
w <- 1/seq(ARI.window, 1)^exponent
# Iterate through every day, starting at the end of the 1st ARI window
for(day in ARI.window:length(files)){
	# Open files within window, including current day
	IMERG <- crop(stack(files[(day - ARI.window + 1):day]), extent(-180, 180, -60, 60))
	# Calculate antecedent rainfall index
	ARI <- calc(w*IMERG,sum)/sum(w)
	# Save to disk
	writeRaster(ARI,filename=gsub('IMERG/','ARI/', files[day]))
}
