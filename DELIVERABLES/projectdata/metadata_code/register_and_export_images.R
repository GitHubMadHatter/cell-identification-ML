library(EBImage)
library(ggplot2)
library(dplyr)

# Load the CSV containing cell boundaries
cell_boundaries <- read.csv("/Users/sgha9047/Downloads/GSE243280_RAW/outs/cell_boundaries.csv")
cell_boundaries$vertex_px <- cell_boundaries$vertex_x / 0.2125
cell_boundaries$vertex_py <- (max(cell_boundaries$vertex_y) - cell_boundaries$vertex_y) / 0.2125


smoothScatter(cell_boundaries[,c("vertex_px", "vertex_py")])
points(cell_boundaries[sample(nrow(cell_boundaries), 200000),c("vertex_px", "vertex_py")], cex = 0.1, pch = 16)

# threepts came from annotating the H&E image in napari
threepts = read.csv("/Users/sgha9047/Downloads/threepts.csv")

points(threepts[,c("axis.1", "axis.0")], cex = 10, pch = 16, col = "red")

# x = locator() # manually selected the three points using the same three points in the scatterplot in R, and saved
# saveRDS(x, file = "/Users/sgha9047/Downloads/threepts_cell_boundaries.Rds")
x = readRDS("/Users/sgha9047/Downloads/threepts_cell_boundaries.Rds")

points(x$x, x$y, cex = 10, pch = 16, col = "yellow")

image_points = as.matrix(threepts[,3:2])
csv_points = matrix(c(x$x, x$y), nrow = 3, byrow = FALSE)

library(Morpho)

tform <- computeTransform(image_points, csv_points, type = "affine")

registered_points <- applyTransform(as.matrix(cell_boundaries[, c("vertex_px", "vertex_py")]), tform)

# Add transformed coordinates to dataframe
cell_boundaries$registered_x <- registered_points[,1]
cell_boundaries$registered_y <- registered_points[,2]

cbr = cell_boundaries[, c("cell_id", "registered_y", "registered_x")]
colnames(cbr) <- c("index", "axis-0","axis-1")  
write.csv(cbr, file = "/Users/sgha9047/Downloads/GSE243280_RAW/cbr.csv", quote = FALSE,
          row.names = FALSE)

meta = read.csv("/Users/sgha9047/Downloads/41467_2023_43458_MOESM4_ESM.csv")

cbr_u = cbr[!duplicated(cbr[,1]),]

plot(cbr_u$`axis-1`, meta$x_centroid)
plot(cbr_u$`axis-0`, meta$y_centroid)

# Load the H&E image
image_path <- "/Users/sgha9047/Downloads/GSE243280_RAW/GSM7780153_Post-Xenium_HE_Rep1.ome.tif"
image <- readImage(image_path)

# subset to the pixels for cell i, and save out in the folders

for (ext in c(50, 100)) {
  for (i in seq_len(nrow(meta))) {
    
    if (i%%100 == 0) print(i)
    
    subpx = apply(apply(cbr[cbr[,1] == i,],2,range), 2,function(x) c(floor(x[1]) - ext, ceiling(x[2]) + ext))[,2:3]
    
    # can't be outside of the boundaries of the image
    subpx[1,1] <- pmax(subpx[1,1], 1) 
    subpx[2,1] <- pmin(subpx[2,1], dim(image)[2])
    subpx[1,2] <- pmax(subpx[1,2], 1) 
    subpx[2,2] <- pmin(subpx[2,2], dim(image)[1])
    
    image_sub = image[seq(subpx[1,2], subpx[2,2], by = 1),
                      seq(subpx[1,1], subpx[2,1], by = 1),]
    
    ct = gsub("&", "and", gsub(" ", "_", meta[i,"Cluster"]))
    
    img_name = paste0("cell_", i,"_", ext, ".png")
    folder_name = paste0("/Users/sgha9047/Downloads/GSE243280_RAW/images/",ext, "/", ct)
    
    if (!file.exists(folder_name)) {
      system(paste0("mkdir ", folder_name))
    }
    
    writeImage(image_sub, paste0(folder_name, "/", img_name))
  }
}

