r <- raster("C:/Users/diamo/Dropbox/Projects/EAB/delineation_code/output/D1_hum_ho_binary.tif")
plot(r)
df <- rasterToPoints(r)
df <- as.data.frame(df)
colnames(df)[3] <- "hu_ho"
mod_binom = glm(hu_ho ~ x+y, data = df, family = binomial)  
summary(mod_binom)

library(spectral)
x <- sort(unique(df$x))
y <- sort(unique(df$y))

n = length(unique(df$x))
m = length(unique(df$y))
mat = matrix(0, nrow=m, ncol=n)
for (i in nrow(df)){
  mat[df[i,"x"], df[i,"y"]] = df[i,"hu_ho"]
}
library(kzfs)
# Raw periodogram (2d spectrum)
rp <- kzp2(mat)
rp_sum <- kzp2
fy <- rp$freq.y
fx <- rp$freq.x
rp <- rp$kzp2d
# smoothing 2D spectrum 2 times
sp <- smooth.kzp2(rp,0.01,k=2)
par(mfrow=c(2,1), cex=0.5)
persp(x=fx, y=fy, z=rp, expand =0.5,
      main = "Raw 2D KZ Periodogram", ltheta=40, shade=0.75,
      theta=-30, phi=15, zlab="",xlab="x", ylab="y",
      ticktype="detailed", col="red")
persp(x=fx, y=fy, z=sp, expand =0.5,
      main = "Smoothed 2D KZ Periodogram", ltheta=40, shade=0.75,
      theta=-30, phi=25, zlab="",xlab="x", ylab="y",
      ticktype="detailed", col="red")
par(mfrow=c(1,1), cex=1)
kzp2.summary(sp) # direction & frequency




FT <- spec.fft(x = round(abs(min(x)) + x, 2), 
               y = round(abs(min(y)) + y, 2),
               z = mat)
# plot
par(mfrow = c(2, 1))
rasterImage2(x = x,
             y = y,
             z = m,
             main = "Propagating Wave")
plot(FT)

# x <- round((sort(unique(df$x)) + abs(min(df$x))) * 1000, 0)
# y <- round((sort(unique(df$y)) + abs(min(df$y))) * 1000, 0)
# m <- matrix(NA, length(x) + 1, length(y) + 1)
# m[cbind(x, y)] <- df$hu_ho
# 
# mat <- matrix(NA, nrow=max(df[[1]])+1, ncol=max(df[[2]])+1 ) 
# mat[data.matrix(df[,1:2])] <- df[,3] 
# mat[cbind(x, y)] <- df$hu_ho
# mat2 <- data.matrix(df)
# 
# mtx <- matrix(NA, nrow=length(unique(df$x)) + 2, ncol=length(unique(df$y)) + 2 )
# mtx[cbind(order(df$x), order(df$y))] <- df$hu_ho




# Test with data from package

x <- seq(0, 1, length.out = 1e2)
y <- seq(0, 1, length.out = 1e2)
# calculate the data
m <- matrix(0, length(x), length(y))
for (i in 1:length(x))
  for (j in 1:length(y))
    m[i, j] <- sin(4 * 2 * pi * x[i] + 10 * 2 * pi * y[j])
# calculate the spectrum
FT <- spec.fft(x = x, y = y, z = m)
# plot
par(mfrow = c(2, 1))
rasterImage2(x = x,
             y = y,
             z = m,
             main = "Propagating Wave")
plot(
  FT,
  main = "2D Spectrum",
  palette = "wb"
  ,
  xlim = c(-20, 20),
  ylim = c(-20, 20),
  zlim = c(0, 0.51)
  ,
  xlab = "fx",
  ylab = "fy",
  zlab = "A",
  ndz = 3,
  z.adj = c(0, 0.5)
  ,
  z.cex = 1
)


# Test with data
x <- sort(unique(df$x))
y <- sort(unique(df$y))

n = length(unique(df$x))
m = length(unique(df$y))
mat = matrix(0, nrow=m, ncol=n)
for (i in nrow(df)){
  mat[df[i,"x"], df[i,"y"]] = df[i,"hu_ho"]
}


x <- sort(unique(round(abs(min(df$x)) + df$x, 2)))
y <- sort(unique(round(abs(min(df$y)) + df$y, 2)))

x2 <- seq(0, 51.3, length.out = 360)
y2 <- seq(0, 54, length.out = 180)

mat2 <- as.matrix(r)
# calculate the data
m2 <- mat
# smooth data for irregular data
library(fields)
mat_smooth <- image.smooth(mat2)
mat3 <- mat_smooth$z
mat3 <- mat3[complete.cases(mat3),complete.cases(t(mat3))]

# calculate the spectrum
FT <- spec.fft(x = x2, y = y2, z = mat3)
# plot
par(mfrow = c(2, 1))
rasterImage2(x = x2,
             y = y2,
             z = mat2,
             main = "Propagating Wave")
plot(
  FT,
  main = "2D Spectrum",
  palette = "wb"
  ,
  # xlim = c(-20, 20),
  # ylim = c(-20, 20),
  # zlim = c(0, 0.51)
  # ,
  xlab = "fx",
  ylab = "fy",
  zlab = "A",
  # ndz = 3,
  # z.adj = c(0, 0.5)
  # ,
  z.cex = 1
)
dev.off()