#*********#
# Sp_T_CV
#********#
install.packages("caret")
Y
library(caret)
library("dplyr")
library("tidyr")

install.packages("spatialrisk") # for radius circle calculate
library(spatialrisk)



str(df_all)
df_all # 136920 = 27384 * 5
df_all_long # 821520 = 136920 * 6 components
tail(df_all)
head(df_all)


df_all_16 <- filter(df_all, Year == 2016)
head(df_all_16)

nfolds <- 10

foldlabels <- rep(1:nfolds, length.out = length(unique(df_all_16[, "ID"])))
str(foldlabels)
tail(foldlabels) # [1]  9 10  1  2  3  4

all(length(unique(df_all_16[, "ID"])) == nrow(unique(df_all_16[, c("Lon", "Lat")])))



head(sample(foldlabels))
head(df_all_16)

set.seed(1)
df_all_16$foldlabs <- sample(foldlabels)
head(df_all_16, 10)



function(data, formulae.list, nfolds = 10, radius = 80000) {
  
  formuae <- lapply(formulae.list, as.formula)
  
  #df_all_16$ID <- seq(nrow(unique(df_all_16[, c("Lon", "Lat")]))) # add ID idx if not exist
  foldlabels <- rep(1:nfolds, length.out = nrow(unique(df_all_16[, "ID"])))
  df_all_16$foldlabs <- sample(foldlabels)
  
  
  MSE.all.folds <- matrix(NA, nrow = nfolds, ncol = length(formulae))
  
  for(fold in nfolds) {
    
    test.rows <- which(df_all_16$foldlabs == fold) # 2739
    test.dat <- df_all_16[test.rows, ]
    train.dat <- df_all_16[-test.rows, ] # 24645
    
    mse.per.test <- matrix(NA, nrow = nrow(test.dat), ncol = length(formulae))
    prediction <- list()
    current.mod <- list()
    new.train.dat <- list()
    circle <- list()
    for(i in seq_along(test.dat[, "Lon"])) {
      circle[[i]] <- points_in_circle(train.dat, lon_center = test.dat[, "Lon"][i],
                                      lat_center = test.dat[, "Lat"][i], 
                                      lon = Lon, lat = Lat,
                                      radius = 80000)
      
      print(length(circle[[i]]$ID)) # the # of points within radius
      
      new.train.dat[[i]] <- train.dat[-circle[[i]]$ID, ]
      
      
      for(j in seq_along(formulae)) {
        
        current.mod[[i]] <- lm(formula = formulae[[j]], data = new.train.dat[[i]])
        
        prediction[[i]] <- predict(current.mod[[i]], data = test.dat[i, ])
        
        mse.per.test[i, j] <- mean((test.dat[i, ]$BC - prediction[[i]]) ^ 2) ## Modify BC
      }
      
      MSE.all.folds[fold, ] <- colMeans(mse.per.test) # each row element in MSE.per.fold is the 
    }
    
    return(colMeans(MSE.all.folds))
    
  }
}

  


test.dat$BC

# tail(seq_along(test.dat[, "Lon"]))
str(circle)
length(circle) # 2739
(new.train.dat)



