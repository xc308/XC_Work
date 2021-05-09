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


#============#
# Review data
#============#

str(df_all)
df_all # 136920 = 27384 * 5
df_all_long # 821520 = 136920 * 6 components
tail(df_all)
head(df_all)


df_all_16 <- filter(df_all, Year == 2016)
head(df_all_16)


foldlabels <- rep(1:nfolds, length.out = length(unique(df_all_16[, "ID"])))
str(foldlabels)
tail(foldlabels) # [1]  9 10  1  2  3  4

all(length(unique(df_all_16[, "ID"])) == nrow(unique(df_all_16[, c("Lon", "Lat")])))

head(sample(foldlabels))
head(df_all_16)

set.seed(1)
df_all_16$foldlabs <- sample(foldlabels)
head(df_all_16, 10)

nfolds <- 10

#========#
# Sp_CV
#========#

formulae.list <- list(BC ~ Lon + Lat + I(Lat ^ 2))
data <- df_all_16

Sp_CV <- function(data, formulae.list, nfolds = 10, radius = 80000) {
  
  formulae <- lapply(formulae.list, as.formula)
  
  set.seed(1)
  #df_all_16$ID <- seq(nrow(unique(df_all_16[, c("Lon", "Lat")]))) # add ID idx if not exist
  foldlabels <- rep(1:nfolds, length.out = length(unique(data[, "ID"])))
  data$foldlabs <- sample(foldlabels)
  
  
  Drop.all.folds <- matrix(NA, nrow = nfolds, ncol = 1)
  MSE.all.folds <- matrix(NA, nrow = nfolds, ncol = length(formulae))
  
  for(fold in 1:nfolds) {
    
    test.rows <- which(data$foldlabs == fold) # 2739
    test.dat <- data[test.rows, ]
    train.dat <- data[-test.rows, ] # 24645
    
    drop.per.test <- vector()
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
      
      if(length(circle[[i]]$ID) != 0) {
        
        print(length(circle[[i]]$ID))
        new.train.dat[[i]] <- train.dat[-circle[[i]]$ID, ]
      } else {
        
        print(length(circle[[i]]$ID))
        new.train.dat[[i]] <- train.dat
      }
      
      #length(circle) # [1] 2739 elements (df)
      #length(new.train.dat) # [1] 2739
      
      drop.per.test[i] <- length(circle[[i]]$ID) # collect # of dropped data
      
      
      for(j in seq_along(formulae)) {
        
        current.mod[[i]] <- lm(formula = formulae[[j]], data = new.train.dat[[i]], na.action=na.exclude)
        prediction[[i]] <- predict(current.mod[[i]], newdata = test.dat[i, ])
        mse.per.test[i, j] <- mean((test.dat[i, ]$BC - prediction[[i]]) ^ 2) ## Modify BC
      }
      
      MSE.all.folds[fold, ] <- colMeans(mse.per.test) # each row element is MSE.per.fold averaging all test per fold 
      
      # within each fold, the percent of dropped data vs total train.dat
      Drop.all.folds[fold, ] <- mean(drop.per.test[i] / nrow(train.dat) * 100)
      
    }
  }
    return(list(MSE.all.folds, colMeans(MSE.all.folds),Drop.all.folds))
}

# 80000
#> MSE.all.folds
#          [,1]
#[1,] 0.7612503
#[2,] 1.0836280
#[3,] 0.7292751
#[4,] 0.9430062
#[5,] 0.8724286
#[6,] 0.9136737
#[7,] 0.8830653
#[8,] 0.8379615
#[9,] 0.9438970
#[10,] 1.0130373
> 

# 4000000
# MSE.all.folds
#          [,1]
#[1,] 0.7812077
#[2,] 1.1009696
#[3,] 0.7483764
#[4,] 0.9638596
#[5,] 0.8925590
#[6,] 0.9313602
#[7,] 0.9001969
#[8,] 0.8604037
#[9,] 0.9610333
#[10,] 1.0287968  

# 0.9168763  
  
  
  
  
#==========#
# Naive CV
#==========#



formulae.list <- list(BC ~ Lon + Lat + I(Lat ^ 2))
data <- df_all_16

function(data, formulae.list, nfolds = 10) {
  
  formulae <- lapply(formulae.list, as.formula)
  
  set.seed(1)
  #df_all_16$ID <- seq(nrow(unique(df_all_16[, c("Lon", "Lat")]))) # add ID idx if not exist
  foldlabels <- rep(1:nfolds, length.out = length(unique(data[, "ID"])))
  data$foldlabs <- sample(foldlabels)
  
  MSE <- matrix(NA, nrow = nfolds, ncol = length(formulae))
  #predictions <- list()
  #current.mod <- list()
  for(fold in 1:nfolds) {
    
    test.rows <- which(data$foldlabs == fold)
    test.dat <- data[test.rows, ]
    train.dat <- data[-test.rows, ]
    
    for(j in 1:length(formulae)) {
      current.mod <- lm(formula = formulae[[1]], data = train.dat)
      predictions <- predict(current.mod, newdata = test.dat)
      MSE[fold, j] <- mean((test.dat$BC - predictions) ^ 2)
    }
  }
  
  return(list(MSE, colMeans(MSE)))
}


# > MSE
#          [,1]
#[1,] 0.7612087
#[2,] 1.0836054
#[3,] 0.7292609
#[4,] 0.9429905
#[5,] 0.8723973
#[6,] 0.9136563
#[7,] 0.8830476
#[8,] 0.8379348
#[9,] 0.9438697
#[10,] 1.0130044

# > colMeans(MSE)
# [1] 0.8980976

test.dat$BC

# tail(seq_along(test.dat[, "Lon"]))
str(circle)
length(circle) # 2739


