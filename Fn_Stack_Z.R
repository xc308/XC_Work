#=============
# 14 Mar. 2024
#=============

# Aim: 
  # Function to stack each Zi in df to become a joint Z

  # Z <- c(df$Z1, df$Z2, ...)

# Arg:
  # df: a data frame has variates name Z1, Z2 


#=========
# Function
#=========
Stack_Z <- function(df){
  
  Z <- c()
  for (i in 1:p) {
    Z <- c(Z, df[[paste0("Z", i)]])
  }
  
  return(Z)
}



#=====
# Test
#=====

#Z_try <- Stack_Z(df = df_2D_TW)
#str(Z_try)  # num [1:1200]



#Z <- c()
#for (i in 1:p) {
 # Z <- c(Z, df[[paste0("Z", i)]])
#}


#try <- df_2D_TW[["Z1"]]
#head(try)
