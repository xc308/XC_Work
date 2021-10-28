#*****************#
# Form_Z function
#*****************#

#' @title Form_Z function
#' @name Form_Z
#' @aliases bisquare_1d
#'
#' @description This function concatenates different variables of a data frame into one vector. 
#' @param h displacement (1d)
#' @details It also allows a reverse version of the order of variables formed into one vector hence create a reverse version of the intended model.
#' @export
#' @examples



#str(df_Cmpt_Res_log_16$BC_Residuals_log) # num [1:27384]
#str(df_Cmpt_Res_log_16[, 4])  # num [1:27384]

#str(matrix(df_Cmpt_Res_log_16$BC_Residuals_log))  # num [1:27384, 1]
#str(df_Cmpt_Res_log_16[4])    # 'data.frame':	27384 obs. of  1 variable:


Form_Z <- function(df, col_num1, col_num2, model_num) {
  
  Z1 <- matrix(df[, col_num1])
  Z2 <- matrix(df[, col_num2])
  
  
  if (model_num > 2 | model_num <= 0) {
    stop("model_num can only be 1 or 2")
  }
  
  
  if (model_num != 1) {
    ## reverse version 
    temp <- Z1
    Z1 <- Z2
    Z2 <- temp
  }
  
  
  Z <- rbind(Z1, Z2)
}



