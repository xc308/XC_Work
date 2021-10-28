#*****************#
# Form_Z function
#*****************#

#' @title Form_Z function
#' @name Form_Z
#' @aliases bisquare_1d
#'
#' @description This function concatenates different variables of a data frame into one joint vector. 
#' @param df data frame 
#' @param col_num1 The 1st col want to extract from df to be the 1st variable in the concatenated joint vector.
#' @param col_num2 The 2nd col want to extract from df to be the 2nd variable in the concatenated joint vector.
#' @param model_num It can only be 1 or 2. If model_num == 2, then concate a reversed version.
#' @details 
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



