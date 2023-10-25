#=============
# 25 Oct. 2023
#=============

# Aim:
  # to generate the corresponding parameter matrices for further use
  # given hierarchy data, produce the parameter matrices


Para_A_mat <- function(p, data){
  A <- matrix(0, p, p)
  for (r in 2 : p) {
    PN = Check_par_node(r, data = data)
    for (t in c(PN)){
      A[r, t] <- paste0("A", r, t)
    }
  }
  A
}


Para_Dlt_mat <- function(p, data){
  Dlt <- matrix(0, p, p)
  for (r in 2 : p) {
    PN = Check_par_node(r, data = data)
    for (t in c(PN)){
      Dlt[r, t] <- paste0("dlt", r, t)
    }
  }
  Dlt
}


#-----
# Test
#-----

hierarchy_data <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4)
)

p <- 5
A_mat <- Para_A_mat(p = 5, data = hierarchy_data)
A_mat[2, 1]
str(A_mat[3, 2]) # chr "A32"




hierarchy_data2 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 5))
)

Para_Dlt_mat(p = 6, data = hierarchy_data2)





