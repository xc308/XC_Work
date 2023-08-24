#=====================
# Checking parent node  (15 Aug 2023)
#=====================

#rm(list = ls())


# Aim: 
  # Select parent nodes for my automation algo

# Steps:
  # 1. Construct the Hierarchical Structure:
      # a data frame with columns like "node_id", "par_id".
  
  # 2. Define a Function:
      # input: desired query input node_id; data with hierarchical structure
      # output: the parent node id corresponds to the input one



#----------------------------------------------------------
# function to check if one node is a parent node of another
#----------------------------------------------------------

# input: 1. query node_id; 
       # 2. data with hierarchy structure.
# output: parent node id for the query node


Check_par_node <- function(Node, data = hierarchy_data) {
  
  if (!all(colnames(data) == c("node_id", "par_id"))) {
    stop("Data must have two columns with the first col name 'node_id' and the second name 'par_id'")
  }
  
  
  # Filter rows where node_id equals desired query node
  rows_matching_condition <- hierarchy_data$node_id == Node
  
  # Extract rows that match the condition
  par_index <- hierarchy_data[rows_matching_condition, ]$par_id

  
  ifelse(is.na(par_index), return("The query node has no parent."), return(par_index))
  
}


#-----
# Test
#------
#Check_par_node(Node = 3, data = hierarchy_data) # [1] 2
#Check_par_node(Node = 4, data = hierarchy_data) # [1] 1


#hierarchy_data2 <- data.frame(
  #node_id = c(1, 2, 3, 3, 4, 4),
  #par_id = c(NA, 1, c(2, 1), c(2, 3))
#)

#Check_par_node(Node = 3, data = hierarchy_data2)
#res <- Check_par_node(Node = 4, data = hierarchy_data2)
#str(res) # num [1:2] 2 3 a vector
#res[1]
#res[2]


#mach_cond <- hierarchy_data2$node_id == 3
#Par_index <- hierarchy_data2[mach_cond, ]$par_id




