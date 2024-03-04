#============
# 4 Mar. 2024
#============

# Aim:
  # investigate the boundary for dlt

# Problem:
  # s = seq(-10 + ds/2, 10 - ds/2, by = ds)
  # displacement H range: -19.9  19.9
  # when dlt set to 0.5 - 20, p.d. are ensured
  # when dlt set to 21, encounter error, see 034c


# conjecture:
  # dlt has upper boundary that is dlt must be smaller than the largest abs(displacement)
  # i.e. dlt < max(abs(H))

# Now want to have a look at the b function when dlt = 21

source("Fn_Wendland_32.R")

B_WL_dlt21 <- WendLd_32(r = H, R = 0.5, dlt = 21, A = 1)

all(B_WL_dlt21 == 0) # [1] TRUE

View(WendLd_32)
range(H) #  -19.9  19.9

# Because R = 0.5, if dlt = 21, then the displacment that will have 
  # a separation with 21 is smaller than 0,5 must be 21 +/- 0.5
  # ie. the displacement must be at least 20.5 to ensure WL function value
  # is not zero,
  # however, the largest H is 19.9, so all the WL value = 0

# conclusion:
  # dlt has upper bound max(abs(H))







