#===============================
# Wave Functions (23 Nov. 2023)
#===============================

# wave_v5: a further modified tri-wave 
  # modifications:
        # 1. smaller compact support range: /h - dlt/ < /dlt/ (not <= 2/dlt/)
        # 2. faster decay of correlation: 2(/h - dlt/ / /dlt/)^2 (not 1/2 (.)^2)
  # advantages:
        # smaller compact support and faster decay of correlation within the spprt
        # result in a contract matrix for p.d.
  # Arg:
        # h: displacement
        # delta: width and skewness
        # A: amplitude
  

# wave_v4: a function of Triangular wave_v4
    # modified on top of wave_v3

# Modification advananges: 
# (read in conjunction with spectral_wave4_vs_tent_report.R)
    # 1. compact support range different
    # 2. values can be +/-/0
    # 3. shapes of asymmetry are different for different delta (shift):
      # more versatile shape of asymmetry won't make b_{lk}(s, s) be
      # constant zero
    # 4. only 1 parameter, more scalable
    # 5. spectral density have more energy than tent 
      # esp. when delta = +/-0.1, when tent spectral has 0 energy
      # across all frequency, fail to capture asymmetry details
      # while wave function are able to capture fine asymmetry details
      # and can capture different shape of asymmetry that changes with locations
        # 1. an intrinsic anisotropic stationary 
        # 2. more pertinent to real cross-correlation



#--------
# wave_v8 (faster decay: phi = 2.5, narrower suppt // < 1/2 //)
#--------
wave_v9 <- function(h, delta, A = 1) {
  
  A * (1 - 2.5 * (abs(h - delta) / abs(delta))^2) * (abs(h-delta) <= 1/2 *abs(delta)) 
}


#--------
# wave_v8 (faster decay: phi = 2, // < 1/2 //)
#--------
wave_v8 <- function(h, delta, A = 1) {
  
  A * (1 - 2 * (abs(h - delta) / abs(delta))^2) * (abs(h-delta) <= 1/2 * abs(delta)) 
}





#---------
# wave_v7 (fast decay phi = 2, larger suppt)
#---------
wave_v7 <- function(h, delta, A = 1) {
  A * (1- 2*(abs(h - delta) / abs(delta))^2)*(abs(h - delta) <= 2 * abs(delta))
}




#---------
# wave_v6 (slow decay phi = 1/2)
#---------
wave_v6 <- function(h, delta, A = 1) {
  A * (1- (1/2)*(abs(h - delta) / abs(delta))^2)*(abs(h - delta) <= abs(delta))
}



#--------
# wave_v5 (fast decay: phi = 2)
#--------
wave_v5 <- function(h, delta, A = 1) {
  
  A * (1 - 2 * (abs(h - delta) / abs(delta))^2) * (abs(h-delta) <= abs(delta)) 
}



#---------
# Wave_v4
#--------
wave_v4 <- function(h, delta, A = 1) {
  A * (1- (1/2)*(abs(h - delta) / abs(delta))^2)*(abs(h - delta) <= 2*abs(delta))
}


#---------
# Wave_v3
#---------
wave_v3 <- function(h, delta, A = 1) {
  A * (1- 0.25*(abs(h - delta) / abs(delta)))*(abs(h - delta) <= 8 * abs(delta))
}



