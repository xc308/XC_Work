# XC_Work

## Contents
- 032c: Tst9c, 1D SG_inv construction, Matern, non-crs-MRF, SpN + Reg, thres = 1e-3, reg_num = 1e-9
- 032d: 1D simu plots functions for non-crs-MRF, in paper only C.I. among p <br><br>

- 034b: SG_inv construction, sparse percentage comparison among crs-MRF and non-crs-MRF for Tri-Wave and Wendland  
- 034c: Tst10c, 1D SG_inv construction, crs-MRF, with SpNorm + Reg for b function, b can be chosen <br><br>


- 046b: generate 1D true processes and noisy data, Tri-Wave and Wendland
- 046c: generate 2D true processes and noisy data, Tri-Wave and Wendland;<br><br>
      - discover consistent relationship between sparsity in uni-SG_inv and joint SG_inv <br><br>

- 047b: optimization, Tri-Wave, Tst10c (crs-MRF)
- 047c: optimization, Wendland, Tst10c (crs-MRF) <br><br>

- 048b: cokrig, Tri-Wave, 1 fold C.V. results
- 048d: cokrig, Wendland, 1 fold C.V. results <br><br>

- 049: neg_logL function of non-crs-MRF, TST9d
- 049b: optimization using 049, Tri-Wave, Wendland <br><br>

- 055: 2D inference (neg_logL_2D, optim) for 6 fields in Fig12, Tri-Wave (converged), Wendland (converged) <br><br>

- 056: 2D cokrig (pure denoising) <br><br>

- 057: Data processing, generate df_Res_log_16_sorted, sorted by Lon (asc), then by Lat (desc); 4 Lon strips <br><br>

- 059: TST12 GPU version <br><br>

- 060: GPU parallel + optim on 1 CPU <br><br>

- 061: GPU parallel + optim on 4 CPUs <br><br>

- 062: pure optim parallel on 51 CPUs, no GPU parallelisation <br><br>

- 063: CAMS data processing <br><br>

- 064a: CAMS data with 060, GPU parallel + optim on 1 CPU <br><br>



## Acknowledgements

- Iain Steison recommended using optimParallel() for parallel L-BFGS-B optimization on the CPU. 
- David Llewellyn-Jones helped set up the HPC resource and answered lots of elementary questions regarding Baskerville HPC. 
- Ryan Chan reminded XC that traditional R code will not automatically utilize GPU resources even when run on HPC. 
