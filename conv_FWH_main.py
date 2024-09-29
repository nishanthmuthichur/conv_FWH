import numpy as np
import conv_FWH_lib as cfwhl
import conv_FWH_io as io

# 1. Write a script to generate an input grid file. 
# 2. Setup the solver functions to read in the input grid files. 
# 3. Write a script to generate input source files. 
# 4. Setup the solver functions to read in the sources. 
# 5. Compute the integrand.
# 6. Compute the far-field pressure. 



# 1. Read the input grid
# 2. Compute the metrics
# 3. Read the sources data
# 4. Compute the integrand
# 5. Compute the far-field pressure

def CFWH_main(cfwh_c):
    
    cfwh_c.M_1   = (cfwh_c.U_1 / cfwh_c.a_inf)
    cfwh_c.beta  = np.sqrt(1 - (cfwh_c.M_1**2))

    cfwh_c.fsurf = io.read_FWH_grid_test(cfwh_c)
    
    cfwh_c.fsurf = cfwhl.comp_surf_metrics(cfwh_c.fsurf)
    
    #fsurf = cfwhl.comp_test_func(fsurf)
    #fsurf = cfwhl.read_cfwh_sources(fsurf)
    
    #fsurf = cfwhl.comp_cfwh_integrand(fsurf)
   
    #fsurf = cfwhl.comp_cfwh_far_field_pres(fsurf)
    
    #print(f'conv_FWH: Numerical surface integral = {fsurf[0].test_ninteg}')
    #print(f'conv_FWH: Analytical surface integral = {fsurf[0].test_ainteg}')    
    
    print('conv_FWH: conv FWH surface successfully created.')