import numpy as np
import cfwh_lib as cfwhl
import cfwh_io as io

PI = np.pi

def cfwh_main(ip):
    
    cfwh_c = cfwhl.init_solver(ip)

    N_surf = ip.N_surf
    N_freq = ip.N_freq
    N_pts  = ip.N_pts

    for surf_idx in range(0, N_surf):

        cfwh_c.fsurf = cfwhl.surf()
        
        cfwh_c = io.read_cfwh_grid(cfwh_c, ip, surf_idx)
       
        cfwh_c.fsurf = cfwhl.comp_surf_metrics(cfwh_c.fsurf)
     
        cfwh_c = io.read_cfwh_sources(cfwh_c, ip, surf_idx)
    
        cfwh_c = cfwhl.comp_fsurf_DFT(cfwh_c)
    
        for freq_idx in range(1, N_freq):
            
            print(f'freq_idx = {freq_idx}')

            for pt_idx  in range(0, N_pts):
    
                cfwh_c = cfwhl.comp_cfwh_integrand(cfwh_c, freq_idx, pt_idx)
   
                ps_Q, ps_F = cfwhl.comp_cfwh_far_field_ps(cfwh_c, freq_idx, pt_idx)

                cfwh_c.ps_Q_list[surf_idx, freq_idx, pt_idx] = ps_Q
                cfwh_c.ps_F_list[surf_idx, freq_idx, pt_idx] = ps_F

        #cfwhl.comp_test_func(cfwh_c.fsurf)        