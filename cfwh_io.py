import numpy as np
import h5py  as h5
import cfwh_lib as cfwhl

PI = np.pi

def read_cfwh_grid(cfwh_c, ip, surf_idx):

    ip_gfname_abs = ip.ip_gpath + ip.ip_gfname

    print(f'CFWH: Reading grid data from file:{ip_gfname_abs}\n')

    fp = h5.File(ip_gfname_abs, 'r')
    
    surf_path = '/surf_' + str(surf_idx) + '/' + 'surf_name'
    surf_name = np.array(fp.get(surf_path))
    print(f'CFWH: Reading surface {surf_name}')    
    
    x_co_path = '/surf_' + str(surf_idx) + '/' + 'x_co'
    x_co = fp.get(x_co_path)        
    x_co = np.array(x_co)

    [N_xi, N_eta] = x_co.shape

    y_co_path = '/surf_' + str(surf_idx) + '/' + 'y_co'
    y_co = fp.get(y_co_path)        
    y_co = np.array(y_co)
        
    z_co_path = '/surf_' + str(surf_idx) + '/' + 'z_co'
    z_co = fp.get(z_co_path)        
    z_co = np.array(z_co)        

    cfwh_c.fsurf.x_co = x_co
    cfwh_c.fsurf.y_co = y_co
    cfwh_c.fsurf.z_co = z_co
    
    cfwh_c.fsurf.N_xi  = N_xi    
    cfwh_c.fsurf.N_eta = N_eta    
       
    fp.close()    

    return cfwh_c
    
def read_cfwh_sources(cfwh_c, ip, surf_idx):

    U1 = cfwh_c.M_1 * cfwh_c.a_inf 

    file_sidx = ip.file_sidx
    file_eidx = ip.file_eidx       

    n1 = cfwh_c.fsurf.n1
    n2 = cfwh_c.fsurf.n2
    n3 = cfwh_c.fsurf.n3    

    [N_xi, N_eta] = n1.shape
    
    N_time = file_eidx - (file_sidx - 1)
    
    cfwh_c.fsurf.Q_src  = np.zeros([N_time, N_xi, N_eta])
    
    cfwh_c.fsurf.F1_src = np.zeros([N_time, N_xi, N_eta])
    cfwh_c.fsurf.F2_src = np.zeros([N_time, N_xi, N_eta])
    cfwh_c.fsurf.F3_src = np.zeros([N_time, N_xi, N_eta])    

    for tidx in range(file_sidx, (file_eidx + 1)):    

        ip_sfname_abs = ip.ip_sfpath + ip.ip_gen_sfname + '_' + str(tidx) + '.h5' 

        print(f'CFWH: Reading solution data from file idx {tidx}')

        fp = h5.File(ip_sfname_abs, 'r')
    
        surf_path = '/surf_' + str(surf_idx) + '/' + 'surf_name'
        surf_name = np.array(fp.get(surf_path))
        print(f'CFWH: Reading surface {surf_name}')    
    
        rho_path = '/surf_' + str(surf_idx) + '/' + 'rho'
        rho = fp.get(rho_path)        
        rho = np.array(rho)

        u_vel_path = '/surf_' + str(surf_idx) + '/' + 'u_vel'
        u_vel = fp.get(u_vel_path)        
        u_vel = np.array(u_vel)        
        
        v_vel_path = '/surf_' + str(surf_idx) + '/' + 'v_vel'
        v_vel = fp.get(v_vel_path)        
        v_vel = np.array(v_vel)                
        
        w_vel_path = '/surf_' + str(surf_idx) + '/' + 'w_vel'
        w_vel = fp.get(w_vel_path)        
        w_vel = np.array(w_vel)        
        
        ps_path = '/surf_' + str(surf_idx) + '/' + 'ps'
        ps = fp.get(ps_path)        
        ps = np.array(ps)        
        
        Q_src = rho * (U1 + u_vel) * n1 + \
                rho * (     v_vel) * n2 + \
                rho * (     w_vel) * n3
 
        F1_src = (rho * u_vel * (U1 + u_vel) + ps) * n1 + \
                 (rho * u_vel * (     v_vel)     ) * n2 + \
                 (rho * u_vel * (     w_vel)     ) * n3

        F2_src = (rho * v_vel * (U1 + u_vel) + ps) * n1 + \
                 (rho * v_vel * (     v_vel)     ) * n2 + \
                 (rho * v_vel * (     w_vel)     ) * n3                 
 
        F3_src = (rho * w_vel * (U1 + u_vel) + ps) * n1 + \
                 (rho * w_vel * (     v_vel)     ) * n2 + \
                 (rho * w_vel * (     w_vel)     ) * n3      
 
        cfwh_c.fsurf.Q_src[tidx, :, :] = Q_src
        
        cfwh_c.fsurf.F1_src[tidx, :, :] = F1_src
        cfwh_c.fsurf.F2_src[tidx, :, :] = F2_src
        cfwh_c.fsurf.F3_src[tidx, :, :] = F3_src        
 
        fp.close()    
  
    return cfwh_c    
