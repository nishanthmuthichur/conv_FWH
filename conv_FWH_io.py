import numpy as np
import h5py  as h5
import conv_FWH_lib as cfwhl

PI = np.pi


def read_FWH_grid_test(ip):

    fsurf = list()    

    ip_grid_fname_abs = ip.ip_grid_path + \
                        ip.ip_grid_fname

    print(f'ponni_solver: Reading grid data from file:{ip_grid_fname_abs}')

    file_dp = h5.File(ip_grid_fname_abs, 'r')

    N_surf = file_dp.get('/N_surf')
    N_surf = np.array(N_surf)

    for surf_idx in range(0, N_surf):
        
        fsurf.append(cfwhl.surf())        
        
        x_coord_path = '/surf_' + str(surf_idx + 1) + '/' + 'x_coord'
        x_coord = file_dp.get(x_coord_path)        
        x_coord = np.array(x_coord)

        [N_xi, N_eta] = x_coord.shape

        y_coord_path = '/surf_' + str(surf_idx + 1) + '/' + 'y_coord'
        y_coord = file_dp.get(y_coord_path)        
        y_coord = np.array(y_coord)
        
        z_coord_path = '/surf_' + str(surf_idx + 1) + '/' + 'z_coord'
        z_coord = file_dp.get(z_coord_path)        
        z_coord = np.array(z_coord)        

        fsurf[surf_idx].x_coord = x_coord
        fsurf[surf_idx].y_coord = y_coord
        fsurf[surf_idx].z_coord = z_coord        

        fsurf[surf_idx].N_xi  = N_xi
        fsurf[surf_idx].N_eta = N_eta

    file_dp.close()    

    return fsurf
    
    
    
    

