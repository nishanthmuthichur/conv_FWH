import numpy as np
import h5py  as h5
import conv_FWH_lib as cfwhl

PI = np.pi
############################INPUT PARAMETERS###################################

fs = 100
f_src_1 = 10

RAD = 1
LEN = 1

N_rad   = 101
N_theta = 101
N_axi   = 101

N_time  = 10

N_surf  = 3

dt = 1/fs

def cons_surf_flow_field(x_co, y_co, z_co, time):
    
    [N_xi, N_eta] = x_co.shape

    rho   = np.ones([N_xi,  N_eta])
    u_vel = np.ones([N_xi,  N_eta])
    v_vel = np.ones([N_xi,  N_eta])
    w_vel = np.ones([N_xi,  N_eta])
    ps    = np.ones([N_xi,  N_eta])

    return rho, u_vel, v_vel, w_vel, ps

#############################CODE EXECUTION####################################

op_gpath  = 'C:\\Users\\Nishanth\\Desktop\\nish_work\\python_proj\\conv_FWH\\test_data\\'
op_gfname = 'cyl_surf_test_grid.h5'

op_sfpath  = op_gpath
op_gen_sfname = 'cyl_surf_test_src_data'

op_gfname_abs = op_gpath + op_gfname

surf_list = []

rad_vec   = np.linspace(0,    RAD, N_rad)
theta_vec = np.linspace(0, 2 * PI, N_theta)
axi_vec   = np.linspace(0,    LEN, N_axi)

#******************************************************************************
#   Surface 1 (Z = 0 surface)
#***************************************

x_co = np.zeros([N_rad, N_theta])
y_co = np.zeros([N_rad, N_theta])    
z_co = np.zeros([N_rad, N_theta])                       

for rad_idx in range(0, N_rad):
  for theta_idx in range(0, N_theta):
    
      x_co[rad_idx, theta_idx] = rad_vec[rad_idx] * np.cos(theta_vec[theta_idx])          
      y_co[rad_idx, theta_idx] = rad_vec[rad_idx] * np.sin(theta_vec[theta_idx])              
      z_co[rad_idx, theta_idx] = 0                       

surf_list.append(cfwhl.surf())           
surf_list[0].x_co = x_co
surf_list[0].y_co = y_co
surf_list[0].z_co = z_co
surf_list[0].name = 'INLET'
            
#***************************************
#   Surface 2 (Curved surface)   
#***************************************

x_co = np.zeros([N_theta, N_axi])
y_co = np.zeros([N_theta, N_axi])    
z_co = np.zeros([N_theta, N_axi])                       

for theta_idx in range(0, N_theta):
  for axi_idx in range(0, N_axi):
    
      x_co[theta_idx, axi_idx] = RAD * np.cos(theta_vec[theta_idx])          
      y_co[theta_idx, axi_idx] = RAD * np.sin(theta_vec[theta_idx])              
      z_co[theta_idx, axi_idx] = axi_vec[axi_idx]                       

surf_list.append(cfwhl.surf())           
surf_list[1].x_co = x_co
surf_list[1].y_co = y_co
surf_list[1].z_co = z_co
surf_list[1].name = 'CSURF'

#***************************************
#   Surface 3 (Z = LEN surface)
#***************************************

x_co = np.zeros([N_rad, N_theta])
y_co = np.zeros([N_rad, N_theta])    
z_co = np.zeros([N_rad, N_theta])                       

for rad_idx in range(0, N_rad):
  for theta_idx in range(0, N_theta):
    
      x_co[rad_idx, theta_idx] = rad_vec[rad_idx] * np.cos(theta_vec[theta_idx])          
      y_co[rad_idx, theta_idx] = rad_vec[rad_idx] * np.sin(theta_vec[theta_idx])              
      z_co[rad_idx, theta_idx] = LEN

surf_list.append(cfwhl.surf())           
surf_list[2].x_co = x_co
surf_list[2].y_co = y_co
surf_list[2].z_co = z_co
surf_list[2].name = 'OUTLET'

#***************************************
# Write out the grid file in HDF5 format
#***************************************

file_dp = h5.File(op_gfname_abs, 'w')

file_dp.create_dataset('/N_surf', data = N_surf)

for idx in range(0, N_surf):
    
    surf_name_path = f'/surf_{idx}/surf_name'
    surf_x_co_path = f'/surf_{idx}/x_co'
    surf_y_co_path = f'/surf_{idx}/y_co'
    surf_z_co_path = f'/surf_{idx}/z_co'
    
    file_dp.create_dataset(surf_name_path, data = surf_list[idx].name)    
    file_dp.create_dataset(surf_x_co_path, data = surf_list[idx].x_co)
    file_dp.create_dataset(surf_y_co_path, data = surf_list[idx].y_co)
    file_dp.create_dataset(surf_z_co_path, data = surf_list[idx].z_co)

file_dp.close()                                          

print(f'Output grid file is printed.')                     

#***************************************
# Write out the test src file in HDF5 format
#***************************************

time_vec = np.linspace(0, dt * (N_time - 1), N_time)

for tidx in range(0, N_time):

    time = time_vec[tidx]    

    op_sfname_abs = op_sfpath + op_gen_sfname + '_' + str(tidx) + '.h5'
    file_dp = h5.File(op_sfname_abs, 'w')

    file_dp.create_dataset('/N_surf', data = N_surf)

    for idx in range(0, N_surf):

        x_co = surf_list[idx].x_co
        y_co = surf_list[idx].y_co        
        z_co = surf_list[idx].z_co                

        rho, u_vel, v_vel, w_vel, ps = cons_surf_flow_field(x_co, y_co, z_co, time)        

        surf_name_path  = f'/surf_{idx}/surf_name'
        surf_rho_path   = f'/surf_{idx}/rho'
        surf_u_vel_path = f'/surf_{idx}/u_vel'
        surf_v_vel_path = f'/surf_{idx}/v_vel'
        surf_w_vel_path = f'/surf_{idx}/w_vel'    
        surf_ps_path    = f'/surf_{idx}/ps'    
    
        file_dp.create_dataset( surf_name_path, data = surf_list[idx].name)    
        file_dp.create_dataset(  surf_rho_path, data = surf_list[idx].rho)
        file_dp.create_dataset(surf_u_vel_path, data = surf_list[idx].u_vel)
        file_dp.create_dataset(surf_v_vel_path, data = surf_list[idx].v_vel)
        file_dp.create_dataset(surf_w_vel_path, data = surf_list[idx].w_vel)    
        file_dp.create_dataset(   surf_ps_path, data = surf_list[idx].ps)

    file_dp.close()                                          

    print(f'Test src file idx: {tidx} is printed.')                     
