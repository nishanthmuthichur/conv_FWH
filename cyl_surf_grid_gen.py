import numpy as np
import h5py  as h5

PI = np.pi
############################INPUT PARAMETERS###################################

RAD = 1
LEN = 1

N_rad   = 101
N_theta = 101
N_axi   = 201

N_surf  = 3

#############################START OF EXECUTION################################

wkdir_grid_write  = 'C:\\Users\\Nishanth\\Desktop\\nish_work\\python_proj\\conv_FWH\\'
op_grid_fname     = 'cyl_surf_grid_test.h5'

op_grid_fname_abs = wkdir_grid_write + \
                       op_grid_fname

#   Surface 1
 
x1_coord = np.zeros([N_rad, N_theta])
y1_coord = np.zeros([N_rad, N_theta])    
z1_coord = np.zeros([N_rad, N_theta])                       

rad_vec   = np.linspace(0,    RAD, N_rad)
theta_vec = np.linspace(0, 2 * PI, N_theta)

for rad_idx in range(0, N_rad):
  for theta_idx in range(0, N_theta):
    
      x1_coord[rad_idx, theta_idx] = RAD * np.cos(theta_vec[theta_idx])          
      y1_coord[rad_idx, theta_idx] = RAD * np.sin(theta_vec[theta_idx])              
      z1_coord[rad_idx, theta_idx] = 0                       
                       
#   Surface 2    
 
x2_coord = np.zeros([N_theta, N_axi])
y2_coord = np.zeros([N_theta, N_axi])    
z2_coord = np.zeros([N_theta, N_axi])                       

theta_vec = np.linspace(0, 2 * PI, N_theta)
axi_vec   = np.linspace(0,      1,   N_axi)

for theta_idx in range(0, N_theta):
  for axi_idx in range(0, N_axi):
    
      x2_coord[theta_idx, axi_idx] = RAD * np.cos(theta_vec[theta_idx])          
      y2_coord[theta_idx, axi_idx] = RAD * np.sin(theta_vec[theta_idx])              
      z2_coord[theta_idx, axi_idx] = axi_vec[axi_idx]                       

#   Surface 3
 
x3_coord = np.zeros([N_rad, N_theta])
y3_coord = np.zeros([N_rad, N_theta])    
z3_coord = np.zeros([N_rad, N_theta])                       

rad_vec   = np.linspace(0,    RAD, N_rad)
theta_vec = np.linspace(0, 2 * PI, N_theta)

for rad_idx in range(0, N_rad):
  for theta_idx in range(0, N_theta):
    
      x3_coord[rad_idx, theta_idx] = RAD * np.cos(theta_vec[theta_idx])          
      y3_coord[rad_idx, theta_idx] = RAD * np.sin(theta_vec[theta_idx])              
      z3_coord[rad_idx, theta_idx] = LEN


file_dp = h5.File(op_grid_fname_abs, 'w')

file_dp.create_dataset('/N_surf'           , data = N_surf)

file_dp.create_dataset('/surf_1/surf_name' , data = 'INLET')    
file_dp.create_dataset('/surf_1/x_coord'   , data = x1_coord)
file_dp.create_dataset('/surf_1/y_coord'   , data = y1_coord)
file_dp.create_dataset('/surf_1/z_coord'   , data = z1_coord)

file_dp.create_dataset('/surf_2/surf_name' , data = 'CSURF')    
file_dp.create_dataset('/surf_2/x_coord'   , data = x2_coord)
file_dp.create_dataset('/surf_2/y_coord'   , data = y2_coord)
file_dp.create_dataset('/surf_2/z_coord'   , data = z2_coord)

file_dp.create_dataset('/surf_3/surf_name' , data = 'OUTLET')    
file_dp.create_dataset('/surf_3/x_coord'   , data = x3_coord)
file_dp.create_dataset('/surf_3/y_coord'   , data = y3_coord)
file_dp.create_dataset('/surf_3/z_coord'   , data = z3_coord)

file_dp.close()                                          

print(f'Output grid file: {op_grid_fname_abs}')                     
print(f'Output grid file is successfully generated.')






