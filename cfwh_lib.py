import numpy as np

import scipy as scp
from scipy import signal

PI = np.pi

DUMMY = -100000

class cfwh_cls():
    
    def __init__(self):

        self.fs            = DUMMY
        self.N_fft         = DUMMY
        self.freq_idx_list = DUMMY        
        self.freq_vec      = DUMMY
    
        self.a_inf = DUMMY
        self.k_inf_vec = DUMMY

        self.gamma = 1.4
        self.R_air = 287

        self.X_list = []

        self.M_1   = DUMMY
        self.beta  = DUMMY

        self.ps_Q_list = DUMMY
        self.ps_F_list = DUMMY

        self.fsurf = []

class surf():
    
    def __init__(self):
     
        self.name = DUMMY   
     
        self.N_xi  = DUMMY
        self.N_eta = DUMMY        
     
        self.x_co = DUMMY
        self.y_co = DUMMY
        self.z_co = DUMMY
   
        self.metric_xE = DUMMY
        self.metric_yE = DUMMY
        self.metric_zE = DUMMY

        self.metric_xN = DUMMY
        self.metric_yN = DUMMY
        self.metric_zN = DUMMY

        self.n1 = DUMMY
        self.n2 = DUMMY
        self.n3 = DUMMY

        self.J = DUMMY

        self.rho   = DUMMY
        self.u_vel = DUMMY
        self.v_vel = DUMMY
        self.w_vel = DUMMY
        self.ps    = DUMMY

        self.Q_src = DUMMY
        self.Q_src_fft = DUMMY
        
        self.F1_src = DUMMY
        self.F2_src = DUMMY
        self.F3_src = DUMMY

        self.F1_src_fft = DUMMY
        self.F2_src_fft = DUMMY
        self.F3_src_fft = DUMMY

        self.test_func  = DUMMY
        self.test_ninteg = DUMMY
        self.test_ainteg = DUMMY


def init_solver(ip):

    cfwh_c = cfwh_cls()

    cfwh_c.fs            = ip.fs
    cfwh_c.N_fft         = ip.N_fft
    cfwh_c.freq_idx_list = ip.freq_idx_list
    
    df               = cfwh_c.fs / cfwh_c.N_fft 
    N_freq           = round(cfwh_c.N_fft / 2)
    cfwh_c.freq_vec  = np.linspace(0, (N_freq - 1) * df,  N_freq)
   
    cfwh_c.a_inf     = ip.a_inf

    cfwh_c.X_list = ip.X_list
    N_pts = len(cfwh_c.X_list)

    cfwh_c.M_1  = ip.M_1
    cfwh_c.beta = np.sqrt(1 - (cfwh_c.M_1**2))

    cfwh_c.ps_Q_list = np.zeros([ip.N_surf, ip.N_freq, ip.N_pts])
    cfwh_c.ps_F_list = np.zeros([ip.N_surf, ip.N_freq, ip.N_pts])    

    print('Conv FWH: Initialization has been complete.')

    return cfwh_c
    

def comp_CD8_deriv(Y):

    N_pts = len(Y);
    dYdX = np.zeros(N_pts);
  
    m = 4;
  
    A = np.array([ \
        [   0    ,    0     ,   0    ,    0   , (-25/12) ,  4    ,   -3    , (4/3)   ,  (-1/4)  ], \
        [   0    ,    0     ,   0    , (-1/4) ,  (-5/6)  , (3/2) , (-1/2)  , (1/12)  ,     0    ], \
        [   0    ,    0     , (1/12) , (-2/3) ,     0    , (2/3) , (-1/12) ,   0     ,     0    ], \
        [   0    , (-1/60)  , (3/20) , (-3/4) ,     0    , (3/4) , (-3/20) , (1/60)  ,     0    ], \
        [(1/280) , (-4/105) , (1/5)  , (-4/5) ,     0    , (4/5) , (-1/5)  , (4/105) , (-1/280) ], \
        [   0    , (-1/60)  , (3/20) , (-3/4) ,     0    , (3/4) , (-3/20) , (1/60)  ,     0    ], \
        [   0    ,    0     , (1/12) , (-2/3) ,     0    , (2/3) , (-1/12) ,   0     ,     0    ], \
        [   0    , (-1/12)  , (1/2)  , (-3/2) ,   (5/6)  , (1/4) ,    0    ,   0     ,     0    ], \
        [ (1/4)  , (-4/3)   ,   3    ,    -4  ,  (25/12) ,   0   ,    0    ,   0     ,     0    ], \
    ])
  
    for idx in range(0, N_pts):
        
        if (idx == 0):
                  
            dYdX[idx] =  A[0,    m   ] * Y[idx    ] + \
                         A[0, (m + 1)] * Y[idx + 1] + \
                         A[0, (m + 2)] * Y[idx + 2] + \
                         A[0, (m + 3)] * Y[idx + 3] + \
                         A[0, (m + 4)] * Y[idx + 4]          
            
        elif (idx == 1):
             
            dYdX[idx] =  A[1, (m - 1)] * Y[idx - 1] + \
                         A[1,    m   ] * Y[idx    ] + \
                         A[1, (m + 1)] * Y[idx + 1] + \
                         A[1, (m + 2)] * Y[idx + 2] + \
                         A[1, (m + 3)] * Y[idx + 3]
                        
        elif (idx == 2):

            dYdX[idx] =  A[2, (m - 2)] * Y[idx - 2] + \
                         A[2, (m - 1)] * Y[idx - 1] + \
                         A[2,    m   ] * Y[idx    ] + \
                         A[2, (m + 1)] * Y[idx + 1] + \
                         A[2, (m + 2)] * Y[idx + 2]
                        
        elif (idx == 3):

            dYdX[idx] =  A[3, (m - 3)] * Y[idx - 3] + \
                         A[3, (m - 2)] * Y[idx - 2] + \
                         A[3, (m - 1)] * Y[idx - 1] + \
                         A[3,    m   ] * Y[idx    ] + \
                         A[3, (m + 1)] * Y[idx + 1] + \
                         A[3, (m + 2)] * Y[idx + 2] + \
                         A[3, (m + 3)] * Y[idx + 3]

        elif (idx == (N_pts - 4)):                       

            dYdX[idx] = A[5, (m - 3)] * Y[idx - 3] + \
                        A[5, (m - 2)] * Y[idx - 2] + \
                        A[5, (m - 1)] * Y[idx - 1] + \
                        A[5,    m   ] * Y[idx    ] + \
                        A[5, (m + 1)] * Y[idx + 1] + \
                        A[5, (m + 2)] * Y[idx + 2] + \
                        A[5, (m + 3)] * Y[idx + 3]
                         
        elif (idx == (N_pts - 3)):
            
            dYdX[idx] = A[6, (m - 2)] * Y[idx - 2] + \
                        A[6, (m - 1)] * Y[idx - 1] + \
                        A[6,    m   ] * Y[idx    ] + \
                        A[6, (m + 1)] * Y[idx + 1] + \
                        A[6, (m + 2)] * Y[idx + 2]
        
        elif (idx == (N_pts - 2)):
            
            dYdX[idx] = A[7, (m - 3)] * Y[idx - 3] + \
                        A[7, (m - 2)] * Y[idx - 2] + \
                        A[7, (m - 1)] * Y[idx - 1] + \
                        A[7,    m   ] * Y[idx    ] + \
                        A[7, (m + 1)] * Y[idx + 1] 
                          
        elif (idx == (N_pts - 1)): 
            
            dYdX[idx] = A[8, (m - 4)] * Y[idx - 4] + \
                        A[8, (m - 3)] * Y[idx - 3] + \
                        A[8, (m - 2)] * Y[idx - 2] + \
                        A[8, (m - 1)] * Y[idx - 1] + \
                        A[8,    m   ] * Y[idx    ]
                          
        else:
            
            dYdX[idx] = A[4, (m - 4)] * Y[idx - 4] + \
                        A[4, (m - 3)] * Y[idx - 3] + \
                        A[4, (m - 2)] * Y[idx - 2] + \
                        A[4, (m - 1)] * Y[idx - 1] + \
                        A[4,    m   ] * Y[idx    ] + \
                        A[4, (m + 1)] * Y[idx + 1] + \
                        A[4, (m + 2)] * Y[idx + 2] + \
                        A[4, (m + 3)] * Y[idx + 3] + \
                        A[4, (m + 4)] * Y[idx + 4]            
            
    return dYdX
        

def comp_surf_metrics(fsurf):

#   i         j        k
#   xe       ye       ze
#   xn       yn       zn

    [N_xi, N_eta] = fsurf.x_co.shape

    print(f'N_xi = {N_xi}; N_eta = {N_eta}')    

    metric_xE = np.zeros([N_xi, N_eta])
    metric_yE = np.zeros([N_xi, N_eta])
    metric_zE = np.zeros([N_xi, N_eta])

    metric_xN = np.zeros([N_xi, N_eta])
    metric_yN = np.zeros([N_xi, N_eta])
    metric_zN = np.zeros([N_xi, N_eta])
    
    J         = np.zeros([N_xi, N_eta])

    for eta_idx in range(0, N_eta):
        
        metric_xE[:, eta_idx] = comp_CD8_deriv(fsurf.x_co[:, eta_idx])
        metric_yE[:, eta_idx] = comp_CD8_deriv(fsurf.y_co[:, eta_idx])
        metric_zE[:, eta_idx] = comp_CD8_deriv(fsurf.z_co[:, eta_idx])

    for xi_idx in range(0, N_xi):
        
        metric_xN[xi_idx, :] = comp_CD8_deriv(fsurf.x_co[xi_idx, :])
        metric_yN[xi_idx, :] = comp_CD8_deriv(fsurf.y_co[xi_idx, :])
        metric_zN[xi_idx, :] = comp_CD8_deriv(fsurf.z_co[xi_idx, :])        

    vec_1 =  (metric_yE * metric_zN - metric_yN * metric_zE)
    vec_2 = -(metric_xE * metric_zN - metric_xN * metric_zE)
    vec_3 =  (metric_xE * metric_yN - metric_xN * metric_yE)
    
    vec_mag = np.sqrt(vec_1**2 + vec_2**2 + vec_3**2)
    
    n1 = vec_1 / vec_mag
    n2 = vec_2 / vec_mag
    n3 = vec_3 / vec_mag    

    J = vec_mag

    fsurf.metric_xE = metric_xE
    fsurf.metric_yE = metric_yE
    fsurf.metric_zE = metric_zE

    fsurf.metric_xN = metric_xN
    fsurf.metric_yN = metric_yN
    fsurf.metric_zN = metric_zN

    fsurf.J         = J    

    fsurf.n1 = n1
    fsurf.n2 = n2 
    fsurf.n3 = n3        

    return fsurf

def comp_fsurf_DFT(cfwh_c):
   
    cfwh_c.fsurf.Q_src_fft  = comp_surf_DFT(cfwh_c.fsurf.Q_src)
    
    cfwh_c.fsurf.F1_src_fft = comp_surf_DFT(cfwh_c.fsurf.F1_src)    
    cfwh_c.fsurf.F2_src_fft = comp_surf_DFT(cfwh_c.fsurf.F2_src)    
    cfwh_c.fsurf.F3_src_fft = comp_surf_DFT(cfwh_c.fsurf.F3_src)    
   
    return cfwh_c

def comp_surf_DFT(g_src):
    
    [N_fft, N_xi, N_eta] = g_src.shape
    
    g_src = g_src - np.mean(g_src, 0)
    
    win = signal.windows.hann(N_fft)
    
    for tidx in range(0, N_fft):
        
        g_src[tidx, :, :] = g_src[tidx, :, :] * win[tidx]
    
    g_src_fft = scp.fft.fft(g_src, \
                            N_fft, \
                                0, \
                 norm = "forward") 
   
    return g_src_fft

def comp_cfwh_integrand(cfwh_c, freq_idx, pt_idx):

    fsurf = cfwh_c.fsurf        
    
    X_pos = cfwh_c.X_list[pt_idx]    

    omega = 2 * PI * cfwh_c.freq_vec[freq_idx]
    
    Q_src_fft  = fsurf.Q_src_fft[freq_idx, :, :]
        
    F1_src_fft = fsurf.F1_src_fft[freq_idx, :, :]
    F2_src_fft = fsurf.F2_src_fft[freq_idx, :, :]
    F3_src_fft = fsurf.F3_src_fft[freq_idx, :, :]

    a_inf = cfwh_c.a_inf
    k_inf = omega / a_inf
    M_1   = cfwh_c.M_1

    beta  = cfwh_c.beta

    x_co = fsurf.x_co
    y_co = fsurf.y_co
    z_co = fsurf.z_co
        
    X = (x_co - X_pos[0])
    Y = (y_co - X_pos[1])
    Z = (z_co - X_pos[2])        
        
    R_star = np.sqrt(X**2 + beta * (Y**2 + Z**2))

    R = (-M_1 * X + R_star) / (beta * beta)
            
    G = np.exp(-1j * k_inf * R) / (4 * PI *R)
        
    DG1 = (np.exp(-1j * k_inf * R) / (4 * PI * R_star))  *       \
              (                                                  \
                  (1j * k_inf * R_star * (M_1 * R_star - X)) +   \
                  (beta**2 * X)                                  \
              )    / (beta**4 * R_star**2) 
                      
    DG2 = (np.exp(-1j * k_inf * R) / (4 * PI * R_star))  * \
              (                                            \
                  (1j * k_inf * R_star * beta * Y)       + \
                  (beta**2 * Y)                            \
              )    / (beta**2 * R_star**2)                  
                      
    DG3 = (np.exp(-1j * k_inf * R) / (4 * PI * R_star))  * \
              (                                            \
                  (1j * k_inf * R_star * beta * Z)       + \
                  (beta**2 * Z)                            \
              )    / (beta**2 * R_star**2)                                        
                      
    fsurf.integ_Q = (1j * omega * G + DG1) * Q_src_fft

    fsurf.integ_F = F1_src_fft * DG1 + \
                    F2_src_fft * DG2 + \
                    F3_src_fft * DG3
   
    cfwh_c.fsurf = fsurf   
   
    return cfwh_c


def comp_cfwh_far_field_ps(cfwh_c, freq_idx, pt_idx):
    
    integ = cfwh_c.fsurf.integ_Q * cfwh_c.fsurf.J
    ps_Q = comp_surf_integ(integ)
    
    integ = cfwh_c.fsurf.integ_F * cfwh_c.fsurf.J
    ps_F = comp_surf_integ(integ)

    return ps_Q, ps_F


def comp_surf_integ(integ):
   
    integ_xi = np.trapz(integ, axis = 1)

    integral = np.trapz(integ_xi, axis = 0)
        
    return integral


def comp_test_func(fsurf):
    
    J = fsurf.J
    
    integral = comp_surf_integ(J)
    
    print(integral)
    
    return 










