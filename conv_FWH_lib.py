import numpy as np

PI = np.pi

DUMMY = -100000

class cfwh_cls():
    
    def __init__(self):
        
        self.a_inf = DUMMY
        self.freq  = DUMMY
        self.omega = DUMMY
        self.k_inf = DUMMY

        self.gamma = 1.4
        self.R_air = 287

        self.X_pos = DUMMY

        self.U_1   = DUMMY
        self.U_2   = DUMMY
        self.U_3   = DUMMY

        self.M_1   = DUMMY
        self.beta  = DUMMY

        self.fsurf = DUMMY




class surf():
    
    def __init__(self):
     
        self.N_xi  = DUMMY
        self.N_eta = DUMMY        
     
        self.x_coord = DUMMY
        self.y_coord = DUMMY
        self.z_coord = DUMMY
   
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

        self.Q = DUMMY
        
        self.F_1 = DUMMY
        self.F_2 = DUMMY
        self.F_3 = DUMMY


        self.test_func  = DUMMY
        self.test_ninteg = DUMMY
        self.test_ainteg = DUMMY


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

    N_surf = len(fsurf)

    for surf_idx in range(0, N_surf):

        [N_xi, N_eta] = fsurf[surf_idx].x_coord.shape

        print(f'N_xi = {N_xi}; N_eta = {N_eta}')    

        metric_xE = np.zeros([N_xi, N_eta])
        metric_yE = np.zeros([N_xi, N_eta])
        metric_zE = np.zeros([N_xi, N_eta])

        metric_xN = np.zeros([N_xi, N_eta])
        metric_yN = np.zeros([N_xi, N_eta])
        metric_zN = np.zeros([N_xi, N_eta])
    
        J         = np.zeros([N_xi, N_eta])

        for eta_idx in range(0, N_eta):
        
            metric_xE[:, eta_idx] = comp_CD8_deriv(fsurf[surf_idx].x_coord[:, eta_idx])
            metric_yE[:, eta_idx] = comp_CD8_deriv(fsurf[surf_idx].y_coord[:, eta_idx])
            metric_zE[:, eta_idx] = comp_CD8_deriv(fsurf[surf_idx].z_coord[:, eta_idx])

        for xi_idx in range(0, N_xi):
        
            metric_xN[xi_idx, :] = comp_CD8_deriv(fsurf[surf_idx].x_coord[xi_idx, :])
            metric_yN[xi_idx, :] = comp_CD8_deriv(fsurf[surf_idx].y_coord[xi_idx, :])
            metric_zN[xi_idx, :] = comp_CD8_deriv(fsurf[surf_idx].z_coord[xi_idx, :])        

        J = np.sqrt( \
                (metric_yE * metric_zN - metric_yN * metric_zE)**2 + \
                (metric_xE * metric_zN - metric_xN * metric_zE)**2 + \
                (metric_xE * metric_yN - metric_xN * metric_yE)**2
            )

        fsurf[surf_idx].metric_xE = metric_xE
        fsurf[surf_idx].metric_yE = metric_yE
        fsurf[surf_idx].metric_zE = metric_zE

        fsurf[surf_idx].metric_xN = metric_xN
        fsurf[surf_idx].metric_yN = metric_yN
        fsurf[surf_idx].metric_zN = metric_zN

        fsurf[surf_idx].J         = J    

    return fsurf






def read_cfwh_sources(fsurf):
    
        
    
    
    return fsurf


def comp_cfwh_integrand(cfwh_c):

    fsurf = cfwh_c.fsurf
    
    a_inf = cfwh_c.a_inf
    k_inf = cfwh_c.k_inf
    U_1   = cfwh_c.U_1
    M_1   = cfwh_c.M_1
    omega = cfwh_c.omega

    beta  = cfwh_c.beta

    X_pos = cfwh_c.X_pos


    
    N_surf = len(fsurf)
    
    for surf_idx in range(0, N_surf):
        
        X_pos = fsurf[surf_idx].X_pos
        
        x_coord = fsurf[surf_idx].x_coord
        y_coord = fsurf[surf_idx].y_coord
        z_coord = fsurf[surf_idx].z_coord        
        
        X = (x_coord - X_pos[0])
        Y = (y_coord - X_pos[1])
        Z = (z_coord - X_pos[2])        
        
        Q = fsurf.Q
        
        F1 = fsurf.F1
        F2 = fsurf.F2
        F3 = fsurf.F3
        
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
                      
        intrnd_Q = (1j * omega * G + DG1) * Q
        
        
                         
               
                   
            
    return fsurf


def comp_cfwh_far_field_pres(fsurf):
    
    return fsurf
 
   
def comp_test_func(fsurf):
    
    x_coord = fsurf[0].x_coord
    y_coord = fsurf[0].y_coord
    z_coord = fsurf[0].z_coord    
    
    theta_coord = np.mod((np.arctan2(y_coord, x_coord) + 2 * PI), 2 * PI)
    
    A = 1
    m = 0.5
    k = 5 * PI / 3
    
    test_func = A * np.exp(1j * (m * theta_coord + k * z_coord))
    
    fsurf[0].test_func = test_func * fsurf[0].J
    fsurf[0].test_ainteg = - (A / (m * k)) * (np.exp(1j * 2 * PI * m) - 1) * (np.exp(1j * k) - 1)
    
    #fsurf[0].test_func   = np.ones([fsurf[0].N_xi, fsurf[0].N_eta]) * fsurf[0].J
    #fsurf[0].test_ainteg = 2 * PI 
    
    return fsurf


def comp_surf_integ(fsurf):
    
    test_func = fsurf[0].test_func
    
    test_ninteg_xi = np.trapz(test_func, axis = 1)

    test_ninteg = np.trapz(test_ninteg_xi, axis = 0)
        
    fsurf[0].test_ninteg = test_ninteg
    
    return fsurf













