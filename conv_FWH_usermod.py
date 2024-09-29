import numpy as np

DUMMY = -100000

class ip_params():
    
    def __init__(self):
        
        self.ip_gpath = DUMMY
        self.ip_gname = DUMMY
        
        self.ip_fpath = DUMMY
        self.ip_fname = DUMMY
        
        self.a_inf = DUMMY
        self.freq  = DUMMY
        self.omega = DUMMY
        self.k_inf = DUMMY
        
        self.X_pos = DUMMY

        self.U_1   = DUMMY
        self.U_2   = DUMMY
        self.U_3   = DUMMY        
