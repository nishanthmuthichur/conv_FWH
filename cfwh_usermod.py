import numpy as np

DUMMY = -100000

class ip_params():
    
    def __init__(self):
        
        self.ip_gpath = DUMMY
        self.ip_gname = DUMMY
        
        self.ip_sfpath = DUMMY        
        self.ip_gen_sfname = DUMMY
 
        self.file_sidx = DUMMY
        self.file_eidx = DUMMY
       
        self.a_inf     = DUMMY
        
        self.fs     = DUMMY
        self.N_fft  = DUMMY        
        self.N_freq = DUMMY
        self.freq_idx_list = DUMMY
        
        self.X_list = DUMMY

        self.M_1   = DUMMY

        self.N_surf = DUMMY