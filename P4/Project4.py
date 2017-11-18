# -*- coding: utf-8 -*-
import os
import numpy as np
#import multiprocessing as mp
#import matplotlib.pyplot as plt
import time as measure
"""
Run as script for simple 2x2 test.
"""
class Ising_Model_2Dlattice_spin:
    def __init__(self,inscript=False,initial_file='noname'):#dimensions,cycles,temp_initial,temp_final,dT,J=1.):
        """
        Note:
            Module created to use a specialized Metropolis algorithm for the Ising model.
            Give initial size of lattice, temperature, MC-cycles and temperature step.
            The order must be [L,T_i,T_f,Tn,MCc], see args below
        Args:
            Note: Given as type = numpy.ndarray or in text file.    
            inscript (numpy.ndarray,semi-optional): This is the preferred option to give model description,
                the included arguments are: 
                    L: Lattice size (corresponding to LxL)
                    T_i: Initial temperature
                    T_f: Final temperature (if T_i = T_f, only one temperature is iterated over)
                    Tn: Number of temperature points in interval [T_i,T_f]
                    MCc: Monte-Carlo cycles
            initial_file (str): IF, for some reason, the initial values are stored in correct order in a file (see above)
                                then give filename and omit 'inscript' argument.
        Additional Note:
            Arguments can be given either as an array or in a text file, while the latter option
            is not necessary. Multiple values as columns is not supported.
        Attributes:
            All arguments.
        """
        if str(initial_file) is not 'noname':
            self.L,self.T_i,self.T_f,self.Tn,self.MCc = np.genfromtxt(str(initial_file))
            self.L = int(self.L)
            self.MCc = int(self.MCc)
        else:
            if inscript is not False:
                self.L,self.T_i,self.T_f,self.Tn,self.MCc = inscript
                self.L = int(self.L)
                self.MCc = int(self.MCc)
            else:
                print 'There are no parameters given\nAborting.'
                sys.exit(1) 
#-------------------------------------------------------------------------------------------          
    def initial_state(self,random_i=True):
        """
        Creates the initial states and calculates the initial total energy -and magnetic moment values.
        Args:
            random_i (bool): If set to 'False', the initial set up will be all s=1, where 's' is the spin value.
                             default is 'True' which means the initial set up is a random state; s in {-1,1}.
        Attributes:
            Makes the total energy, magnetization and spin matrix class attributes.
            spin_mat: a (Lx2)x(Lx2) numpy.ndarray containing the system configuration.
            E (float): Total energy for inital state.
            M (float): Total magnetization for initial state.
        Returns:
            None-type
        """
        a = np.array([-1.,1.])
        L = self.L
        if random_i is True:
            spin_mat = np.random.choice(a,size=(L+2,L+2))
        else:
            spin_mat = np.ones((L+2,L+2))
        # Setting periodic BC      
        spin_mat[:,0] = spin_mat[:,L]
        spin_mat[:,L+1] = spin_mat[:,1]
        spin_mat[0,:] = spin_mat[L,:]
        spin_mat[L+1,:] = spin_mat[1,:]
        
        # Setting initial energy and magnetic mom
        for i in xrange(1,L+1):
            for j in xrange(1,L+1):
                self.E -= spin_mat[i,j]*(spin_mat[i+1,j]+spin_mat[i,j+1])
                self.M += spin_mat[i,j]
        self.spin_mat = spin_mat
        return None
#------------------------------------------------------------------------------------------- 
    def energy(self): # This is added in ordert to extract a PDF from the results
        """
        Note:
            This is not implemented to be a part of the class, alterations in the code must be done,
            see script: 'supplementary'
        """
        E_st = 0
        for i in xrange(1,int(self.L)+1):
            for j in xrange(1,int(self.L)+1):
                E_st -= self.spin_mat[i,j]*(self.spin_mat[i+1,j]+self.spin_mat[i,j+1])
        self.E_count.append(E_st)        
        return None
#-------------------------------------------------------------------------------------------         
    def flipping(self,i,j):
        """
        Flipping one spin, correct BC's if needed
        """
        L = self.L
        self.spin_mat[i,j] *= -1.
        # Check and update PBC's if i and/or j is the boundary  
        if i == L:
            self.spin_mat[0,j] = self.spin_mat[i,j]
        elif i == 1:
            self.spin_mat[L+1,j] = self.spin_mat[i,j]
        if j == L:
            self.spin_mat[i,0] = self.spin_mat[i,j]
        elif j == 1:
            self.spin_mat[i,L+1] = self.spin_mat[i,j]    
        
        return None
        
#-------------------------------------------------------------------------------------------          
    def Metropolis(self,w):
        """
        Note:
            Calculates a complete MC-cycle optimized to suit the Ising model.
        """
        L = int(self.L)
        for iterable in xrange(L*L):
            i,j = np.random.randint(1,L+1,2)
            dE_t = 2.*self.spin_mat[i,j]*(self.spin_mat[i,j+1] + \
                self.spin_mat[i,j-1] + self.spin_mat[i+1,j] + self.spin_mat[i-1,j])
            if dE_t <= 0:
                self.flipping(i,j)
                self.E += dE_t            
                self.M += 2.*self.spin_mat[i,j]
                self.accepted += 1
            else:
                if np.random.random() <= w[8+int(dE_t)]:
                    self.flipping(i,j)
                    self.E += dE_t
                    self.M += 2.*self.spin_mat[i,j]
                    self.accepted += 1              
        return None      

#-------------------------------------------------------------------------------------------            
            
    def Solve(self,random_i=True,store_values=False,filename='noname',when_save=1):
        """
        Solves the system for given MC-cycles.
        Args:
            random_i (bool,default=True): Set if the initial state should be random or s=1 for all.
            store_values (bool,optional): Set 'True' if values are to be stored then give:
            filename (str,semi-optional): if store_values is set True, give a filename.
            when_save (int,optional): If the above two are given, then this will set the number of
                                      MC-cycles should be between values are written. 
        Attributes:
            The structure follows the example C++ code found in the Computational Physics lecture notes (2015) at UiO.
            The different aspects are assumed known for the intended reader.
        Returns:
            Default is a None-type, but if total energy is wanted, then change this below.
        """
        # Make them local var, define T interval
        MCc = self.MCc; T_i = self.T_i
        T_f = self.T_f; Tn = self.Tn
        if T_i == T_f:
            T_ = [T_i]
        else:            
            T_ = np.arange(T_i,T_f+1E-10,Tn)       
        start = measure.time()
        
        for T in T_: 
            self.E_count = []
            self.accepted=0
            np.random.RandomState()
            print "Running MC with T = %g, lattice: %gX%g\n"%(T,self.L,self.L)
            # Create initial state
            self.E = 0; self.M = 0
            self.initial_state(random_i)
            avg = np.zeros(5)
            beta = 1./T
            w = np.zeros(17) # Larger than 5 in ored to acount for the correct index later
            dE = np.array([-8.,-4.,0,4.,8.])
            for i in range(5):
                w[int(dE[i])+8] = np.exp(-beta*dE[i])
            
            cycle = 0
            # Runs trough the MC without storing values. Done in order to avoid in-loop if test
            if store_values is False:
                while cycle < MCc:
                    self.Metropolis(w)
                    avg[0] += self.E; avg[1] += self.E**2
                    avg[2] += self.M; avg[3] += self.M**2
                    avg[4] += abs(self.M)
                    cycle += 1
                    #self.energy() #  Use to get Gaussian, but check '!!!!' below
                self.average(avg,T,cycle)
            # If values are to be stored
            else:
                if len(T_) == 1: # Check if single temperature is to be used
                    while cycle < MCc:
                        self.Metropolis(w)
                        avg[0] += self.E; avg[1] += self.E**2
                        avg[2] += self.M; avg[3] += self.M**2
                        avg[4] += abs(self.M)
                        cycle += 1
                        if cycle % when_save == 0:
                            self.average(avg,T,cycle,filename)
                        self.accepted = 0
                else: 
                    while cycle < MCc:
                        self.Metropolis(w)
                        avg[0] += self.E; avg[1] += self.E**2
                        avg[2] += self.M; avg[3] += self.M**2
                        avg[4] += abs(self.M)
                        cycle += 1
                    self.average(avg,T,cycle,filename)
                
        end = measure.time()
        print '\nTotal time: %g s'%(end-start)
        return None #self.E_count,self.E_var                       #!!!!!!
#-------------------------------------------------------------------------------------------                 
             
    def average(self,avg,T,cycle,filename='noname'): # Filename to be implemented
        """
        Note:
            Calculates physical properties relating to the system state.
            The structure is again assumed known for the intended reader, but with
            changes in how values are to be stored. Due to the different sets, and 
            memory usage, there are no 'preset' on which values are stored.
            All presets used are shown below, there are three of them, two are 
            always commented out.
        
        Additional comments:
            The statement: 'if hasattr(self,'fil')' is set in order to know if the file is to be created
                            and store values for the first time. It will then assign a filename as a class attribute
                            and the method will then know that values are to be 'appended'.
        Returns:
            None-type
        """
        
        Mc =  float(cycle) # 
        T = float(T)
        L_tot = float(self.L*self.L)
        norm = 1./Mc     
        E_bar, E2_bar, M_bar, M2_bar, M_bar_abs = avg*norm 
        E_var = (E2_bar - E_bar*E_bar)/L_tot # = (<E^2> - <E>^2)/numberspin
        M_abs_var = (M2_bar - M_bar_abs*M_bar_abs)/L_tot # = <M^2> - <|M|>^2 = <M>^2 - <|M|>^2
        Cv = E_var/(T*T)             # Heat capacity per spin
        Chi = M_abs_var/T            # Susceptibility per spin
        M_bar_abs *= 1./L_tot        # <|M|> per spin
        E_bar *= 1./L_tot            # <E> per spin
        
        # Some checks and tweaks if the values are to be stored.
        # This is mostly a result of me being tired of overwriting data...
        if filename is not 'noname': # Check if values are to be stored
            #values = np.array([E_bar,E_var,M_bar_abs,M_abs_var,T,Mc,self.accepted/L_tot])
            #values = np.array([E_bar,E_var,M_bar_abs,M_abs_var,Cv,Chi,T,Mc]) # This is original (initial)
            values = np.array([E_bar,M_bar_abs,Cv,Chi,T]) # For multiple temp run
            if hasattr(self,'fil') is False:    # Check if class has initialized a filename
                if not '.txt' in str(filename): # Check if filename is given in full
                    filename += '.txt'
                if os.path.isfile(filename):    # Check if file exists, respond accordingly
                    print 'Your chosen filename already exists:'
                    answer = str(raw_input("""type 'y' to go ahead anyway, this will overwrite data,
or:
type in new filename (end with .txt) to avoid overwriting  """))
                    if answer is 'y':
                        os.remove(filename)
                        print 'New values will be stored with existing filename.'
                    else:
                        filename = answer
                print 'Data are stored in columns (see script, they change alot).'
                with open(filename,'w') as f_name:
                    np.savetxt(f_name,values.reshape(1,values.shape[0]))
                self.fil = filename # Make the class aware that a filename now exist
            else:
                with open(self.fil,'a') as f_name:
                    np.savetxt(f_name,values.reshape(1,values.shape[0]))
        else:
            print """All values are per spin:
<E> = %g\n<|M|> = %g\nCv = %g\nChi = %g"""%(E_bar,M_bar_abs,Cv,Chi)
        self.E_var = E_var
        self.Mabsbar = M_bar_abs
        self.ebar = E_bar
        return None
    
    def test_run(self,eps=1E-2):
        line = '--------------------------------------'
        print line
        print 'Testing expectation values for energy and\nand abs magnetization for a 2x2-lattice'
        print 'Accept eps=%g'%eps
        print 'MC-cycles = 1E4\n'
        E_bar = -1.99
        M_bar = 0.998
        print 'Solver output:'
        print line
        self.Solve()
        print line
        print 'Test result:'
        if abs(self.ebar - E_bar) < eps and abs(self.Mabsbar - M_bar) < eps:
            print 'Success!'
        else:
            print 'Failed!\nCheck again to be sure, it may be a fluke due\nto test speed>MC-cycles.'
        print line
        return None
            
        
#------------------------------------------------------------------------------------------- 
if __name__ == '__main__':
    
    # test
    L = 2; T_i = 1.0; T_f = 1.; Tn = 1; MCc = 1E4
    ins = np.array([L,T_i,T_f,Tn,MCc])
    test = Ising_Model_2Dlattice_spin(inscript = ins)
    test.test_run()
    del test
    
    
