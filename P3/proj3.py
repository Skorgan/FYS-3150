import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from Solvers import Solvers as S
"""
Note:
This module is more specialized than generelized to project 3 in Computational physics at UIO.
Some formalities in structure (and unwritten 'rules' (norms)) are therefore disregarded, but one could extend the main body to suit more general scenarios.
I would not recommend using this in any other way under the restrictions written in the doc strings. Importing this class will however
import multiple modules (see above)

Also, many methods require the module 'Solvers' which is in a seperate py script in order to shorten this script. The module
'Solvers' was created especially for this project and contains integrator methods and 2D, 3D cross product in order to speed up certain
numerical parts (for instance; avoid Numpy's 'cross' method), but is general.
"""
class SolarSystem: 
#------------------------------------------------------------------------------------------------------------------    
    def __init__(self,dimensions,num_bodies,pos,V0,mass,beta=2): # pos is given as nested array [number_of_bodies,dimension]
        """
        Initializes the solar system in question by creating the instance, currently supports only the Sun as a star, 
        and weight in solar mass and distances in astronomical units: 1AU = 1.5E11 meters for further usage.
         
        Args:
            dimensions (int): Currently supports 2 or 3 dimensions, sets dimensionality of system.
            num_bodies (int): number of planets + 1, theoretical unlimited, practically, well...
            pos (numpy.ndarray): [[vec_sun],[vec_planet1],...], the positions of the sun and planets (pos[0] = sun position x,y (,z))
                                 dimensions should relate as: pos.size = dimension*num_planets, pos[i].size = dimension.
            V0 (numpy.ndarray): same set-up as 'pos'
            mass (numpy.ndarray): array containing the masses IN solar units; array([M_sun,M_planet1,M_planet2,...]), (M_sun = 1)
            beta (int): This does not need to be here, but it is a fudge factor in the force calculation (seen througout this script),
                        for beta=2, the physics are correct (since this is special to an assignment, I assume this is known throughout)
        Attributes:
            G: Newtons gravitational constant in units [(AU)^3(yr)^-2/M_sun].
            All from args, but shorter names; num_bodies=n_b, dimensions=dim..
            
        Example:
            2 celestial bodies (including the Sun), 3-dimensions, where Sun position = [0,0,0] and planet position = [1,0,0], then
            pos = np.array([[0,0,0],[1,0,0]]) (same type for V0), M_planet = M_planet_in_kg/1.99E30kg such that 
            mass = np.array([1.,M_planet]). The instance is created by: arg = SolarSystem(3,2,pos,V0,mass), but as is, this serves
            no function. 
        """
        self.lines = '-------------------------------------------------------------------------------'
        n_b = int(num_bodies); dim = int(dimensions) 
        if hasattr(mass,'__len__'):
            if len(mass) == n_b:    
                self.M = mass
            else:
                print 'ERROR:\n',self.lines,'Number of celectial bodies do not correspond to number of masses given.\nAborted program!\n',self.lines
                sys.exit(1)
        else:
            print 'ERROR:\n',self.lines,'Data: mass, does not seem to be an array.. Check and run again\nAborted program!\n',self.lines
            sys.exit(2)
        if len(pos) == n_b and len(V0) == n_b and len(pos[0]) == dim and len(V0[0]) == dim:
            self.R0 = pos 
            self.V0 = V0
        else:
            print 'ERROR:\n',self.lines,'You have not given correct dimensions for positions and/or velocities of celectial bodies.\nPlease check that you have given correct number dimensions in the given values\nAborted program!\n',self.lines
            sys.exit(3)
        self.G = 4.*np.pi*np.pi      #(AU)^3(yr)^-2/M_sun
        self.n_b = n_b
        self.dim = dim
        self.M_sun = 1. 
        self.beta = beta
        print self.lines
        print 'System with %g celestial bodies in %g dimensions is initialized.'%(self.n_b,self.dim)
        print self.lines
#------------------------------------------------------------------------------------------------------------------        
class Planets(SolarSystem): # Properties to the celestial bodies  
    """
    Contains different properties relating to celestial bodies, parent class: SolarSystem, inherits __init__.
    Current methods:
        distance, Force, acceleration_single, acceleration_all, acceleration_stat_sun_N, Kinetic_E, Potential_E, Momentum, AngMomentum
    """
    def distance(self,R): # R is distance vector
        """
        Note:
            REQUIRED dimensionality of given vector, see below. Supports N-dimensions.      
        Args:
            R (numpy.ndarray): dimensionality; R must be as: np.array([vec]) where 'vec' = x,y(,z) such that R.size = dimension.
        Returns:
            Lengths of given vector
        """        
        dim = self.dim
        temp = 0
        for i in range(dim):
            temp += R[i]**2
        return np.sqrt(temp)   
#-----------------------------------------------------------------------------------------------------------------------------------        
    def Force(self,R_rel,M1,M2,beta=2): # R_rel is positional vector from the body with mass M1 to body of mass M2
        """
        Caution:
            This method calls 'distance' in class 'Planets'. In the case there is an error in 'distance', it is inherited here.
            However, this was done in order to shorten the overall coding
        Note:
            Calculates gravitational force vector; Newtonian
        Args:
            R_rel (numpy.ndarray): relative positional/coordinate vector from body (1) to body (2), units: [AU] 
            M1 and M2 (float): Mass of body (1) (M1), mass of body (2) (M2) units: solar_mass
        
        Returns:
            Force vector in correct direction given arg: 'R_rel', is 'relatively' correct.      
        """
        r = self.distance(R_rel)
        return self.G*M1*M2/(r**(beta+1))*R_rel  
#-----------------------------------------------------------------------------------------------------------------------------------       
    def acceleration_single(self,R1,M2,beta=2): # acceleration of body:1 caused by body:2
        """
        Caution:
            Calls method 'distance' in class 'Planets'. Created for stationary Sun.
        Note:
            Calculates acceleration of body (1) caused by body (2). One could use method 'Force' and derive the acceleration,
            but this has special use that shortens future code lines for the creator.
        Args:
            R1 (numpy.ndarray): positional/coordinate vector, R.size = dimension = len(R).
            M2 (float): mass of body (2) at origo (for all intent and purpose; the Sun).
        Returns:
            Acceleration (vector) for body (1).
        """
        G = self.G
        beta = self.beta
        R2 = np.zeros(self.dim)
        R_rel = R2 - R1 # This is not really necessary for its use, could return the negative value instead (with R_rel = R1)
        r = self.distance(R_rel)
        return G*M2/(r**(beta+1))*R_rel
#-----------------------------------------------------------------------------------------------------------------------------------            
    def acceleration_all(self,R,M,beta=2):  
        """
        Caution:
            Calls method 'distance' in class 'Planets'.
        Note:
            Calculates acceleration vectors for all bodies (Newtonian).
        Args:
            R (numpy.ndarray): nested positional vector array, dimensionality: R.size = number_of_bodies*dimension=len(R)*len(R[i]), 
                R[:,i].size = number_of_bodies.
            M (numpy.ndarray) (or list): Contains all body masses, M.size = number of bodies
        Returns:
            Acceleration, same type (and similar order) as arg: R. (a[0] = acceleration of body at R[0])
        """
        G = self.G
        n_b = self.n_b; dim = self.dim
        a = np.zeros((n_b,dim))
        for i in xrange(n_b):
            for j in xrange(i+1,n_b):
                R_rel = R[j] - R[i] # Relative vector from current body [i] to body [j]
                r = self.distance(R_rel)
                a[i] += G*M[j]/(r**(beta+1))*R_rel
                a[j] -= G*M[i]/(r**(beta+1))*R_rel  # Opposite acceleration, switched mass
        return a
#-----------------------------------------------------------------------------------------------------------------------------------        
    def acceleration_stat_sun_N(self,R,M,beta=2):
        """
        Caution:
            Calls method 'distanse' from class 'Planets'.
        Note:
            Calculates the acceleration of N-planets given a stationary Sun.
        Args:
            R (numpy.ndarray): nested positional/coordinate vector array, R[0] = vector_pos_sun. Units: AU
            M (numpy.ndarray): Contains the masses, M[0] = 1 = sun_solar_mass, index (i) corresponds to body (i) at position R[i].
        Returns: 
            acceleration for all bodies - the sun, a[0] = origo.
        """
        G = self.G
        n_b = self.n_b; dim = self.dim
        a = np.zeros((n_b,dim))
        for i in xrange(1,n_b):
            r = self.distance(R[i])
            a[i] = -G*M[0]/(r**(beta+1))*R[i]
            for j in xrange(i+1,n_b):
                R_rel = R[j] - R[i]
                r = self.distance(R_rel)
                temp = G/(r**(beta+1))
                a[i] += M[j]*temp*R_rel
                a[j] -= M[i]*temp*R_rel # opposite acceleration      
        return a
#-----------------------------------------------------------------------------------------------------------------------------------                      
    def Kinetic_E(self,V,M): 
        """
        Note:
            Calculates the kinetic energy. It uses numpys module 'dot' (for reasons.., laziness?), 
            which may be slow for large arrays. No restriction in units.
        Args:
            V (numpy.ndarray): nested velocity array, dimensionality: V.size = number_of_velocity_vectors*dimension, V[i].size = dimension,
            M (float, int): objects mass
        Returns:
            K: Objects kinetic energy as type: (numpy.ndarray).
        """
        n = len(V)
        K = .5*M*np.array([np.dot(V[i],V[i]) for i in xrange(n)])
        return K
#-----------------------------------------------------------------------------------------------------------------------------------           
    def Potential_E(self,R,M,number_of_points):
        """
        Caution:
            Calls method 'distance' in class 'Planets'. Restricted for planet-sun system,
              (not generalized yet, but few alterations required to get that done).
        Note:
            Calculates the gravitational potential energy for a single body/planet relative to the Sun.
        Args:
            R (numpy.ndarray): nested positional/coordinate array (optionally: single positional array).
            M (float): Mass of planet in solar masses.
            number_of_points (int): number of positions/points in array: R.
        Returns:
            T (numpy.ndarray) (or float,int): potential energy of planet-sun system.   
        """
        r = np.array([self.distance(R[i]) for i in xrange(int(number_of_points))]) 
        T = - self.G*M*self.M_sun/r
        return T
#-----------------------------------------------------------------------------------------------------------------------------------           
    def Momentum(self,M,V):
        """
        Note:
            Calculates translational momentum. Not unit binding.
        Args:
            M (float,((int)): Objects mass.
            V (numpy.ndarray) or (float): Velocity vector (if type(V) is 'int', type(M) must be a float and vice versa). 
        Returns:
            Translational momentum, either as an array, or type: float, depended on type(V).
        """
        return M*V
#-----------------------------------------------------------------------------------------------------------------------------------            
    def AngMomentum(self,R,p,abs_val = False): # p = momentum, R rel to center or center of mass, R = n_points_calculatedXdim = p
        """
        Note:
            Calculates the orbital angular momentum.
        Caution/remark:
            Even though the parallell module 'Solvers' has 2D and 3D vector cross product solver, Numpy's method 'cross' is used here.
            Reason is chronological in nature (this was created first).
        Args:
            R (numpy.ndarray): 2D or 3D positional vector array.
            p (numpy.ndarray): 2D or 3D translational momentum vector array.
            abs_val (bool): if 'True'; calculates abs_value of the orbital angular momentum (either as array or float),
                 if 'False' does not calculate abs val.
        Returns:
            Either abs_val, value, abs_val array or value array for orbital angular momentum. 
        """
        ang_mom = np.cross(R,p)
        if abs_val is False:
            return ang_mom
        else:
            if hasattr(ang_mom,'__len__'): 
                abs_ang = np.array([np.linalg.norm(ang_mom[i]) for i in xrange(len(ang_mom))])
                return abs_ang
            else:
                return ang_mom
#-----------------------------------------SOME TESTS-------------------------------------------------------------------------   
    def test_distance(self):
        line = self.lines
        print line
        print "Testing: 'distance' in 'Planets':"
        correct = 1.
        R_test = self.R0[1]
        num = self.distance(R_test)
        if (correct==num) == False:
            print 'ERROR:\nSomething is altered, test failed....'
            print line
            return False
        else:
            print 'Success!'
            print line
            return True

    def test_Force(self):
        line = self.lines
        print line
        print "Testing: 'Force' in 'Planets':" 
        correct = np.array([-2*np.pi**2,0])
        R_rel = -self.R0[1]
        num = self.Force(R_rel,self.M[1],self.M[0])
        if False in (correct==num):
            print 'ERROR:\nSomething is altered, test failed....'
            print line
            return False
        else:
            print 'Success!'
            print line
            return True       
    
    def test_Kinetic_E(self):
        line = self.lines
        print line
        print "Testing: 'Kinetic_E' in 'Planets':"
        M = self.M[1]
        V = np.array([self.V0[1]])
        correct = np.array([.25])
        num = self.Kinetic_E(V,M)
        if False in (correct==num):
            print 'ERROR:\nSomething is altered, test failed....'
            print line
            return False
        else:
            print 'Success!'
            print line
            return True
            
#------------------------------------------------------------------------------------------------------------------    
class Orbit(Planets,SolarSystem):
    """
    A more specialized sub class, inherits __init__ from class 'SolarSys', can use all methods in class 'Planets'.
    Current methods:
        Escape_V, Orbit_2_body_Earth_Sun, Orbit_3_body, acceleration_Mercury, Mercury_precession.
    """
#-----------------------------------------------------------------------------------------------------------------------------------        
    def Escape_V(self,time,dt,V0x,eps):
        """
        Caution:
            Calls method 'acceleration_single' from class 'Planets', imports module 'Solvers' and use Numpy's module 'linalg' 
            method 'norm' (like 'distance' in class 'Planets', but since it is called few times, it does not matter) 
            and will not work without. Currently supports 2D (no reason for 3D as any line in 3D can be projected as a 2D plane).
        Note: 
            Aims to calculate the escape velocity using conservation of energy as requirement. Planet escapes the gravitational
            well of the Sun if the kinetic energy is non-zero for a potential energy value = 0 (with margin eps numerically).
            The user must give the time and time step to use as well as single or multiple initial velocity candidate(s)
            to test.
        Args:
            time (float): Time to simulate in years
            dt (float): Time step to use, in years. (time/dt) then corresponds to maximum number of iterations.
            V0x (numpy.ndarray): Containing initial velocities to test, Unit: AU/year, NB: 
                                 make sure that the values increase from L, to R.
            eps (float): Error margin corresponding to a numerical zero value.
        Returns:
            0, if no values were found, or the first passing escape velocity.
        """
        from Solvers import Solvers as S
        # Initial set up, pos,-vec,-time ets -arrays
        N = int(round(time/float(dt)))   
        R = np.zeros((N+1,2))           
        V = np.zeros((N+1,2))           
        M = self.M                      
        a = np.zeros((N+1,2))
        t = np.zeros(N+1)
        R[0] = self.R0[1]
        a[0] = self.acceleration_single(R[0],M[0]) 
        value = np.zeros(len(V0x))    # value = 0 would save memory, an array is not needed, but good to have if changes are made here
        o_v = 0 # iterative index value to index array 'value' if a values is found, 
        if self.dim != 2:
            print 'Function: Escape_V, works only in 2D\nAborting'
            sys.exit()
        for Vx in V0x:
            V[0] = np.array([Vx,0]) # New initial test val to be iterated
            for i in xrange(N):
                R[i+1],V[i+1],a[i+1],t[i+1] = S().Velocity_Verlet(R[i],V[i],a[i],t[i],dt,self.M_sun,self.acceleration_single)
                
                if V[i,0]-V[i+1,0] > V[i,0]+V[i+1,0]: #Check if V[i+1,0] is neg, breaks for loop for current Vx val, saves time,num error
                    break
                if np.linalg.norm(V[i+1]) > eps: # Condition that the velocity is not zero
                    y = np.array([R[i+1]])       # y is a bad name, this is basically the position, but 'potential_E' takes special type
                    T_trial = self.Potential_E(y,M[1],1) # Potential energy
                    if abs(T_trial) < eps: # Condition if potential energy is zero within chosen numerical zero val: eps
                        value[o_v] = Vx
                if value[o_v] != 0:
                    o_v +=1
                    break
            if True in (value != 0): # Here comes the 'value=0 save memory part', by commenting out this bool, iteration will continue
                break
        print 'Escape vel: Vx0 = %g AU/yr'%value[0]
        return value[0]          
#-----------------------------------------------------------------------------------------------------------------------------------       
    def Orbit_2_body_Earth_Sun(self,time,dt,method='V',save_dat=False,filename='noname'): 
        """
        Caution:
            Calls method 'acceleration_singe' from class 'Planets', imports module 'Solvers', will not work without.
            Highly specialised method, however, even though the name implies 'Earth', any other planet, or made up body can be used.
            The sun is stationary (at origo) in this method.
        Note:
            Solves a 'Sun-planet' system with the option of two iterative methods (see module 'Solvers'; methods: 'Euler_Chromer' and
            'Velocity_Verlet'. If a filename is given (see Args below), then the postions, velocities and time values are also written
            to a file (stored as columns [R V t] ex.  x[i], y[i], vx[i], vy[i], t[i]). (see Args below for dimensionalities etc.). 
            Default method is the 'Velocity_Verlet' method in 'Solvers'. All methods in 'Solvers' has a boolian test function, 
            a test of the chosen function before proceeding with calculations is more safe, but this is not implemented, always run 
            'Solvers.py' to see if it works as intended.
        Args:
            time (float (or int)): Time to calculate in years.
            dt (float): Time differential in years.
            method ('str',optional): allowed; 'V','E' (no failsafe), 'V' uses 'Velocity_Verlet'- while 'E' uses 'Euler_Chromer' 
                   -method in module 'Solvers'.
            save_dat (bool): True; saves data/values [position,velocity,time] to given filename, if filename is not provided as
                an argument, the user is prompt to provide one in the terminal.
            filename ('str',optional): IF 'save_dat' is set to 'True', data is stored, if 'False', no data is written (even though
                user has provided a filename).
        Returns:
            Position, -velocity and time array data; R,V,t. Both; R and V are for the planet alone, R.size=tot_points*dim,
            R[:,i].size = tot_points, and similarly for V.
                      
        """
        from Solvers import Solvers as S # Import the solver class; see 'Solvers.py'
        # Initial set up
        dim = self.dim; n_b = self.n_b
        M = self.M
        N = int(round(time/float(dt)))
        R = np.zeros((N+1,dim))
        V = np.zeros((N+1,dim))
        t = np.zeros(N+1)       
        R[0] = self.R0[1]
        V[0] = self.V0[1]
        # Which method:
        if method == 'E': # Euler Chromer
            # Calculate...
            for i in xrange(N):
                a = self.acceleration_single(R[i],self.M_sun)
                R[i+1],V[i+1],t[i+1] = S().Euler_Chromer(R[i],V[i],a,t[i],dt)
            # Save data...
            if save_dat is True:
                if str(filename) != 'noname':
                    if not '.txt' in str(filename):
                        filename += '.txt'
                    values = np.array([R,V,t])
                    np.savetxt(str(filename),np.c_[R,V,t])
                    print 'Saved data in file: %s\n'%(str(filename))
                else:
                    filename = str(raw_input("Please provide a filename that ends with '.txt' (less coding): ")) 
                    np.savetxt(str(filename),np.c_[R,V,t])
                    print 'Saved data in file: %s\n'%(str(filename))         
            else:
                pass
            return R,V,t
        else:  # Velocity Verlet    
            a = np.zeros((N+1,dim)) 
            a[0] = self.acceleration_single(R[0],M[0]) 
            for i in xrange(N):
                R[i+1],V[i+1],a[i+1],t[i+1] = S().Velocity_Verlet(R[i],V[i],a[i],t[i],dt,self.M_sun,self.acceleration_single)
            if save_dat is True:
                if str(filename) != 'noname':
                    if not '.txt' in str(filename):
                        filename += '.txt'
                    np.savetxt(str(filename),np.c_[R,V,t])
                    print 'Saved data in file: %s\n'%(str(filename))
                else:
                    filename = str(raw_input("Please provide a filename that ends with '.txt' (less coding): ")) 
                    np.savetxt(str(filename),np.c_[R,V,t])
                    print 'Saved data in file: %s\n'%(str(filename))                    
            else:
                pass
            return R, V, t
#-----------------------------------------------------------------------------------------------------------------------------------            
    def Orbit_3_body(self,time,dt,V0_other=None,filename = 'noname',stationary = False,when_save=1,fold=None):  
        """
        Note/Caution:
            The name is misleading. This method can calculate N-body system, or N-1 system (stationary sun), but
            should only be used with dynamic systems (this is also default). 
            It import 'Solvers' (use Velocity_Verlet) and uses 'acceleration_all' -and 'acceleration_stat_sun_N' from class 'Planets'.
        Args:
            time (float,int): Total time to calculate in years.
            dt (float): Time differential in years.
            V0_other (numpy.ndarray, optional): provide different initial velocities if one does not want to use __init__ values.
            stationary (bool): True; stationary Sun, False; dynamical.
            filename (numpy.ndarray ['str'],semi-optional):  IF 'stationary' is set to 'True', then user must provide an array
                containing string text names (at least 2), program aborts if this is not the case, if 'stationary' is false, no filename
                is needed.
            when_save (float,int): Sets when values are written while they are calculated, ex: when_save=2 means every other
                second calculated values is written in the files (save memory for large calculations).
            fold ('str',optional): Provide a folder where values will be stored, either existing or not (if directory does not
                exist, it will create it, use with caution). This is optional, if None-type, files are stored in running directory.
                Example: fold = 'Folder', or fold='Folder/'.
        Returns:
            None-type if data gets written (N-body sys), or R,V,t (position,Velocity,time) arrays if stationart sun. Dimensions
            follow that of the class formality. NB: NOT recommended to use stationary Sun for large calculations, this conter-
            intuitive fact is due to memory usage, and a fault of the specialized creation.       
        """
        from Solvers import Solvers as S
        dim = self.dim; n_b = self.n_b
        M = self.M
        N = int(round(time/float(dt)))
        if stationary is False: 
            if not hasattr(filename,'__len__'): # Check if filenames are given
                print 'ERROR:\nYou must provide filenames as array for all objects to be calculated (at least 2)\nAborting program.'
                sys.exit()
            if V0_other is None: # Check if one will use input vel or given vel initialized in instance
                V_prev = self.V0
            else:
                V_prev = V0_other
            R_prev = self.R0          
            n_files = self.n_b
            a_prev = self.acceleration_all(R_prev,M)
            t_prev = 0
            
            if fold is not None: # Use/create folder if this is True
                if not os.path.exists(fold):
                    os.makedirs(fold)
                fold = fold+'/'
                filename = np.array([fold + filename[i] for i in range(n_files)]) # Apply folder to name
            for j in xrange(n_files): # Create files and store first vals.
                with open(str(filename[j]),'w') as f_name:
                    data = np.zeros(2*dim+1)
                    data[0:dim] = R_prev[j]; data[dim:2*dim] = V_prev[j]; data[-1] = t_prev
                    np.savetxt(f_name,data.reshape(1,data.shape[0]))
            
            for i in xrange(N): # Calculate and store values
                R_next, V_next, a_next, t_next = S().Velocity_Verlet(R_prev,V_prev,a_prev,t_prev,dt,M,self.acceleration_all)
                R_prev = R_next;V_prev = V_next; a_prev = a_next; t_prev = t_next
                if (i+1)%when_save == 0:
                    for j in xrange(n_files):
                        with open(str(filename[j]),'a') as f_name:
                            data = np.zeros(2*dim+1)
                            data[0:dim] = R_next[j]; data[dim:2*dim] = V_next[j]; data[-1] = t_next
                            np.savetxt(f_name,data.reshape(1,data.shape[0]))
            return None       
        else: # Stationary sun
            R = np.zeros((n_b,N+1,dim))
            V = np.zeros((n_b,N+1,dim))
            a = np.zeros((n_b,N+1,dim))
            t = np.zeros(N+1)       
            R[:,0] = self.R0 
            if V0_other is None:
                V[:,0] = self.V0
            else:  
                V[:,0] = V0_other
            a[:,0] = self.acceleration_stat_sun_N(R[:,0],M)
            for i in xrange(N):
                R[:,i+1], V[:,i+1], a[:,i+1], t[i+1] = S().Velocity_Verlet(R[:,i],V[:,i],a[:,i],\
                             t[i],dt,M,self.acceleration_stat_sun_N)
            return R,V,t
#-----------------------------------------------------------------------------------------------------------------------------------   

# Tweeks in order to make Mercury precision calculation faster, not for individual use, not logically connected to system class
# however, unecessary instance in 'Orbit' is 'necessary' in order to use class attributes.
# Explenation of methods is reduced.
    def acceleration_Mercury(self,R,M):
        r = self.distance(R)
        a = self.acceleration_single(R,M[0])*(1. + self.cons/((r**2))) # GR corrected acceleration
        return a
        
    def Mercury_precession(self,time,dt,filename='default_mercury.txt',when_save=1,fold=None,prec=True):
        """
        In order to make a Newtonian and GR calculation of Mercurys orbit around the sun, one can calculate both here (in turn).
        Changing argument 'prec' to 'False' switches to Newtonian calculation. All values are stored in a file. Default filename
        is given, but no failsafe if one does not give a filename when calculating both ways in one turn from a script, it would overwrite.
        """
        from Solvers import Solvers as S
        dim = self.dim; n_b = self.n_b
        M = self.M
        yr = 365.*24*60*60 # sec
        AU = 1.5E11 # m
        c = 3E8*yr/AU # speed of light in AU/yr       
        N = int(round(time/float(dt)))
        
        if fold is not None:
            if not os.path.exists(fold):
                os.makedirs(fold)
            fold = fold+'/'
            filename = fold + filename
        else:
            pass
        # Initial vals      
        V_prev = self.V0[1]
        R_prev = self.R0[1]
        t_prev = 0
        if prec is True: # calculate with relativistic correction
            l = abs(S().cross_product2D(R_prev,V_prev))
            self.cons = 3.*(l/c)**2 # now a class attrib
            a_prev = self.acceleration_Mercury(R_prev,M)
            # Create file and store first vals
            with open(str(filename),'w') as f_name:
                data = np.zeros(dim+1)
                data[0:dim] = R_prev; data[-1] = t_prev
                np.savetxt(f_name,data.reshape(1,data.shape[0]))
        
            # Calculate...
            for i in xrange(N):
                R_next, V_next, a_next, t_next = \
                    S().Velocity_Verlet(R_prev,V_prev,a_prev,t_prev,dt,M,self.acceleration_Mercury)
                R_prev = R_next;V_prev = V_next; a_prev = a_next; t_prev = t_next
                if (i+1)%when_save == 0:
                    with open(str(filename),'a') as f_name:
                        data = np.zeros(dim+1)
                        data[0:dim] = R_next; data[-1] = t_next
                        np.savetxt(f_name,data.reshape(1,data.shape[0]))
        else: # Calculate without relativistic correction
            a_prev = self.acceleration_single(R_prev,M[0])
            with open(str(filename),'w') as f_name:
                data = np.zeros(dim+1)
                data[0:dim] = R_prev; data[-1] = t_prev
                np.savetxt(f_name,data.reshape(1,data.shape[0]))
            # Calculate...    
            for i in xrange(N):
                R_next, V_next, a_next, t_next = S().Velocity_Verlet(R_prev,V_prev,a_prev,t_prev,dt,M[0],self.acceleration_single)
                R_prev = R_next;V_prev = V_next; a_prev = a_next; t_prev = t_next
                if (i+1)%when_save == 0:
                    with open(str(filename),'a') as f_name:
                        data = np.zeros(dim+1)
                        data[0:dim] = R_next; data[-1] = t_next
                        np.savetxt(f_name,data.reshape(1,data.shape[0]))
        return None      
#----------------------------------TEST RUNS IN PLANETS--------------------------------------------------
if __name__ == '__main__':
    dim = 2
    n_b = 2
    M = np.array([1.,.5])
    R0 = np.array([[0,0],[1.,0]])
    V0 = np.array([[0,0],[0,1.]])
    test = Planets(dim,n_b,R0,V0,M)
    all_test = np.array([test.test_distance(),test.test_Force(),test.test_Kinetic_E()])
    if False in all_test:
        print 'Check error for failed test(s) before using this module.'
    else:
        print 'The few tested implemented were all good.'
    
# Last test 25.10.17
    """
-------------------------------------------------------------------------------
System with 2 celestial bodies in 2 dimensions is initialized.
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
Testing: 'distance' in 'Planets':
Success!
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
Testing: 'Force' in 'Planets':
Success!
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
Testing: 'Kinetic_E' in 'Planets':
Success!
-------------------------------------------------------------------------------
The few tested implemented were all good.
    """
