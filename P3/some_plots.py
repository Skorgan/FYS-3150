from proj3 import *
import time as measure
from mpl_toolkits.mplot3d import Axes3D
Ms = 1.
Me = 6E24/1.99E30
MJ = 1.898E27/1.99E30
G = 4.*np.pi*np.pi


def initial_vel2D(R,M_s): # Calculated from v = v0(-sin,cos) in orbit
    r = np.linalg.norm(R)
    k = np.sqrt(G*M_s/float(r))
    R_ = R[::-1]
    R_[0] *= -1.
    return k*R_/r

#-----------------------------------------------------------------------------------------------------------------------------------    
def escape_vel():
    dim = 2
    n_b = 2
    time = 5.; dt = 1E-4; eps = 1E-5
    R0 = np.array([[0,0],[1.,0]])
    M = np.array([Ms,Me])
    V0x = np.linspace(8.5,9,5)
    V0 = np.zeros((2,2)) # Just to initialize instance
    test_v = Orbit(dim,n_b,R0,V0,M)
    print 'Calculating escape velocity for values in linspace(8.5,9,5)\nTime span allotted = 5 years, dt = 1E-4 years'
    value = test_v.Escape_V(time,dt,V0x,eps)
    del test_v
    return None
#-----------------------------------------------------------------------------------------------------------------------------------    
def different_beta():
    dim = 2
    n_b = 2
    time = 5.; dt = 1E-3
    R0 = np.array([[0,0],[1.,0]])
    M = np.array([Ms,Me])
    Ve = initial_vel2D(R0[1],M_s=1.)
    V0 = np.array([[0,0],Ve])
    #-------------------------------------
    test = Orbit(dim,n_b,R0,V0,M,beta=2)
    R,V,t = test.Orbit_2_body_Earth_Sun(time,dt,method='V')
    r_o = np.array([np.linalg.norm(R[i]) for i in range(len(R[:,0]))])
    del test
    #-------------------------------------
    test = Orbit(dim,n_b,R0,V0,M,beta=2.5)
    R,V,t = test.Orbit_2_body_Earth_Sun(time,dt,method='V')
    r = np.array([np.linalg.norm(R[i]) for i in range(len(R[:,0]))])
    del test
    #-------------------------------------
    test = Orbit(dim,n_b,R0,V0,M,beta=2.9)
    Rv,Vv,t = test.Orbit_2_body_Earth_Sun(time,dt,method='V')
    rv = np.array([np.linalg.norm(Rv[i]) for i in range(len(Rv[:,0]))])
    del test
    #-------------------------------------
    plt.plot(t,1E3*r_o,label='$\\beta=2$')
    plt.plot(t,1E3*r,label='$\\beta=2.5$')
    plt.plot(t,1E3*rv,label='$\\beta=2.9$')
    plt.xlabel('$t,[\mathrm{yrs}]$',size = 15)
    plt.ylabel('$r(t),[10^{3}\mathrm{AU}]$',size=15)
    plt.legend(loc='upper left')
    plt.xlim(0,5)
    plt.grid()
    #plt.savefig('diff_beta.png')
    plt.show()
    return None
#different_beta()

#----------------------------------------------------------------------------------------------------------------------------------- 
# Plots Earths orbit for given factor of Jupiter mass Mj_fac (makes: Mj = Mj_fac*Mj) default is Mj_fac = 1.   
def E_J(time=5.,dt=1E-4,plot_orbit=True,Mj_fac = 1.):
    M = np.array([Ms,Me,Mj_fac*MJ])
    dim = 2
    n_b = 3
    R0 = np.array([[0,0],[1.,0],[5.2,0]])
    Ve = initial_vel2D(R0[1],M_s=1.)
    Vj = initial_vel2D(R0[2],M_s=1.)
    V0 = np.array([[0,0],Ve,Vj])
    test = Orbit(dim,n_b,R0,V0,M)
    fargs = ['-k','-b']
    R,V,t = test.Orbit_3_body(time,dt,stationary = True)
    if plot_orbit is True:
        plt.figure()
        plt.plot(0,0,'ro')      
        for i in range(1):
            plt.plot(R[i+1,:,0],R[i+1,:,1],fargs[i])
        plt.legend(['$\mathrm{Sun}$','$\mathrm{Earth}$'])#,'$\mathrm{Jupiter}$'])
        plt.grid()
        plt.axis('equal')
        plt.xlabel('$x, [\mathrm{AU}]$',size=15)
        plt.ylabel('$y, [\mathrm{AU}]$',size=15)
        plt.show()
        return None
    else:
        return R,V,t
#E_J(time=5.,dt=1E-4,plot_orbit=True,Mj_fac = 1.) 

#-----------------------------------------------------------------------------------------------------------------------------------    
def calc_3b_sys(time=10.,dt=1E-4,fold=None,when_s=1):
      # initial set up and calc v0y_sun from ang_mom = L about CM is zero
      R_e=np.array([1.,0])
      R_J=np.array([5.2,0])
      Ve_y = initial_vel2D(R_e,Ms)[1]
      Vj_y = initial_vel2D(R_J,Ms)[1]
      xe = R_e[0]; xJ = R_J[0]
      M_tot = Ms+Me+MJ
      X_cm = (1./M_tot)*(Me*xe + MJ*xJ)
      xe_cm = xe - X_cm
      xJ_cm = xJ - X_cm
      xs_cm = -X_cm
      R_e[0] = xe_cm
      R_J[0] = xJ_cm
      Vs_y = -(1./xs_cm)*(Me*xe_cm*Ve_y + MJ*xJ_cm*Vj_y)  
      V0 = np.array([[0,-Vs_y*.96E-3],[0,Ve_y],[0,Vj_y]]) # fudge factor here : .96E-3, num corr
      R0 = np.array([[xs_cm,0],R_e,R_J])
      M = np.array([Ms,Me,MJ])   
      dim = 2
      n_b = 3    
      filenames = np.array(['Sun_3b.txt','Earth_3b.txt','Jupiter_3b.txt'])
      test = Orbit(dim,n_b,R0,V0,M)
      print 'Calculating orbits and storing values, please wait..'
      start = measure.time()
      test.Orbit_3_body(time,dt,filename = filenames,fold=fold,when_save=when_s)  
      tot = measure.time() - start
      print "Finish!\nValues are stored in files: 'Sun_3b.txt', 'Earth_3b.txt' and 'Jupiter_3b.txt'\nin directory: 'SEJ_3b'"
      print 'Total time: ',tot
      return None
#calc_3b_sys(time=10.,dt=1E-4,fold='SEJ_3b',when_s=1)

#-----------------------------------------------------------------------------------------------------------------------------------    
def plot_SEJ(fold=None):
    filenames = np.array(['Sun_3b.txt','Earth_3b.txt','Jupiter_3b.txt'])
    if fold is not None:
        fold = fold+'/'
    filenames = np.array([fold + filenames[i] for i in range(len(filenames))])     
    for file_n in filenames:
        r = np.genfromtxt(file_n)
        plt.plot(r[:,0],r[:,1])
    plt.grid()
    plt.axis('equal')
    plt.xlabel('$x, [\mathrm{AU}]$',size=15)
    plt.ylabel('$y,[\mathrm{AU}]$',size=15)
    plt.legend(['$\mathrm{Sun}$','$\mathrm{Earth}$','$\mathrm{Jupiter}$'])
    plt.show()
    return None
#plot_SEJ(fold='SEJ_3b')

#-----------------------------------------------------------------------------------------------------------------------------------    
def compare_rad_SEJ(time=10.,dt=1E-4,fold=None,filename='Earth_3b.txt'):
    R_stat,V, t = E_J(time,dt,plot_orbit=False)
    R_stat = R_stat[1] # Earth
    if fold is not None:
        fold = fold+'/'
    filename = fold + filename
    data = np.genfromtxt(filename)
    X_dyn = data[:,0]
    Y_dyn = data[:,1]
    r_s = np.array([np.linalg.norm(R_stat[i]) for i in range(len(R_stat[:,0]))])
    r_d = np.array([np.sqrt(X_dyn[i]**2+Y_dyn[i]**2) for i in range(len(X_dyn))])
    plt.plot(t,r_s,label='$\mathrm{Stat}$')
    plt.plot(t,r_d,label='$\mathrm{Dyn}$')
    plt.legend(loc='best')
    plt.xlabel('$t, [\mathrm{yrs}]$',size=15)
    plt.ylabel('$r(t), [\mathrm{AU}]$',size=15)
    plt.grid()
    #plt.savefig('compare_dyn_stat_earth_t_10_dt_4.png')
    plt.show()
    return None
    
#compare_rad_SEJ(time=10.,dt=1E-4,fold='SEJ_3b',filename='Earth_3b.txt')    

#-----------------------------------------------------------------------------------------------------------------------------------    
if __name__=='__main__':
    if len(sys.argv) == 2:
        if str(sys.argv[1]) is 's':
            print "Running all:"
            print "Escape velocity:\n"
            escape_vel()
            print "\nDifferent beta:\n"
            different_beta()
            print "\nPlotting Earth orbit due to different Jupiter mass:"
            mass_f = float(raw_input('Type Jupiter mass factor (regular is = 1): '))
            E_J(time=5.,dt=1E-4,plot_orbit=True,Mj_fac = mass_f) 
            print "\nPlotting 3body Sun, Earth and Jupiter:\n"
            plot_SEJ(fold='SEJ_3b')
            print "\nPlotting Earth's radius comparison 3-body static and non -static Sun:\n"
            compare_rad_SEJ(time=10.,dt=1E-4,fold='SEJ_3b',filename='Earth_3b.txt')  
            print "Finish."
            sys.exit()
        else:
            print "Wrong system argument, 's' to run all."
            sys.exit()
    if len(sys.argv) > 2:
        print 'Aborting, too many system arguments.'
        sys.exit()
    print "\nYou can here recreate:\n-escape velocity value,\n-different beta plot,\n-Earth Sun plot for different factor >=1 for Jupiter mass,\n-3_body Sun, Earth and Jupiter calculation and plot,\n-Plot Earth radius comparison between static and dynamic Sun."
    print ' '
    print "If you have the files 'Sun_3b.txt', 'Earth_3b.txt' and 'Jupiter_3b.txt' in\na folder 'SEJ_3b' you can skip this part of the program (see GitHub) and\nrun the whole thing by adding system argument 's'.\nIf not, we need to create these files, even though it is not needed for\n'escape velocity' and 'different beta' part of the program."
    print ' '
    print 'The files will not be created just yet, you will be warned before the files are generated.'
    cont = str(raw_input('Do you wish to continue? (y/n): '))
    if cont is not 'y':
        print 'Program aborted.'
        sys.exit()
    print """\nTo save time, and just recreate, the time stamp, dt etc is chosen beforehand.
These values are:
Time = 10 years,
dt = 1E-4 years,
and it will write every value calculated; three files with 500k+1 numbers."""
    print 'I need a directory to store these values (in order to save the creator some\ntime).\n'
    f = str(raw_input("Is the folder/directory name 'SEJ_3b' free to be used? (y/n)"))
    if str(f) is not 'y':
        print 'Make sure it is, and run this program again'
        sys.exit()
    else:
        print 'Okay!'
        f = str(raw_input("Type 's' to start generating files, anything else aborts program. "))
        if f is not 's':
            print "Aborted."
            sys.exit()
        else:
            calc_3b_sys(time=10.,dt=1E-4,fold='SEJ_3b',when_s=1)
            print "You can now run the entire script with the system argument 's'."
            sys.exit()
            
    
    
