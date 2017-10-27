from proj3 import *
import time as measure
from mpl_toolkits.mplot3d import Axes3D

global filenames
global Planets

filenames = np.array(['Sun.txt','Mercury.txt','Venus.txt','Earth.txt','Mars.txt','Jupiter.txt','Saturn.txt',\
        'Uranus.txt','Neptun.txt','Pluto.txt'])
Planets = np.array(['$\mathrm{Sun}$','$\mathrm{Mercury}$','$\mathrm{Venus}$','$\mathrm{Earth}$',\
        '$\mathrm{Mars}$','$\mathrm{Jupiter}$','$\mathrm{Saturn}$','$\mathrm{Uranus}$','$\mathrm{Neptun}$','$\mathrm{Pluto}$'])

def solarsystem(time=20.,dt=1E-4,fold='SolSys',when_s=50.):
    yr = 365. # days
    dim=3
    n_b=10
    # initial pos and vel for the Sun, Mercury ,..., Pluto
    R_sun = np.array([2.192279039707197E-3,5.762703037438353E-3,-1.295390859661947E-4])
    V_sun = np.array([-5.271633313328310E-6,5.466209510902422E-6,1.241657817440020E-7])
    
    R_mercury = np.array([-2.327991908789227E-1,-3.904847890125498E-1,-1.095001531240173E-2])
    V_mercury = np.array([1.851753087417236E-2,-1.299992586068514E-2,-2.761863418814871E-3])
    
    R_venus = np.array([-6.859478142184278E-1,2.103240483071906E-1,4.238650573048407E-2])
    V_venus = np.array([-5.867791574672555E-3,-1.947450171655608E-2,7.121929284791517E-5])   
    
    R_earth = np.array([8.679004239479287E-1,4.962296089500106E-1,-1.1566030887927330E-4])
    V_earth = np.array([-8.771827063871312E-3,1.491359022662084E-2,-3.259453001542365E-7])  
    
    R_mars = np.array([-1.586900999129547,5.001555256353886E-1,4.922989415831799E-2])
    V_mars = np.array([-3.638385221590926E-3,-1.216093048374219E-2,-1.656655640290856E-4]) 
    
    R_jupiter = np.array([-4.560765931270251,-2.9577037362943823,1.142760349921668E-1])
    V_jupiter = np.array([4.015936257121555E-3,-5.973351131084072E-3,-6.505975036158468E-5]) 
    
    R_saturn = np.array([-3.211170438350750E-1,-1.005045599494378E1,1.875288301719218E-1])
    V_saturn = np.array([5.270389045341049E-3,-1.958598975300664E-4,-2.065248375899736E-4]) 
    
    R_uranus = np.array([1.784901934033485E1,8.829883531158330,-1.984425267310287E-1])
    V_uranus = np.array([-1.772851373238405E-3,3.341974951808675E-3,3.527445872207531E-5]) 
    
    R_neptun = np.array([2.861925865229266E1,-8.803228391659726,-4.782740103805204E-1])
    V_neptun = np.array([9.022088893823943E-4,3.018763798642282E-3,-8.336640811463031E-5]) 
    
    R_pluto = np.array([1.056531600375861E1,-3.171056197967332E1,3.370971929509492E-1])
    V_pluto = np.array([3.035458138931598E-3,3.257531361158181E-4,-9.164683322714200E-4]) 
    
    R0 = np.array([R_sun,R_mercury,R_venus,R_earth,R_mars,R_jupiter,R_saturn,R_uranus,R_neptun,R_pluto])
    V0 = np.array([V_sun,V_mercury,V_venus,V_earth,V_mars,V_jupiter,V_saturn,V_uranus,V_neptun,V_pluto])*yr
    M = np.array([1.99E30,3.3E23,4.9E24,6E24,6.6E23,1.9E27,5.5E26,8.8E25,1.03E26,1.3E22])/1.99E30

    print 'Calculating orbits and storing values, please wait...'
    test = Orbit(dim,n_b,R0,V0,M)
    test.Orbit_3_body(time,dt,filename = filenames,when_save=when_s,fold=fold)  
    print 'Finish! Values are written to individual files.'
    return None

#-----------------------------------------------------------------------------------------------------------------------------------  
def plot_solarsystem2D(fold=None):
    plt.figure(figsize=(10,6),dpi=80)
    Plan = Planets
    if fold is not None:
        if not '/' in fold:
            fold = fold+'/' 
        files = np.array([fold + filenames[i] for i in range(len(filenames))])    
    for i in range(len(files)):
        r = np.genfromtxt(files[i])
        plt.plot(r[:,0],r[:,1],label=Plan[i])
    plt.legend(loc='best')
    plt.xlabel('$x, [\mathrm{AU}]$',size=15)
    plt.ylabel('$y,[\mathrm{AU}]$',size=15)
    plt.xlim(-25,40)
    plt.ylim(-41,30)
    plt.grid() 
    #plt.savefig('solsys_t_20_dt_4.png') 
    plt.show()
    return None

#-----------------------------------------------------------------------------------------------------------------------------------    
def plot_solarsystem3D(fold=None):
    Plan = Planets
    if fold is not None:
        if not '/' in fold:
            fold = fold+'/'
        files = np.array([fold + filenames[i] for i in range(len(filenames))])   
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    for i in range(len(files)-5):
        r = np.genfromtxt(files[i])
        ax.plot(r[:,0],r[:,1],r[:,2],label=Plan[i])
    plt.legend(loc='upper left')
    ax.set_xlabel('$x, [\mathrm{AU}]$',size=15)
    ax.set_ylabel('$y, [\mathrm{AU}]$',size=15)
    ax.set_zlabel('$z, [\mathrm{AU}]$',size=15)
    #plt.savefig('solsys_t_20_dt_4_3d.png')
    plt.show()
    return None
#-----------------------------------------------------------------------------------------------------------------------------------  
if __name__=='__main__':
    def plots(fold_):
        d2 = str(raw_input('Plot 2D? (y/n): '))
        d3 = str(raw_input('Plot 3D? (y/n): '))
        if str(d2) is 'y' and str(d3) is 'y':
            print "You must close the first plot in order to produce the other."
        if str(d2) is 'y':
            plot_solarsystem2D(fold=fold_)
        if str(d3) is 'y':
            plot_solarsystem3D(fold=fold_)
        return None
    
    if len(sys.argv)> 2:
        print "Aborted, too many sys args."
        sys.exit()
            
    if len(sys.argv) == 2:
        fold_ = str(sys.argv[1])
        if fold_ is '1':
            fold = None
        plots(fold_) # Calls above plot function
        sys.exit()

    else:
        print 'Here you can plot the solar system in 2D or in 3D (3D up to mars).\nFirst of:'    
        print """If you have the planet files (all), and stored them in a directory (or current) you can run this script as:\nterminal> python solsys.py foldername\nSet foldername as '1' if the files is in current directory.\nThe files must either be created here, or downloaded from GitHub. No figures\nwill be saved."""
        print "Contiue ONLY if you want to create the files, and be sure that no other\nfiles of names 'Planetname.txt' (with capital letter, includes 'Sun.txt') is in\nyour current directory or in a chosen folder."
        print "\nFiles will be created with following parameters:\ntime=20 years\ndt=1E-4 years\nsaving every 50 iteration (every 5E-3 year)"
        print '(you will be asked for a folder before the calculation starts)'
        cont = str(raw_input('Do you wish to continue? (y/n): '))
        if str(cont) is 'n':
            print "Aborted."
            sys.exit()
        if str(cont) is 'y':
            fold_ = str(raw_input('Type in a directory/folder name you wish to store the text files\ntype 1 for current directory: '))
            if str(fold_) is '1':
                sure = str(raw_input('You have chosen current working directory,\ncalculation will start if you continue (y/n): '))
                if str(sure) is 'y':
                    fold_ = None
                else:
                    print "Aborted."
                    sys.exit()
            else:
                check = str(raw_input('Use folder name: %s ? Calculation will start if you continue. (y/n): '%fold_))
                if str(check) is not 'y':
                    print 'Aborted.'
                    sys.exit()
                else:
                    pass
        else:
            print "Neither 'y' or 'n' was typed, aborting due to danger of overwriting files."
            sys.exit()
        print "Starting calculation"
        solarsystem(time=20,dt=1E-4,fold=fold_,when_s=50) # This line promts the calculation and creation of files (and optional folder)
        now = str(raw_input('Wish to choose plots now? (y/n): '))
        if str(now) is 'y':
            plots(fold_)
        else:
            print "We are finish!"
            sys.exit()
        
        
        
