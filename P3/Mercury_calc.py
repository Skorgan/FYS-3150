from proj3 import *
import time as measure
Ms = 1.
Me = 6E24/1.99E30
MJ = 1.898E27/1.99E30
G = 4.*np.pi*np.pi

# Calculate and store Mercury GR corr data file
def GR_Merc(time=100.,dt=1E-4,file_='Mercury_gr.txt',fold='Merc',when_s = 1):
    n_b = 2; dim = 2
    R0 = np.array([[0,0],[.3075,0]])  # AU
    V0 = np.array([[0,0],[0,12.44]])  # AU/yr
    Mm = 3.3E23/1.99E30       # mercury mass in solar mass
    M = np.array([1.,Mm])     
    
    test = Orbit(dim,n_b,R0,V0,M)
    print 'Calculating GR Mercury, data will be stored as columns: [x y t]'
    start = measure.time()
    test.Mercury_precession(time,dt,filename=file_,when_save=when_s,fold=fold)
    tot_t = measure.time()-start
    print 'Finish! Values are stored in file: %s\nTotal time: %g sec'%(str(file_),tot_t)
    return None

#GR_Merc(time=4.,dt=1E-4,file_='Mercury_gr_t_4_dt4.txt',fold='Merc',when_s=1) 

# Calculate and store Mercury Newtonian
def Newton_Merc(time=100.,dt=1E-4,file_='Mercury_r.txt',fold='Merc',when_s = 1):
    n_b = 2; dim = 2
    R0 = np.array([[0,0],[.3075,0]])  # AU
    V0 = np.array([[0,0],[0,12.44]])  # Au/yr
    Mm = 3.3E23/1.99E30       # mercury mass in solar mass
    M = np.array([1.,Mm])     
    
    test = Orbit(dim,n_b,R0,V0,M)
    print 'Calculating Newtonian Mercury, data will be stored as columns: [x y t]'
    start = measure.time()
    test.Mercury_precession(time,dt,filename=file_,when_save=when_s,fold=fold,prec=False)
    tot_t = measure.time() - start
    print 'Finish! Values are stored in file: %s\nTotal time: %g sec'%(str(file_),tot_t)
    return None
    
#Newton_Merc(time=4.,dt=1E-4,file_='Mercury_n_t_4_dt4.txt',fold='Merc',when_s=1)  
#-----------------------------------------------------------------------------------------------------------------------------------    
def plot_Mercury(file_='Mercury_gr.txt',fold=None):
    if fold is not None:
        if not '/' in fold:
            fold = fold+'/'  
        file_ = fold+file_
    data = np.genfromtxt(str(file_))
    #plt.figure(figsize=(10,6),dpi=80) 
    plt.plot(0,0,'ro',label='$\mathrm{Sun}$')
    plt.plot(data[:,0],data[:,1],label='$\mathrm{Mercury}$')
    plt.legend()
    plt.xlabel('$x, [\mathrm{AU}]$',size=15)
    plt.ylabel("$y, [\mathrm{AU}]$",size=15)
    plt.grid()
    plt.show()
    return None
#plot_Mercury(file_='Mercury_gr_t_4_dt4.txt',fold='Merc')    

def Check_Mercury(file_1='Mercury_gr.txt',file_2='Mercury_ne.txt',fold=None,fig=True,extrapolate=True): 
    if fold is not None:
        if not '/' in fold:
            fold = fold+'/'  
        file_1 = fold+file_1
        file_2 = fold+file_2
    rad_to_arcsec = 180.*60*60/np.pi
    # Get data
    data1 = np.genfromtxt(str(file_1))
    data2 = np.genfromtxt(str(file_2))
    x_gr = data1[:,0]; y_gr = data1[:,1] 
    r_gr = np.array([np.sqrt(x_gr[i]**2+y_gr[i]**2) for i in range(len(x_gr))])

    t = data1[:,-1]
    x_n = data2[:,0]; y_n = data2[:,1]
    r_n = np.array([np.sqrt(x_n[i]**2+y_n[i]**2) for i in range(len(x_n))])

    # NB: np.diff(here) gives approx 48,  ->(diff*step_length*dt*yr) 48*50*1E-4*365 days = 87.6 (not quite right, fix tomorrow)
    R_gr_per = np.zeros((x_gr.size/4,2))
    R_n_per = np.zeros((x_n.size/4,2))
    j = -1
    for i in xrange(1,x_gr.size-1):
        if r_gr[i-1] > r_gr[i] and r_gr[i] < r_gr[i+1]:
            j += 1
            R_gr_per[j] = x_gr[i],y_gr[i]
    o = -1
    for i in xrange(1,x_n.size-1):
        if r_n[i-1] > r_n[i] and r_n[i] < r_n[i+1]:
            o += 1
            R_n_per[o] = x_n[i],y_n[i]
    R_gr_per = R_gr_per[:j+1]
    R_n_per = R_n_per[:o+1]
    # due to computation heavy, use file that is t=4 years, dt = 1E-4 years, however extrapolation not recommended
    if extrapolate is True: 
        fac = 25. # = 100./4
        y = 0    
        for i in range(R_gr_per[:,1].size):
            y += R_gr_per[i,1]
        y_yr = y/R_gr_per[:,1].size # average increase in y per 4 years
        # In hundred years:      
        y_100 = fac*y_yr
        x_100 = np.sqrt(0.3075**2 - y_100**2)
        ang = np.arctan(y_100/x_100)*rad_to_arcsec  # linear estimate of ang size in arcsec with GR corr
        # Same extrapolation for Newtonian case:
        y = 0
        for i in range(R_n_per[:,1].size):
            y += R_n_per[i,1]
        y_yr = y/R_n_per[:,1].size
        y_100 = fac*y_yr
        x_100 = np.sqrt(0.3075**2 - y_100**2)
        ang1 = np.arctan(y_100/x_100)*rad_to_arcsec # Lineat estimate of ang size in arcsec classical
        print 'Extrapolated results (in arcseconds):'
        print 'ang GR: ',ang,'\nang Newton: ',ang1,'\nabs diff: ',abs(ang-ang1)
      
    if fig is True:   
        ang_gr = np.arctan(R_gr_per[:,1]/R_gr_per[:,0])*rad_to_arcsec
        ang_n = np.arctan(R_n_per[:,1]/R_n_per[:,0])*rad_to_arcsec
        
        time = t[-1]
        t_gr = np.linspace(0,time,ang_gr.size)     
        t_n = np.linspace(0,time,ang_n.size)
        
        plt.plot(t_gr,ang_gr,'o',label='$\mathrm{GR}\,\mathrm{correction}$')
        plt.plot(t_n,ang_n,label='$\mathrm{Newtonian}$')
        plt.grid()
        plt.xlabel('$t, [\mathrm{yr}]$',size=15)
        plt.ylabel("$\\theta'', [\mathrm{arcseconds}]$",size=15)
        plt.legend(loc='lower right')
        sav = raw_input('Save plot? (you will be asked for a name) (y/n): ')
        if str(sav) is 'y':
            filename = raw_input("Enter filename (end with '.png', saves me time) ")
            plt.savefig(str(filename))
            print 'Figure saved.'
        plt.show()       
    return None  

#Check_Mercury(file_1='Mercury_gr_t_4_dt4.txt',file_2='Mercury_n_t_4_dt4.txt',fold='Merc',fig=True) 
if __name__ == '__main__':
    if len(sys.argv) == 2:
        tt = str(sys.argv[1])
        if tt is 's':
            Check_Mercury(file_1='Mercury_gr_t_4_dt4.txt',file_2='Mercury_n_t_4_dt4.txt',fold='Merc',fig=True) 
        else:
            print 'Invalid system argument.'
            sys.exit()
    else:
    # Because it takes a long time:      
        print "Make sure you have the files 'Mercury_gr_t_r_dt4.txt' and\n'Mercury_n_t_4_dt4.txt' in your current directory or in dir: 'Merc'.\nIf this is the case, run the script by adding 's' as system argument,\nIt will then calculate the extrapolated precession value and generate one plot\nused in report.\nCheck script and comment out this part to test some plots."
        sys.exit()
    
