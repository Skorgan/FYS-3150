from proj3 import * # includes
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
        
def sys_2body(time=1.,dt=1E-4,plot2D=False,plot_radius_comp=False,check_e = False,check_ang=False):
    dim = 2
    n_b = 2
    R0 = np.array([[0,0],[1.,0]])
    M = np.array([Ms,Me])
    Ve = initial_vel2D(R0[1],M_s=1.)
    V0 = np.array([[0,0],Ve])
    if plot2D is True:
        print '-Earth-Sun'
        plt.figure()
        test = Orbit(dim,n_b,R0,V0,M) # Gives solar system data, creates instance
        R,V,t = test.Orbit_2_body_Earth_Sun(time,dt,method='E')
        plt.plot(R[:,0],R[:,1],label='$\mathrm{Earth}$')
        plt.plot(0,0,'or',label='$\mathrm{Sun}$')
        plt.legend()
        plt.axis('equal')
        plt.xlabel('$x, [\mathrm{AU}]$',size=15)
        plt.ylabel('$y, [\mathrm{AU}]$',size=15)
        plt.grid()
        plt.show()
        del test
    
    if plot_radius_comp is True:
        print '-Radius comparison'
        plt.figure()
        test_E = Orbit(dim,n_b,R0,V0,M)
        R,V,t = test_E.Orbit_2_body_Earth_Sun(time,dt,method='E')
        del test_E
        
        test_V = Orbit(dim,n_b,R0,V0,M)
        Rv,Vv,t = test_V.Orbit_2_body_Earth_Sun(time,dt,method='V')
        del test_V
        
        r_Euler = np.array([np.linalg.norm(R[i]) for i in range(len(R[:,0]))])
        r_Vel = np.array([np.linalg.norm(Rv[i]) for i in range(len(Rv[:,0]))])
        
        plt.plot(t,r_Euler,'-b',label='$\mathrm{EC},$')
        plt.plot(t,r_Vel,'-r',label='$\mathrm{VV}$')
        plt.grid()
        plt.ylabel('$r, [\mathrm{AU}]$',size=15)
        plt.xlabel('$t, [\mathrm{yrs}]$',size=15)
        plt.text(.1,1.32,'$dt=%g$'%dt,size=20)
        plt.legend(loc='best')
        plt.show()  
    
    if check_e is True:
        print '-Energy conservation'
        # Calculate with Euler Chromer method
        test_EC = Orbit(dim,n_b,R0,V0,M)
        Re,Ve,t = test_EC.Orbit_2_body_Earth_Sun(time,dt,method='E')
        Ke = test_EC.Kinetic_E(Ve,M[1])
        Te = test_EC.Potential_E(Re,M[1],len(Re[:,0]))
        Tot_e = Ke+Te
        Rel_tot_e_EC = Tot_e/Tot_e.max()
        del test_EC
        # Same for Velocity Verlet solver
        test_VV = Orbit(dim,n_b,R0,V0,M)
        Rv,Vv,t = test_VV.Orbit_2_body_Earth_Sun(time,dt,method='V')
        Kv = test_VV.Kinetic_E(Vv,M[1])
        Tv = test_VV.Potential_E(Rv,M[1],len(Rv[:,0]))
        Tot_e = Kv+Tv
        Rel_tot_e_VV = Tot_e/Tot_e.max()
        del test_VV
        # Plot
        plt.subplot(2,1,1)
        plt.plot(t,Rel_tot_e_EC,label='$\mathrm{EC}$')
        plt.ylabel('$|E_{\mathrm{tot}}|^*$',size=20)
        plt.legend(loc='best')
        plt.subplot(2,1,2)
        plt.plot(t,Rel_tot_e_VV,label='$\mathrm{VV}$')
        plt.ylabel('$|E_{\mathrm{tot}}|^*$',size=20)
        plt.xlabel('$t, [\mathrm{yrs}]$',size=20)
        plt.legend()
        #plt.savefig('rel_energycons_EC_VV_one_year_N_1E4.png')
        plt.show()
    if check_ang is True:
        print '-Conservation of angular momentum'
        plt.figure()
        dt = np.array([1E-1,1E-2,1E-3,1E-4])
        test_VV = Orbit(dim,n_b,R0,V0,M)
        for i in range((len(dt))):
            test_VV = Orbit(dim,n_b,R0,V0,M)
            R,V,t = test_VV.Orbit_2_body_Earth_Sun(time,dt[i],method='V')
            pv = test_VV.Momentum(M[1],V)
            Lv = test_VV.AngMomentum(R,pv,abs_val=True)
            plt.plot(t,Lv,label='$dt=%g$'%dt[i])          
        del test_VV
        plt.legend()
        plt.ylabel('$L_z(t),[M_{\odot}\mathrm{AU}^2/\mathrm{yrs}^2]$',size=15)
        plt.xlabel('$t, [\mathrm{yrs}]$',size=15)
        plt.show()
    return None

if __name__ == '__main__':
     print "Here, a preset of plots can be created which was used in a rapport.\nThese are the 'comparison' and the angular momentum plots.\nI recommend one year with time differential 1E-4 or 1E-3.\nNo figures will be saved."
     print 'So:'
     time_ = float(raw_input('Enter total time (in years) calculations will run: '))
     dt_ = float(raw_input('Choose time differential in years (example: 1E-4): '))
     n = int(round(time_/dt_))
     print 'This corresponds to N-1 = %g points to calculate whenever it is required,'%n
     okay = str(raw_input('Is that okay? (y/n): '))
     if okay is not 'y':
         print 'Okay, aborting.. Run again'
         sys.exit() 
     else:
         pass
     print 'Now choose plots:\n'
     plot_2D = str(raw_input('Plot Earth-Sun? (y/n): '))
     print ' '
     plot_rad_comp = str(raw_input('Plot radius comparison between EC and VV method? (y/n): '))
     print ' '
     chec_tot_e = str(raw_input('Plot comparison of energy conservation between EC and VV method? (y/n): '))
     print ' '
     chec_tot_ang = str(raw_input('Plot ang_mom conservation for VV method (this takes time, runs four\ntimes with different dt)? (y/n): '))
     print ' '
     if plot_2D is 'n' and plot_rad_comp is 'n':
         if chec_tot_e is 'n' and chec_tot_ang is 'n':
             print 'There is nothing more to plot, so nothing will happen.'
             sys.exit()
         else:
             pass
     else:
         pass       
     print 'Okay, working..'
     if plot_2D is 'y':
         plot_2D = True
     else:
         plot_2D = False
     
     if plot_rad_comp is 'y':
         plot_rad_comp = True
     else:
         plot_rad_comp = False
     
     if chec_tot_e is 'y':
         chec_tot_e = True
     else:
         chec_tot_e = False   
     
     if chec_tot_ang is 'y':
         chec_tot_ang = True
     else:
         chec_tot_ang = False 
     sys_2body(time=time_,dt=dt_,plot2D=plot_2D,plot_radius_comp=plot_rad_comp,\
             check_e = chec_tot_e,check_ang=chec_tot_ang)
