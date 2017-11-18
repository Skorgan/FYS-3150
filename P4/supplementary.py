# -*- coding: utf-8 -*-
from Project4 import *
import matplotlib.pyplot as plt
"""
NB! If function 'P' is run with argument 'get_vals=True', you need to uncomment the two expressions
at lines 186 and 211 in the class 'Project4'. If you copy all text files
in my GitHub, no 'writing' runs are necessary. The one exception is the last two functions, 
the files were (despite of efforts against it) overwritten and lost, so in order run, 
one must run 'several_runs()' (which I highly advice against due to the time needed = 52h) 
but check the arguments regarding folder and names.

The structure should be chronological following the assignment 4. in Computational Physics at UiO.
"""
def two_two_lattice(MCc=1E6): 
    # Run parameters
    L = 2; T_i = 1.0; T_f = 1.; Tn = 1
    ins = np.array([L,T_i,T_f,Tn,MCc])
    # Run MC 
    test = Ising_Model_2Dlattice_spin(inscript = ins)
    test.Solve(random_i = True)
    del test
    return None
#two_two_lattice()

def write_expectation_vals_cycle(): # Only use if you want to write data
    L = 20; T_i = 2.4; T_f = 2.4; Tn = 1; MCc = 5E3
    ins = np.array([L,T_i,T_f,Tn,MCc])
    file_name = 'expec_per_M5E3_T2'
    test = Ising_Model_2Dlattice_spin(inscript = ins)
    test.Solve(random_i=True,store_values=True,filename=file_name,when_save=1)
    return None
#write_expectation_vals_cycle()

def plotting_expectation_vals_cycle():
    # T=1, unordered initial    
    data = np.genfromtxt('expec_per_M5E3_T1.txt')
    E_bar = data[:,0]
    M_bar_abs = data[:,2]
    MCc = data[:,-2]
    accepted_per = data[:,-1]

    # 2) T=1, ordered initial
    data2 = np.genfromtxt('expec_per_M5E3_T1_ordered.txt')
    E_bar2 = data2[:,0]
    M_bar_abs2 = data2[:,2]
    
    accepted_per2 = data2[:,-1]

    # 3) T=2.4, unordered
    data3 = np.genfromtxt('expec_per_M5E3_T2.txt')
    E_bar3 = data3[:,0]
    M_bar_abs3 = data3[:,2]
    
    accepted_per3 = data3[:,-1]    
    
    # 4) T= 2.4, ordered initial
    data4 = np.genfromtxt('expec_per_M5E3_T2_ordered.txt')
    E_bar4 = data4[:,0]
    M_bar_abs4 = data4[:,2]
    
    accepted_per4 = data4[:,-1]    
    
    # Plot <E> and <|M|> for T=1 and T=2.4 unordered:
    plt.subplot(2,1,1)
    plt.plot(MCc,E_bar,label='$T=1.0$')
    plt.plot(MCc,E_bar3,label='$T=2.4$')
    plt.ylabel('$\\langle E \\rangle$',size=15)
    plt.legend(loc='center')
    plt.grid()
    plt.ylim(-2.1,-1)
    
    plt.subplot(2,1,2)
    plt.plot(MCc,M_bar_abs,label='$T=1.0$')
    plt.plot(MCc,M_bar_abs3,label='$T=2.4$')
    plt.ylabel('$\\langle |M| \\rangle$',size=15)
    plt.legend(loc='upper right')
    plt.grid()
    plt.ylim(0,2)
    plt.xlabel('$\mathrm{time;\,}t,\,[\mathrm{cycles}]$',size=15)
    #plt.savefig('EM_unordered.png')
    plt.show()
    
    # Plot <E> and <|M|> for T=1 and T=2.4 ordered
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(MCc,E_bar2,label='$T=1.0$')
    plt.plot(MCc,E_bar4,label='$T=2.4$')
    plt.ylabel('$\\langle E \\rangle$',size=15)
    plt.legend(loc='center')
    plt.grid()
    plt.ylim(-2.1,-1)

    plt.subplot(2,1,2)
    plt.plot(MCc,M_bar_abs2,label='$T=1.0$')
    plt.plot(MCc,M_bar_abs4,label='$T=2.4$')
    plt.ylabel('$\\langle |M| \\rangle$',size=15)
    plt.legend(loc='upper right')
    plt.ylim(0,2)
    plt.grid()
    plt.xlabel('$\mathrm{time;\,}t,\,[\mathrm{cycles}]$',size=15)
    #plt.savefig('EM_ordered.png')
    plt.show()
    # Plot acceptance probability per cycle sweep as function of cycle, T=1 and T=2.4 unordered
    plt.figure()
    plt.subplot(2,1,1)
    plt.semilogx(MCc,accepted_per,label='$T=1.0$')
    plt.semilogx(MCc,accepted_per3,label='$T=2.4$')
    plt.ylabel('$A$',size=15)
    plt.ylim(0,1)
    plt.legend(loc='upper center')
    plt.grid()

    plt.subplot(2,1,2)
    plt.semilogx(MCc,accepted_per2,label='$T=1.0$')
    plt.semilogx(MCc,accepted_per4,label='$T=2.4$')
    plt.ylabel('$A$',size=15)
    plt.xlabel('$\mathrm{Cycle}$',size=15)
    plt.legend(loc='upper left')
    plt.grid()
    #plt.savefig('accept_ordered_unordered.png')
    plt.show()
    return None
#plotting_expectation_vals_cycle()


def P(get_values=False,MCc=1E6,L=20,T=1.0): # 
    ins = np.array([L,T,T,1,MCc])
    if T == 1:
        filename = 'Energies_L%g_MCc1E%g.txt'%(L,np.log10(MCc))
    else:
        filename = 'Energies_L%g_MCc1E%gT2.txt'%(L,np.log10(MCc))
    
    if get_values is True:
        test = Ising_Model_2Dlattice_spin(inscript = ins)
        E_count, E_var = test.Solve(random_i=True)
        vals = np.append(E_var,E_count)
        if os.path.isfile(filename) is True:
            filename = 'Edefaulttttt.txt'
            print 'Filename has been changed to; Edefaulttttt.txt'
        np.savetxt(filename,np.c_[vals]) # NB first value is the variance
    else:
        vals = np.genfromtxt(filename)
        E_var = vals[0]
        E_count = vals[1:]
        ss = 1000 # When approximately steady state occurs
        E_count = np.array(E_count[ss:])/400.
        plt.hist(E_count,weights=np.ones_like(E_count)/float(len(E_count)),bins=25)
        plt.ylabel('$\mathrm{PDF}$',size=15)
        plt.xlabel('$E,\,[\mathrm{J}]$',size=15)
        print np.sqrt(E_var)/20.  # E_var is given per spin, so in order to find variance, sqrt(400)/400 = 1/20
        plt.show()
    return None
#P(T=2.4)  # Remove T=2.4 to get T=1. histo



def several_runs(): # Writes data for L=40,60,80,100 over 11 temperatures in the interval [2.2,2.3]
    
    T_i = 2.2; T_f = 2.3; Tn = 0.01
    MCc = 1E5
    L_ = np.array([40,60,80,100])
    for L in L_:
        file_name = 'data/L%sMCc%s'%(str(L),str(int(MCc)))
        ins = np.array([L,T_i,T_f,Tn,MCc])
        test = Ising_Model_2Dlattice_spin(inscript = ins)
        test.Solve(store_values=True,filename=file_name)
        del test
    return None
#several_runs()
"""
These functions do not have the files needed in my GitHub.
Comments on this is given in a report.
"""
def plot_several_runs():
    data = np.genfromtxt('data/L40MCc100000.txt')   # Ordered as  [<E>,<|M|>,Cv,chi,T]
    data2 = np.genfromtxt('data/L60MCc100000.txt')
    data3 = np.genfromtxt('data/L80MCc100000.txt')
    data4 = np.genfromtxt('data/L100MCc100000.txt')
    
    T = data[:,-1]
    
    # Energies
    E_bar = data[:,0]
    E_bar2 = data2[:,0]
    E_bar3 = data3[:,0]
    E_bar4 = data4[:,0]
    plt.plot(T,E_bar,label='$L=40$')   
    plt.plot(T,E_bar2,label='$L=60$')
    plt.plot(T,E_bar3,label='$L=80$')    
    plt.plot(T,E_bar4,label='$L=100$')    
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('$T,\,[\mathrm{J}]$',size=15)
    plt.ylabel('$\\langle E \\rangle,\,[1]$',size=15)
    #plt.savefig('EbarL40_100.png')
    plt.show()
    # <|M|>
    M_bar_abs = data[:,1]
    M_bar_abs2 = data2[:,1]
    M_bar_abs3 = data3[:,1]
    M_bar_abs4 = data4[:,1]
    plt.plot(T,M_bar_abs,label='$L=40$')   
    plt.plot(T,M_bar_abs2,label='$L=60$')
    plt.plot(T,M_bar_abs3,label='$L=80$')    
    plt.plot(T,M_bar_abs4,label='$L=100$')    
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('$T,\,[\mathrm{J}]$',size=15)
    plt.ylabel('$\\langle |M| \\rangle,\,[1]$',size=15)
    #plt.savefig('Mbar_absL40_100.png')
    plt.show()
    
    # Specific heat capacity; Cv
    Cv = data[:,2]
    Cv2 = data2[:,2]
    Cv3 = data3[:,2]
    Cv4 = data4[:,2]
    plt.plot(T,Cv,label='$L=40$')   
    plt.plot(T,Cv2,label='$L=60$')
    plt.plot(T,Cv3,label='$L=80$')    
    plt.plot(T,Cv4,label='$L=100$')    
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('$T,\,[\mathrm{J}]$',size=15)
    plt.ylabel('$C_v,\,[1]$',size=15)
    #plt.savefig('CvL40_100.png')
    plt.show()
    
    # Magnetic susceptibility; chi
    chi = data[:,3]
    chi2 = data2[:,3]
    chi3 = data3[:,3]
    chi4 = data4[:,3]
    plt.plot(T,chi,label='$L=40$')   
    plt.plot(T,chi2,label='$L=60$')
    plt.plot(T,chi3,label='$L=80$')    
    plt.plot(T,chi4,label='$L=100$')    
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('$T,\,[\mathrm{J}]$',size=15)
    plt.ylabel('$\chi,\,[1/J]$',size=15)
    #plt.savefig('chiL40_100.png')
    plt.show()
    return None
#plot_several_runs()    

def find_TC():
    data = np.genfromtxt('data/L40MCc100000.txt')   # Ordered as  [<E>,<|M|>,Cv,chi,T]
    data2 = np.genfromtxt('data/L60MCc100000.txt')
    data3 = np.genfromtxt('data/L80MCc100000.txt')
    data4 = np.genfromtxt('data/L100MCc100000.txt')
    T = data[:,-1]
    
    M_bar_abs = data[:,1]
    M_bar_abs2 = data2[:,1]
    M_bar_abs3 = data3[:,1]
    M_bar_abs4 = data4[:,1]
    
    Cv_ = data[:,2]
    Cv2 = data2[:,2]
    Cv3 = data3[:,2]
    Cv4 = data4[:,2]
    beta = 1./8
    
    Cv = np.array([Cv_,Cv2,Cv3,Cv4])
    Mbar = np.array([M_bar_abs,M_bar_abs2,M_bar_abs3,M_bar_abs4])
    L = np.array([40.,60.,80.,100.])
    Tc = np.zeros(len(L))
    o = 5 # start index in order to avoid wrong values from file
    a = 0
    for i in range(len(L)):
        Tc[i] = T[np.argmax(Cv[i,o:]) + o]
        ind = np.argmax(Cv[i,o:]) + o
        a += L[i]*(Mbar[i,ind]**(1./beta))
    a /= float(len(L))
    print a      
    
    Tc_L_lim = 0.
    for i in range(len(L)):
        Tc_L_lim += Tc[i] - a/L[i]
    Tc_L_lim /= float(len(L))
    print Tc_L_lim
    return None
#find_TC()
    
    
