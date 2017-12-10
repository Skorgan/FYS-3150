import numpy as np
from scipy.optimize import curve_fit as fitting
import matplotlib.pyplot as plt
from scipy.special import gamma as G
"""
Follows order in report, creates figures that are not direct plot of data.
Name of the functions does not reflect figure numbers.
"""
global fits,Patriarca

def fits(m,nu):
    return m**(-(1+nu)) 
    
def Patriarca(x,lam):
    n = 1. + 3.*lam/(1.-lam)
    return (n**n/G(n))*x**(n-1.)*np.exp(-n*x)
# Saving
def first(w):
    H1 = np.genfromtxt('textfiles/HistModel2_lam0_25.txt')
    H2 = np.genfromtxt('textfiles/HistModel2_lam0_5.txt')
    H3 = np.genfromtxt('textfiles/HistModel2_lam0_9.txt')
    m = H1[:,1]
    dm = m[1]-m[0]
    m = np.append(0,m[:-1])
    lam025 = H1[:,0]/float(sum(H1[:,0]))
    lam05 = H2[:,0]/float(sum(H2[:,0]))
    lam09 = H3[:,0]/float(sum(H3[:,0]))
    cs = ['ro','bo','go']
    cs2 = ['r','b','g']
    #---------------------------------------------------
    
    def end_powerlaw():
        here = int(np.where(m==1.5)[0])
        here2 = int(np.where(m==4.)[0])
        m_ = m[here:here2+1]
        lam025_ = Patriarca(m_,0.25)*dm
        p = np.zeros(len(lam025_))
        nu_trial = np.linspace(1.9,2.,1000)
        delta = 1E6
        nu_use = 0
        for i in xrange(len(nu_trial)):
           p_trial = np.array([fits(mm,nu_trial[i])*dm for mm in m_])
           delta_ = sum(abs(p_trial-lam025_))
           if delta_ < delta:
               p = p_trial
               delta = delta_
               nu_use = nu_trial[i]
               del delta_
        
        plt.semilogy(m_[::10],p[::10],cs[0],label='$\mathrm{Pareto\,power;\,}\\nu=%g$'%nu_use)
        plt.semilogy(m_,lam025_,cs2[0],label='$\mathrm{Param;\,}\lambda=0.25$')
        plt.ylabel('$P_m\,[1]$',size=15)  
        plt.xlabel('$\mathrm{Income;\,} m\, [1]$',size=15)
        plt.legend()
        plt.grid()
        plt.show()
        
        lam05_ = Patriarca(m_,0.5)*dm
        nu_trial = np.linspace(1,3,1000)
        delta = 1E6
        for i in xrange(len(nu_trial)):
            p_t = np.array([fits(mm,nu_trial[i])*dm for mm in m_])
            delta_ = sum(abs(p_t-lam05_))
            if delta_ < delta:
                p = p_t
                delta = delta_
                nu_use = nu_trial[i]
                del delta_
        plt.subplot(2,1,1)
        plt.semilogy(m_[::10],p[::10],cs[1],label='$\mathrm{Pareto\,power;\,}\\nu=%.2f$'%nu_use)
        plt.semilogy(m_,lam05_,cs2[1],label='$\mathrm{Param;\,}\lambda=0.5$')
        plt.legend(loc='best')
        plt.ylabel('$P_m$',size=15)
        plt.ylim(1E-6,1E-1)
        plt.grid()
        lam09_ = Patriarca(m_,0.9)*dm
        delta = 1E6
        nu_trial = np.linspace(6,8,1000)
        for i in xrange(len(nu_trial)):
            p_t = np.array([fits(mm,nu_trial[i])*dm for mm in m_])
            delta_ = sum(abs(p_t-lam09_))
            if delta_ < delta:
                p = p_t
                delta = delta_
                nu_use = nu_trial[i]
                del delta_
        plt.subplot(2,1,2)
        plt.semilogy(m_[::10],p[::10],cs[2],label='$\mathrm{Pareto\,power;\,}\\nu=%.2f$'%nu_use)
        plt.semilogy(m_,lam09_,cs2[2],label='$\mathrm{Param;\,}\lambda=0.9$')
        plt.legend(loc='best')
        plt.grid()
        plt.xlabel('$\mathrm{Income;\,} m\, [1]$',size=15)
        plt.ylabel('$P_m$',size=15)
        plt.show()    
        return None
    
    def check_gamma_distr():
        plt.plot(m,lam025,cs[0])
        plt.plot(m,lam05,cs[1])
        plt.plot(m,lam09,cs[2])
        plt.xlim(0,3)
        
        x = np.linspace(0,5,1000)
        j = 0
        for lam in [0.25,.5,.9]:
            p = np.array([Patriarca(x[i],lam) for i in xrange(len(x))])      
            plt.plot(x,p*dm,cs2[j],label='$\lambda = %g$'%lam)
            j += 1
        plt.grid()
        plt.legend()
        plt.xlabel('$\mathrm{Income;\,} m\, [1]$',size=15)
        plt.ylabel('$P_m\,[1]$',size=15)    
        plt.show()
        return None
    if w == 1:
        check_gamma_distr()
    elif w == 2:
        end_powerlaw()
        
    return None

#------------ With nearest Ne----------------
def second():
    cs = ['ro','bo','go','ko']
    cs2 = ['r','b','g','k']
    H1 = np.genfromtxt('textfiles/HistModel3N2_alpha0_5.txt')
    H2 = np.genfromtxt('textfiles/HistModel3N2_alpha_1.txt')
    H3 = np.genfromtxt('textfiles/HistModel3N2_alpha1_5.txt')
    H4 = np.genfromtxt('textfiles/HistModel3N2_alpha_2.txt')
    m = H1[:,1]
    dm = m[1]-m[0]
    m = np.append(0,m[:-1])
    alph0_5 = H1[:,0]/float(sum(H1[:,0]))
    alph1 = H2[:,0]/float(sum(H2[:,0]))
    alph1_5 = H3[:,0]/float(sum(H3[:,0]))
    alph2 = H4[:,0]/float(sum(H4[:,0]))
    
    def nearest_analysis():
        here = int(np.where(m==4.)[0])
        here2 = int(np.where(m==20)[0])
        m_ = m[here:here2+1]
        alph0_5_ = alph0_5[here:here2+1]
        alph1_ = alph1[here:here2+1]
        alph1_5_ = alph1_5[here:here2+1]
        alph2_ = alph2[here:here2+1]
        alph0_5_ = alph0_5_[::20]
        alph1_ = alph1_[::20]
        alph1_5_ = alph1_5_[::20]
        alph2_ = alph2_[::20]
        m_ = m_[::20]
        alphas = np.array([alph0_5_,alph1_,alph1_5_,alph2_])
        nu_fit = np.zeros(len(alphas))
        al = np.array([0.5,1.0,1.5,2.0])
        for i in xrange(nu_fit.size):
            nu_fit[i] = fitting(fits,m_,alphas[i]/dm)[0]
            plt.semilogy(m_[::3],(fits(m_,nu_fit[i])*dm)[::3],cs[i],label='$\\nu=%.2f$'%nu_fit[i])
            plt.semilogy(m_,alphas[i],cs2[i],label='$\\alpha=%.2f$'%al[i])
        
        plt.grid()
        plt.legend(loc='best')
        plt.xlim(m_.min(),m_.max())
        plt.ylim(4E-6,1E-3)
        plt.xlabel('$\mathrm{Income;\,} m\, [1]$',size=15)
        plt.ylabel('$P_m\,[1]$',size=15)   
        plt.show()
        
        return None
    nearest_analysis()
    return None
    

#------- Nearest and saving--------
def third():
    H1 = np.genfromtxt('textfiles/HistModel4N2_alpha0_5_lam_0_4.txt')
    H2 = np.genfromtxt('textfiles/HistModel4N2_alpha_1_lam_0_4.txt')
    H3 = np.genfromtxt('textfiles/HistModel4N2_alpha1_5_lam_0_4.txt')
    H4 = np.genfromtxt('textfiles/HistModel4N2_alpha_2_lam_0_4.txt')
    alph0_5 = H1[:,0]/float(sum(H1[:,0]))
    alph1 = H2[:,0]/float(sum(H2[:,0]))
    alph1_5 = H3[:,0]/float(sum(H3[:,0]))
    alph2 = H4[:,0]/float(sum(H4[:,0]))
    m = H1[:,1]
    m = np.append(0,m[:-1])
    cs2 = ['r','b','g','k']
    
    def nearest_save_analysis():
        here = int(np.where(m==1)[0])
        here2 = int(np.where(m==10)[0])
        m_ = m[here:here2+1]
        alph0_5_ = alph0_5[here:here2+1]
        alph1_ = alph1[here:here2+1]
        alph1_5_ = alph1_5[here:here2+1]
        alph2_ = alph2[here:here2+1]
        alph0_5_ = alph0_5_[::10]
        alph1_ = alph1_[::10]
        alph1_5_ = alph1_5_[::10]
        alph2_ = alph2_[::10]
        m_ = m_[::10]
        alphas = np.array([alph0_5_,alph1_,alph1_5_,alph2_])
        al = np.array([.5,1.,1.5,2.])
        for i in xrange(4):
            plt.loglog(m_,alphas[i],cs2[i],label='$\\alpha=%.2f$'%al[i])
        plt.grid()
        plt.legend(loc='best')
        plt.xlim(m_.min(),m_.max())
        plt.xlabel('$\mathrm{Income;\,} m\, [1]$',size=15)
        plt.ylabel('$P_m\,[1]$',size=15)   
        plt.show()
        return None
    nearest_save_analysis()
    return None

#---- Nearest, former alpha = 1,2-------------
def fourth(w):
    cs2 = ['r','b','g','k']
    cs = ['ro','bo','go','ko']
    if w == 1:
        H1 = np.genfromtxt('textfiles/HistModel5_gamma_1_alpha_1_0.txt')
        H2 = np.genfromtxt('textfiles/HistModel5_gamma_2_alpha_1_0.txt')
        m = H1[:,1]
        m = np.append(0,m[:-1])
        
        here = int(np.where(m==3)[0])
        here2 = int(np.where(m==15)[0])
        m = m[here:here2+1]
        dm = m[1]-m[0]
        gam1 = H1[here:here2+1,0]/float(sum(H1[:,0]))
        gam2 = H2[here:here2+1,0]/float(sum(H2[:,0]))
        gammas = np.array([gam1[::20],gam2[::20]])
        m = m[::20]
        nu_fit = np.zeros(len(gammas))
        gams = np.array([1,2,3,4])
        for i in xrange(nu_fit.size):
            nu_fit[i] = fitting(fits,m,gammas[i]/dm)[0]
            plt.loglog(m[::3],(fits(m,nu_fit[i])*dm)[::3],cs[i],label='$\\nu=%.2f$'%nu_fit[i])
            plt.loglog(m,gammas[i],cs2[i],label='$\gamma=%.2f$'%gams[i])
        plt.legend(loc='best')
        plt.xlim(m.min(),m.max())
        plt.ylim(4E-6,1E-3)
        plt.grid()
        plt.xlabel('$\mathrm{Income;\,} m\, [1]$',size=15)
        plt.ylabel('$P_m\,[1]$',size=15)   
        plt.show()
    elif w==2:
        H1 = np.genfromtxt('textfiles/HistModel5_gamma_3_alpha_2_0.txt')
        H2 = np.genfromtxt('textfiles/HistModel5_gamma_4_alpha_2_0.txt')
        m = H1[:,1]
        m = np.append(0,m[:-1])
        
        here = int(np.where(m==4.4)[0])
        here2 = int(np.where(m==25)[0])
        m = m[here:here2+1]
        dm = m[1]-m[0]
        gam1 = H1[here:here2+1,0]/float(sum(H1[:,0]))
        gam2 = H2[here:here2+1,0]/float(sum(H2[:,0]))
        gammas = np.array([gam1[::20],gam2[::20]])
        m = m[::20]
        
        nu_fit = np.zeros(len(gammas))
        gams = np.array([3,4])
        for i in xrange(nu_fit.size):
            nu_fit[i] = fitting(fits,m,gammas[i]/dm)[0]
            plt.loglog(m[::3],(fits(m,nu_fit[i])*dm)[::3],cs[i],label='$\\nu=%.2f$'%nu_fit[i])
            plt.loglog(m,gammas[i],cs2[i],label='$\gamma=%.2f$'%gams[i])
        plt.legend(loc='best')
        plt.xlim(m.min(),m.max())
        plt.ylim(4E-6,5E-4)
        plt.grid()
        plt.xlabel('$\mathrm{Income;\,} m\, [1]$',size=15)
        plt.ylabel('$P_m\,[1]$',size=15)   
        plt.show()
    return None

# Nearest, former, savings
def fifth(w):
    cs2 = ['r','b','g','k']
    cs = ['ro','bo','go','ko']
    if w == 1:
        H1 = np.genfromtxt('textfiles/HistModel6lam0_4_gamma_1_alpha_1_0.txt')
        m = H1[:,1]
        here = int(np.where(m==.9)[0])
        here2 = int(np.where(m==8)[0])
        m = m[here:here2+1]
        dm = m[1]-m[0]
        gam1 = H1[here:here2+1,0]/float(sum(H1[:,0]))
        gammas = np.array([gam1[::20]])
        m = m[::20]
        nu_fit = np.zeros(len(gammas))
        gams = np.array([1])
        for i in xrange(nu_fit.size):
            Par = np.array([Patriarca(m_,.3) for m_ in m])*dm
            plt.semilogy(m[::3],Par[::3],cs[3],label='$\mathrm{Param;}\,\lambda=0.3$')
            nu_fit[i] = fitting(fits,m,gammas[i]/dm)[0]
            plt.semilogy(m[::3],(fits(m,nu_fit[i])*dm)[::3],cs[0],label=\
                '$\mathrm{Power;\,}\\nu=%.2f$'%nu_fit[i])
            plt.semilogy(m,gammas[i],cs2[1],label=\
                '$\gamma=%.1f,\,\lambda=0.4,\,\\alpha=1.0$'%gams[i])
        plt.legend(loc='lower left')
        plt.xlim(m.min(),6)
        plt.ylim(1E-6,1E-2)
        plt.grid()
        plt.xlabel('$\mathrm{Income;\,} m\, [1]$',size=15)
        plt.ylabel('$P_m\,[1]$',size=15)   
        plt.show()
    elif w==2:
        H2 = np.genfromtxt('textfiles/HistModel6lam0_4_gamma_4_alpha_2_0.txt')
        m = H2[:,1]
        
        here = int(np.where(m==0.8)[0])
        here2 = int(np.where(m==7)[0])
        m = m[here:here2+1]
        
        dm = m[1]-m[0]
        gam2 = H2[here:here2+1,0]/float(sum(H2[:,0]))
        gammas = np.array([gam2[::20]])
        m = m[::20]
        
        nu_fit = np.zeros(len(gammas))
        gams = np.array([4])
        for i in xrange(nu_fit.size):
            nu_fit[i] = fitting(fits,m,gammas[i]/dm)[0]
            Par = np.array([Patriarca(m_,.1) for m_ in m])*dm
            plt.semilogy(m[::3],Par[::3],cs[3],label='$\mathrm{Param;}\,\lambda=0.1$')
            plt.semilogy(m[::3],(fits(m,nu_fit[i])*dm)[::3],cs[0],label=\
                '$\mathrm{Power;}\,\\nu=%.2f$'%nu_fit[i])
            plt.semilogy(m,gammas[i],cs2[1],label=\
                '$\gamma=%.1f,\,\lambda=0.4,\,\\alpha=2.0$'%gams[i])
        plt.legend(loc='lower left')
        plt.xlim(m.min(),m.max())
        plt.grid()
        plt.xlabel('$\mathrm{Income;\,} m\, [1]$',size=15)
        plt.ylabel('$P_m\,[1]$',size=15)   
        plt.show()
    
    return None
    