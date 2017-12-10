from Market import *
import matplotlib.pyplot as plt


def write_Model1(tot_money=500,N=500,Nt=1E3,tot_time=1E5,dm=0.05,filename='noname',fold=False):
    test = Market(tot_money,N)
    test.Deals(Nt,tot_time,dm,filename,fold,model=1)
    del test
    return None

def plot_Model1(filename,fold=False):
    if fold is False:
        fil = filename
        filb = 'Hist' + filename
    else:
        fil = fold + '/' + filename
        filb = fold + '/' + 'Hist' + filename
    avg = np.genfromtxt(fil)
    Histogram = np.genfromtxt(filb)
    H = Histogram[:,0]; bins = Histogram[:,1]
    dm = bins[1] - bins[0]
    beta = len(avg)/sum(avg) # beta = 1/(sum(m)/N) = 1/<m>
    # Plot logarithm :
    # Extract w from numerical results
    m_num = bins
    P_num = H/float(sum(H))
    plt.semilogy(m_num,P_num,label='$\mathrm{Numeric}$') # m spacing is in the bins created

    # Add analytic
    m = np.linspace(0,6,10)
    P_an = beta*np.exp(-beta*m)*dm
    #------
    plt.semilogy(m,P_an,label='$\mathrm{Gibbs}$')
    plt.xlim(0,6)
    plt.ylim(1E-4,.1)
    plt.xlabel('$\mathrm{Income;}\, m\, [1]$',size=15)
    plt.ylabel('$P_m\,[1]$',size=15)
    plt.legend(loc='best')
    plt.grid()
    plt.show()
    return None

#------------Model 1------------------
#write_Model1(Nt=1E3,tot_time=1E5,dm=.05,filename='Model1.txt',fold='Model1')

#----------------------------------------------------------------------------------
def write_Model2(tot_money=500,N=500,Nt=1E3,tot_time=1E6,dm=0.05,filename='noname',fold=False):
    lams = np.array([.25,.5,.9])
    for lamb in lams:
        lam = '%g'%lamb
        lam = lam[:lam.find('.')]+ '_' + lam[lam.find('.')+1:]
        f_name = filename + '_' + 'lam' + lam 
        inst = Market(tot_money,N)
        inst.Deals(Nt,tot_time,dm,f_name,fold,model=2,lam=lamb)
        del inst
    return None

def plot_Model2(filename,fold):
    lams = np.array([.25,.5,.9])
    names = []
    namesb = []
    i = 0
    for lamb in lams:
        lam = '%g'%lamb
        lam = lam[:lam.find('.')]+ '_' + lam[lam.find('.')+1:]
        names.append(filename + '_' + 'lam' + lam) 
        namesb.append('Hist' + names[i])
        if fold is not False:
            namesb[i] = fold + '/' + namesb[i]
        names[i] += '.txt'
        namesb[i] += '.txt'
        i += 1
    del names
    for i in xrange(len(namesb)):
        Histogram = np.genfromtxt(namesb[i])
        H = Histogram[:,0]; m_num = Histogram[:,1]
        P_num = H/sum(H)
        plt.plot(m_num,P_num,label='$\lambda = %g$'%lams[i])
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel('$\mathrm{Income;}\, m\, [1]$',size=15)
    plt.ylabel('$P_m\,[1]$',size=15)
    plt.xlim(0,3)
    plt.show()
    return None
    
#------------Model 2---------------
#write_Model2(Nt=1E3,tot_time=1E5,dm=0.05,filename='Model2',fold='Model2') # Used; t=12.35 min
#-------------------------------------------------------------------------------------
def write_Model3(tot_money=500,N=500,Nt=1E3,tot_time=1E7,dm=0.02,filename='noname',fold=False):
    alphas = np.array([.5,1.,1.5,2.])
    for al_ in alphas:
        al_s = '%g'%al_
        al_s = al_s[:al_s.find('.')]+ '_' + al_s[al_s.find('.')+1:]
        f_name = filename + '_' + 'alpha' + al_s 
        inst = Market(tot_money,N)
        inst.Deals(Nt,tot_time,dm,f_name,fold,model=3,alpha=al_)
        del inst
    return None

def plot_Model3(filename,fold=False):
    alphas = np.array([.5,1.,1.5,2.])
    names = []
    namesb = []
    i = 0
    
    for alpha in alphas:
        al = '%g'%alpha
        al = al[:al.find('.')]+ '_' + al[al.find('.')+1:]
        names.append(filename + '_' + 'alpha' + al) 
        namesb.append('Hist' + names[i])
        if fold is not False:
            names[i] = fold + '/' + names[i]
            namesb[i] = fold + '/' + namesb[i]
        names[i] += '.txt'
        namesb[i] += '.txt'
        i += 1
    for i in xrange(len(namesb)):
        Histogram = np.genfromtxt(namesb[i])
        H = Histogram[:,0]; m_num = Histogram[:,1]
        P_num = H/sum(H)
        plt.loglog(m_num,P_num,label='$\\alpha = %g$'%alphas[i])
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel('$\mathrm{Income;}\, m\, [1]$',size=15)
    plt.ylabel('$P_m\,[1]$',size=15)
    plt.xlim(0.02,100)
    plt.ylim(1E-6,.1)
    plt.show()
    return None

#-----------Model 3 alpha .5,1,1.5,2, N=500 and N=1000-------------
#write_Model3(Nt=1E3,filename='Model3N1',fold='Model3') #
#write_Model3(tot_money=1000,N=1000,filename='Model3N2',fold='Model3')    


#-------------------------------------------------------------------------------------
def write_Model4(tot_money=500,N=500,Nt=1E3,tot_time=1E5,dm=0.02,lamb=0.4,filename='noname',fold=False):
    alphas = np.array([.5,1.,1.5,2.])
    lamb=float(lamb)
    ll = str(lamb)
    ll = ll[:ll.find('.')]+ '_' + ll[ll.find('.')+1:]
    for al_ in alphas:
        al_s = '%g'%al_
        al_s = al_s[:al_s.find('.')]+ '_' + al_s[al_s.find('.')+1:]
        f_name = filename + '_' + 'alpha' + al_s + '_' + 'lam' + '_' + ll
        inst = Market(tot_money,N)
        inst.Deals(Nt,tot_time,dm,f_name,fold,model=4,lam=lamb,alpha=al_)
        del inst
    return None

def plot_Model4(filename,fold=False,lamb=0.4):
    alphas = np.array([.5,1.,1.5,2.])
    names = []
    namesb = []
    lamb=float(lamb)
    ll = str(lamb)
    ll = ll[:ll.find('.')]+ '_' + ll[ll.find('.')+1:]
    i = 0
    for alpha in alphas:
        al = '%g'%alpha
        al = al[:al.find('.')]+ '_' + al[al.find('.')+1:]
        names.append(filename + '_' + 'alpha' + al + '_' + 'lam' + '_' + ll) 
        namesb.append('Hist' + names[i])
        if fold is not False:
            names[i] = fold + '/' + names[i]
            namesb[i] = fold + '/' + namesb[i]
        names[i] += '.txt'
        namesb[i] += '.txt'
        i += 1
    del names
    for i in xrange(len(namesb)):
        Histogram = np.genfromtxt(namesb[i])
        H = Histogram[:,0]; m_num = Histogram[:,1]
        P_num = H/sum(H)
        plt.loglog(m_num,P_num,label='$\\alpha = %g$'%alphas[i])
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel('$\mathrm{Income;}\, m\, [1]$',size=15)
    plt.ylabel('$P_m\,[1]$',size=15)
    plt.xlim(2E-2,20)
    plt.ylim(1E-6,.1)
    plt.show()
    return None
    
# -- Model 4, nearest with savings----
#write_Model4(tot_money=1000,N=1000,filename='Model4N2',fold='Model4')    

#-------------------------------------------------------------------------------------
def write_Model5_6(tot_money=1000,N=1000,Nt=1E3,tot_time=1E7,dm=0.01,alpha=1.,\
                  filename='noname',fold=False,mod=5,lamb=0.0,er='Default',Iwantp=False):
    gammas = np.arange(0,4.1,1)
    alph = float(alpha)
    al = str(alph)
    al = al[:al.find('.')]+ '_' + al[al.find('.')+1:]
    for g_ in gammas:
        g_s = '%g'%g_
        g_s = g_s[:g_s.find('.')]+ '_' + g_s[g_s.find('.')+1:]
        f_name = filename + '_' + 'gamma' + g_s + '_' + 'alpha' + '_' + al
        inst = Market(tot_money,N)
        inst.Deals(Nt,tot_time,dm,f_name,fold,model=mod,lam=lamb,alpha=alph,gamma=g_,err=er,\
                   Iwantprint=Iwantp)
        del inst
    return None
    
def plot_Model5_6(filename,fold=False,alpha=1.0,Nt=1E3,lam=0.0,mod=5,i_w=2,i_c=1,Plot_C=False,\
                  Plot_H=False,w=1,double=False):
    gammas = np.arange(0,4.1,1)
    i_w = int(i_w)
    i_c = int(i_c)
    names = []
    namesb = []
    alph = float(alpha)
    al = str(alph)
    lam = str(lam)
    al = al[:al.find('.')]+ '_' + al[al.find('.')+1:]
    lam = lam[:lam.find('.')] + '_' + lam[lam.find('.')+1:]
    i = 0
    for g_ in gammas:
        g_s = '%g'%g_
        g_s = g_s[:g_s.find('.')]+ '_' + g_s[g_s.find('.')+1:]
        names.append(filename + '_' + 'gamma' + g_s + '_' + 'alpha' + '_' + al) 
        namesb.append('Hist' + names[i])
        if fold is not False:
            names[i] = fold + '/' + names[i]
            namesb[i] = fold + '/' + namesb[i]
        names[i] += '.txt'
        namesb[i] += '.txt'
        i += 1
    del names
    C = np.zeros((len(namesb),Nt))
    if Plot_H is True:
        #plt.figure(figsize=(10,6),dpi=100)
        for i in xrange(len(namesb)):
            Histogram = np.genfromtxt(namesb[i])
            H = Histogram[:,0]; m_num = Histogram[:,1]
            
            P_num = H/sum(H)
            plt.loglog(m_num[::i_w],P_num[::i_w],label='$\gamma = %g$'%gammas[i])
        plt.legend(loc='best')
        plt.grid()
        #m_num2 = np.array([0.02,1,100])
        #w_num2 = np.array([1E-6,1E-2,0.1])
        #plt.xticks(m_num2,m_num)
        #plt.yticks(w_num2,w_num)
        plt.xlabel('$\mathrm{Income;}\, m\, [1]$',size=15)
        plt.ylabel('$P_m\,[1]$',size=15)
        plt.xlim(1E-2,100)
        plt.ylim(1E-6,.05)
        plt.show()
    
    if Plot_C is True:
        markers = ['--o','--v','--d','--s','--^']
        plt.subplot(2,1,w)
        for i in xrange(len(namesb)):
            if fold is not False:
                C[i,:] = np.genfromtxt(fold+'/C_Model%i'%mod+'_'+'gamma' +str(i)+'_'+'0'+'_'+ 'lam' +\
            lam + '_' + 'alpha' + al + '.txt')
            else:
                C[i,:] = np.genfromtxt('C_Model%i'%mod+'_'+'gamma' +str(i)+'_'+'0'+'_'+ 'lam' +\
            lam + '_' + 'alpha' + al + '.txt')
        for i in xrange(len(C)):
            binn = np.arange(0,C[i].max()+1.1,1)
            CH = np.histogram(C[i],binn)[0]/float(sum(C[i]))
            maxim = binn[np.argmax(CH)]
            k = binn[:-1]
            plt.semilogy(k[::i_c],CH[::i_c],markers[i],label='$\gamma=%g\,\\ \mathrm{peak:\,} %g$'%(gammas[i],maxim))
        plt.grid()
        plt.xlim(0,12)
        plt.text(1,0.003,'$\\alpha = %g$'%alpha,size=15)
        plt.legend(loc='best')
        plt.ylabel('$D(k)\,[1]$',size=15)
        if w is 1:
            if double is False:
                plt.xlabel('$k\,[1]$',size=15)
                plt.show()
        else:
            plt.xlabel('$k\,[1]$',size=15)
        
    return None



# ----Model 5------ 
#write_Model5_6(filename='Model5',fold='Model5',mod=5,alpha=1.0) 


#write_Model5_6(filename='Model5',fold='Model5',mod=5,alpha=2.0) 


# -----For model 6, lam = 0.4----------

#write_Model5_6(filename='Model6lam0_4',fold='Model6',alpha=1.0,lamb=0.4,er=0.2,mod=6,\
#        Iwantp=False) 


#write_Model5_6(filename='Model6lam0_4',fold='Model6',alpha=2.0,lamb=0.4,er=0.2,mod=6,\
#              Iwantp=False) 


# -----The C Vals------
# Separate run for Nt=2E3, N=tot_money and tot_time=N_agents

if __name__ == '__main__':
    line = '--------------------------------------'
    ufold = 'textfiles'
    print line
    if len(sys.argv) == 1:
        print 'For specific figure; give figure number as sys.arg.'
        print 'To show all figures; sys arg = a,'
        print "all files must be present from GitHub in folder named 'textfile',"
        print "and script 'Sup2'"
        print line
        sys.exit(1)
    else:
        if len(sys.argv) != 2:
            print 'To many system arguments.'
            print line
            sys.exit(2)
        else:
            pass
    what = str(sys.argv[1])
    if what is '1':
        print 'Figure 1.'
        plot_Model1('Model1.txt',ufold)
    
    elif what is '2':
        print 'Figure 2.'
        from Sup2 import first
        first(1)
    elif what is '3':
        print 'Figure 3.'
        from Sup2 import first
        first(2)
    elif what is '4':
        print 'Figure 4.'
        print 'first: N=500'
        plot_Model3(filename='Model3N1',fold=ufold)
        print 'second: N=1000'
        plot_Model3(filename='Model3N2',fold=ufold)
    elif what is '5':
        print 'Figure 5.'
        from Sup2 import second
        second()
    elif what is '6':
        print 'Figure 6.'
        from Sup2 import third
        plot_Model4(filename='Model4N2',fold=ufold)
        third()
    elif what is '7':
        print 'Figure 7.'
        from Sup2 import fourth
        print 'First: alpha = 1.0,'
        plot_Model5_6(filename='Model5',fold=ufold,alpha=1.0,i_w=1,i_c=1,Plot_H=True)
        fourth(1)
        print 'second: alpha = 2.0'
        plot_Model5_6(filename='Model5',fold=ufold,alpha=2.0,i_w=1,i_c=1,Plot_H=True)
        fourth(2)
    elif what is '8':
        print 'Figure 8.'
        from Sup2 import fifth
        print 'First: alpha=1.0, lam=0.4,'
        plot_Model5_6(filename='Model6lam0_4',fold=ufold,alpha=1.0,lam=0.4,mod=6,i_w=2,Plot_H=True,\
           Plot_C=False)
        fifth(1)
        print 'second: alpha=2.0, lam=0.4'
        plot_Model5_6(filename='Model6lam0_4',fold=ufold,alpha=2.0,mod=6,lam=0.4,i_w=2,Plot_H=True,\
            Plot_C=False)
        fifth(2)
    elif what is '9':
        print 'Figure 9'
        plot_Model5_6(filename='Model5',fold=ufold,alpha=1.0,i_w=1,i_c=1,double=True,Plot_C=True)
        plot_Model5_6(filename='Model5',fold=ufold,alpha=2.0,i_w=1,i_c=1,w=2,double=True,Plot_C=True)
        plt.show()
    elif what is 'a':
        from Sup2 import first,second,third,fourth,fifth
        print 'Going from Fig.1-9\n'
        print 'Figure 1.'
        plot_Model1('Model1.txt',ufold)
        print 'Figure 2.'
        first(1)
        print 'Figure 3.'
        first(2)
        print 'Figure 4.'
        print 'first: N=500'
        plot_Model3(filename='Model3N1',fold=ufold)
        print 'second: N=1000'
        print 'Figure 5.'
        second()
        print 'Figure 6.'
        plot_Model4(filename='Model4N2',fold=ufold)
        third()
        print 'Figure 7.'
        print 'First: alpha = 1.0,'
        plot_Model5_6(filename='Model5',fold=ufold,alpha=1.0,i_w=1,i_c=1,Plot_H=True)
        fourth(1)
        print 'second: alpha = 2.0'
        plot_Model5_6(filename='Model5',fold=ufold,alpha=2.0,i_w=1,i_c=1,Plot_H=True)
        fourth(2)
        print 'Figure 8.'
        print 'First: alpha=1.0, lam=0.4,'
        plot_Model5_6(filename='Model6lam0_4',fold=ufold,alpha=1.0,lam=0.4,mod=6,i_w=2,Plot_H=True,\
           Plot_C=False)
        fifth(1)
        print 'second: alpha=2.0, lam=0.4'
        plot_Model5_6(filename='Model6lam0_4',fold=ufold,alpha=2.0,mod=6,lam=0.4,i_w=2,Plot_H=True,\
            Plot_C=False)
        fifth(2)
        print 'Figure 9.'
        plot_Model5_6(filename='Model5',fold=ufold,alpha=1.0,i_w=1,i_c=1,double=True,Plot_C=True)
        plot_Model5_6(filename='Model5',fold=ufold,alpha=2.0,i_w=1,i_c=1,w=2,double=True,Plot_C=True)
        plt.show()
        print 'Done.'
    print line
    