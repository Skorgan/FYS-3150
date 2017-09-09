import matplotlib.pyplot as plt
import sys
import numpy as np
import matplotlib as mpl
mpl.rcParams['axes.labelsize'] = 'x-large'

# Analytical solution:
# ----------------------------------------------------------
def u(x):
    return 1 + (np.exp(-10) - 1.)*x - np.exp(-10*x)
    
xx = np.linspace(0,1,1E5)
uu = u(xx) # Containes analytical solution for interval n=1E5 points

# ----------------------------------------------------------
# Read filenames, extract values, plot and save analytical and numerical solution:

if len(sys.argv) == 1:  # Check if command line arguments are given, if not exit program
    print "ERROR:\nInclude filename(s) to be read in order to execute program: %s"%sys.argv[0]
    sys.exit(1)

else:     
    filenames = np.array([str(name) for name in sys.argv[1:]])  
    save_fig = False      # Boolean expression, set 'True' to save figures, just
    for filename in filenames:
        values = np.genfromtxt(filename)  # get values from given filename in iteration
        
        # Extract inteval x and numerical solution v
        x = values[:,0]       
        v = values[:,2]
        n = len(x)
        
        # Plot numerical and analytical solution
        plt.figure() 
        plt.plot(xx[0::1500],uu[0::1500],'ro')
        plt.plot(x,v,'-b')
        plt.legend(["Analytical: $u(x)$","Numerical: $v(x)$, $n = %g$"%n])
        plt.xlabel('$x$, $[x] = 1$')
        plt.ylabel('$[v(x)] \wedge [u(x)]=1$')
        if save_fig == True:
            plt.savefig("P1b_n_%s.png"%str(n))
        else:
            pass
        plt.show()
    



