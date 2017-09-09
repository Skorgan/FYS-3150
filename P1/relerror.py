import matplotlib.pyplot as plt
import sys
import glob, os
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
# Read filenames, extract values, plot and save analytical and numerical solution
if len(sys.argv) == 1:
    print "ERROR:\nInclude filename(s) to be read in order to execute program: %s"%sys.argv[0]
    sys.exit(1)

else:     
    filenames = np.array([str(name) for name in sys.argv[1:]])     # Create array containing filenames given as com.arg.
    log_eps_max = np.zeros(len(filenames))                         # To store log10(max rel error)
    log_h = np.zeros(len(filenames))                               # To store log10(h) value (step length in x)

    save_fig = True             # Boolean expression, set 'True' to save figures
    save_err_h = True           # Boolean expression, set 'True' to save a text file containing errors 
    o = 0
    # Check if it is a test run
    if filenames[0] == 'Test.txt':
        print "Extracting values from test file"
        True_Test_Values = np.genfromtxt('TrueTest.txt')
        Test_values = np.genfromtxt(filenames[0])
        x_T = True_Test_Values[:,0]; u_T = True_Test_Values[:,1]; v_T = True_Test_Values[:,2]
        x = Test_values[:,0]; u = Test_values[:,1]; v = Test_values[:,2]
        if len(x) != len(x_T):
            print "Dimensions in test file does not fit dimensions in solution file, set int n = %g"%len(x_T)
            sys.exit(2)
        else:
            print "Comparing values"
            epsilon = 1E-5
            test_x = abs(x_T - x) < epsilon 
            test_u = abs(u_T - u) < epsilon 
            test_v = abs(v_T - v) < epsilon
            if False in test_x or False in test_u or False in test_v:
                print "###FAILED!###\naccuracy not kept in the interfall eps: %g "%epsilon
                
            else:
                print "###SUCCESS!###"
            os.remove("Test.txt")
            print "Test file removed"
            sys.exit(3)
    else:
        pass
                
    for filename in filenames:
        values = np.genfromtxt(filename)
        #x = values[:,0]
        u = values[:,1]  # Extract analytic values
        v = values[:,2]  # Extract numerical values
        n = len(u)      
        log_h[o] = np.log10(1./(n+1))  # Step length
        log_eps_max[o] = np.log10(abs((v[1:-1] - u[1:-1])/u[1:-1])).max() # Find maximum relative error
        o += 1
    # Plotting figure    
    plt.plot(log_h,log_eps_max,'o')
    plt.ylabel('$\epsilon_{max}$, $[1]$')
    plt.xlabel('$\log_{10}(h)$, $[1]$')
    if save_fig == True:
        plt.savefig('relerr.png')
    plt.show()    
    if save_err_h == True:
        np.savetxt('err_h.txt', np.c_[log_h,log_eps_max])
        print "Values stored in 'err_h.txt' as columns: log_10(h), eps"

