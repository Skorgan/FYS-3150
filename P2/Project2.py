import numpy as np
import matplotlib.pyplot as plt
import sys
import time
class Jacobi:
    def __init__(self,N=100,rho_min = 0, rho_max = 6.,tol=1E-12,omega_r = 0.01): 
        """
        N: Number of steps in the interval, also corresponds to dimensionality, default = 100.
        rho_min and rho:max: The coordinate interval, default: 0,6 respectably.
        tol: The tolerance for off diagonal elements is zero, default = 1E-12.
        omega_r: Strength of the potential in the case for two electrons, constant (=1) for single electron, default=0.001.
        """
        # Create intervall; rho, max iteration
        N = int(N)     
        iter_max = N*N*N    # Maximum iterations the while loop in 'Solve' is allowed
        h = (rho_max - rho_min)/float(N+1)  
        rho = np.array([rho_min + i*h for i in range(1,N+1)]) # Interval [rho_1,rho_N]  'original' interval: [rho_0,rho_N+1]
        
        # Set dummy eigenvector matrix, To contain the eigenvectors, R[:,j] = eigenvector with the eigenvalue lamb[j] 
        R = np.eye(N)  # now dot(R[:,i],R[:,j]) = delta_i,j to start, as the solution space is hermitian --> orthogonal        
        # Global variables within class
        self.R = R 
        self.rho = rho
        self.N = N
        self.iter_max = iter_max 
        self.tol = tol
        self.h = h
        self.omega_r = omega_r
#----------------------------------------------------------------------------------------     
    def Create_A(self,interacting=False): # Create the matrix for non -and interacting case in HO potential
        t_start = time.time()
        h = self.h; N = self.N; rho = self.rho
        Diagonal = 2./(h**2) + self.potential(rho,interacting)
        offDiagonal = -1./(h**2)
        A = np.diag(Diagonal)
        A[0,1] = offDiagonal
        A[-1,-2] = offDiagonal
        for i in xrange(1,N-1):
            A[i,i+1] = offDiagonal
            A[i,i-1] = offDiagonal
        t_end = time.time()
        #print "Elapsed time %f s "%(t_end-t_start)
        return A      
#---------------------------------------------------------------------------   
    def potential(self,rho,interacting=False):  # Potential used, rho is scaled
        if interacting==False:
            return rho**2
        else:
            return self.omega_r**2*rho**2 + 1./rho
#----------------------------------------------------------------------------        
    def Solve(self,A):
        #print "---------------------------------------------------------"
        print"Finding eigenvalues and eigenvectors"
        t_start = time.time()
        R = self.R
        N = self.N; tol = self.tol
        maxOffDiag,k,l = self.maximum_off_diag(A,N)
        iter_max = self.iter_max
        iteration = 0
        while (iteration <= iter_max and maxOffDiag > tol):
            maxOffDiag,k,l = self.maximum_off_diag(A,N)
            A,R = self.Jacobi_rotate(A,R,k,l)
            iteration += 1
        # Extract eigenvalues, sort lowest to highest eigenvalues and corresponding eigenvectors
        eigenvals = np.array([A[i,i] for i in range(N)])
        i_sort = eigenvals.argsort()
        eigenvals = eigenvals[i_sort]
        eigenvecs = R[:,i_sort]
        t_end = time.time()
        print "Done\nelapsed time: %f s"%(t_end-t_start)
        #print "---------------------------------------------------------"
        return eigenvals,eigenvecs,self.rho
#---------------------------------------------------------------------------    
    def maximum_off_diag(self,A,N): # Finds maximum off diagonal value, returns value and indexes
        maxi = 0
        for i in range(N):
            for j in range(i+1,N):
                if abs(A[i,j]) > maxi:
                    maxi = abs(A[i,j])
                    k = i; l = j
        return maxi, k, l
#---------------------------------------------------------------------------    
    def Jacobi_rotate(self,A,R,k,l): # rotates matrix and returns the sim trans
        # Calculate angles for min, ensures A[k,l] = A[l,k] = 0
        if A[k,l] != 0.0:
            tau = (A[l,l] - A[k,k])/(2.*A[k,l])
            if tau >= 0:
                t = -tau + np.sqrt(1. + tau**2)   
            else:
                t = -tau - np.sqrt(1. + tau**2)   
            c = 1./np.sqrt(1. + t**2)             
            s = c*t                               
        else:
            c = 1.; s = 0 
        a_kk = A[k,k]; a_ll = A[l,l] 
        A[k,k] = a_kk*c**2 - 2.*A[k,l]*c*s + a_ll*s**2 # New diag element in sim_trans matrix
        A[l,l] = a_ll*c**2 + 2.*A[k,l]*c*s + a_kk*s**2 # --""--
        A[k,l] = 0.; A[l,k] = 0.  # Set to zero in order to avoid any float errors leaving a small non zero number
        # Calculate rest of the sim_trans mat indexes
        for i in range(self.N):
            if (i != k and i != l):
                a_ik = A[i,k]; a_il = A[i,l]
                A[i,k] = a_ik*c - a_il*s
                A[k,i] = A[i,k]
                A[i,l] = a_il*c + a_ik*s
                A[l,i] = A[i,l]
            # Update the corresponding eigenvectors in the sim_trans, corresponds to change in A above
            r_ik = R[i,k]; r_il = R[i,l]
            R[i,k] = c*r_ik - s*r_il
            R[i,l] = c*r_il + s*r_ik
        return A,R
#------------------------------------------------------------------------------------------------
    def test_eigenvalues(self):
        N = 5
        R = np.eye(N)
        A = np.matrix([[2.,1.,1.,0,-1.],[1.,2.,0,1.,1.],[1.,0,2.,1.,-1.],[0,1.,1.,2.,1.],[-1.,1.,-1.,1.,2.]])
        exact = np.zeros(5)
        exact[0:2] = 0. 
        exact[2] = 2. 
        exact[3:] = 4.
        self.tol = 1E-15
        self.max_iter = 10000
        self.A = A
        self.R = R
        self.N = N
        lambs,R,rho = self.Solve(A)
        testing = abs(exact - lambs) < self.tol
        print "---------------------------------------------------------"
        if False in testing:
            number_fails = sum(testing==False)
            print "Eigenvalue test:\nFailed!\n%g/5 eigenvalues not within tolerance %g"%(number_fails,self.tol)
        else:
            print "Eigenvalue test within error eps = %g:\nSuccess!"%self.tol
        print "---------------------------------------------------------"
        return None
#---------------------------------------------------------------------------------------------   
    def test_largest_offdiag(self):
        N = 4
        A = np.matrix([[2.,1.,4.,6],[1.,9.,8.,2.],[4.,8.,2.,9.1],[6.,2.,9.1,10]])
        maxOffDiag,k,l = self.maximum_off_diag(A,N)
        maxi = 9.1
        print "---------------------------------------------------------"
        if maxi == maxOffDiag and A[k,l] == maxi:
            print "Maximum offdiag element test:\nSuccess!"
        else:
            print "Maximum offdiag element test:\nFailed!"
        print "---------------------------------------------------------"
        return None
#----------------------------------------------------------------------------------------------  
    def test_orthogonality(self):
        N = 5
        R = np.eye(N)
        self.tol = 1E-10
        A = np.matrix([[2.,1.,1.,0,-1.],[1.,2.,0,1.,1.],[1.,0,2.,1.,-1.],[0,1.,1.,2.,1.],[-1.,1.,-1.,1.,2.]])
        self.R = R
        self.N = N
        A,R,rho = self.Solve(A)
        test_mat = np.inner(R,R)  # Should be the identity matrix 
        criteria = np.allclose(test_mat,np.eye(N),atol=self.tol) 
        print "---------------------------------------------------------"
        if criteria == True:
            print "Orthogonality test:\nSuccess!"
        else:
            print "Orthogonaloty test:\nFailed!"
        print "---------------------------------------------------------"
        return None
#----------------------------------------------------------------------------------------   
    def plot_different_omega_groundstate(self,omega_r = [5,1,.5,.01],fig_name=False):
        N = self.N
        eigenvectors = np.zeros((4,N))
        lams = np.zeros(len(omega_r))
        for i in range(len(omega_r)):
            v = Jacobi(N,omega_r = omega_r[i])
            A = v.Create_A(True)
            lam,vec,rho = v.Solve(A)
            eigenvectors[i] = vec[:,0]
            lams[i] = lam[0]
        for i in range(len(omega_r)):
            plt.plot(rho,eigenvectors[i]**2,label="$\omega_r = %g$"%omega_r[i])
        plt.ylabel('$|u(\\rho)|^2,$ $[1]$',size=20)
        plt.xlabel('$\\rho,$ $[1]$',size=20)
        # Compute no interaction simeq to omega_r = 1 and omiting interaction terms
        vv = Jacobi(N)
        A = vv.Create_A()
        lam,vec,rho = vv.Solve(A)
        plt.plot(rho,vec[:,0]**2,label = 'non-interacting')
        plt.legend()
        plt.grid()
        if fig_name != False:
            plt.savefig(str(fig_name))
        plt.show()
        return None
        
#----------------------------------------------------------------------------------------------
if __name__ == '__main__':
    test = Jacobi()
    test.test_orthogonality()
    test.test_largest_offdiag()
    test.test_eigenvalues()  

# Last test: 30.09.17    
"""
Finding eigenvalues and eigenvectors
Done
elapsed time: 0.001820 s
---------------------------------------------------------
Orthogonality test:
Success!
---------------------------------------------------------
---------------------------------------------------------
Maximum offdiag element test:
Success!
---------------------------------------------------------
Finding eigenvalues and eigenvectors
Done
elapsed time: 0.002056 s
---------------------------------------------------------
Eigenvalue test within error eps = 1e-15:
Success!
---------------------------------------------------------
""" 
