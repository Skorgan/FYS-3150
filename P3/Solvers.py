import numpy as np
import sys
# Run module script 'python Solvers.py' in order to test all methods in class.
class Solvers:
    def Euler_Chromer(self,R,V,a,t,dt): # Euler Chromer/Forward Euler method.
        V_next = V + a*dt
        R_next = R + V_next*dt
        t += dt
        return R_next, V_next, t
    
    def Velocity_Verlet(self,R,V,a,t,dt,M,a_func): # NB: a_func is a function argument, must be callable, no failsafe here.
        R_next = R + V*dt + .5*a*dt**2
        a_next = a_func(R_next,M)
        V_next = V + .5*dt*(a + a_next)
        t += dt
        return R_next, V_next, a_next, t         
    
    def cross_product3D(self,vec1,vec2):    # Used to avoid Numpy's own cross method.
        a = vec1[1]*vec2[2]-vec2[1]*vec1[2]
        b = vec1[0]*vec2[2]-vec2[0]*vec1[2]
        c = vec1[0]*vec2[1]-vec2[0]*vec1[1]
        return np.array([a,b,c])
    
    def cross_product2D(self,vec1,vec2):    # ditto above.
        return vec1[0]*vec2[1]-vec2[0]*vec1[1]

#---------------------------------TESTS---------------------------------------------  
    def test_Euler_Chromer(self,check_error=False):
        line = '-------------------------------------------------'
        R = 3; V = 2; t = 0; dt = 2; a = -4
        R_test,V_test,t_test = self.Euler_Chromer(R,V,a,t,dt)
        R_c = -9.; V_c = -6.; t_c = 2      
        num = np.array([R_test,V_test,t_test])
        correct = np.array([R_c,V_c,t_c]) 
        if check_error is False:
            print line
            print "Testing: 'Euler_Chromer' in 'Solvers':"   
            if False in (num==correct):
                print 'ERROR:\nSomething is altered, test failed....'
                print line
                return False
            else:
                print 'Success!'
                print line
                return True
        else:
            print line
            print "Data type: [R,V,t]\nCorrect: [-9,-6,2]\nCalculated: [%g,%g,%g]"%(num[0],num[1],num[2])
            print line
            return None
            
    def test_Velocity_Verlet(self,check_error=False):
        line = '-------------------------------------------------'
        def f(R,M): return M*R  # ambiguous test acceleration func without care for units
        R = 2.; V = -1.; t = 1.; dt = 1.; a = 2.; M = 2.         
        R_test,V_test,a_test,t_test = self.Velocity_Verlet(R,V,a,t,dt,M,f)
        R_c = 2.; V_c = 2.;a_c=4.; t_c = 2.      
        num = np.array([R_test,V_test,a_test,t_test])
        correct = np.array([R_c,V_c,a_c,t_c]) 
        if check_error is False:
            print line
            print "Testing: 'Velocity_Verlet' in 'Solvers':"      
            if False in (num==correct):
                print 'ERROR:\nSomething is altered, test failed....'
                print line
                return False
            else:
                print 'Success!'
                print line
                return True   
        else:
            print line
            print "Data type: [R,V,a,t]\nCorrect: [2,2,4,2]\nCalculated: [%g,%g,%g,%g]"%(num[0],num[1],num[2],num[3])
            print line
            return None
                   
    def test_cross_product2D(self,check_error=False):
        line = '-------------------------------------------------'
        vec1 = np.array([2.,1.])
        vec2 = np.array([3.,4.])
        num = self.cross_product2D(vec1,vec2)
        correct = 5.
        if check_error is False:
            print line
            print "Testing: 'cross_product2D' in 'Solvers':"      
            if (num == correct) == False:
                print 'ERROR:\nSomething is altered, test failed....'
                print line
                return False
            else:
                print 'Success!'
                print line
                return True
        else:
            print line
            print "Correct: %g\nCalculated: %g"%(correct,num)
            print line
            return None
            
    def test_cross_product3D(self,check_error=False):
        line = '-------------------------------------------------'
        vec1 = np.array([2,-1,3])
        vec2 = np.array([3,0,4])
        correct = np.array([-4,-1,3])
        num = self.cross_product3D(vec1,vec2)
        if check_error is False:
            print line
            print "Testing: 'cross_product3D' in 'Solvers':"
            if False in (num == correct):
                print 'ERROR:\nSomething is altered, test failed....'
                print line
                return False
            else:
                print 'Success!'
                print line
                return True
        else:
            print line
            print """Correct: [%g,%g,%g]\nCalculated: [%g,%g,%g]"""%(correct[0],correct[1],correct[2],num[0],\
                num[1],num[2])
            print line
            return None
                  
if __name__=='__main__':
    args = len(sys.argv)
    test = Solvers()
    if args == 1:
        all_test = np.array([test.test_Velocity_Verlet(),test.test_Euler_Chromer(),\
            test.test_cross_product2D(),test.test_cross_product3D()])
        if False in all_test:
            print """Check error before using the faulty method(s) in 'Solvers'\nYou can check error by giving a sys.arg ex:
terminal> python Solvers.py 'module_name'"""
        else:
            print "Everything in 'Solvers' is working as intended,\nthe only question that remains is:\nDo you trust me?"
    if args == 2:
        methods = np.array(['cross_product2D','cross_product3D','Euler_Chromer','Velocity_Verlet'])
        arg = str(sys.argv[1])
        index = np.where(methods==arg)[0]
        if index.size == 0:
            print """You have tried to check an error in a method, 
but the name entered does not match any module, please try again.
Ex: terminal> python Solvers.py Euler_Chromer
or give the name as a string: 'Euler_Chromer' """
        if index == 0:
            test.test_cross_product2D(True)
        elif index == 1:
            test.test_cross_product3D(True)
        elif index == 2:
            test.test_Euler_Chromer(True)
        elif index == 3:
            test.test_Velocity_Verlet(True)

    # Last test: 25.10.17        
    """
-------------------------------------------------
Testing: 'Velocity_Verlet' in 'Solvers':
Success!
-------------------------------------------------
-------------------------------------------------
Testing: 'Euler_Chromer' in 'Solvers':
Success!
-------------------------------------------------
-------------------------------------------------
Testing: 'cross_product2D' in 'Solvers':
Success!
-------------------------------------------------
-------------------------------------------------
Testing: 'cross_product3D' in 'Solvers':
Success!
-------------------------------------------------
Everything in 'Solvers' is working as intended,
the only question that remains is:
Do you trust me?
    """    
           
            
    
    
      
      
