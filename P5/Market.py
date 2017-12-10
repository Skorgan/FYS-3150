import numpy as np
import time as measure
import os
import sys
"""
Note:
    This module is more specialized than generelized to the structure of project 5 in 
    Computational physics at UIO. Some formalities in structure (and unwritten 'rules' (norms))
    are therefore disregarded. I would not recommend using this in any other way under the 
    restrictions written in the doc strings. Importing this class will however import multiple 
    modules (see above). Main method 'Deals' contain several 'models' in order to reduce FLOPs,
    all 'models' can be run with model=6 (see doc_string), but this is the most computational heavy
    and is made to produce what was needed.  
"""
class Market:
    def __init__(self,tot_money,agents): 
        """
        Initializes the market, the total number of agents and total money.
        
        Note:
            Multiple methods are an integrated part of other method(s) and will not work
            on their own i.e they are to be though of as 'private' or a part of the 
            call function. All functions can be called by name, '__'-formality is not used.  
        
        Args:
            tot_money (float,int): total money (>0) in the market, to be distributed.
            agents (float,int): total number of agents (>2) in the market.
        
        Attributes:
            N (int): number of agents.
            agents (numpy.ndarray): function; to contain the income of the agents.
            indexes (numpy.ndarray): contains integer numbers from [0,agents.size-1],
                                     will be used to draw random indexes.
            lines,shortlines,nl (str): to increase readability in the terminal.
            From args; tot_money
        
        Current methods (see details in their respective doc string):
            errors, random_state, Transaction, Transaction_and_Savings, Transaction_Nearest_Ne,
            Transaction_Nearest_Ne_Savings, Transaction_Nearest_Former, Deals, check_filename,
            save_values, check_eq, Check_Eq_ErrorUsage.
        """
        
        self.lines = "---------------------------------------------------"
        self.shortlines = "-----------------------"
        self.nl = '\n'
        if tot_money < 0:
            self.errors('tot_money')
        agents = int(agents)
        if agents is 0 or agents < 0:
            self.errors('agents')
        self.tot_money = float(tot_money)
        self.agents = np.zeros(agents)
        self.N = len(self.agents)
        self.indexes = np.arange(0,agents,1)
#-----------------------------------------------------------------------    
    def errors(self,a):
        """
        Note:
            Arguments that fail certain test will prompt this function to display the error
            message and aborting the program.     
        Args:
            a (str): cause of error
        """
        if a is 'tot_money':
            print "Error:\nWe do not force agents to start with dept at this time,",\
                "\nargument: 'tot_money' must be >= 0"
            sys.exit(1)
        elif a is 'agents':
            print "Error:\nargument 'agents' must be >= 2"
            sys.exit(2)
        elif a is 'tot_time':
            print "Error:\nArgument 'tot_time' must be int > 0"
            sys.exit(3)
        elif a is 'C':
            print 'MemoryError:\nYour chosen number of agents is too high'
            print self.lines
            sys.exit(4)
        return None
#-----------------------------------------------------------------------    
    def random_state(self):
            """
            Note:
                Used for 'model >= 3' see method; 'Deals' below.
                No arguments are required as it will alter the global attr.
            
            Attribute change:
                agents (numpy.ndarray): change from equal valued indexes agents[i], to
                randomized within [1/2,1/5]*initial_equal_values (agents[i]). 
            
            Returns:
                agn (numpy.ndarray): randomized money distribution, could set 'self.agents=agn'
                and avoid a return value.
            """
            agn = np.zeros(self.N)
            agn[:] = self.tot_money/len(agn)
            equal = agn[0]
            ints = np.random.randint(2,5,len(agn))  #1,4
            agn /= ints
            missing = equal - agn
            np.random.shuffle(missing)
            agn += missing 
            del ints,missing,equal
            return agn
#-----------------------------------------------------------------------               
    def Transaction(self):
        """ 
        Note:
            Two agents are chosen at random for potential wealth exchange.
            This model follows that of Patriarca et. al, without saving. This is for internal
            usage within the class. Money is conserved in the transaction.
        
        Attribute change:
            agents (numpy.ndarray): change index [i] and [j] for i != j.
        
        Returns:
            None-type
        """
        i,j = np.random.choice(self.indexes,2)
        if i == j:
            self.Transaction()
        else:
            eps = np.random.uniform()
            mimj = self.agents[i] + self.agents[j]
            mi = eps*mimj
            self.agents[i] = mi
            self.agents[j] = (1.-eps)*mimj     
            return None
    #--------------------        
    def Transaction_and_Savings(self,lam):
        """ 
        Note:
            Private method:
            Two agents are chosen at random for potential wealth exchange.
            This model follows that of Patriarca et. al, with saving. This is for internal
            usage within the class. Money is conserved in the transaction.
 
        Args:
            lam (float): corresponds to saving amount, see method 'Deals', typically within [0,1>
            where (lam = 0) correspond to the above method 'Transaction'.
            
        Attribute change:
            agents (numpy.ndarray): change index [i] and [j] for i != j.
        
        Returns:
            None-type
        """
        i,j = np.random.choice(self.indexes,2)
        if i == j: 
            self.Transaction_and_Savings(lam)
            return None
        else:
            eps = np.random.uniform()
            delta_m = (1.-lam)*(eps*self.agents[j] - (1.-eps)*self.agents[i])
            self.agents[i] += delta_m
            self.agents[j] -= delta_m     
            return None
    #-------------------- 
    def Transaction_Nearest_Ne(self,alpha):
        """ 
        Note:
            Two agents are chosen at random for potential wealth exchange.
            This model follows that of Patriarca et. al, without saving for
            alpha = 0, and Goswami and Sen for gamma >=1 where wealth exchange
            will also depend on similar wealth status. This is for internal
            usage within the class. Money is conserved in the transaction.
 
        Args:
            alpha (int,float): determines the exponential strength of agents wealth/income
                               difference, see method 'Deals'.
            
        Attribute change:
            agents (numpy.ndarray): change index [i] and [j] for i != j.
        
        Returns:
            None-type
        """
        np.seterr(divide='ignore')
        i,j = np.random.choice(self.indexes,2)
        o = (abs(float(self.agents[i]) - self.agents[j]))
        try:
            p = 1./float(o)**alpha
        except ZeroDivisionError:
            p = 1. 
        
        if p > np.random.uniform():
            eps = np.random.uniform()
            mimj = self.agents[i] + self.agents[j]
            mi = eps*mimj
            self.agents[i] = mi
            self.agents[j] = (1.-eps)*mimj
        else:
            self.Transaction_Nearest_Ne(alpha)
        return None
    #--------------------    
    def Transaction_Nearest_Ne_Savings(self,alpha,lam): 
        """ 
        Note:
            Two agents are chosen at random for potential wealth exchange.
            This model is a combination of Patriarca et. al, with saving
            and Goswami and Sen for gamma >=1 where wealth exchange
            will also depend on similar wealth status (combined methods: Transaction and
            Transaction_and_Savings, done in order to reduce FLOPs). 
            This is for internal usage within the class. Money is conserved in the transaction.
            
        Args:
            alpha (int,float): determines the exponential strength of agents wealth/income
                               difference, see method 'Deals'.
            lam (float): corresponds to saving amount, see method 'Deals',
                         usage for lam within <0,1>.
              
        Attribute change:
            agents (numpy.ndarray): change index [i] and [j] for i != j.
        
        Returns:
            None-type
        """
        i,j = np.random.choice(self.indexes,2)
        o = (abs(float(self.agents[i]) - self.agents[j]))
        try:
            p = 1./float(o)**alpha
        except ZeroDivisionError:
            p = 1
            
        if p > np.random.uniform():
            eps = np.random.uniform()
            delta_m = (1.-lam)*(eps*self.agents[j] - (1.-eps)*self.agents[i])
            self.agents[i] += delta_m
            self.agents[j] -= delta_m 
            return None
        else:
            self.Transaction_Nearest_Ne_Savings(alpha,lam) 
            return None    
    #-------------------- 
    def Transaction_Nearest_Former(self,alpha,gamma,lam=0.): # With or with saving
        """ 
        Note:
            Two agents are chosen at random for potential wealth exchange.
            This model is a generalisation of Patriarca et. al,
            and Goswami and Sen, For general purpose, this method replicates
            previous transaction methods, but will contain additional operations 
            which makes it slower. This is for internal usage within the class. 
            Money is conserved in the transaction.
            
        Args:
            alpha (float): determines the exponential strength of agents wealth/income
                           difference, see method 'Deals'.
            gamma (int,float): corresponds to the previous transaction weight exponent,
                               usage for gamma >= 0
            lam (float): default set to 0.0, corresponds to saving amount, see method 'Deals',
                         general usage; lam within [0,1>. lam=1 would prevent transactions. 
              
        Attribute change:
            agents (numpy.ndarray): change index [i] and [j] for i != j.
            C (numpy.ndarray): nested array with dimensions [agents.size,agents.size],
                               keep track of previous transactions between agent (i) and
                               (j) (vv.). 
        Returns:
            None-type
        """
        np.seterr(divide='ignore')
        i,j = np.random.choice(self.indexes,2)
        
        o = (abs(float(self.agents[i]) - self.agents[j]))
        try:
            p = 1./float(o)**alpha
        except ZeroDivisionError:
            p = 1
        if (self.C[i,j]**gamma)*p > np.random.uniform(): ############
            eps = np.random.uniform()
            delta_m = (1.-lam)*(eps*self.agents[j] - (1.-eps)*self.agents[i])
            self.agents[i] += delta_m
            self.agents[j] -= delta_m
            self.C[i,j] += 1
            self.C[j,i] += 1
            return None
        else:
            self.Transaction_Nearest_Former(alpha,gamma,lam)
            return None
            
#-----------------------------------------------------------------------            
    def Deals(self,Nt=1E3,tot_time=1E7,dm=0.02,filename='noname',fold=False,model=6,lam=0.,\
               alpha=0.,gamma=0.,err='Default',Iwantprint=False):
        """
        Note:
            Main method; 
            Sets up the Monte Carlo simulation for wealth exchange distribution. For general use,
            model=6 is default which in turn use the generalized transcation method; 'Transaction_Nearest_Former'.
            See description for argument 'model' below. The sorted average income distribution
            and the averaged histogram (number,bins) are stored in separate files, see method
            'store_values' and 'check_filename'. 
            
        Args:
            Nt (float,int): default set as 1E3; number of total runs (i.e how many times the approx stationary
                            state is acheived and count towards eq distribution).
            tot_time (float,int): default set as 1E7; total iterations per run, in model 5; tot_time >= N, where
                                  N is the number of agents. 
            dm (float,int): default set as 2E-2; set the bin spacing in histogram/function.
            filename (str): default set as 'noname', if filename is not given, the user will
                            be asked to give one, see method: 'check_filename'.
            fold (bool, optional; str): default set as 'False', give folder name, either existing
                                        or new, see method 'check_filename'.
            lam (float): default set to 0.0, corresponds to saving rate during transaction,
                         general use within [0,1>, lam=1 would prevent money exchange.
            alpha (int,float): determines the exponential strength of agents wealth/income
                               difference, general usage alpha >= 0
            gamma (int,float): corresponds to the previous transaction weight exponent,
                               usage for gamma >= 0.
            err (float,semi-optional): the max difference in distribution after one MC cycle for which
                                       stationary state is determined, see method 'Check_Eq_ErrorUsage'.
                                       One should provide a value, hence 'semi-optional'.                 
            Iwantprint (bool): default; 'False', if True and model >= 4, then; current run, time
                               of previous run and iterations to eq will be printed, additionally for model >=5,
                               the averaged degree of a random agent at run (i) in Nt.
                               Negative effect on simulation time.
            
            model (int): default is 5, this is for generalized wealth distribution model. Can
                         take intiger numbers 1-5, where;
                         model=1: 
                             - use method 'Transaction',
                         model=2:
                             - --"-- 'Transaction_and_Savings',
                         model=3:
                             - --"-- 'Transaction_Nearest_Ne',
                         model=4:
                             - --"-- 'Transaction_Nearest_Ne_Savings',
                         model=5 = model=6:
                             - --"-- 'Transaction_Nearest_Former'
                             difference between 5 and 6 is the name of the files, 
                             model 6 includes lam values in filename, 5 does not.
            
            NB:
                - total_time/number_of_agents should be an integer.
                
                - if arg 'err' is 'default' and total money > number of agents,
                  danger is set 'True', default 'err'-value will be scaled up
                  (in order to maximize efficiency at the cost of accuracy)
                  and user is warned. Including an 'err' value will avoid this lock.
                
                - Uses Numpy's 'histogram' module, slow for high number_of_agents/money_step
                  values (i.e; len(agents)/dm )
                
            Returns:
                None-type
       
        """         
        # --------Verify and prepare for simulation---------
        if tot_time <= 0:
            self.errors('tot_time')
        N = self.N
        if self.tot_money > N:
            print '\nNB: Model(s) work best with tot_money <=',\
'agents.\nTweeks are done in eq handling for optimalization, error scale: 10*tot_money/agents.' 
            self.danger = True
        else:
            self.danger = False
        
        self.binss = np.linspace(0,self.tot_money,self.tot_money/dm + 1)
        H = np.zeros(len(self.binss)-1)
        average_income = np.zeros(N)
       
        self.filename = filename
        self.fold = fold
        self.check_filename()
        
        tot_time = int(tot_time)
        Nt = int(Nt)
        print self.nl + self.lines
        
        if err is 'Default':
            err = self.Check_Eq_ErrorUsage(model,lam)
            
        if model == 5 or model == 6:
            try:
                self.C = np.ones((len(self.agents),len(self.agents))) 
            except MemoryError:
                self.errors('C')
        
        avg_it = 0   
        # ------------Simulations------------
        if model == 1: 
            print 'Started simulations; model 1, please wait..'
            start = measure.time()
            np.random.seed()
            for runs in xrange(Nt):
                self.agents[:] = self.tot_money/float(N) 
                agents_was = self.agents.copy()
                agents_was.sort()
                for time in xrange(tot_time):
                    self.Transaction()
                    if (time+1)%(int(2*N)) is 0: # Checks after two total MC trans
                        agents_now = self.agents.copy()
                        agents_now.sort()
                        if bool(self.check_eq(agents_was,agents_now,eps=err)) is True:
                            break
                        else:
                            agents_was = agents_now.copy()
                self.agents.sort()
                average_income += self.agents  
                H += self.Histogram(self.agents,self.binss)
        #------------           
        elif model == 2:
            print 'Started simulations; model 2, with lambda=%g, please wait..'%lam
            start = measure.time() 
            np.random.seed()
            for runs in xrange(Nt):
                self.agents[:] = self.tot_money/len(self.agents)
                agents_was = self.agents.copy()
                agents_was.sort()
                for time in xrange(tot_time):
                    self.Transaction_and_Savings(lam)
                    if (time+1)%(int(2*N)) is 0:
                        agents_now = self.agents.copy()
                        agents_now.sort()
                        if bool(self.check_eq(agents_was,agents_now,eps=err)) is True:
                            break
                        else:
                            agents_was = agents_now.copy()
                self.agents.sort()
                average_income += self.agents
                H += self.Histogram(self.agents,self.binss)
        #------------ 
        elif model == 3: # Nearest_Ne transactions
            print 'Started simulations; model 3 with alpha=%g, please wait..'%alpha
            start = measure.time()
            np.random.seed()
            for runs in xrange(Nt):
                self.agents = self.random_state()               
                agents_was = self.agents.copy()
                agents_was.sort()
                times = measure.time()
                for time in xrange(tot_time):
                    self.Transaction_Nearest_Ne(alpha)
                    if (time+1)%(int(N)) is 0:
                        agents_now = self.agents.copy()
                        agents_now.sort()
                        if bool(self.check_eq(agents_was,agents_now,eps=err)) is True:
                            break
                        else:
                            agents_was = agents_now.copy()
                self.agents.sort()
                average_income += self.agents
                avg_it += time + 1.
                H += self.Histogram(self.agents,self.binss)
                if Iwantprint is True:
                    print measure.time() - times
                    print 'Run = ',runs+1
                    print 'Iterations = ',time+1
                    print self.shortlines
            print 'Average iterations = %g'%(avg_it/float(Nt))
        #------------ 
        elif model == 4:
            print 'Started simulations; model 4 with:',\
    'alpha=%g and lambda=%g, please wait..'%(alpha,lam)
            start = measure.time()
            for runs in xrange(Nt):
                np.random.seed()
                self.agents = self.random_state()
                agents_was = self.agents.copy()
                agents_was.sort()
                for time in xrange(tot_time):
                    self.Transaction_Nearest_Ne_Savings(alpha,lam)
                    if (time+1)%(int(N)) is 0:
                        agents_now = self.agents.copy()
                        agents_now.sort()
                        if bool(self.check_eq(agents_was,agents_now,eps=err)) is True:
                            break
                        else:
                            agents_was = agents_now.copy()
                avg_it += time + 1
                self.agents.sort()
                average_income += self.agents
                H += np.histogram(self.agents,self.binss)[0]
            print 'Average iterations = %g'%(avg_it/float(Nt))
        #------------ 
        elif model == 5 or model == 6: 
            print 'Started simulations; model %g with:'%model,\
    'alpha=%g, lam=%g and gamma=%g, please wait..'%(alpha,lam,gamma)
            start = measure.time()
            C_agent = np.zeros(Nt)       
            C_i = np.random.randint(0,N) # choose a random agent to evaluate its degree
            for runs in xrange(Nt):
                np.random.seed()
                self.agents = self.random_state()
                agents_was = self.agents.copy()
                agents_was.sort()
                times = measure.time()
                MC_cycles = 0
                for time in xrange(tot_time/N):
                    for MC_step in xrange(N):
                        self.Transaction_Nearest_Former(alpha,gamma,lam)
                    agents_now = self.agents.copy()
                    agents_now.sort()
                    MC_cycles += 1.
                    if False in (abs(agents_was - agents_now) < err): 
                        agents_was = agents_now.copy()                        
                    else:
                        break
                        
                C_agent[runs] = sum(self.C[C_i,:]-1)/MC_cycles 
                if Iwantprint is True:
                    print 'C agent: ', C_agent[runs]
                    print 'time: %g s'%(measure.time() - times)
                    print 'Run = ',runs+1
                    print 'Iterations = ',(time+1)*N
                    print self.shortlines
                avg_it += (time + 1)*N    
                self.agents.sort()
                average_income += self.agents
                H += np.histogram(self.agents,self.binss)[0]
                self.C[:,:] = 1 # Remove the initial numb
            self.alpha = alpha
            self.lam = lam
            self.gamma = gamma
            self.C_agent = C_agent
            print 'Average iterations = %g'%(avg_it/float(Nt))
        #------------ 
            
        # Finishing touches
        average_income /= float(Nt)   
        self.H = H/float(Nt)
        self.average_in = average_income
        self.model = model
        self.save_values()
        print 'Done!\nTotal time = %g min'%((measure.time() - start)/60.)
        print self.lines + self.nl
        
        return None
#-----------------------------------------------------------------------
    def save_values(self):
        """
        Writes 2-3 files (value(s) stored in column(s));
             average income distribution (not wealth distribution),
             histogram/distribution values with bins in columns and
             if model >= 5, attribute 'C_agent' (average degree).      
 
        """
        print 'Storing values: average_income and histogram (latter in '\
'separate\ncolumns [val,bin] and name Hist+filename'
        avg = self.average_in
        binss = self.binss[:-1]
        H = self.H 
        np.savetxt(self.filename,np.c_[avg])
        np.savetxt(self.filename_b,np.c_[H,binss])
        if hasattr(self,'C_agent'):
            print 'Storing the average number agent (i) has interacted with\nall j in column,'
            print 'filename = C_Model%g_gamma%g_lam%g_alpha%g'%(self.model,self.gamma,self.lam,self.alpha)
            print "with '.' changed to '_'.\n" 
            fold = self.fold
            C_agent = self.C_agent
            gam = str(float(self.gamma))
            l = str(float(self.lam))
            al = str(self.alpha)
            gam = gam[:gam.find('.')] + '_' + gam[gam.find('.')+1:]
            l = l[:l.find('.')] + '_' + l[l.find('.')+1:]
            al = al[:al.find('.')] + '_' + al[al.find('.')+1:]
            fil = 'C_Model%g_'%(self.model) + 'gamma' + gam + '_'+ 'lam' + l + \
                   '_' + 'alpha' + al + '.txt'
            if fold is not False:
                if not '/' in fold[-1]:
                    fold = fold + '/'
                fil = fold + fil
            if os.path.isfile(fil):
                pass # Implement C check
            np.savetxt(fil,np.c_[C_agent])
        print 'Values stored'
        return None

#-----------------------------------------------------------------------
    def check_filename(self):
        """
        Checks filenames created from class.
        
        NB: Does not check 'C' file -> will not remove this file before a potential
            overwrite.
        """
        
        filename = self.filename
        del self.filename 
        fold = self.fold
        
        if filename is 'noname':
            filename = str(raw_input("You must provide a filename: "))
        
        if not '.txt' in filename:
                filename += '.txt'
        filename_b = 'Hist' + filename
        if fold is not False:
            if not os.path.exists(fold):
                os.makedirs(fold)
            if not '/' in fold[-1]:
                fold += '/'
            filename = fold + filename
            filename_b = fold + filename_b
        else:
            pass # Structural
        print '\nChecking path/filename: %s'%filename
        if os.path.isfile(filename):
            print "Warning: Filename (path): '%s' already exists"%filename
            ans = str(raw_input("Type 'y' to overwrite, or type a new name: " ))
            
            if ans is not 'y':
                if not '.txt' in ans:
                    ans += '.txt'
                print "Using filename: '%s'"%ans
                filename = ans # NB, not folder, add later, danger in overwriting now
                filename_b = 'Hist' + filename
            else:
                if str(ans) is 'y':
                    print "Using existing name"
                    os.remove(filename)
                    if os.path.isfile(filename_b):
                        os.remove(filename_b)
                else:
                    pass # To be rewritten
        else:
            print 'All clear!' 
        self.filename= filename
        self.filename_b = filename_b
        return None
#-----------------------------------------------------------------------        
    def check_eq(self,agents_was,agents_now,eps=.2):
        """
        Checks stability of current and prior state, both arguments are therefore
        sorted prior to call (done with Numpy's 'sort' in 'Deals',
        unfortunatelly, rather slow.). 
        """
        if False in (abs(agents_was - agents_now) < eps):
            return False
        else:
            return True
#----------------------------------------------------------------------- 
    def Check_Eq_ErrorUsage(self,model,lam):
        if model < 3:
            if lam <= 0.5:
                if self.danger is False:
                    return .3 # Set dynamic error estimate for equilibrium
                    
                else:
                    return 3.*self.tot_money/len(self.agents) # NB too large
                    
            else:
                if self.danger is False:
                    return .05
                    
                else:
                    return .5*self.tot_money/len(self.agents)                
        else:
            return 0.2
#----------------------------Some Tests---------------------------------------------
    def test_random_state(self):
        v = np.ones(10)
        vv = self.random_state()
        if False in (v-vv==0):
            return True
        else:
            return False
    
    def test_Check_Eq_ErrorUsage(self):
        model1 = 1
        model4 = 4
        sol = np.array([0.3,0.2])
        lam = 0.5
        self.danger = False
        test = np.array([self.Check_Eq_ErrorUsage(model1,lam),self.Check_Eq_ErrorUsage(model4,lam)])
        if not False in (test == sol):
            return True
        else:
            return False
#---------------------------------------------------------------------------------------------
if __name__ == '__main__':
    ins = Market(10,10)
    line = '----------------------------'
    print line
    print "Testing: 'random_state'"
    if ins.test_random_state() is True:
        print 'Success!'
    else:
        print 'Failed!'
    print line
    print "Testing: 'Check_Eq_ErrorUsage'"
    if ins.test_Check_Eq_ErrorUsage() is True:
        print 'Success!'
    else:
        print 'Failed!'
    print line
    del ins
    # No more tests included
    