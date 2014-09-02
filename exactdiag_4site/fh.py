import numpy as np

class lattice():
    """Contains functions to help out calculate matrices
       in the Fermi-Hubbard model"""
    def __init__(self, xs,ys,zs):
        '''The dimensions of the grid are given to initialize the lattice.
           Recommended max of 4 sites, otherwise it can take too long to 
           complete.'''
      
        # x, y, and z have the shape of the grid, and contain the
        # respective (x,y,z) coordinates of the latttic sites:
        self.x, self.y, self.z = np.mgrid[ 0:xs, 0:ys, 0:zs] 

        self.xs = xs
        self.ys = ys
        self.zs = zs
   
    def show(self,spins):
        ''' This prints a particular state to the terminal'''
        for i in np.ravel(spins):
            print "%d "%i,
        print

    def state(self,m):
        '''
        # Each site can have 4 possible configurations, we have 
        # labeled them as follows:
        # 
        #  0  = vacuum
        #  1  = spin up
        #  2  = spin down
        #  3  = doubly occupied
        #
        #  All possible states are numbered with an index m.  This function
        #  constructs the m_th state in the lattice.  The spin configuration of
        #  the m_th state is stored in the 'spins' matrix and returned.  
        #
        #  Since there are 4 possible states per site (see above) the 
        #  convention is that m be represented in base-4 (quaternary) and 
        #  each digit can be assigned using the 0,1,2,3 convention above.
        #  
        '''
        spins = np.zeros_like( self.x)
        i = 0
        end = False
        while m > 0:
            if  i>=spins.size:
                end =True
                break 
            spins.flat[i] =  (m%4)
            m = m /4
            i = i +1
        if end:
            return None
        else:
            return spins
 
    def sector(self):
        # Finds the spin sector for the current state
        s = 0
        for i in self.spins.flat:
            if i == 0 : s = s+0
            elif i == 1 : s = s+1
            elif i == 2 : s = s-1
            elif i == 3 : s = s+0
        return s

    def filling(self):
        # Finds the filling for the current state
        f = 0
        for i in self.spins.flat:
            if i == 0 : f = f+0
            elif i == 1 : f = f+1
            elif i == 2 : f = f+1
            elif i == 3 : f = f+2
        return f
        
  
    def defstates(self):
        '''This function defines the half filling states of the 
           Fermi-Hubbard model in a 3D lattice. 

           It creates a dictionary where the keys correspond to the 
           different spin sectors available, and the values are a list
           of the states in the spin sector. 
 
           For a balanced spin mixture one only needs to consider the
           spin=0 sector. 
        '''
        end = False
        n = 0
        self.states = {}
        while n < 300:
            self.spins = self.state(n)
            
            # The condition on this specifies only HALF-FILLING states 
            if self.spins is not None and self.filling() == self.spins.size:
                sec = self.sector()
                if sec in self.states.keys():
                    self.states[ sec].append(self.spins) 
                else:
                    self.states[ sec]=[self.spins] 
            n = n+1
        for k in self.states.keys():
            print "Sector %d, %d states:"%(k,len(self.states[k]))
            for spins in self.states[k]:
                self.show(spins)

    def nearest(self):
        '''This function makes a list of the nearest neighbor 
           pairs in the lattice'''
        print "\nNearest neighbors:"
       
        # First we create a flat list of all the lattice sites. 
        # each element in the list is (x[i], y[i], z[i], i) 
        sites = [] 
        for i in range(self.x.size):
            sites.append( (self.x.flat[i], self.y.flat[i], self.z.flat[i], i))
 
        # We do a nested iteration over the lists and create a list
        # of pairs which are nearest neighbors.
        neighbors = []
        for i,s1 in enumerate(sites):
           for j,s2 in enumerate(sites): 
               if j > i: 
                   d2 = (s1[0]-s2[0])**2 + (s1[1]-s2[1])**2 + (s1[2]-s2[2])**2 
                   print s1,"--",s2," = ",d2
                   if d2 == 1: 
                       neighbors.append( (s1[3],s2[3]))
        print 
        print "Final neighbor list: "
        print neighbors
        self.neighbors = neighbors
    

    def kinetic0(self):
        r'''This function calculates the kinetic energy matrix
           in the spin=0 sector. 

        The matrix is constructed by iterating over the nearest neighbors.
        As a reminder, the kinertic enrgy is given by 

         K = -t \sum_{\langle i j \rangle} a_{i\sigma}^{\dagger} a_{j\sigma}

        So in order to find it's matrix elements we need to apply first an 
        annihilation operator and then a creation operator.  The tricky part 
        is keeping track of the signs.   
 
        '''
        print
        msize = len(self.states[0])
        kinetic = np.zeros((msize,msize))

        for i,s1 in enumerate(self.states[0]):
            for j,s2 in enumerate(self.states[0]):

                # We will calculate the matrix element 
                #  < s1 | K | s2 >
                # This matrix element involves a sum over nearest neighbors
                # and sum over spins, so we go ahead and iterate: 

                t = 0.
                for n in self.neighbors:
                    PRINT = False
                    for spin in ['up','down']:
                        if PRINT:
                            print 
                            print "<", np.ravel(s1)," | K | ", np.ravel(s2),">"

               
                        # Annihilates 'spin' at site n[0] 
                        signA, stateA = annihilate( n[0], spin, s2) 
                        # Create 'spin' at site n[1]   
                        signC, stateC = create(n[1], spin, stateA)
                        if PRINT: 
                            print "annihilate %d,%5s"%(n[0],spin)," -->",stateA
                            print "    create %d,%5s"%(n[1],spin)," -->",stateC

                        #  If  K|s2> has  a projecton on <s1| then we add it to
                        #  t 
                        if np.array_equal(stateC,np.ravel(s1)):
                            if PRINT: print " tmatrix --> % d" % (signA*signC )
                            t+= signA*signC

                        r'''
                        Notice that sometimes people write the kinetic energy as 
                         
                           K = -t \sum_{\langle i j \rangle} 
                              a_{i\sigma}^{\dagger} a_{j\sigma}  + c.c.
                        
                        where the letters c.c. refer to the complex conjugate.  
                        If they do that, then it means that the sum over nearest 
                        neighbors must only occur for one ordering of the 
                        neighbor pair,  for instance just 1-2 whereas the sum 
                        over both orderings includes 1-2 and 2-1.  
                        
                        Here we just run the sum over both orderings.  
                        '''
                        #  We repeat the process with the different neighbor 
                        #  ordering:
                        signA, stateA = annihilate( n[1], spin, s2) 
                        signC, stateC = create(n[0], spin, stateA)
                        if PRINT: 
                            print "annihilate %d,%5s"%(n[1],spin)," -->",stateA
                            print "    create %d,%5s"%(n[0],spin)," -->",stateC
                        if np.array_equal(stateC,np.ravel(s1)):
                            if PRINT: print " tmatrix --> % d" % (signA*signC )
                            t+= signA*signC
 
                kinetic[i,j] = t
 
        print "\nKinetic energy matrix: ",kinetic.shape
        print kinetic 
        self.kinetic = kinetic
 
    def interaction0(self):
        '''This fuction calculates the interaction energy matrix
           in the spin=0 sector'''
        print
        msize = len(self.states[0])
        inter = np.zeros((msize,msize))
        # The basis we have chose is of number states,
        # so the interaction energy is diagonal
        for i,s1 in enumerate(self.states[0]):
            for site in s1.flat:
                if site == 3: # 3=double occupancy
                    inter[i,i] = inter[i,i] + 1
        print "\nInteraction energy matrix:i ",inter.shape
        print inter
        self.inter = inter

    def diagonal0(self):
        '''This fuction calculates a diagonal matrix
           in the spin=0 sector'''
        print
        msize = len(self.states[0])
        diag = np.zeros((msize,msize))
        # The basis we have chose is of number states,
        # so the interaction energy is diagonal
        for i,s1 in enumerate(self.states[0]):
            for site in s1.flat:
                diag[i,i] = 1.0
        self.diag = diag

def annihilate( i, spin, state):
    # The order for the creation operators is lower site number
    # to the left,  and then spin-up to the left
    s = np.ravel(state)
    out = np.copy(s)
    samespin = {'up':1, 'down':2} 
    flipspin = {'up':2, 'down':1} 
    
    ncommute = 0. 
    for j in range(i):
        if s[j] == 3: ncommute +=2
        if s[j] == 1 or s[j] == 2: ncommute+=1 
    sign = (-1)**ncommute
    
    if s[i] == 0: 
        out = np.zeros_like(s) 

    if s[i] == flipspin[spin]:
        out = np.zeros_like(s) 

    if s[i] == 3: 
        out[i] = flipspin[spin] 
        if spin == 'up': sign*= 1
        if spin == 'down': sign*=-1

    if s[i] == samespin[spin]:
        out[i] = 0 

    #print s, ", annihilate %d,%5s"%(i,spin)," --> %+d"%sign, out
    return sign, out
    
def create( i, spin, state):
    # The order for the creation operators is lower site number
    # to the left,  and then spin-up to the left
    s = np.ravel(state)
    out = np.copy(s)
    samespin = {'up':1, 'down':2} 
    flipspin = {'up':2, 'down':1} 
    
    ncommute = 0. 
    for j in range(i):
        if s[j] == 3: ncommute +=2
        if s[j] == 1 or s[j] == 2: ncommute+=1 
    sign = (-1)**ncommute
    
    if s[i] == 0: 
        out[i] = samespin[spin]

    if s[i] == flipspin[spin]:
        out[i] = 3 
        if spin == 'up': sign*=1
        if spin == 'down': sign*=-1

    if s[i] == 3: 
        out = np.zeros_like(s) 

    if s[i] == samespin[spin]:
        out = np.zeros_like(s) 

    #print s, ", create %d,%5s"%(i,spin)," --> %+d"%sign, out
    return sign, out

def puretext(state):
    out = r'|'
    for j,i in enumerate(np.ravel(state)):
        if i == 0 : out+='0'
        elif i == 1 : out+=r'1'
        elif i == 2 : out+=r'2'
        elif i == 3 : out+=r'D'
        if  j+1< state.size:
            out+=','
        else:
            out+= r'>'
    return out

def latex(state):
    MATRIX = True
    if MATRIX:
        # one of the dimensions needs to be 1 to do matrix output
        dims = [] 
        idx =  []
        one = None
        for ss,s in enumerate(state.shape):
            if s > 1 : 
                dims.append( s ) 
                idx.append( ss ) 
            if s == 1 : one = ss
        assert( len(dims) == 2 ) 
        assert( one is not None ) 

        just=''
        for i in range( dims[0] ) :
            just = just+ 'c'
            if i < dims[0]-1:
                just = just +'|' 
        #print state        
        out = r"$ \begin{array}{"+just+"} "
        for mm in range(dims[0]):
            for nn in range(dims[1]):
                tup = np.empty_like( state.shape ) 
                tup[one] = 0 
                tup[idx[0]] = mm 
                tup[idx[1]] = nn
                i  = state[ tuple( tup.tolist() ) ] 
                
                if i == 0 : out+='0'
                elif i == 1 : out+=r'\uparrow'
                elif i == 2 : out+=r'\downarrow'
                elif i == 3 : out+=r'\uparrow\! \downarrow'
                if nn < dims[1] - 1 :
                    out += ' & ' 
                else:
                    out += r' \\ '
            if mm <  dims[0] - 1 : 
                out += r'\hline'
        out += r"\end{array}$" 
        return out  
             
                
    else:
        out = r"$|"
        for j,i in enumerate(np.ravel(state)):
            if i == 0 : out+='0'
            elif i == 1 : out+=r'\uparrow'
            elif i == 2 : out+=r'\downarrow'
            elif i == 3 : out+=r'\uparrow\! \downarrow'
            if  j+1< state.size:
                out+=','
            else:
                out+= r'\rangle'
        out+=r'$'
        return out

        
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif'})
      
      
if __name__=="__main__":
    #a = lattice(2,2,1)
    #a.defstates()
    #a.nearest()
    #a.kinetic0()
    #a.interaction0()
    #np.savetxt('221_t.dat', a.kinetic)
    #np.savetxt('221_U.dat', a.inter)
    

    SITES = 4
    if SITES == 4 : 
        b = lattice(1,2,2)
        b.defstates()
        b.nearest()
        b.kinetic0()
        b.interaction0()
        b.diagonal0()
        np.savetxt('221_t.dat', b.kinetic, fmt='% 01d')
        np.savetxt('221_U.dat', b.inter, fmt='% 01d')
        outfile = 'Ut_eigenvalues_4site.png'
 
    elif SITES == 2:    
        b = lattice(2,1,1)
        b.defstates()
        b.nearest()
        b.kinetic0()
        b.interaction0()
        b.diagonal0()
        np.savetxt('211_t.dat', b.kinetic, fmt='%01d')
        np.savetxt('211_U.dat', b.inter, fmt='%01d')
        outfile = 'Ut_eigenvalues_2site.png'


    # SOLUTION IS CALCULATED FOR A SET OF U VALUES
    t = 1. 
    U = np.linspace(0.1,18.,32)
    eva = []
    eve = []
    for u in U: 
        H = t*b.kinetic + u*b.inter 
        ##print H
        evals,evecs = np.linalg.eigh(H)
        ##print "U = ",u
        ##print evals 
        ##print evecs
        SORT = True
        if SORT:
            # Sort the eigenvals and eigenvecs 
            index = np.argsort(evals) 
            eva.append(evals[index]) 
            # Ensure the eigenvecs have correct phase
            vecs=[]
            for i in index:
                vec = evecs[:,index[i]] 
                #Find first entry that is non-zero
                i = list(np.abs(vec) > 1e-5).index(True)
                vec = vec / np.sign(vec[i])
                vecs.append(vec) 
            vecs =  np.transpose( np.array(vecs) )
            eve.append(vecs)
        else:
            eva.append(evals) 
            eve.append(evecs)

        if False:
            print 
            print evals[index]
            print evecs[index]
            print "#################"
            print index
            print 
            for i in index:
                print "Eigenvalue  %d = "%i, evals[index[i]]
                print "Eigenvector %d = "%i, evecs[:,index[i]]
                print "H*ev        %d = "%i, np.dot(H, evecs[:,index[i]]) 
                #print np.dot(H, evecs[index[i]]) / evecs[index[i]]
                print  
        
    eva = np.array(eva)
    eve = np.array(eve)
    print "Eigenvalues", eva.shape
    print "Eigenvectors", eve.shape



    # SOLUTIONS ARE PLOTTED   
    # Start matplotlib
    from matplotlib import rc
    rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = [
           r'\usepackage{bm}',        # for bold math
    ]  
    plt.rcParams['axes.linewidth'] = 0.6
    plt.rcParams['patch.linewidth'] = 0.4 

    nstates = len(b.states[0])
    print "Number of States in Sector 0 = ", nstates
    #This number should be a square:
    if np.abs( np.sqrt(nstates) % 1. ) > 1e-4:
        "Error, number of states in Sector 0 is not a square."

    SQUARE = False
    if SQUARE:
        plotrows = int(np.sqrt(nstates))
        plotcols = plotrows
    else:
        plotrows = 4 
        plotcols = 9
    
    figure = plt.figure(figsize=(2.6*plotrows,2.0*plotrows))
    print "Making %d x %d figure" % (plotrows, plotcols)

    gs0 = matplotlib.gridspec.GridSpec( 1,1, left=0.3, right=0.7,\
               bottom=0.62, top=0.98) 

    gs = matplotlib.gridspec.GridSpec( plotrows, plotcols, \
               left=0.03, right=0.98, bottom=0.05, top=0.55, \
                wspace=0.14, hspace=0.05) 
    figure.suptitle('')
    ax = plt.subplot( gs0[0] ) 
    #ax = plt.subplot( gs[0:plotrows,0:plotcols0] ) 
    axvs = []
    for i in range(plotrows):
        for j in range(plotcols):
            axvs.append( plt.subplot( gs[i,j]))

    # Find indices for the ground state, and other relevant states
    ground = 0 
    high = nstates-1 
    if SITES == 4:
        important = [ground, high, 6] 
    if SITES == 2:
        important = [ground, high] 

    #  Find if there is a state with energy U 
    Uindex = -1 
    for nn in range(nstates):
        if np.abs( eva[Uindex,nn] - U[Uindex] ) < 1e-4:
            #important.append(nn)
            break
    print "Importaant states = ", important 

    cc = 0 
    c=['blue','green','red','black','purple','limegreen','orange','brown']
    for col in range(eva.shape[1]):
        labeltxt = '%d'%col
        if col in important:
            color = c[cc % len(c)]
            ax.plot( U, eva[:,col], '-', c=color,lw=1.5,\
                label=labeltxt)

            for i,axv in enumerate(axvs):
                 if i >= len(eve[0,:,0]):
                     continue
                 if col == 6:
                     subset = U > 4 
                     axv.plot( U[subset], eve[:,i,col][subset],\
                           '-',c=color,lw=1.1,alpha=1.0)

                 else:
                     axv.plot( U, eve[:,i,col],\
                           '-',c=color,lw=1.1,alpha=1.0 ) 
            cc = cc + 1 
        else: 
            ax.plot( U, eva[:,col], '-', c='0.5',lw=0.8, alpha=0.4)

  
    # Print out the ground state for various Us
    Uindex=U.size-1
    Eindex =0
    print
    print "Ground state U=",U[Uindex],"  E=",eva[Uindex,Eindex], ":"

    # Organize the basis states by the magnitude of their projection
    # onto the ground state 
    order =  np.argsort(np.abs(eve[Uindex,:,Eindex]))[::-1]
    for i in order: 
        print "%02d --> % 02.6f   %s" %(i,eve[Uindex,i,Eindex], \
                                        puretext(b.states[0][i]))
        

    frame_coding = { \
        14: 'blue',\
        21: 'blue',\

        34: 'green',\
        25: 'green',\
        32: 'green',\
        19: 'green',\
        16: 'green',\
        10: 'green',\
        01: 'green',\
        03: 'green',\
 
         8: 'red',\
        27: 'red',\
       }

    for i,axv in enumerate(axvs):
        if i in frame_coding.keys():
            if False:
                for spine in axv.spines.values():
                    spine.set_edgecolor( frame_coding[i] )

            axv.text( 0.26,0.20,latex( b.states[0][i]), rotation=0 ,\
                      ha='center',va='center', fontsize=6,\
                      color = frame_coding[i],\
                      bbox=dict(facecolor='white', lw=0., pad=1.),\
                      transform=axv.transAxes)
        else:
            axv.text( 0.26,0.20,latex( b.states[0][i]), rotation=0 ,\
                      ha='center',va='center', fontsize=6,\
                      bbox=dict(facecolor='white', lw=0., pad=1.),\
                      transform=axv.transAxes)
          
        axv.yaxis.grid(which='both', alpha=0.3)
        axv.xaxis.grid(which='major', alpha=0.3)
        axv.set_ylim(-1.1,1.1)
        axv.set_xlim(0., 18.2)
        axv.xaxis.set_major_locator( matplotlib.ticker.MultipleLocator(6.) ) 
        axv.xaxis.set_minor_locator( matplotlib.ticker.MultipleLocator(3.) ) 
        axv.yaxis.set_major_locator( matplotlib.ticker.MultipleLocator(1.) ) 
        axv.yaxis.set_minor_locator( matplotlib.ticker.MultipleLocator(0.5) )
        axv.tick_params(axis='both', which='major', labelsize=9.,  length=1.5) 
        axv.tick_params(axis='both', which='minor', labelsize=9.,  length=1.0) 
        if i // plotcols < plotrows-1:
            axv.xaxis.set_ticklabels([]) 
        if i % plotcols > 0 :
            axv.yaxis.set_ticklabels([]) 
      


    #ax.grid()
    ax.set_xlabel('$U/t$')
    ax.set_ylabel('$E/t$')
    #ax.legend(loc='best',numpoints=1,ncol=int(nstates)//8,\
    #     prop={'size':10}, \
    #     handlelength=1.2,handletextpad=0.5)
    #gs.tight_layout(  figure, rect=[0.0, 0.0, 1.00, 0.7])
    figure.savefig(outfile, dpi=250)
  
        
        
  
    
