import numpy as np

class lattice():
    """Contains functions to help out calculate matrices
       in the Fermi-Hubbard model"""
    def __init__(self, xs,ys,zs):
        self.x, self.y, self.z = np.mgrid[ 0:xs, 0:ys, 0:zs] 
        self.xs = xs
        self.ys = ys
        self.zs = zs
   
    def show(self,spins):
        for i in np.ravel(spins):
            print "%d "%i,
        print

    def state(self,m):
        # Each site can have 4 possible configurations, we have 
        # labeled them as follows:
        # 
        #  0  = vacuum
        #  1  = spin up
        #  2  = spin down
        #  3  = doubly occupied
        # 
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
        # Calculates the spin sector for the current state
        s = 0
        for i in self.spins.flat:
            if i == 0 : s = s+0
            elif i == 1 : s = s+1
            elif i == 2 : s = s-1
            elif i == 3 : s = s+0
        return s

    def filling(self):
        # Calculates the fillign for the current state
        f = 0
        for i in self.spins.flat:
            if i == 0 : f = f+0
            elif i == 1 : f = f+1
            elif i == 2 : f = f+1
            elif i == 3 : f = f+2
        return f
        
  
    def defstates(self):
        '''This function calculates the half filling states of the 
           Fermi-Hubbard model in a 3D lattice'''
        end = False
        n = 0
        self.states = {}
        while n < 300:
            self.spins = self.state(n)
            
            # The condition onf this if specifies only HALF-FILLING states 
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
       
        sites = [] 
        for i in range(self.x.size):
            sites.append( (self.x.flat[i], self.y.flat[i], self.z.flat[i], i))
        neighbors = []
        for i,s1 in enumerate(sites):
           for j,s2 in enumerate(sites): 
               if j > i: 
                   d2 = (s1[0]-s2[0])**2 + (s1[1]-s2[1])**2 + (s1[2]-s2[2])**2 
                   print s1,"--",s2," = ",d2
                   if d2 == 1: 
                       neighbors.append( (s1[3],s2[3]))
        print "Neighbor list: "
        print neighbors
        self.neighbors = neighbors
    
    def kinetic0(self):
        '''This function calculates the kinetic energy matrix
           in the spin=0 sector'''
        connected = [(0,3,1,2),\
                     (0,3,2,1),\
                     (3,0,1,2),\
                     (3,0,2,1),\
                     (1,2,0,3),\
                     (1,2,3,0),\
                     (2,1,0,3),\
                     (2,1,3,0)]
        tsign = [ 1, -1, 1, -1, 1, 1, -1, -1]
        print
        msize = len(self.states[0])
        kinetic = np.zeros((msize,msize))
        for i,s1 in enumerate(self.states[0]):
            for j,s2 in enumerate(self.states[0]):
                for n in self.neighbors:
                    # Here we have two sites with two states, we will write them
                    # in a tuple as:
                    c = (s1.flat[n[0]], s1.flat[n[1]], s2.flat[n[0]], s2.flat[n[1]])
                     
                    print "States %d,%d"%(i,j),"Neighbor pair ",n,\
                          "  --> %d,%d and %d,%d"%c,
		    if c in connected:
                        kinetic[i,j] = kinetic[i,j] - tsign[ connected.index(c) ] 
                        print " -t" 
                    else:
                        print 
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

def latex(state):
    out = r"$|"
    for j,i in enumerate(np.ravel(state)):
        if i == 0 : out+='0'
        elif i == 1 : out+=r'\!\uparrow'
        elif i == 2 : out+=r'\!\downarrow'
        elif i == 3 : out+=r'\!\uparrow\! \downarrow'
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
   
    b = lattice(2,1,1)
    b.defstates()
    b.nearest()
    b.kinetic0()
    b.interaction0()
    b.diagonal0()
    np.savetxt('211_t.dat', b.kinetic, fmt='%01d')
    np.savetxt('211_U.dat', b.inter, fmt='%01d')

    t = 1. 
    U = np.concatenate( ( np.linspace(0.1,2.,12), np.linspace(2.0,10.,8)))
    eva = []
    eve = []
    for u in U: 
        H = t*b.kinetic + u*b.inter 
        ##print H
        evals,evecs = np.linalg.eigh(H)
        ##print "U = ",u
        ##print evals 
        ##print evecs
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
        #eve.append(evecs[index])

        ##print 
        ##print evals[index]
        ##print evecs[index]
        ##print "#################"
        ##print index
        ##print 
        ##for i in index:
        ##    print "Eigenvalue  %d = "%i, evals[index[i]]
        ##    print "Eigenvector %d = "%i, evecs[:,index[i]]
        ##    print "H*ev        %d = "%i, np.dot(H, evecs[:,index[i]]) 
        ##    #print np.dot(H, evecs[index[i]]) / evecs[index[i]]
        ##    print  
        
    eva = np.array(eva)
    eve = np.array(eve)

    from matplotlib import rc 
    rc('font', **{'family':'serif'})
    rc('text', usetex=True)
    
    figure = plt.figure(figsize=(8.,4))
    gs = matplotlib.gridspec.GridSpec( 2,11) 
    figure.suptitle('')
    ax = plt.subplot( gs[0:2,0:5] )
    ax0 = plt.subplot( gs[0,5:8] )
    ax1 = plt.subplot( gs[1,5:8] )
    ax2 = plt.subplot( gs[0,8:11] )
    ax3 = plt.subplot( gs[1,8:11] )
    axvs = [ax0,ax1,ax2,ax3] 

    c=['blue','green','red','black']
    for col in range(eva.shape[1]):
        ax.plot( U, eva[:,col], '-', c=c[col],lw=2.,\
            label='%d'%col)
        for i,axv in enumerate(axvs):
             axv.plot( U, eve[:,i,col],\
                       '-',c=c[col],lw=1.5,\
            label='%d'%col)
    txts = ['a', 'b', 'c', 'd']
    for i,axv in enumerate(axvs):
        axv.set_ylabel( latex( b.states[0][i]), rotation=0 , labelpad=12)
        axv.grid()
        axv.set_ylim(-1,1.)
        axv.text( 0.05, 0.02, txts[i], ha='left', va='bottom', transform=axv.transAxes)

    statestext = [latex( b.states[0][i] )[1:-1] for i in range(4) ]
    figure.text( 0.74, 0.94, r'$\psi = a\,%s + b\,%s + c\,%s + d\,%s$'% tuple(statestext ),\
                 ha='center', fontsize=16)
    ax.grid()
    ax.set_xlabel('$U/t$', fontsize=14)
    ax.set_ylabel('$E/t$', fontsize=14)
    ax.legend(loc='best',numpoints=1,\
         prop={'size':10}, \
         handlelength=1.1,handletextpad=0.5)
    gs.tight_layout(figure, rect=[0,0.0,1.0,0.92])
    outfile = 'Ut_eigenvalues_2site.png'
    figure.savefig(outfile, dpi=250)
  
        
        
  
    
