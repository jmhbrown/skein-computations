import logging


logging.basicConfig(level='INFO') # options: DEBUG, INFO, WARN, ERROR, CRITICAL
logger = logging.getLogger(__name__)


def W_rep(sigma): 
    """
    Makes a matrix representation of the Weyl group on the basis lattice.
    
    sage: W = SymmetricGroup(3)
    sage: sigma = W([1,3,2])

    sage : W_rep(sigma)

    Parameters:
        sigma SymmetricGroup.element -- Weyl group element, acts on the basis lattice.

    Returns:
        matrix -- represents the action of sigma on the basis lattice.
    """

    n = sigma.parent().degree()
    index = Matrix.identity(n).rows()

    A = matrix.block_diagonal(Matrix(sigma(index)), Matrix(sigma(index)))
    return A

def MCG_rep(gamma, n):
    """
    Makes a matrix representation of the action of the MCG on the basis lattice.

    sage: gamma = SL2Z.0*SL2Z.1^4

    sage: MCG_rep(SL2Z.0*SL2Z.1^4,3)
    [ 0  0  0|-1  0  0]
    [ 0  0  0| 0 -1  0]
    [ 0  0  0| 0  0 -1]
    [--------+--------]
    [ 1  0  0| 4  0  0]
    [ 0  1  0| 0  4  0]
    [ 0  0  1| 0  0  4]


    Parameters:
        gamma SL2Z.element_class -- MCG element
        n Int -- rank of G, on which our skein theory is defined. Half the rank of the lattice.

    Returns:
        matrix -- represents the action of gamma on the basis lattice
    """
    
    return gamma.matrix().tensor_product(Matrix.identity(n))


def get_cokernel(gamma,sigma,take_coinvariants=True): #OUTDATED
    """
    sage: gamma = SL2Z.0
    sage: sigma = SymmetricGroup(5).an_element()
    sage: sigma = SymmetricGroup(5)((1,3,4))
    sage: get_cokernel(gamma,sigma)
    Finitely generated module V/W over Integer Ring with invariants (2, 2)
    """
    n = sigma.parent().degree()
    omega = matrix.block([[0,-Matrix.identity(n)],[Matrix.identity(n),0]])

    lattice = FreeModule(ZZ,2*n)

    M = (Matrix.identity(2*n) - W_rep(sigma)*MCG_rep(gamma,n))
    U = M.right_kernel()
    List = U.gens()
    if List == ():
        List=lattice.zero()
    
    omega_on_U = omega*Matrix(List).transpose()
    
    Uperp = omega_on_U.left_kernel()
    
    sublattice_gens = M.columns()
    
    if take_coinvariants:
        for tau in SymmetricGroup(n).centralizer(sigma):
            sublattice_gens += [(Matrix.identity(2*n) - W_rep(tau))*v for v in Uperp.gens()]

    submodule = lattice.submodule(sublattice_gens)
                 
    return Uperp.quotient(submodule)


def try_all_conjugacy_classes(gamma,n,take_coinvariants=True): #OUTDATED
    """
    sage: gamma = SL2Z.1^4
    sage: try_all_conjugacy_classes(gamma,5)

    """
    conjugacy_classes = [cc.an_element() for cc in SymmetricGroup(n).conjugacy_classes()]

    cokernels = [
            {'partition':s, 'homology': get_cokernel(gamma,s,take_coinvariants)}
            for s in conjugacy_classes
            ]
    for co in cokernels:
        logger.debug("{0} - {1}".format(co['partition'],co['homology']))

    # compute dimensions
    for co in cokernels:
        co['dim'] = co['homology'].cardinality()
        co['torsion_dim'] = product(co['homology'].invariants())
    return cokernels

def dim_skeinmod(gamma,n,take_coinvariants=True): #OUTDATED
    cokernels = try_all_conjugacy_classes(gamma,n,take_coinvariants)
    sum = 0
    for co in cokernels:
        sum += co['torsion_dim']
    return sum

def change_ring(module,R=QQ):
    """
    change_ring isn't implemented, so this does it by hand.
    
    Parameters:
        module Finitely generated module V/W - we want to run module.change_ring(R) on this
        R - want want to change the ring to this

    Returns:
        module over R.
    """

    return module.V().change_ring(R).quotient(module.W().change_ring(R))

def Monica_table(L,n,take_coinvariants=True): #OUTDATED
    """
    sage: gamma = SL2Z.1^4
    sage: try_all_conjugacy_classes(gamma,5,take_coinvariants)

    """
    conjugacy_classes = [cc.an_element() for cc in SymmetricGroup(n).conjugacy_classes()]
    
    return table([[conjugacy_classes[0]] + conjugacy_classes] + [[l] + [get_cokernel(l,sigma,take_coinvariants).invariants() for sigma in conjugacy_classes] for l in L])


def Blurange_table(L,n=8,take_coinvariants=True): #OUTDATED
    """
    Parameters:
        L a vector of gammas in SL2
        n: table for GL_i from 1 to n
    Returns:
        table nxdim(L) of the dimension of the gamma-twisted skein modules in GL_i
    """
    return table([[l] + [dim_skeinmod(l,i,take_coinvariants) for i in IntegerRange(1,n+1)] for l in L], header_row=[i for i in IntegerRange(0,n+1)], frame = True)


def get_dim_empty_skein ( gamma ):
    I = matrix (ZZ , 2 , [1 , 0 , 0 , 1])
    D_plus , U_plus , V_plus = (I - gamma ). smith_form ()
    a_plus = [a for a in D_plus . diagonal () if a != 0]
    D_minus , U_minus , V_minus = (I + gamma ). smith_form ()
    a_minus = [a for a in D_minus. diagonal() if a != 0]
    p_plus = len ([ a for a in a_plus if a %2 == 0])
    p_minus = len ([ a for a in a_minus if a %2 == 0])
    return (prod(a_plus)+2**(p_plus))/2+(prod(a_minus) + 2**p_minus)/2

def Patrick_dim_GL2_skein_table(L):
    """
    Parameters:
        L a vector of gammas in SL2
    Returns:
        table dim(L) of the dimension of the gamma-twisted skein modules in GL_2 according to Patrick's paper
    """
    return table([[l] + [get_dim_empty_skein(l)] for l in L], frame = True)

def Blue_table(L,n,take_coinvariants=True): #OUTDATED
    """
    Exactly same as Monica_table + extra column containing |trace(gamma-id)|

    sage: gamma = SL2Z.1^4
    sage: try_all_conjugacy_classes(gamma,5,take_coinvariants)

    """
    conjugacy_classes = [cc.an_element() for cc in SymmetricGroup(n).conjugacy_classes()]
    
    return table([[conjugacy_classes[0]] + conjugacy_classes+ ["|trace(gamma-id)|"]] + [[l] + [get_cokernel(l,sigma,take_coinvariants).invariants() for sigma in conjugacy_classes]+ [abs(l[0,0]+l[1,1]-2)] for l in L])

def get_cokernel_reduced(gamma,n):    
    #get cokernel for a cycle of order n, using MOnica's Smith reduction.
    lattice = FreeModule(ZZ,2)
    omega = Matrix(ZZ,[[0,-1],[1,0]])
    M = (Matrix.identity(2) - gamma^n)
    U = M.right_kernel()
    List = U.gens()
    if List == ():
        List=lattice.zero()
    omega_on_U = omega*Matrix(List).transpose()
    Uperp = omega_on_U.left_kernel()
    sublattice_gens = M.columns()
    submodule = lattice.submodule(sublattice_gens)
    coker = Uperp.quotient(submodule) #up until here everything same as before
    elts = list(coker)                #putting all elements of coker in a list
    counter = 0                       #setting counter to zero
    while elts != []:                 #we're gradually removing orbits from coker
        counter = counter+1           #if we haven't removed everything yet, there's a new orbit to consider so +1
        c = elts[0]                   #take first (remaining) element in coker
        orbit = [ ((gamma.matrix())^i)*(c.lift()) for i in range(0,n)] #calculate orbit in Uperp, so first lift it
        orbit_in_coker = [coker(i) for i in orbit] #project down to coker, (coker(i) projects from Uperp to quotient)
        elts = [x for x in elts if x not in orbit_in_coker]#remove all elements in the orbit
        #we cannot remove elements that are already gone, e.g. for zero element, orbit is a list of repeated entries of the zero element
    return counter

def get_cokernel_v2(gamma,sigma):
    
    n = sigma.parent().degree()
    omega = matrix.block([[0,-Matrix.identity(n)],[Matrix.identity(n),0]])
    
    lattice = FreeModule(ZZ,2*n)
    
    M = (Matrix.identity(2*n) - W_rep(sigma)*MCG_rep(gamma,n))
    U = M.right_kernel()
    List = U.gens()
    if List == ():
        List=lattice.zero()
    
    omega_on_U = omega*Matrix(List).transpose()
    
    Uperp = omega_on_U.left_kernel()
    
    sublattice_gens = M.columns()
    submodule = lattice.submodule(sublattice_gens)
    
    coker = Uperp.quotient(submodule) #up until here everything same as before
    
    elts = list(coker)                #putting all elements of coker in a list
    counter = 0                       #setting counter to zero
    while elts != []:                 #we're gradually removing orbits from coker
        counter = counter+1           #if we haven't removed everything yet, there's a new orbit to consider so +1
        c = elts[0]                   #take first (remaining) element in coker
        orbit = [ W_rep(s)*(c.lift()) for s in SymmetricGroup(n).centralizer(sigma)] #calculate orbit in Uperp, so first lift it
        orbit_in_coker = [coker(i) for i in orbit] #project down to coker, (coker(i) projects from Uperp to quotient)
        elts = [x for x in elts if x not in orbit_in_coker]#remove all elements in the orbit
        #we cannot remove elements that are already gone, e.g. for zero element, orbit is a list of repeated entries of the zero element
                
    
    return counter

def Monica_table_v2(L,n): #Monica table with coker with coinvar calculated via brute forcing, counting orbits
    """
    sage: gamma = SL2Z.1^4
    sage: try_all_conjugacy_classes(gamma,5,take_coinvariants)

    """
    conjugacy_classes = [cc.an_element() for cc in SymmetricGroup(n).conjugacy_classes()]
    
    return table([[conjugacy_classes[0]] + conjugacy_classes] + [[l] + [get_cokernel_v2(l,sigma) for sigma in conjugacy_classes] for l in L])

def Monica_conjecture(gamma,m):
    #Monica's conjecture for dim HH_0 with take coinvar True for a single m-cycle. IMPORTANT: only works if gamma^m-I is invertible
    result = 0
    I = identity_matrix(2)
    for d in divisors(m):
        result = result + euler_phi(m/d)*abs(det(gamma**d-I))
    result= result/m
    return result

def skein_mod_dim(gamma,n):
    #gives  GLn skein module dimension, using Monica's conjecture #only works for regular gamma!!!
    single_cycle = [0]*(n+1)                         #calculate result for single m-cycle (more efficient to do it once in the beginning)
    for i in range(1,n+1):                           #m ranging from 1 till n
        single_cycle[i] = Monica_conjecture(gamma,i) #using Monica's conjecture
        
    result_per_conjugacyclass =[]                    #calculate HH_0 with coinvar for each conjugacyclass, to later sum them
    conjugacy_classes = [cc.an_element() for cc in SymmetricGroup(n).conjugacy_classes()] #list all conjugacyclasses
    for cc in conjugacy_classes:
        cycle_type = list(cc.cycle_type())           #cycle type will give list of cycle lengths appearing. e.g. (12)(34)(5) in S5 will give [2,2,1]
        count = [0]*(n+1)                            #count is a list where count[i] gives the number of cycles of length i
        for i in range(1,n+1):
            count[i] = cycle_type.count(i)
        result= 1                                    #multiply all contributions of different cycle lengths together
        for i in range(1,n+1):
                result= result*binomial(single_cycle[i] + count[i]-1,count[i]) #using the stars and bar formula if there are multiple cycles of the same length
        result_per_conjugacyclass.append(result)     #add result to list
    endresult = sum(result_per_conjugacyclass)       #sum each contribution HH_0 with coinvar together to get end result
    return endresult

def get_cokernel_v3(gamma,sigma): #calculates HH_0 with coinvar for given gamma and sigma. Note this one is not exactly used in skein_mod_dim for effciency reasons. But it does exactly the same thing as what is happening there. 
    n = sigma.parent().degree()
    cycle_type = list(sigma.cycle_type())
    count = [0]*(n+1)
    for i in range(1,n+1):
        count[i] = cycle_type.count(i)
    result= 1
    for i in range(1,n+1):
        if count[i] != 0:
            result= result*binomial(Monica_conjecture(gamma,i) + count[i]-1,count[i])
    return result

def Monica_table_v3(L,n): #Monica table with coker with coinvar calculated via Monica's conjecture
    """
    sage: gamma = SL2Z.1^4
    sage: try_all_conjugacy_classes(gamma,5,take_coinvariants)

    """
    conjugacy_classes = [cc.an_element() for cc in SymmetricGroup(n).conjugacy_classes()]
    
    return table([[conjugacy_classes[0]] + conjugacy_classes] + [[l] + [get_cokernel_v3(l,sigma) for sigma in conjugacy_classes] for l in L])

def gen_fct(gamma,n): 
    #generating function up to and including t^n, brute forcing
    Q.<t> = LazyPowerSeriesRing(SR)
    f=1
    for i in range(1,n+1):
        f= f + skein_mod_dim(gamma,i)*t**i
    return f

def solve_for_ansatz_v3(gamma, N): # input N is the highest power of t we want to consider
    c = var(['c{}'.format(i-1) for i in range(1, N+2)]) # create the variables c_i (such that c_i = c[i])
    h = gen_fct(gamma,N) # our generating function brute force calculated up to t^N
    G = Q.prod(lambda n: ((1 - t**n)**(-1))**(c[n]), range(1, N+1)) # the ansatz, product up to and including N
    eqs = [h[k] == G[k] for k in range(0,N+1)] # equating the coefficients
    return solve(eqs,c) #solve for c_i
