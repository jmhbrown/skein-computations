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


def get_cokernel(gamma,sigma,take_coinvariants=True):
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


def try_all_conjugacy_classes(gamma,n,take_coinvariants=True):
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

def dim_skeinmod(gamma,n):
    cokernels = try_all_conjugacy_classes(gamma,n,take_coinvariants=True)
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

def Monica_table(L,n):
    """
    sage: gamma = SL2Z.1^4
    sage: try_all_conjugacy_classes(gamma,5)

    """
    conjugacy_classes = [cc.an_element() for cc in SymmetricGroup(n).conjugacy_classes()]
    
    return table([[conjugacy_classes[0]] + conjugacy_classes] + [[l] + [get_cokernel(l,sigma).invariants() for sigma in conjugacy_classes] for l in L])


def Blurange_table(L,n=8):
    """
    Parameters:
        L a vector of gammas in SL2
        n: table for GL_i from 1 to n
    Returns:
        table nxdim(L) of the dimension of the gamma-twisted skein modules in GL_i
    """
    return table([[l] + [dim_skeinmod(l,i) for i in IntegerRange(1,n+1)] for l in L], header_row=[i for i in IntegerRange(0,n+1)], frame = True)
