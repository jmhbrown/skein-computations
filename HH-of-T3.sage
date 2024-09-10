import logging


logger = logging.getLogger(__name__)
logging.basicConfig(level='INFO') # options: DEBUG, INFO, WARN, ERROR, CRITICAL


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

    lattice = FreeModule(ZZ,2*n)

    sublattice_gens = (Matrix.identity(2*n) - W_rep(sigma)*MCG_rep(gamma,n)).columns()
    if take_coinvariants:
        for tau in SymmetricGroup(n).centralizer(sigma):
            sublattice_gens += (Matrix.identity(2*n) - W_rep(tau)).columns()

    return lattice.quotient(lattice.submodule(sublattice_gens))


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
        co['torsion_dim'] = sum(co['homology'].invariants())
    return cokernels

