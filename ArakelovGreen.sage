# A collection of mathematical tools related to metrized graphs. Most
# importantly, this includes the MetrizedGraph class and a method to 
# compute the Arakelov-Green function g_muD for given a divisor D.

var('x','y')

### METRIZED GRAPH CLASS

class MetrizedGraph(DiGraph):

    """ Metrized Graph.

    This implementation equips a metrized graph with a fixed (user-specified) 
    vertex set, edge set and orientations.

    INPUT:

    A metrized graph is initialized with MetrizedGraph(V,E), where
        
        V   is an array of enumerated vertices;
        
        E   is an array of enumerated edges, which are encoded in the form
            (p,q,l). Here, p is the first point of the edge, q the second,
            and l the length.

    ATTRIBUTES:

    Upon initialization, several attributes are computed:

        L       The discrete Laplacian matrix of the metrized graph

        Lplus   The Moore-Penrose generalized inverse of L

        C       The connectivity matrix of the metrized graph

    METHODS:

    The following methods are available:

        v(i)            Return the ith vertex

        e(i)            Return the ith edge

        p(i)            Return the first vertex of the ith edge

        q(i)            Return the second vertex of the ith edge

        lth(i)          Return the length of the ith edge

        tot_lth()       Compute the total length of the metrized graph

        val(p)          Compute the valence of a vertex p

        res(p,q)        Compute the resistance between vertices p and q
        
        vol(p,q,s)      Compute the voltage function j_s(p,q) on vertices

        tau()           Compute the tau constant of the metrized graph

        can_div()       Compute the canonical divisor of the metrized graph

        r(D)            Compute the function r(D,.) on each edge

        val_mat(D)      Compute the value matrix of g_muD

        eps_inv(D)      Compute the epsilon invariant wrt the divisor D """

    def __init__(self, V, E):
        # Store the vertices, edges, and the order they were given in
        self._order_vertices = V
        self._order_edges = E

        # Inherit everything from the DiGraph class
        super().__init__([V, E], weighted=True)

        # Compute discrete Laplacian matrix and its Moore-Penrose generalized inverse
        self.L = discrete_laplacian_matrix(self).simplify_rational()
        self.Lplus = moore_penrose(self.L).simplify_rational()

        # Compute connectivity matrix
        self.C = connectivity(self)

    def __repr__(self):
        n = self.order()
        m = self.size()
        return "Metrized graph with %s vertices and %s edges" %(n, m)

    # Basic attributes that use the order in which vertices and edges were given
    def v(self, i):
        return self._order_vertices[i]

    def e(self, i):
        return self._order_edges[i]

    def p(self, i):
        return self._order_edges[i][0]

    def q(self, i):
        return self._order_edges[i][1]

    def lth(self, i):
        return self._order_edges[i][2]

    def tot_lth(self):
        return sum(self.lth(i) for i in range(self.size()))

    # Valence of a vertex
    def val(self, p):
        val = 0

        for i in range(self.size()):
            if p in (self.p(i),self.q(i)):
                val += 1

        return val

    # Resistance function on vertices
    def res(self, p, q):
        i = self._order_vertices.index(p)
        j = self._order_vertices.index(q)

        resist = self.Lplus[i][i] - 2*self.Lplus[i][j] + self.Lplus[j][j]
        return resist.full_simplify()

    # Voltage function on vertices
    def vol(self, p, q, s):
        i = self._order_vertices.index(p)
        j = self._order_vertices.index(q)
        k = self._order_vertices.index(s)

        volt = self.Lplus[k][k] - self.Lplus[i][k] - self.Lplus[j][k] + self.Lplus[i][j]
        return volt.full_simplify()

    # Tau constant
    def tau(self):
        L = self.L
        Lplus = self.Lplus

        n = self.order()

        tau = - 1/12 * sum(L[i][j]*(1/L[i][j] + Lplus[i][i] - 2*Lplus[i][j] + Lplus[j][j])^2 for i in range(n) for j in range(n) if self.has_edge((self.v(i),self.v(j))))\
              + 1/4 * sum(L[i][j]*Lplus[i][i]*Lplus[j][j] for i in range(n) for j in range(n)) \
              + 1/self.order()*Lplus.trace()

        return tau.full_simplify()

    # Canonical divisor
    def can_div(self, q=None):
        """Computes the canonical divisor of the metrized graph given a polarization 
        on each vertex.

        The polarization is encoded as an array where the ith entry corresponds to
        the value on the ith vertex. If no polarization is specified, all values
        are assumed to be 0."""

        n = self.order()

        if q == None:
            q = [0 for i in range(n)]
        else:
            assert len(q) == n

        D = [self.val(self.v(i)) - 2 + 2*q[i] for i in range(n)]

        return D

    # r(D)-array per edge; uses parallelization
    def r(self, D=None):
        """Computes r(D,.) on each edge.

        If the divisor D is not specified, D = 0 is used."""

        n = self.order()
        m = self.size()

        if D is None:
            D = [0 for i in range(n)]
        else:
            assert len(D) == n

        R = [0 for i in range(m)]
        generator = r([(self,i,D) for i in range(m)])
        for g in generator:
            R[g[0][0][1]] = g[1]

        return R

    # Value matrix for g_muD; uses parallelization
    def val_mat(self, D=None):
        """Computes the value matrix for g_muD.

        The (i,j)th entry of the value matrix is a formula for g_muD(x,y) when x lies 
        on the ith and y lies on the jth edge. If the divisor D is not specified, D = 0 
        is used, yielding the value matrix for the canonical Arakelov-Green function."""

        n = self.order()
        m = self.size()

        if D is None:
            D = [0 for i in range(n)]
        else:
            assert len(D) == n

        # Determine r-array
        R = self.r(D)

        # Compute tau constant and deg(D)
        tau = self.tau()
        deg = sum(D)

        # Compute constant cD
        cD = 8*tau*(deg+1) + sum(D[i]*D[j]*self.res(self.v(i),self.v(j)) for i in range(n) for j in range(n))
        cD *= 1/(2*(deg+2)^2)

        # Compute value matrix
        ZD = matrix(SR,m)
        generator = Arakelov_Green([(self,i,j,D,R,tau,cD) for i in range(m) for j in range(m)])
        for g in generator:
            (i,j) = (g[0][0][1],g[0][0][2])
            ZD[i,j] = g[1]

        return ZD

    # Epsilon invariant
    def eps_inv(self, D=None, method='thm'):
        """Computes the epsilon invariant with respect to a given divisor.

        This function can use two different methods: lem, which uses the lemma
        by Moriwaki; or thm, which uses the theorem by Cinkir. If no method is 
        specified, thm is used."""

        n = self.order()
        m = self.size()

        assert method == 'thm' or method == 'lem'
        assert len(D) == n

        if method == 'thm':
            eps = 4*??.tau()*sum(D) + sum(D[i]*D[j]*??.res(??.v(i),??.v(j)) for i in range(n) for j in range(n))
            eps /= sum(D) + 2
            return eps.full_simplify()

        else:
            ZD = ??.val_mat(D)
            eps = 0

            for k in range(n):
                v = ??.v(k)

                for l in range(m):
                    if v is ??.p(l):
                        gD = ZD[0,l](x=0,y=0)
                        break
                    elif v is ??.q(l):
                        gD = ZD[0,l](x=0,y=??.lth(l))
                        break

                eps += D[k] * gD

            eps *= sum(D) + 2
            eps += sum(D[k] * ??.res(??.p(0),??.v(k)) for k in range(n))

            return eps.full_simplify()


### DISCRETE LAPLACIAN MATRIX

# Function to compute the discrete Laplacian matrix of a metrized graph ??.
def discrete_laplacian_matrix(??):
    n = ??.order()

    L = matrix(SR,n)

    for i in range(n):
        for j in range(i+1,n):
            if Graph(??).has_edge((??.v(i),??.v(j))):
                L[i,j] = -1/Graph(??).edge_label(??.v(i),??.v(j))
                L[j,i] = L[i,j]

        L[i,i] = -sum(L[i,k] for k in range(n) if k != i)

    return L

# Function to compute the Moore-Penrose generalized inverse of a discrete Laplacian matrix.
def moore_penrose(L):
    n = L.nrows()

    Jn = matrix.ones(n)
    Lplus = (L - (1/n)*Jn).inverse() + (1/n)*Jn

    return Lplus


### CONNECTIVITY MATRIX

# Helper function alpha. 
def alpha(??, i, j):
    # Delete e_i to obtain ?? - e_i
    ??.delete_edge(??._order_edges[i])
    
    # Determine which connected subgraph p_{e_j} belongs to; this is the
    # subgraph that e_j belongs to.
    if ??.p(i) in ??.connected_component_containing_vertex(??.p(j)):
        alpha = 0
    else:
        alpha = 1

    # Restore e_i
    ??.add_edge(??._order_edges[i])

    return alpha

# Helper function beta. 
def beta(??, i, j):
    # Temporarily set the lengths of e_i and e_j to 1, and the lengths of other 
    # edges to 0. The closest neighbours of e_i and e_j is then points on the 
    # respective edges that are distance 0 apart from each other.

    m = ??.size()

    for k in range(m):
        if k == i or k == j:
            ??.set_edge_label(??.p(k),??.q(k),1)
        else:
            ??.set_edge_label(??.p(k),??.q(k),0)

    # Determine the closest neighbours
    for index in range(2):
        for jndex in range(2):
            dist = Graph(??).shortest_path_length(??.e(i)[index],??.e(j)[jndex],by_weight=True) 
            if dist == 0:
                # Restore the edge lengths
                for k in range(m):
                    ??.set_edge_label(??.p(k),??.q(k),??.lth(k))

                return 10*index + jndex

# Function to compute the connectivity matrix
def connectivity(??):
    m = ??.size()

    # Initialize C
    C = matrix(m)

    # Set diagonal entries
    for i in range(m):
        if ??.is_cut_edge(??.e(i)):
            C[i,i] = 1

    # Set non-diagonal entries
    for i in range(m):
        if C[i,i] == 1:                     # If e_i is a bridge,
            for j in range(i+1,m):          # Loop over the next edges e_j
                C[i,j] = alpha(??,i,j)       # Assume that e_j is not a bridge
                if C[j,j] == 1:             # If e_j is a bridge, adjust
                    C[i,j] = 100*C[i,j] + beta(??,i,j)
                    C[j,i] = 100*alpha(??,j,i) + beta(??,j,i)
                else:
                    C[j,i] = C[i,j]         # If e_j is not a bridge, set c_ji=c_ij

        else:                               # If e_i is not a bridge,
            for j in range(i+1,m):          # Loop over the next edges e_j
                if C[j,j] == 1:             # If e_j is a bridge, set 
                    C[i,j] = alpha(??,j,i)   # nonzero entry
                    C[j,i] = C[i,j]

    return C


### ARAKELOV-GREEN FUNCTION

# Given a divisor D, this function computes r(D,.) on the edge e_i of the 
# metrized graph ??. With the @parallel decorator, this function can make 
# use of multiple processors when it is applied to multiple edges.
@parallel
def r(??, i, D):
    n = ??.order()
    m = ??.size()

    # Initialize r
    rD = 0

    # If e_i is not a bridge
    if ??.C[i,i] == 0:
        for k in range(n):
            tmp = - x^2 * (??.lth(i) - ??.res(??.p(i),??.q(i)))/??.lth(i)^2 \
                  + x * (??.lth(i) - ??.res(??.p(i),??.q(i)) + ??.res(??.v(k),??.q(i)) - ??.res(??.v(k),??.p(i)))/??.lth(i) \
                  + ??.res(??.v(k),??.p(i))
            rD += D[k] * tmp

    # If e_i is a bridge
    else:
        # Iterate over all vertices, add the right expression to r for each vertex
        for k in range(n):
            if ??.v(k) == ??.p(i):
                rD += D[k] * x
            elif ??.v(k) == ??.q(i):
                rD += D[k] * (-x + ??.lth(i))
            else:
                # If vertex k is not an end point of e_i, search for an edge
                # that the vertex is an end point of. Then use the
                # connectivity matrix to determine the suitable expression.
                for j in range(m):
                    if j != i and (??.v(k) == ??.p(j) or ??.v(k) == ??.q(j)):
                        if (??.C[j,j] == 0 and ??.C[i,j] == 0) or (??.C[j,j] == 1 and ??.C[i,j] < 100):
                            rD += D[k] * (x + ??.res(??.v(k),??.p(i)))
                        else:
                            rD += D[k] * (- x + ??.lth(i) + ??.res(??.v(k),??.q(i)))
                        break

    return rD

# For a given divisor D, this function computes the Arakelov-Green function g_muD on 
# the pair of edges (e_i,e_j) of the metrized graph ??. The input consists several values
# that should be computed beforehand. With the @parallel decorator, this function can 
# make use of multiple processors when it is applied to multiple pairs of edges.
@parallel
def Arakelov_Green(??, i, j, D, R, tau, cD):
    deg = sum(D)

    # Short-hand for clunky recurring expression
    tauD = lambda i,j: 1/(deg+2) * (4*tau + 1/2 * (R[i] + R[j].substitute(x=y))) - cD

    # Same edge
    if i == j:
        if ??.C[i,i] == 0:   # Not a bridge
            gD = tauD(i,i) - 1/2*abs(x-y) + (x-y)^2 * (??.lth(i) - ??.res(??.p(i),??.q(i)))/(2*??.lth(i)^2)
        else:               # Bridge
            gD = tauD(i,i) - 1/2*abs(x-y)

    # Different edges
    else:
        # Neither a bridge
        if ??.C[i,i] == 0 and ??.C[j,j] == 0:     
            gD = tauD(i,j) \
                + x^2 * (??.lth(i) - ??.res(??.p(i),??.q(i)))/(2*??.lth(i)^2) \
                + y^2 * (??.lth(j) - ??.res(??.p(j),??.q(j)))/(2*??.lth(j)^2) \
                - (x*y)/(??.lth(i)*??.lth(j)) * (??.vol(??.p(i),??.q(j),??.p(j)) - ??.vol(??.q(i),??.q(j),??.p(j))) \
                - x/(2*??.lth(i)) * (??.lth(i) - 2*??.vol(??.q(i),??.p(j),??.p(i))) \
                - y/(2*??.lth(j)) * (??.lth(j) - 2*??.vol(??.p(i),??.q(j),??.p(j))) \
                - 1/2*??.res(??.p(i),??.p(j))
            
        # Edge e_i a bridge, e_j not
        elif ??.C[i,i] == 1 and ??.C[j,j] == 0:
            if ??.C[i,j] == 0: # If e_j in ??_{p_{e_i}}
                T = y * (??.lth(j) - ??.res(??.p(j),??.q(j)) + ??.res(??.p(i),??.q(j)) - ??.res(??.p(i),??.p(j)))/??.lth(j) + x + ??.res(??.p(i),??.p(j))
                gD = tauD(i,j) + y^2 * (??.lth(j) - ??.res(??.p(j),??.q(j)))/(2*??.lth(j)^2) - 1/2*T

            else:             # If e_j in ??_{q_{e_i}}
                T = y * (??.lth(j) - ??.res(??.p(j),??.q(j)) + ??.res(??.q(i),??.q(j)) - ??.res(??.q(i),??.p(j)))/??.lth(j) - x + ??.lth(i) + ??.res(??.q(i),??.p(j))
                gD = tauD(i,j) + y^2 * (??.lth(j) - ??.res(??.p(j),??.q(j)))/(2*??.lth(j)^2) - 1/2*T

        # Edge e_j a bridge, e_i not
        elif ??.C[i,i] == 0 and ??.C[j,j] == 1: 
            if ??.C[i,j] == 0: # If e_i in ??_{p_{e_j}}
                T = x * (??.lth(i) - ??.res(??.p(i),??.q(i)) + ??.res(??.p(j),??.q(i)) - ??.res(??.p(j),??.p(i)))/??.lth(i) + y + ??.res(??.p(j),??.p(i))
                gD = tauD(i,j) + x^2 * (??.lth(i) - ??.res(??.p(i),??.q(i)))/(2*??.lth(i)^2) - 1/2*T

            else:             # If e_i in ??_{q_{e_j}}
                T = x * (??.lth(i) - ??.res(??.p(i),??.q(i)) + ??.res(??.q(j),??.q(i)) - ??.res(??.q(j),??.p(i)))/??.lth(i) - y + ??.lth(j) + ??.res(??.q(j),??.p(i))
                gD = tauD(i,j) + x^2 * (??.lth(i) - ??.res(??.p(i),??.q(i)))/(2*??.lth(i)^2) - 1/2*T
            
        # Both bridges
        else:
            if ??.C[i,j] % 100 == 0:
                gD = tauD(i,j) - 1/2*(x + y + ??.res(??.p(i),??.p(j)))
                     
            elif ??.C[i,j] % 100 == 1:
                gD = tauD(i,j) - 1/2*(x - y + ??.lth(j) + ??.res(??.p(i),??.q(j)))
                    
            elif ??.C[i,j] % 100 == 10:
                gD = tauD(i,j) - 1/2*(- x + y + ??.lth(i) + ??.res(??.q(i),??.p(j))) 
                    
            else:
                gD = tauD(i,j) - 1/2*(- x - y + ??.lth(i) + ??.lth(j) + ??.res(??.q(i),??.q(j)))

    return gD.full_simplify()