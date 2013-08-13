import numpy as np

def getDegrees(A):
    """Return an numpy ndarray containing the 
    degrees of the adjacency matrix

    >>>>getDegrees(np.identity(10))
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
    """
    return np.array([A[i,:].sum() for i in range(A.shape[0])])
def getAdjEigVals(A):
    """Return a numpy ndarray containing the eigenvalues
    of the adjacency matrix, sorted in ascending order
    """
    return np.sort(np.linalg.eigvals(A))
def getAdjEigVects(A):
    """Return a numpy ndarray containing the eigenvectors
    of the adjacency matrix, sorted so that the first
    column corresponds to the smallest eigenvalue
    """
    eigvals, eigvects = np.linalg.eig(A)
    sortedIndices = np.argsort(eigvals)
    return [eigvects[:,i] for i in sortedIndices]

def getAdjLeadingEigVect(A):
    n = A.shape[0]
    eigvals, eigvects = np.linalg.eig(A)
    sortedIndices = np.argsort(eigvals)
    return np.sort(eigvects[:,sortedIndices[n-1]])

def getLaplEigVals(A):
    """Return a numpy ndarray containing the eigenvalues
    of the Laplacian matrix, sorted in ascending order
    """
    return np.sort(np.linalg.eigvals(np.diag(getDegrees(A, A.shape[0]))-A))

def getLaplEigVects(A):
    """Return a numpy ndarray containing the eigenvectors
    of the Laplacian matrix, sorted so that the first
    column corresponds to the smallest eigenvalue
    """
    eigvals, eigvects = np.linalg.eig(np.diag(getDegrees(A, A.shape[0]))-A)
    sortedIndices = np.argsort(eigvals)
    #return [eigvects[:,i] for i in sortedIndices]
    return eigvects[:, sortedIndices]

def fitXYFunction(X, Y, Z, fns):
    """returns lambda function that represents the linear combination
    of fns with appropriate least-squares-fitted coefficients. assumes
    X.shape = Y.shape"""
    m = X.shape[0]
    n = X.shape[1]
    nPoints = n*m
    nCoeff = len(fns)
    Zvect = np.reshape(Z, nPoints)
    A = np.zeros((nPoints,nCoeff))
    count = 0
    for i in range(m):
        for j in range(n):
            for k in range(nCoeff):
                #real poor handling of division by zero
                try:
                    A[count, k] = fns[k](X[i,j], Y[i,j])
                except ZeroDivisionError:
                    A[count, k] = 20
            count = count + 1
    coeffs = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(A), A)), np.transpose(A)), Zvect)
    return (lambda x,y: np.sum([coeffs[i]*fns[i](x,y) for i in range(nCoeff)], 0))

def getSVD(A):
    return np.linalg.svd(A)

def getSVs(A):
    u, s, v = np.linalg.svd(A)
    return s

def getSVLeadingEigVect(A):
    u, s, v = np.linalg.svd(A)
    return u[:,0]

def getEigenReconstruction(A):
    n = A.shape[0]
    vals = getAdjEigVals(A)
    vects = np.array(getAdjEigVects(A))
    u = vects[:,n-1]
    u.shape = (1,n)
    ANew = vals[n-1]*np.dot(np.transpose(u), np.conj(u))
    degs = getDegrees(ANew)
    i = np.argsort(degs)
    cpy = ANew
    ANew[:,(n-1) - np.arange(n)] = cpy[:,i]
    cpy = ANew
    ANew[(n-1) - np.arange(n),:] = cpy[i,:]
    return np.transpose(np.rint(ANew))
