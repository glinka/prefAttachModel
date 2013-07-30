import numpy as np

def getDegrees(A, n):
    """Return an numpy ndarray containing the 
    degrees of the adjacency matrix

    >>>>getDegrees(np.identity(10))
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
    """
    return np.array([A[i,:].sum() for i in range(n)])
def getAdjEigVals(A, n):
    """Return a numpy ndarray containing the eigenvalues
    of the adjacency matrix, sorted in ascending order
    """
    return np.sort(np.linalg.eigvals(A))
def getAdjEigVects(A, n):
    """Return a numpy ndarray containing the eigenvectors
    of the adjacency matrix, sorted so that the first
    column corresponds to the smallest eigenvalue
    """
    eigvals, eigvects = np.linalg.eig(A)
    sortedIndices = np.argsort(eigvals)
    return [eigvects[:,i] for i in sortedIndices]
def getLaplEigVals(A, n):
    """Return a numpy ndarray containing the eigenvalues
    of the Laplacian matrix, sorted in ascending order
    """
    return np.sort(np.linalg.eigvals(np.diag(getDegrees(A, n))-A))
def getLaplEigVects(A, n):
    """Return a numpy ndarray containing the eigenvectors
    of the Laplacian matrix, sorted so that the first
    column corresponds to the smallest eigenvalue
    """
    eigvals, eigvects = np.linalg.eig(np.diag(getDegrees(A, n))-A)
    sortedIndices = np.argsort(eigvals)
    return [eigvects[:,i] for i in sortedIndices]

def fitXYFunction(X, Y, Z, fns):
    """returns lambda function that represents the linear combination
    of fns with appropriate least-squares-fitted coefficients. assumes
    X.shape = Y.shape"""
    m = X.shape[0]
    n = X.shape[0]
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
                    A[count, k] = fns[k](i,j)
                except ZeroDivisionError:
                    A[count, k] = 20
            count = count + 1
    coeffs = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(A), A)), np.transpose(A)), Zvect)
    return (lambda x,y: np.sum([coeffs[i]*fns[i](x,y) for i in range(nCoeff)], 0))
