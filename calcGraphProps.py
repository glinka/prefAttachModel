import numpy as np

class calcGraphProps:
    def __init__(self, AdjMtrx):
        self.A = AdjMtrx
        self.n = AdjMtrx.shape[0]
    def getDegrees(self):
        """Return an numpy ndarray containing the 
        degrees of the adjacency matrix

        >>>>getDegrees(np.identity(10))
        array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
        """
        return np.array([self.A[i,:].sum() for i in range(self.n)])
    def getAdjEigVals(self):
        """Return a numpy ndarray containing the eigenvalues
        of the adjacency matrix, sorted in ascending order
        """
        return np.sort(np.linalg.eigvals(self.A))
    def getAdjEigVects(self):
        """Return a numpy ndarray containing the eigenvectors
        of the adjacency matrix, sorted so that the first
        column corresponds to the smallest eigenvalue
        """
        eigvals, eigvects = np.linalg.eig(self.A)
        sortedIndices = np.argsort(eigvals)
        return [eigvects[:,i] for i in sortedIndices]
    def getLaplEigVals(self):
        """Return a numpy ndarray containing the eigenvalues
        of the Laplacian matrix, sorted in ascending order
        """
        return np.sort(np.linalg.eigvals(np.diag(getDegrees())-self.A))
    def getLaplEigVects(self):
        """Return a numpy ndarray containing the eigenvectors
        of the Laplacian matrix, sorted so that the first
        column corresponds to the smallest eigenvalue
        """
        eigvals, eigvects = np.linalg.eig(np.diag(self.getDegrees())-self.A)
        sortedIndices = np.argsort(eigvals)
        return [eigvects[:,i] for i in sortedIndices]
