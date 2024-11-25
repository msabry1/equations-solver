from .Iterations import iterations

class jacobi (iterations):

    def __init__(self, arrayOFInitials, MatrixA, MatrixB):
        super().__init__(arrayOFInitials, MatrixA, MatrixB)

    def getResult(self, iteration):
        currentRes = self.Initials 
        for i in range(iteration):
            finalRes = [
            for j in range(len(self.Initials)):
                sum_ = sum(self.MatrixA[j][k] * currentRes[k] for k in range(len(self.Initials)) if k != j)
                finalRes.append((self.MatrixB[j] - sum_) / self.MatrixA[j][j]) #B is matrix of values
            
            currentRes = finalRes

        return currentRes

    def nextStep(self):
        currentRes = []
        for j in range(len(self.Initials)):
            sum_ = sum(self.MatrixA[j][k] * self.res[k] for k in range(len(self.Initials)) if k != j)
            currentRes.append((self.MatrixB[j] - sum_) / self.MatrixA[j][j]) #B is matrix of values
                
        self.res = currentRes
        return self.res