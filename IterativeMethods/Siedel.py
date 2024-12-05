from .Iterations import iterations


class siedel (iterations):
    
    def __init__(self, arrayOFInitials, MatrixA, MatrixB):
        super().__init__(arrayOFInitials, MatrixA, MatrixB)
        
    def getResult(self, iteration = 0):
        currentRes = self.Initials
        for i in range(0, iteration):
            for j in range(0, len(self.Initials)):
                sum_ = sum(self.MatrixA[j][k] * currentRes[k] for k in range(len(self.Initials)) if k != j)
                currentRes[j] = (self.MatrixB[j] - sum_) / self.MatrixA[j][j]
        
        return currentRes
    
    def nextStep(self):
        for j in range(0, len(self.Initials)):
            sum_ = sum(self.MatrixA[j][k] * self.res[k] for k in range(len(self.Initials)) if k != j)
            self.res[j] = (self.MatrixB[j] - sum_) / self.MatrixA[j][j]
            
        return self.res