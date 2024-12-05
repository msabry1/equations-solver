from abc import ABC, abstractmethod

class iterations (ABC):
    
    def __init__ (self, array, MatrixA, MatrixB):
        self.Initials = array.copy()
        self.MatrixA = MatrixA.copy()
        self.MatrixB = MatrixB.copy()
        self.res = array.copy()

    @abstractmethod
    def getResult(self, iteration = 0):
        pass
    
    @abstractmethod
    def nextStep(self):
        pass