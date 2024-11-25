import numpy as np

class gauss:
    array=np.ones((3,3))
    answer=np.zeros(3)
    scalers=np.zeros(3)
    tol=1e-9
    def __init__(self,n):
        self.array=n
        self.answer=np.zeros(len(self.array))
        self.scalers=np.zeros(len(self.array))
    
    def scaling(self):
        arr = self.array
        self.scalers = [0]*len(arr)

        for i in range(len(arr)): 
            max_val = abs(arr[i][0]) 

            for j in range(1, len(arr)): 
                if abs(arr[i][j]) > max_val:
                    max_val = abs(arr[i][j])
        
            self.scalers[i] = max_val  
        return self.scalers   
    def pivoting(self, row):
        index = row
        mx = abs(self.array[row][row] / self.scalers[row]) 
        for i in range(row + 1, len(self.array)):
            tmp = abs(self.array[i][row] / self.scalers[i])
            if tmp > mx:
                mx = tmp
                index = i
        if index != row:
            self.array[row], self.array[index] = self.array[index], self.array[row]
            self.scalers[row], self.scalers[index] = self.scalers[index], self.scalers[row]
        return self.array
    def Forward_Elimination(self):
        arr = self.array
        rows = len(arr)
        cols = len(arr[0])
        for k in range(rows):
            self.pivoting(k)
            if(abs(arr[k][k])/self.scalers[k]<self.tol):
                return -1

            for i in range(k + 1, rows):
                factor = arr[i][k] / arr[k][k]
                for j in range(k, cols):
                    arr[i][j] -= factor * arr[k][j]
            if(abs(arr[rows-1][rows-1]/self.scalers[rows-1])<self.tol):
                return -1
        self.array=arr

    def Backward_Substitution(self):
        res=self.answer
        arr=self.array
        col=len(arr[0])
        rows=len(arr)
        res[rows-1]=arr[rows-1][col-1]/arr[rows-1][rows-1]
        for i in range(rows-2,-1,-1):
            sum=0
            for j in range (i+1,rows):
                sum+=arr[i][j]*res[j]
            res[i]=(arr[i][col-1]-sum)/arr[i][i]
        return res        
    def Gauss(self):
        self.scalers=self.scaling()
        if(self.Forward_Elimination()!=-1):
            self.answer=self.Backward_Substitution()
        else:
            return -1    
        return self.answer                
m1=gauss([[1,5,1,4],[1,-1,5,7],[1,-2,5,12]])

print(m1.Gauss())
    

