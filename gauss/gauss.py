import numpy as np
from collections import deque
import copy
import Forward_elimination

class gauss:
    def __init__(self,n):
        
        self.array=n
        self.answer=np.zeros(len(self.array),dtype=float)
        self.factors=np.zeros(len(self.array),dtype=float)
        self.scalers=self.scaling()
        self.eliminator=Forward_elimination.forward_eliminator(self.array,self.scalers)

    
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
    def Backward_Substitution(self):
       
        arr=self.array
        col=len(arr[0])
        rows=len(arr)
        self.answer[rows-1]=arr[rows-1][col-1]/arr[rows-1][rows-1]
        for i in range(rows-2,-1,-1):
            sum=0
            for j in range (i+1,rows):
                sum+=arr[i][j]*self.answer[j]
            self.answer[i]=(arr[i][col-1]-sum)/arr[i][i]
         
        return self.answer        
    def Gauss(self):
        if(self.eliminator.Forward_Elimination()!=-1):
            self.answer=self.Backward_Substitution()
        else:
            return -1    
        return np.array(self.answer)    
    def next_Step(self):
        ind=self.eliminator.Forward_Elimination_nextStep()
        if ind==-1:
            return -1
        flag,arr=ind
        if flag==True:
            return self.Backward_Substitution()
        else:
            self.array=arr  
            return arr
m1=gauss(np.array([[1,2,2,3],[5,5,2,5],[3,4,6,7]],dtype=float))
print(m1.Gauss())





  










    

