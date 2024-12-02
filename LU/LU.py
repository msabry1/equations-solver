import numpy as np

class LU:
    def __init__(self, matrix, B):
       self.factors = []  # Changed from zeros to list
       self.tol = 1e-9
       self.array = matrix.copy()
       self.coefficient_matrix = matrix.copy()
       self.constant_vector = B.copy()
       
       # Consume the generator and get the scalers
       scaling_generator = self.scaling()
       scaling_steps = list(scaling_generator)
       self.scalers = self.scaling()
       self.results = B.copy()
       self.F_results = np.zeros(len(B))
       self.U = np.zeros((len(self.array), len(self.array)))
       self.L = np.zeros((len(self.array), len(self.array)))
       self.solution_type = self.detect_system_type()
    
    def detect_system_type(self):
        A = self.coefficient_matrix
        Aug = np.column_stack((A, self.constant_vector))
        rank_A = np.linalg.matrix_rank(A)
        rank_Aug = np.linalg.matrix_rank(Aug)
        n = A.shape[1]
        
        if rank_A == rank_Aug:
            if rank_A == n:
                self.solution_type = 'unique'
                return 'unique'
            else:
                self.solution_type = 'infinite'
                return 'infinite'
        else:
            self.solution_type = 'no solution'
            return 'no solution'
           
    def scaling(self):
        arr = self.array
        scalers = [0] * len(arr)
       
        for i in range(len(arr)):
            max_val = abs(arr[i][0])
            for j in range(1, len(arr)):
                if abs(arr[i][j]) > max_val:
                   max_val = abs(arr[i][j])
        
            scalers[i] = max_val
        
        self.scalers = scalers
       
        return scalers
    def pivoting(self, row):
        index = row
      
        max_ratio = abs(self.array[row][row] / self.scalers[row])
        for i in range(row+1,len(self.array)):
            current_ratio = abs(self.array[i][row] / self.scalers[i])
            if current_ratio > max_ratio:
                max_ratio = current_ratio
                index = i
        if index != row:
            self.array[[row, index]] = self.array[[index, row]]
            self.scalers[row], self.scalers[index] = self.scalers[index], self.scalers[row]
       
        return self.array
    def get_U_generator(self):
        if self.solution_type == 'no solution':
            raise Exception("There is no solution")
      
        if self.solution_type == 'infinite':
            raise Exception("There is an infinite number of solutions")
        
        arr = self.array.copy()
        rows = len(arr)
        
        yield  arr.copy()
        
        for k in range(rows):
   
            arr=self.pivoting(k)
           
            
 
            for i in range(k + 1, rows):
                factor = arr[i][k] / arr[k][k]
                self.factors.append(factor)
                
                for j in range(k, rows):
                    arr[i][j] -= factor * arr[k][j]
                
                yield arr.copy()
        

        self.U = arr.copy()

        yield self.U.copy()
    
    def forward_substitution_generator(self):
        if self.solution_type == 'no solution':
            raise Exception("There is no solution")
      
        if self.solution_type == 'infinite':
            raise Exception("There is an infinite number of solutions")
        
        L = self.get_L()
        n = np.zeros((len(self.array), len(self.array) + 1))
        
        # Prepare augmented matrix
        for i in range(len(n)):
            for j in range(len(n)):
                n[i][j] = L[i][j]
        for i in range(len(n)):
            n[i][len(n)] = self.results[i]
        
        yield n.copy()
        
        rows = len(n)
        col = len(n[0])
        
  
        self.F_results[0] = n[0][col-1] / n[0][0]
    
        for i in range(1, len(n)):
            sum_val = 0
            for j in range(0, rows):
                if j != i:
                    sum_val += n[i][j] * self.F_results[j]
            
            self.F_results[i] = (n[i][col-1] - sum_val) / n[i][i]
            
        yield self.F_results
        
        return self.F_results
    
    def backward_substitution_generator(self):
        if self.solution_type == 'no solution':
            raise Exception("There is no solution")
      
        if self.solution_type == 'infinite':
            raise Exception("There is an infinite number of solutions")
        
        U = self.U
        F_results = self.F_results
        
        n = np.zeros((len(self.array), len(self.array) + 1))
        for i in range(len(n)):
            for j in range(len(n)):
                n[i][j] = U[i][j]
        for i in range(len(n)):
            n[i][len(n)] = F_results[i]
        
        yield n.copy()
        
        rows = len(n)
        col = len(n[0])
        
    
        self.results[rows-1] = n[rows-1][col-1] / n[rows-1][rows-1]
        
        
        
        for i in range(rows-2, -1, -1):
            sum_val = 0
            for j in range(i+1, rows):
                sum_val += n[i][j] * self.results[j]
            
            self.results[i] = (n[i][col-1] - sum_val) / n[i][i]
            
        yield self.results
        
        return self.results
    
    def get_L(self):

        ptr = 0
        for i in range(len(self.array)):
            for j in range(len(self.array)):
                if i == j:
                    self.L[i][j] = 1
                elif i < j:
                    self.L[j][i] = self.factors[ptr]
                    ptr += 1
        return self.L
    def getfinal(self):
        U=0
        F_res=0
        results=0
       
        for step in self.get_U_generator():
            print(step)
            U=step
        for step in self.forward_substitution_generator():
            print(step)
            F_res=step
        for step in self.backward_substitution_generator():
            print(step)
            results=step
            

        return U,self.L,F_res,results       
        

def main():
    m1=LU(np.array([[0,2,5],[2,1,1],[3,1,0]],dtype=float),np.array([1,1,2],dtype=float))
    print(m1.solution_type)
    print(m1.getfinal())
main()  