import numpy as np
class LU:
    def __init__(self, matrix, B):
       self.factors = []  # Changed from zeros to list
       self.tol = 1e-9
       self.array = matrix.copy()
       self.coefficient_matrix=matrix.copy()
       self.constant_vector=B.copy()
    # Consume the generator and get the scalers
       scaling_generator = self.scaling()
       scaling_steps = list(scaling_generator)
       self.scalers = [step['max_value'] for step in scaling_steps if 'max_value' in step]
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
        yield {
        'step': 'Scaling Initialization', 
        'description': 'Starting scaling process'
        }
        for i in range(len(arr)):
            max_val = abs(arr[i][0])
            for j in range(1, len(arr)):
                if abs(arr[i][j]) > max_val:
                   max_val = abs(arr[i][j])
        
            scalers[i] = max_val
        
            yield {
            'step': f'Scaling Row {i}', 
            'max_value': max_val,
           'description': f'Maximum value for row {i}'
           }
        self.scalers = scalers
        return scalers
    def get_U_generator(self):
        
        if(self.solution_type=='no solution'):
            yield{
                "error":'there is no solution'
            }
            raise Exception("there is no solution")
      
        if(self.solution_type=='infinite'):
            yield{
                "error":'there is an infinite number of solution'
                
            }    
            raise Exception("there is an infinite number of solution")
        arr = self.array.copy()
        rows = len(arr)
        
        # Yield initial matrix
        yield {
            'step': 'Initial Matrix', 
            'matrix': arr.copy(), 
            'description': 'Starting matrix before decomposition'
        }
        for k in range(rows):
            # Check pivot
            if abs(arr[k][k] / self.scalers[k]) < self.tol:
                yield {
                    'step': 'Pivot Check Failed', 
                    'matrix': arr.copy(), 
                    'description': f'Pivot at row {k} is too small'
                }
                return
            # Perform elimination
            for i in range(k + 1, rows):
                factor = arr[i][k] / arr[k][k]
                self.factors.append(factor)
                
                yield {
                    'step': f'Elimination Step (Row {i}, Column {k})', 
                    'factor': factor,
                    'matrix_before': arr.copy(),
                    'description': f'Eliminating element at row {i}, column {k}'
                }
                for j in range(k, rows):
                    arr[i][j] -= factor * arr[k][j]
                
                yield {
                    'step': f'Matrix After Elimination (Row {i}, Column {k})', 
                    'matrix': arr.copy(),
                    'description': 'Matrix after performing elimination'
                }
        # Final U matrix
        self.U = arr.copy()
        yield {
            'step': 'Final U Matrix', 
            'matrix': self.U.copy(), 
            'description': 'Completed Upper Triangular Matrix'
        }
    def scaling(self):
        """Scaling method with generator"""
        arr = self.array
        scalers = [0] * len(arr)
        
        yield {
            'step': 'Scaling Initialization', 
            'description': 'Starting scaling process'
        }
        for i in range(len(arr)):
            max_val = abs(arr[i][0])
            for j in range(1, len(arr)):
                if abs(arr[i][j]) > max_val:
                    max_val = abs(arr[i][j])
            
            scalers[i] = max_val
            
            yield {
                'step': f'Scaling Row {i}', 
                'max_value': max_val,
                'description': f'Maximum value for row {i}'
            }
        self.scalers = scalers
        return scalers
    def forward_substitution_generator(self):
        if(self.solution_type=='no solution'):
            yield{
                "error":'there is no solution'
            }
            raise Exception("there is no solution")
      
        if(self.solution_type=='infinite'):
            yield{
                "error":'there is an infinite number of solution'
                
            }    
            raise Exception("there is an infinite number of solution")
        L = self.get_L()
        n = np.zeros((len(self.array), len(self.array) + 1))
        # Prepare augmented matrix
        for i in range(len(n)):
            for j in range(len(n)):
                n[i][j] = L[i][j]
        for i in range(len(n)):
            n[i][len(n)] = self.results[i]
        yield {
            'step': 'Initial Augmented Matrix', 
            'matrix': n.copy(),
            'description': 'Augmented matrix for forward substitution'
        }
        rows = len(n)
        col = len(n[0])
        # Initialize results
        self.F_results[0] = n[0][col-1] / n[0][0]
        
        yield {
            'step': 'First Solution', 
            'value': self.F_results[0],
            'description': 'First value in forward substitution'
        }
        # Compute remaining solutions
        for i in range(1, len(n)):
            sum_val = 0
            for j in range(0, rows):
                if j != i:
                    sum_val += n[i][j] * self.F_results[j]
            
            self.F_results[i] = (n[i][col-1] - sum_val) / n[i][i]
            
            yield {
                'step': f'Solution for Row {i}', 
                'value': self.F_results[i],
                'sum': sum_val,
                'description': 'Computed solution for current row'
            }
        return self.F_results
    def backward_substitution_generator(self):
        
        if(self.solution_type=='no solution'):
            yield{
                "error":'there is no solution'
            }
            raise Exception("there is no solution")
      
        if(self.solution_type=='infinite'):
            yield{
                "error":'there is an infinite number of solution'
                
            }    
            raise Exception("there is an infinite number of solution")
        U = self.U
        F_results = self.F_results
        
        n = np.zeros((len(self.array), len(self.array) + 1))
        for i in range(len(n)):
            for j in range(len(n)):
                n[i][j] = U[i][j]
        for i in range(len(n)):
            n[i][len(n)] = F_results[i]
        yield {
            'step': 'Initial Augmented Matrix', 
            'matrix': n.copy(),
            'description': 'Augmented matrix for backward substitution'
        }
        rows = len(n)
        col = len(n[0])
        # Compute last solution
        self.results[rows-1] = n[rows-1][col-1] / n[rows-1][rows-1]
        
        yield {
            'step': 'Last Solution', 
            'value': self.results[rows-1],
            'description': 'First solution from the end'
        }
        # Compute remaining solutions
        for i in range(rows-2, -1, -1):
            sum_val = 0
            for j in range(i+1, rows):
                sum_val += n[i][j] * self.results[j]
            
            self.results[i] = (n[i][col-1] - sum_val) / n[i][i]
            
            yield {
                'step': f'Solution for Row {i}', 
                'value': self.results[i],
                'sum': sum_val,
                'description': 'Computed solution for current row'
            }
        return self.results
    def get_L(self):
        """Compute Lower Triangular Matrix"""
        ptr = 0
        for i in range(len(self.array)):
            for j in range(len(self.array)):
                if i == j:
                    self.L[i][j] = 1
                elif i < j:
                    self.L[j][i] = self.factors[ptr]
                    ptr += 1
        return self.L
def main():

    m1 = LU(np.array([[25,5,1],[25,5,1],[144,12,1]], dtype=float), [1,1,3])
    print(m1.solution_type)
    print("U Matrix Decomposition Steps:")
    for step in m1.get_U_generator():
        print(step)
    
    print("\nForward Substitution Steps:")
    for step in m1.forward_substitution_generator():
        print(step)
    
    print("\nBackward Substitution Steps:")
    for step in m1.backward_substitution_generator():
        print(step)
if __name__ == "__main__":
    main()