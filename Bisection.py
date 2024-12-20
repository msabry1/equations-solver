import sympy as sp
def relativeError(x1,x2):
    if x1==0 or x2==0:
        return 17.1717
    return abs((x2-x1)/x2)
def check_validity(x):
    return  x == sp.zoo or not x.is_real or not x.is_finite
def Bisection(str,a,b,percision=5,tol=1e-10,mx_iter=15):
    a=sp.N(a,percision)
    b=sp.N(b,percision)
    
    x=sp.symbols('x')
    equ=sp.N(sp.sympify(str),percision)
    if(equ.subs(x,a)>0):
        tmp=sp.N(b,percision)
        b=sp.N(a,percision)
        a=tmp
    prev=0
    print(a,b)
    if(check_validity(equ.subs(x,a)) or check_validity(equ.subs(x,b))):
        raise Exception("the function is not continuous at the given range")
    while mx_iter>0:
        mx_iter-=1
        mid=(b+a)/2
        if(check_validity(equ.subs(x,mid))):
            raise Exception("the function is not continuous at the given range")
        print(mid)
        if(equ.subs(x,a)*equ.subs(x,b)>1):
            raise Exception("this equation cannot be solved by the Bisection method")
        if(equ.subs(x,mid)>0):
            b=mid
        else :
            a=mid
        yield [mid,relativeError(prev,mid)]
        if(relativeError(prev,mid)<tol):
            break
        prev =mid
def main():   
    # equation = "x**4 - 2*x**3 - 4*x**2 + 4*x + 4"
    # root, steps, iterations = bisection(equation, -2, -1, precision=5, eps=1e-3)
    equation = "1.5*x - 6 - 1/2 * sin(2*x)"
    # error, steps, roots = bisection(equation, 4, 5, precision=5, eps=1e-5)
    # equation = "x**3 - 0.165*x**2 + 3.993*10**-4"
    # error, steps, roots = bisection(equation, 0, 0.11, precision=4, eps=2e-3)
    # equation = "x**4 + 3*x - 4"
    # error, steps, roots = bisection(equation, 0, 3, precision=5, eps=1e-5, max_iterations=100)
    #equation = "x-4"
    gene=Bisection(equation, 4, 5,percision=7,tol=1e-4)
    for i in gene:
        print(i)
main()            