import sympy as sp
def relativeError(x1,x2):
    if x1==0 or x2==0:
        return 17.1717
    return abs((x2-x1)/x2)
def check_validity(x):
    return  x == sp.zoo or not x.is_real or not x.is_finite
def FixedPointIteration(equ,startPoint,precision=5,tol=1e-4,mx_iter=50):
    iter=mx_iter
    startPoint=sp.N(startPoint,precision)
    root=startPoint
    x=sp.Symbol('x')
    diff=sp.N(sp.diff(equ),precision)
    while iter>=0:
        if(diff.subs(x,root)>1):
            raise Exception("method doesn't converge")
        iter-=1
        new_root=equ.subs(x,root)
        if(check_validity(new_root)):
            raise Exception("cannot be solved")
        relative=relativeError(root,new_root)
        if(relative<tol):
            yield {
            "iteration":mx_iter-iter,
            "x_i":root,
            "x_i+1":new_root,
            "RelativeError":relative
            }
            break
        if(relative==17.1717):
            relative="no relative error"
        yield {
            "iteration":mx_iter-iter,
            "x_i":root,
            "x_i+1":new_root,
            "RelativeError":relative
            }
            
        root=new_root
    if mx_iter==0:
        raise Exception("method didn't converge")    

def final_result(equ,startPoint,precision=5,tol=1e-4,mx_iter=50):
    generator=FixedPointIteration(equ,startPoint,precision,tol,mx_iter)
    finalResult=0
    for step in generator:
        finalResult=step["x_i+1"]
    return finalResult         
def main():
    equation="(2*x+3)**0.5"
    equation=sp.N(sp.sympify(equation),6)
    generator = FixedPointIteration(equation,4)
    for i in generator:
       print(i)
    print (final_result(equation,4) )  

main()
