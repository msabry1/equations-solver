import sympy as sp
def relativeError(x1,x2):
    if x1==0 or x2==0:
        return 17
    return abs((x2-x1)/x2)
def check_validity(x):
    return  x == sp.zoo or not x.is_real or not x.is_finite
def FormulaCalc(equ,a,b,x):
   fa=equ.subs(x,a)
   fb=equ.subs(x,b)
   return (a*fb-b*fa)/(fb-fa)
def FalsePosMethod(equ,a,b,precision=5,tol=1e-10,mx_iter=15):
    a=sp.N(a,precision)
    b=sp.N(b,precision)
    
    x=sp.symbols('x')
    if(equ.subs(x,a)==0):
        yield[a,"no relative error"]
        return
    if(equ.subs(x,b)==0):
        yield[b,"no relative error"]  
        return  
    if(equ.subs(x,a)>0):
        tmp=sp.N(b,precision)
        b=sp.N(a,precision)
        a=tmp
        
    prev=0

    if(check_validity(equ.subs(x,a)) or check_validity(equ.subs(x,b))):
        raise Exception("the function is not continuous at the given range")
    while mx_iter>0:
        mx_iter-=1
        mid=FormulaCalc(equ,a,b,x)
        if(check_validity(equ.subs(x,mid))):
            raise Exception("the function is not continuous at the given range")
        print(mid)
        if(equ.subs(x,a)*equ.subs(x,b)>1):
            raise Exception("this equation cannot be solved by the False Position method")
        if(equ.subs(x,mid)==0):
            relative=relativeError(mid,prev)
            if relative==17:
                relative="no relative error"
            yield [mid,relative]
            break
        if(equ.subs(x,mid)>0):
            b=mid
        else:
            a=mid
        if(relativeError(prev,mid)<tol):
            break    
        relative=relativeError(mid,prev)
        if relative==17:
            relative="no relative error"
        yield [mid,relative]
        prev =mid
def main():
    x=sp.symbols('x')
    precision=5
    equation=sp.N(sp.sympify("x**2+6*x+4"),precision)
    g=equation-x
    generator=FalsePosMethod(g,-2,5)
    for i in generator:
        print(i)
main()  