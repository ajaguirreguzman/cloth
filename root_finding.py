
def fxANDdfx(x,f,xg):
    x1=x[x<xg][-1]
    x2=x[x>xg][0]
    f1=f[x<xg][-1]
    f2=f[x>xg][0]
    m=(f2-f1)/(x2-x1)
    b=f1-m*x1
    return m*xg+b,m

def newton_raphson(x,f,xg): # xg is initial guess
    tol=1e-4
    i=1
    while i<10000:
        fxg,dfxg=fxANDdfx(x,f,xg)
        print('%i \t %.8f \t %.8f' % (i,xg,fxg))
        #if abs(fxg)<tol: break
        if abs(fxg/dfxg)<tol: break
        xg=xg-fxg/dfxg
        i+=1
    return xg
