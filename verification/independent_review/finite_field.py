"""Quadratic field F_p[omega]/(omega^2+omega+1), for p=2 mod 3."""
class CyclotomicField:
    def __init__(self,p):
        assert p%3==2
        self.p=p;self.zero=(0,0);self.one=(1,0)
    def c(self,a):
        return (a%self.p,0) if isinstance(a,int) else (int(a[0])%self.p,int(a[1])%self.p)
    def add(self,a,b):
        a=self.c(a);b=self.c(b);return ((a[0]+b[0])%self.p,(a[1]+b[1])%self.p)
    def neg(self,a):a=self.c(a);return (-a[0]%self.p,-a[1]%self.p)
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):
        a=self.c(a);b=self.c(b)
        return ((a[0]*b[0]-a[1]*b[1])%self.p,(a[0]*b[1]+a[1]*b[0]-a[1]*b[1])%self.p)
    def inv(self,a):
        a=self.c(a);n=(a[0]*a[0]-a[0]*a[1]+a[1]*a[1])%self.p
        q=pow(n,-1,self.p);return ((a[0]-a[1])*q%self.p,-a[1]*q%self.p)
    def div(self,a,b):return self.mul(a,self.inv(b))
    def pow(self,a,n):
        if n<0:return self.pow(self.inv(a),-n)
        r=self.one;a=self.c(a)
        while n:
            if n&1:r=self.mul(r,a)
            a=self.mul(a,a);n>>=1
        return r
    def sum(self,x):
        r=self.zero
        for a in x:r=self.add(r,a)
        return r
    def det(self,A):
        a=[[self.c(q) for q in r] for r in A];n=len(a);d=self.one
        for j in range(n):
            k=next((i for i in range(j,n) if a[i][j]!=self.zero),None)
            if k is None:return self.zero
            if k!=j:a[k],a[j]=a[j],a[k];d=self.neg(d)
            pivot=a[j][j];d=self.mul(d,pivot)
            for i in range(j+1,n):
                c=self.div(a[i][j],pivot)
                for k in range(j+1,n):a[i][k]=self.sub(a[i][k],self.mul(c,a[j][k]))
        return d
    def sqrt(self,a):
        a=self.c(a)
        for i in range(self.p):
            for j in range(self.p):
                if self.mul((i,j),(i,j))==a:return (i,j)
        raise ArithmeticError('not a square in this quadratic field')
    def eval(self,poly,values):
        total=self.zero
        for monom,coefficient in poly.terms():
            term=self.c(int(coefficient))
            for a,n in zip(values,monom):term=self.mul(term,self.pow(a,n))
            total=self.add(total,term)
        return total
