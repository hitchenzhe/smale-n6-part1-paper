"""Exact F17(x,y) arithmetic using FLINT multivariate polynomial gcds."""
from flint import fmpz_mod_mpoly_ctx

CTX=fmpz_mod_mpoly_ctx.get(('x','y'),17,'lex')


class Rat:
    def __init__(self,n=0,d=None,reduce=True):
        if isinstance(n,Rat) and d is None:
            self.numer=n.numer;self.denom=n.denom;return
        n=CTX.constant(n) if isinstance(n,int) else n
        d=CTX.constant(1) if d is None else (CTX.constant(d) if isinstance(d,int) else d)
        if not d:raise ZeroDivisionError
        if not n:self.numer=CTX.constant(0);self.denom=CTX.constant(1);return
        if reduce:
            g=n.gcd(d);n=n//g;d=d//g
        self.numer=n;self.denom=d
    def __add__(self,other):
        other=Rat(other)
        g=self.denom.gcd(other.denom)
        a=self.denom//g;b=other.denom//g
        return Rat(self.numer*b+other.numer*a,a*other.denom)
    __radd__=__add__
    def __neg__(self):return Rat(-self.numer,self.denom,reduce=False)
    def __sub__(self,other):return self+-Rat(other)
    def __rsub__(self,other):return Rat(other)+-self
    def __mul__(self,other):
        other=Rat(other)
        g=self.numer.gcd(other.denom);h=other.numer.gcd(self.denom)
        return Rat((self.numer//g)*(other.numer//h),(self.denom//h)*(other.denom//g),reduce=False)
    __rmul__=__mul__
    def __truediv__(self,other):
        other=Rat(other)
        if not other.numer:raise ZeroDivisionError
        return self*Rat(other.denom,other.numer,reduce=False)
    def __rtruediv__(self,other):return Rat(other)/self
    def __pow__(self,n):
        if n<0:return Rat(self.denom,self.numer,reduce=False)**(-n)
        return Rat(self.numer**n,self.denom**n,reduce=False)
    def __bool__(self):return bool(self.numer)


def generators():return tuple(Rat(a) for a in CTX.gens())
