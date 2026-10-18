"""Independent No. 98 root-disk and 1,925,120 torque-label replay.

The generator reconstructs Euler quintics and Lagrange phase quadratics from
the manuscript, certifies root disks by exact rational Rouché bounds, then
streams rational disks for every crossed root and signed torque sum.
Requires SymPy 1.14 for high-precision candidate centers; every acceptance
inequality itself is checked with fractions.Fraction.
"""
from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path

import sympy as sp


Q = 10**18
ROOT_RADIUS = F(1, 10**9)
OMEGA_SCALE = 10**24
CROSS_RADIUS = F(1, 10**4)
T = sp.symbols("T")
R = sp.symbols("R")


def cadd(a,b): return (a[0]+b[0],a[1]+b[1])
def csub(a,b): return (a[0]-b[0],a[1]-b[1])
def cmul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cscale(a,k): return (a[0]*k,a[1]*k)
def cabs1(a): return abs(a[0])+abs(a[1])
def cinf(a): return max(abs(a[0]),abs(a[1]))
def cinv(a):
    d=a[0]*a[0]+a[1]*a[1]
    if not d: raise ZeroDivisionError
    return (a[0]/d,-a[1]/d)


def qcenter(z):
    z=complex(z)
    return (F(round(z.real*Q),Q),F(round(z.imag*Q),Q))


def disk(center,radius=F(0)):
    return (center,radius)


def dadd(a,b): return (cadd(a[0],b[0]),a[1]+b[1])
def dsub(a,b): return (csub(a[0],b[0]),a[1]+b[1])
def dscale(a,k): return (cscale(a[0],F(k)),a[1]*abs(F(k)))
def dmul(a,b):
    return (cmul(a[0],b[0]),cabs1(a[0])*b[1]+cabs1(b[0])*a[1]+a[1]*b[1])
def drecip(a):
    lower=F(cinf(a[0]))
    if lower<=a[1]: raise ArithmeticError("reciprocal disk meets zero")
    return (cinv(a[0]),a[1]/(lower*(lower-a[1])))
def ddiv(a,b): return dmul(a,drecip(b))


def quantized_disk(z,radius=F(0)):
    c=qcenter(z)
    return (c,radius+F(1,Q))


def rouche_integer(poly, center, radius):
    """Certify exactly one root in a disk for integer coefficients."""
    c=center
    def eval_derivative(order):
        out=(F(0),F(0))
        for j in range(len(poly)-1,order-1,-1):
            coeff=math.factorial(j)//math.factorial(j-order)
            out=cadd(cmul(out,c),(F(poly[j]*coeff),F(0)))
        return out
    p0=eval_derivative(0); p1=eval_derivative(1)
    lhs=cabs1(p0)
    for k in range(2,len(poly)):
        lhs+=cabs1(eval_derivative(k))*radius**k/F(math.factorial(k))
    rhs=radius*cinf(p1)
    if not lhs<rhs:
        raise ArithmeticError(("Euler Rouché failed",float(lhs),float(rhs),center))
    return lhs,rhs


def euler_polynomial(masses,signs):
    a,b,c=masses; M=a+b+c; e12,e13,e23=signs
    p1=-(b+c*(1+T))/M
    p2=1-(b+c*(1+T))/M
    A1=-b*e12-c*e13/(1+T)**2
    A2=a*e12-c*e23/T**2
    poly=sp.Poly(sp.cancel(M*T**2*(1+T)**2*(p2*A1-p1*A2)),T)
    den,ints=sp.polys.polytools.primitive(poly.clear_denoms()[1])
    coeff=[int(ints.nth(i)) for i in range(ints.degree()+1)]
    while len(coeff)>1 and coeff[-1]==0: coeff.pop()
    return coeff


def euler_roots(masses):
    a,b,c=masses;M=a+b+c;rows=[]
    for signs in itertools.product((-1,1),repeat=3):
        coeff=euler_polynomial(masses,signs)
        if len(coeff)!=6: raise ArithmeticError(("Euler degree",masses,signs,coeff))
        poly=sp.Poly.from_list(list(reversed(coeff)),gens=T)
        if poly.discriminant()==0: raise ArithmeticError(("multiple Euler root",masses,signs))
        roots=sp.nroots(poly,n=45,maxsteps=300)
        rec=[]
        for root in roots:
            d=quantized_disk(complex(root),ROOT_RADIUS)
            lhs,rhs=rouche_integer(coeff,d[0],ROOT_RADIUS)
            rec.append({"signs":signs,"coefficients_ascending":coeff,
                        "center":d[0],"radius":ROOT_RADIUS,
                        "rouche_lhs_upper":lhs,"rouche_rhs_lower":rhs})
        # The Rouché disks for a single quintic must be disjoint.
        for i in range(len(rec)):
            for j in range(i):
                if cinf(csub(rec[i]["center"],rec[j]["center"]))<=2*ROOT_RADIUS:
                    raise ArithmeticError("Euler root disks overlap")
        rows.extend(rec)
    if len(rows)!=40: raise AssertionError(len(rows))
    return rows


def sqrt3_omega_disks():
    scale=OMEGA_SCALE
    n=math.isqrt(3*scale*scale)
    # sqrt(3) lies in [n/S,(n+1)/S); omega has half this interval width.
    imag=F(2*n+1,4*scale)
    rad=F(1,4*scale)
    omega=((F(-1,2),imag),rad)
    omega2=((F(-1,2),-imag),rad)
    return ((F(1),F(0)),F(0)),omega,omega2


def lagrange_root_disks(masses):
    a,b,c=masses; M=a+b+c
    one,omega,omega2=sqrt3_omega_disks()
    phase_disks=(one,omega,omega2)
    phase_expr=(sp.Integer(1),(-1+sp.sqrt(3)*sp.I)/2,(-1-sp.sqrt(3)*sp.I)/2)
    rows=[]
    for ids in itertools.product(range(3),repeat=3):
        qd=[phase_disks[i] for i in ids]
        qe=[phase_expr[i] for i in ids]
        q12,q13,q23=qd
        aa=q12
        bb=dsub(dadd(q12,q23),q13)
        cc=q23
        discr=sp.simplify((qe[0]+qe[2]-qe[1])**2-4*qe[0]*qe[2])
        exact_double=(discr==0)
        if exact_double:
            root_disk=dscale(ddiv(dscale(bb,-1),dscale(aa,2)),1)
            # Keep the two algebraic-multiplicity labels used in the paper.
            # Both labels denote the same exact shape and prove the same sums.
            roots=[root_disk,root_disk]
        else:
            approx_omega=complex(-.5,math.sqrt(3)/2)
            qnum=[1+0j,approx_omega,approx_omega.conjugate()]
            q1,q2,q3=(qnum[i] for i in ids)
            coeff_a=q1; coeff_b=q1+q3-q2; coeff_c=q3
            disc_num=coeff_b*coeff_b-4*coeff_a*coeff_c
            sr=complex(sp.N(sp.sqrt(sp.expand_complex(discr)),35))
            z1=(-complex(sp.N(qe[0]+qe[2]-qe[1],35))+sr)/(2*complex(sp.N(qe[0],35)))
            z2=(-complex(sp.N(qe[0]+qe[2]-qe[1],35))-sr)/(2*complex(sp.N(qe[0],35)))
            roots=[]
            for approx in (z1,z2):
                rd=quantized_disk(approx,ROOT_RADIUS)
                # Coefficient disks include the certified omega enclosure.
                point=(rd[0],F(0))
                val=dadd(dadd(dmul(aa,dmul(point,point)),dmul(bb,point)),cc)
                der=dadd(dscale(dmul(aa,point),2),bb)
                remainder=(cabs1(aa[0])+aa[1])*ROOT_RADIUS**2
                lhs=cabs1(val[0])+val[1]+remainder
                rhs=ROOT_RADIUS*(cinf(der[0])-der[1])
                if rhs<=0 or not lhs<rhs:
                    raise ArithmeticError(("Lagrange Rouché failed",ids,float(lhs),float(rhs),approx))
                roots.append((rd[0],ROOT_RADIUS))
            if cinf(csub(roots[0][0],roots[1][0]))<=2*ROOT_RADIUS:
                raise ArithmeticError(("Lagrange root disks overlap",ids))
        # Root polynomial is q12*r^2+(q12+q23-q13)r+q23.
        for root_disk in roots:
            rr=root_disk
            # W1,W2,W3 are mass-centred from W12=q12, W13=q13/(1+r).
            q12d,q13d,q23d=qd
            one_d=disk(((F(1),F(0))),F(0))
            center=dscale(dadd(dscale(q12d,b),dscale(ddiv(q13d,dadd(one_d,rr)),c)),F(1,M))
            # W1= -center; W2=q12-center; W3=q13/(1+r)-center.
            W=(dscale(center,-1),dsub(q12d,center),dsub(ddiv(q13d,dadd(one_d,rr)),center))
            rows.append({"phase_ids":ids,"double_root":exact_double,
                         "r_center":rr[0],"r_radius":rr[1],"W_disks":W})
    if len(rows)!=54: raise AssertionError(len(rows))
    # Repeated quadratic roots retain their two multiplicity labels.
    return rows


def positions_euler(root_record,masses):
    a,b,c=masses;M=a+b+c
    td=(root_record["center"],root_record["radius"])
    one=disk((F(1),F(0)),F(0))
    center=dscale(dadd(disk((F(b),F(0)),F(0)),dscale(dadd(one,td),c)),F(1,M))
    return (dscale(center,-1),dsub(one,center),dsub(dadd(one,td),center))


def cross_root_disk(zdisk,wdisk):
    product=dmul(zdisk,wdisk)
    u=dscale(drecip(product),-1)
    approx=complex(float(u[0][0]),float(u[0][1]))
    kcenter=qcenter(complex(math.sqrt(abs(approx))*math.cos(math.atan2(approx.imag,approx.real)/2),
                            math.sqrt(abs(approx))*math.sin(math.atan2(approx.imag,approx.real)/2)))
    rad=CROSS_RADIUS
    residual=cabs1(csub(cmul(kcenter,kcenter),u[0]))+u[1]+rad*rad
    linear=2*rad*cinf(kcenter)
    if not residual<linear:
        raise ArithmeticError(("cross-root Rouché failed",float(residual),float(linear),approx))
    return (kcenter,rad)


def torque_sweep(upper_positions,lower_positions,upper_masses,lower_masses):
    signs=tuple(itertools.product((-1,1),repeat=9))
    count=0; least_num=None; least_exact=None; least_where=None
    for ui,up in enumerate(upper_positions):
        for li,lo in enumerate(lower_positions):
            terms=[]
            radius=F(0)
            for a in range(3):
                for b in range(3):
                    root=cross_root_disk(up[a],lo[b])
                    weight=upper_masses[a]*lower_masses[b]
                    terms.append(cscale(root[0],weight))
                    radius+=root[1]*weight
            for phase in signs:
                sr=F(0);si=F(0)
                for eps,(re,im) in zip(phase,terms):
                    sr+=eps*re;si+=eps*im
                lower=max(abs(sr),abs(si))-radius
                if lower<=0:
                    raise ArithmeticError(("torque interval contains zero",ui,li,phase,float(lower)))
                scaled=int(lower*10**24)
                if least_exact is None or lower<least_exact:
                    least_exact=lower;least_num=scaled;least_where=(ui,li,phase)
                count+=1
    return {"labels":count,"proved_nonzero":count==1925120,
            "minimum_lower_bound_times_1e24":least_num,
            "minimum_lower_bound_exact":str(least_exact),
            "least_label_indices":least_where}


def main():
    upper_masses=(30,50,330);lower_masses=(1,2,4)
    upper_roots=euler_roots(upper_masses)
    lower_roots=euler_roots(lower_masses)
    lag_roots=lagrange_root_disks(lower_masses)
    upper_positions=[positions_euler(r,upper_masses) for r in upper_roots]
    lower_positions=[positions_euler(r,lower_masses) for r in lower_roots]
    lower_positions.extend(r["W_disks"] for r in lag_roots)
    torque=torque_sweep(upper_positions,lower_positions,upper_masses,lower_masses)
    if torque["labels"]!=1925120: raise AssertionError(torque)
    expected_margin=1052550837275735549440000000
    if int(F(torque["minimum_lower_bound_exact"])*10**24)!=expected_margin:
        raise AssertionError((torque["minimum_lower_bound_exact"],expected_margin))
    record={"masses_upper":upper_masses,"masses_lower":lower_masses,
            "body_order_m1_to_m6":lower_masses+upper_masses,
            "witness_role":"upper block 456 is Euler; lower block 123 is Euler or Lagrange",
            "upper_euler_root_disks":len(upper_roots),
            "lower_euler_root_disks":len(lower_roots),
            "lower_lagrange_root_disks":len(lag_roots),
            "root_radius":str(ROOT_RADIUS),"omega_radius":str(F(1,4*OMEGA_SCALE)),
            "cross_root_radius":str(CROSS_RADIUS),"torque_sweep":torque,
            "margin_scale":10**24,
            "minimum_lower_bound_exact":torque["minimum_lower_bound_exact"],
            "manuscript_printed_margin":1052550837275735549440000000,
            "manuscript_margin_normalization_defined_in_source":True,
            "margin_comparison_status":"matches floor of the exact minimum certified lower bound after multiplication by 10^24",
            "euler_disk_records":[{"signs":r["signs"],
                                  "coefficients_ascending":r["coefficients_ascending"],
                                  "center":[str(x) for x in r["center"]],
                                  "radius":str(r["radius"]),
                                  "rouche_lhs_upper":str(r["rouche_lhs_upper"]),
                                  "rouche_rhs_lower":str(r["rouche_rhs_lower"])}
                                 for r in upper_roots+lower_roots],
            "lagrange_shape_disks":[{"phase_ids":r["phase_ids"],"double_root":r["double_root"],
                                      "r_center":[str(x) for x in r["r_center"]],
                                      "r_radius":str(r["r_radius"])} for r in lag_roots]}
    out=Path(__file__).with_name("no98_interval_certificate.json")
    out.write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:record[k] for k in ("upper_euler_root_disks","lower_euler_root_disks",
                                              "lower_lagrange_root_disks","torque_sweep")},indent=2))
    print(f"saved {out}")


if __name__=="__main__":main()
