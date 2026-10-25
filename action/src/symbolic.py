"""Bianchi identity and propagated constraint, independently with formal symbols."""
from common import *
import sympy as s
def run():
    H,K,hp,xi,S,Sd,v2,u2,r,qUg,Hd=s.symbols('H K hp xi S Sd v2 u2 r qUg Hd')
    R=6*(Hd+2*H**2+K)
    Sdd=2*u2-3*H*Sd-2*qUg-2*xi*R*S
    rhod=-3*H*(v2+u2)-4*H*r-xi*R*Sd/2
    psid=2*H*(Hd-K)
    Cdot=3*hp*psid-rhod-3*xi*(psid*S+(H**2+K)*Sd+Hd*Sd+H*Sdd)
    den=2*(hp-xi*S)+12*xi**2*S
    num=2*(hp-xi*S)*K-v2+(2*xi-1)*u2-4*r/3-2*xi*qUg-4*xi*H*Sd-12*xi**2*S*(2*H**2+K)
    result=s.factor(Cdot.subs(Hd,num/den))
    rho_plus_p=v2+u2+4*r/3+xi*(2*(K-Hd)*S+H*Sd-Sdd)
    accel=s.factor((2*hp*(Hd-K)+rho_plus_p).subs(Hd,num/den))
    pressure_identity=s.factor((rhod+3*xi*(psid*S+(H**2+K)*Sd+Hd*Sd+H*Sdd)+3*H*rho_plus_p).subs(Hd,num/den))
    save('symbolic.json',dict(constraint_derivative=str(result),acceleration_residual=str(accel),
        stress_continuity_residual=str(pressure_identity),pass_exact=bool(result==0 and accel==0 and pressure_identity==0)))
    print('SYMBOLIC',result,accel,pressure_identity)
if __name__=='__main__':run()
