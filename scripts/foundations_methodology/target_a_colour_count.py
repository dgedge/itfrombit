from fractions import Fraction as F
# Per-generation left-Weyl content with N_c colours, SM hypercharges:
# Q_L doublet (2 states, Y=1/6) xN_c ; u^c (Y=-2/3) xN_c ; d^c (Y=1/3) xN_c ;
# L_L doublet (2, Y=-1/2) ; e^c (Y=+1).  (no nu_R -> chiral)
def conds(Nc):
    Nc=F(Nc)
    doub_Y=[(F(1,6),2*Nc),(F(-1,2),2)]                 # (Y, multiplicity) for SU(2) doublets' members
    allY =[(F(1,6),2*Nc),(F(-2,3),Nc),(F(1,3),Nc),(F(-1,2),2),(F(1,1),1)]
    col  =[(F(1,6),2*Nc),(F(-2,3),Nc),(F(1,3),Nc)]     # coloured
    A1=sum(y*m for y,m in col)                           # [SU(Nc)]^2 U(1)
    A2=sum(y*(m//2 if isinstance(m,int) else m/2) for y,m in doub_Y)  # per-doublet
    A2=F(1,2)*sum(y*m for y,m in doub_Y)                 # = sum over doublets of Y
    A3=sum((y**3)*m for y,m in allY)                     # [U(1)]^3
    A4=sum(y*m for y,m in allY)                          # grav-U(1)
    return A1,A2,A3,A4
print("Nc :  A1      A2        A3        A4     anomaly-free?")
for Nc in range(1,9):
    a1,a2,a3,a4=conds(Nc)
    print(f"{Nc:2d} : {str(a1):6s} {str(a2):8s} {str(a3):8s} {str(a4):6s}   {(a1,a2,a3,a4)==(0,0,0,0)}")
print("\nclosed forms:  A2 = Nc/6 - 1/2  (->0 at Nc=3);  A3 = (3-Nc)/4  (->0 at Nc=3)")
print("=> anomaly cancellation FORCES Nc=3 (given SM hypercharges).")
