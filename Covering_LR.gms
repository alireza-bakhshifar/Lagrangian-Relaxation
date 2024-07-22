********************************************************
************** Covering ********************************
********************************************************

Sets
I /i1*i20/
J /j1*j100/

*******

Parameters
B
cap(i)
c(i)
dem(j)
w(j)
dis(i,j)
r
;




B        = 1000;
cap(i)   = uniform(100,200);
c(i)     = uniform(100,200);
dem(j)   = uniform(20,50);
w(j)     = uniform(0,1);
dis(i,j) = uniform(5,50);
r        = 15;

Display
B
cap
c
dem
w
dis
r


Binary Variable
x(i)
u(i,j)
;

Positive Variable
y(i,j)
;
y.up(i,j)=dem(j);

Free Variable
ZL
;


Equations
obj
cons1
cons2
cons3
cons4
cons5
;

Scalar Lamda;

obj.. ZL =e= sum({i,j},w(j)*y(i,j)) - Lamda*( sum({i},c(i)*x(i)) - B);

cons2(i)..  sum({j},y(i,j)) =l= cap(i)*x(i);

cons3(i,j)..  u(i,j) =l= r/dis(i,j)*x(i);

cons4(i,j)..  y(i,j) =l= dem(j)*u(i,j);

cons5(j)..   sum({i},y(i,j)) =l= dem(j);

*********

Model Covering_LR
/
obj
cons2
cons3
cons4
cons5
/
;

Options
RESLIM=100
OPTCR=0
MIP=CPLEX
;


************** Setting

Set iter /1*10/

Parameters
MRE /0.01/
RE
Phi /0.0001/
Lamda /0/
Result(iter,*)
FC 'feasibility checker'
convergency
UB /5000/
LB /0/
;
convergency = NO;
**************


LOOP(iter$(NOT(convergency)),

Solve  Covering_LR us MIP max ZL;

Result(iter,'L') = Lamda;
UB=ZL.l;
Result(iter,'UB = ZLR(L)') = ZL.l;

** FC and LB
IF( sum({i},c(i)*x.l(i)) - B > 0,

FC=NO;

ELSE
FC=YES;
LB = sum({i,j},w(j)*y.l(i,j)) ;
)
;
**

Result(iter,'FC') = FC;
Result(iter,'LB') = LB;

RE = (UB - LB)/UB;

Result(iter,'RE') = RE;

*

IF( RE <= MRE ,
convergency = YES;
);


** update L

Lamda = max { 0 , Lamda  + Phi*(UB-LB)*(sum({i},c(i)*x.l(i)) - B)/abs((sum({i},c(i)*x.l(i)) - B)) }


)
*End of Loop
;

Display
Result
;



Execute_Unload 'LR_IP401'




















