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
Z
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

obj.. Z =e= sum({i,j},w(j)*y(i,j)) - Lamda*( sum({i},c(i)*x(i)) - B);

cons1..     sum({i},c(i)*x(i)) =l= B ;

cons2(i)..  sum({j},y(i,j)) =l= cap(i)*x(i);

cons3(i,j)..  u(i,j) =l= r/dis(i,j)*x(i);

cons4(i,j)..  y(i,j) =l= dem(j)*u(i,j);

cons5(j)..   sum({i},y(i,j)) =l= dem(j);

*********

Model Covering
/
obj
*cons1
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

Set Iter /1*50/;
Parameter
Lm(iter)
ZLR(iter,*);
Lm(iter)= (ord(iter)-1)/card(iter)*2;


LOOP(iter,

Lamda=Lm(iter);

Solve  Covering us MIP max Z;

ZLR(iter,'L')= Lamda;
ZLR(iter,'Z(L)')= Z.l;

)
;

Display
ZLR

;

Execute_Unload 'LR_IP401'




















