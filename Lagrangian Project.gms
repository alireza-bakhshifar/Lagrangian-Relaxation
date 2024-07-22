Set
i           /1*4/
j          /1*4/
;

Parameter
kmin
kmax
;
kmin=1;
kmax=10;

Table hb(i,j)
    1   2   3   4
1   6   24  40  12
2   4.5 18  30  9
3   6   24  40  12
4   7.5 30  50  15;

Table R(i,j)
    1   2   3   4
1   18  70  4   14
2   75  8   45  22
3   32  10  55  30
4   28  9   5   7;

Table AB(i,j)
    1       2        3        4
1   0.0003  0.00012  0.0002   0.0006
2   0.0045  0.00018  0.0030   0.0090
3   0.0030  0.00012  0.0002   0.0006
4   0.0060  0.00240  0.0040   0.0120
;

Parameter As(i) /1 2100,2 8400,3 14000,4  4200/;
Parameter hs(i) /1 4.2,2 16.8,3 28.0,4  8.4/;
Parameter d(j) /1 30,2 150,3 300,4  250/
cob
cos
chb
chs
csb
css
hb(i,j)
II(i,j)
vi(i)
vj(j)
;
cob=uniform(5,50);
cos=uniform(5,50);
chb=uniform(5,50);
chs=uniform(5,50);
csb=uniform(5,50);
css=uniform(5,50);
hb(i,j) =uniform(5,50);
II(i,j) =uniform(200,10000);
vi(i)  =uniform(2,5);
vj(j) =uniform(100000,500000);

Free variable
Z1
Z2
;

Positive variable
K(I,J)
KB(J)
KS
T(I,J)
TJ(J)
Q(I,J)
FC
VC
;
k.lo(i,j)=1;

EQUATION
Obj1
Obj2
const1
const2
const3
const4
const5
const6
const7
const8
const9
const10
const11
const12
const13
const14
;
*******************************************************************************
Set iter /1*5/
;
Parameter
Iteration
Feas
Result(iter,*)
Lambda /0/
Convergency
g
Feasibilty(iter)
MaxRE /0.001/
phi /2/
;
Convergency = NO;
Scalar
UB  /+INF/
LB  /-INF/
;


Obj1  ..  Z1=e=-(sum(j,kb(j))-Lambda*(sum(j,Tj(j)*sum(i,Vi(i)*K(i,j)*R(i,j))-vj(j))));
Obj2  ..  Z2=e=ks;

const1(i,j)  .. T(i,j)=e=K(i,j)*Tj(j);
const2(i,j)  .. Q(i,j)=e=T(i,j)*R(i,j);
const3       .. cob=e=sum((i,j),AB(i,j)/0.01+T(i,j));
const4       .. sum((i,j),(hb(i,j)*Q(i,j))/2)=g=chb;
const5       .. csb=e=sum((i,j),T(i,j)*R(i,j)*d(j)*vc);


const8  ..cos=e=sum((i,j),As(i)/0.1+Tj(j));

const9  ..chs=e=sum((i,j),(hs(i)*Q(i,j)/2));


const6(j)  .. kb(j) =e= sum(i,((Ab(i,j)/K(i,j)*Tj(j))
                   +(hb(i,j)*K(i,j)*Tj(j)/2) + Tj(j)*K(i,j)*R(i,j)*d(j)*vc )) ;
const7(j)  .. Tj(j) =g= sqrt( 2*sum(i,AB(i,j))/sum(i,(K(i,j)**2)
                   *R(i,j)*(hb(i,j)+2*d(j)*vc) ));
const10    .. ks    =e= sum((i,j),(As(i)/K(i,j)*Tj(j))
                   + FC/K(i,j)*Tj(j)+(hs(i)*K(i,j)*R(i,j)*Tj(j)/2)) ;


const11(i,j)  .. Tj(j)*K(i,j)*R(i,j)=l=II(i,j);
const12(j)    .. Tj(j)*sum(i,Vi(i)*K(i,j)*R(i,j))=l=vj(j);
const13(i,j)  .. k(i,j)=g=1;
const14(i,j)  .. k(i,j)=l=10;





model scm
/
Obj1
Obj2
const1
const2
const3
const4
const5
const6
const7
*const8
*const9
const10
const11
*const12
const13
const14
/
;

option
optcr=0
minlp=bonmin
;

solve scm us minlp min Z1;

display
Z1.l
Z2.l
KB.l
KS.l
T.l
TJ.l
Q.l
FC.l
VC.l;


********************************************

Scalars
counter     /0/
Max_Counter /3/
;

Loop(iter$(Not(Convergency)),


if(counter=Max_Counter,
phi=phi/2;
counter=0;
);


Solve scm us minlp MAX Z1 ;

Result(iter,'UB')= Z1.l;

if ( Z1.l < UB ,
UB= Z1.l;

else
counter=counter+1;
)
;

g = sum(j,Tj.l(j)*sum(i,Vi(i)*K.l(i,j)*R(i,j))-vj(j)) ;

if ( g <= 0 ,
Feasibilty(iter) = 1;
)
;

Result(iter,'Feasibility')= Feasibilty(iter);

if ( Feasibilty(iter)=1,

Result(iter,'LB')= -(sum(j,kb.l(j)));

if ( Result(iter,'LB') > LB,

LB= Result(iter,'LB');

);

);


**********************Stop Critertia********************************************
if( (UB-LB)/UB < MaxRE ,
Convergency=YES
)
;

Iteration=ord(iter);

feas = Feasibilty(iter) ;

Display
"Iteration"
Iteration
"-------------------"
phi
counter
"-------------------------"
feas

Lambda
g

;


*************************LAMBDA UPDATE******************************************
Lambda = max ( 0 , Lambda + [phi*(Z1.l-LB)/abs(g)]*(sum(j,Tj.l(j)*sum(i,Vi(i)*K.l(i,j)*R(i,j))-vj(j)))) ;


)
;

UB=-UB;
LB=-LB;

Display
UB
LB
;

Execute_Unload "Project"