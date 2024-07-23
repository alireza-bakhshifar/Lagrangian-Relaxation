Sets
N                                        /1*9/
v set of vehicles /1*3/
;

Alias(i,j,k,N);
ALIAS(V,VP);
*************************************************************************************
PARAMETER

Q(v) max bare kamion k /1 20000, 2 20000,3 20000 /
;

************************************************************************************************
Parameters
Num / 3/
distance(i,j)

demand(i) /1 0,2 10,3 11,4 18,5 20,6 23,7 13,8 21,9 12 /


Cap(V)

;

Cap(V) =52;
table distance(i,j)

         1                  2                  3                  4                  5                  6                 7                  8                  9
1        0                  23.34523506        82.41965785        40.36087214        31.01612484        74.40430095       85.80209788        23                 34.71310992
2        23.34523506        0                  62.22539674        24.04163056        54.12023651        60.6052803        107.1867529        29.83286778        48.41487375
3        82.41965785        62.22539674        0                  42.54409477        112.6809656        92.26592004       168.2171216        71.51223672        110.4536102
4        40.36087214        24.04163056        42.54409477        0                  70.17834424        77.07788269       126.0357092        30.2654919         71.34423593
5        31.01612484        54.12023651        112.6809656        70.17834424        0                  94.9210198        56.56854249        44.68780594        31.82766093
6        74.40430095        60.6052803         92.26592004        77.07788269        94.9210198         0                 129.1897829        89.93886813        67.74215822
7        85.80209788        107.1867529        168.2171216        126.0357092        56.56854249        129.1897829        0                 101.1780609        64.28841264
8        23                 29.83286778        71.51223672        30.2654919         44.68780594        89.93886813       101.1780609        0                  57.42821606
9        34.71310992        48.41487375        110.4536102        71.34423593        31.82766093        67.74215822       64.28841264        57.42821606        0



parameter cost(v) /1 2,2 3,3 4/ ;

distance(i,i) = 1000;



display distance;
**********************************************************************************************

******************

Variables
Z1

;

Binary Variables
X(i,j,v) agar az i be j beravim
Zik(i,v)

Positive Variables
S1(j)
T(v,i) zamane residane kamion k be gere i
F1(v)
L(n,n) bare vasile dar kamane i be j
number(i,v)
;

Equations
Objective1

EQ1
EQ2
EQ3
EQ4
EQ5
EQ6
EQ7
EQ8
EQ9


;
PARAMETER Lamda(V);

positive variable ss;
Objective1                 .. Z1 =e=-(sum((i,j,v),X(i,j,v)*cost(v))
                              - SUM(V,Lamda(V)*(sum((i,j),demand(i)*X(i,j,v))- Cap(V))));

EQ1(i) $(ord(i) = 1)       .. Sum((v,j), X(i,j,v)) =e= Num;
EQ2(v)                     .. Sum(j, X('1',j,v)) =l= 1;
EQ3(v)                     .. Sum(I, X(I,'1',v)) =l= 1;
EQ4(i) $(ord(i) <> 1)      .. Sum((j,v), X(i,j,v)) =e= 1;
EQ5(k,v)                   .. Sum(i $(ord(k)<> ord(i)), X(i,k,v)) =e= Sum(j $(ord(k)<>ord(j)), X(k,j,v));
EQ6(i,j,v) $( ord(j)<> 1)  .. T(v,j) =g= T(v,i) + distance(i,j) - 10000 * (1 - X(i,j,v));
EQ7(i,j,v) $( ord(j)<> 1)  .. T(v,j) =l= T(v,i) + distance(i,j) + 10000 * (1 - X(i,j,v));
EQ8(i,j,v)                 .. F1(v)  =g= T(v,i) + distance(i,'1') - 10000 * (1 - X(i,'1',v));
EQ9(i,j,v)                 .. F1(v)  =l= T(v,i) + distance(i,'1') + 10000 * (1 - X(i,'1',v));

***********************************************
Model EXP /
ALL/ ;


************** Setting

Set iter /1*30/

Parameters
MRE /0.01/
RE
Phi /0.0001/
Lamda(V)
Result(iter,*)
FC 'feasibility checker'
convergency
UB /5000/
LB /0/
;
Lamda(V) =0;
convergency = NO;
**************


LOOP(iter$(NOT(convergency)),

Solve  EXP us MIP max Z1;
LOOP(V,
Result(iter,'L') = Lamda(V);
UB=Z1.l;
Result(iter,'UB = ZLR(L)') = Z1.l;
);
** FC and LB
LOOP(V,
IF( (sum((i,j),demand(i)*X.L(i,j,v)- Cap(V)))> 0,
FC=NO;


ELSE

FC=YES;
);
LB = -(sum((i,j),X.L(i,j,v)*cost(v))) ;
) ;
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

Lamda(V) = max { 0 , Lamda(V)  + Phi*(UB-LB)*((sum((i,j),demand(i)*X.L(i,j,v))- Cap(V)))/abs((sum((i,j),demand(i)*X.L(i,j,v))- Cap(V))) }


)
*End of Loop
;

Display
Result
;



Execute_Unload 'LR_IP401'




















