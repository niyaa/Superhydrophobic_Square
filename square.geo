// Gmsh project created on Tue Apr 28 16:30:39 2015

nx=8; //number of nodes in the x dir
ny=18; //number of nodes in the y dir

S=0.75;
a=1;




Lam = (2*Pi)/a;
 
Point(0)={0,1-S,0,1};
Point(1)={Lam/4,1-S,0,1};
Point(2)={Lam/4,1+S,0,1};
Point(3)={3*Lam/4,1+S,0,1};
Point(4)={3*Lam/4,1-S,0,1};
Point(5)={Lam,1-S,0,1};

Point(6)={0,-(1-S),0,1};
Point(7)={Lam/4,-(1-S),0,1};
Point(8)={Lam/4,-(1+S),0,1};
Point(9)={3*Lam/4,-(1+S),0,1};
Point(10)={3*Lam/4,-(1-S),0,1};
Point(11)={Lam,-(1-S),0,1};



Line(1) = {0, 1};
Line(2) = {1, 4};
Line(3) = {4, 5};
Line(4) = {5, 11};
Line(5) = {11, 10};
Line(6) = {10, 7};
Line(7) = {7, 6};
Line(8)=  {0,6};

Line(9) = {2, 1};
Line(10) = {2, 3};
Line(11) = {3, 4};
Line(12) = {10, 9};
Line(13) = {7, 8};
Line(14) = {8, 9};


Line(15) = {1, 7};
Line(16) = {4, 10};


//+
Line Loop(17) = {1, 15, 7, -8};
//+
Plane Surface(18) = {17};
//+
Line Loop(19) = {10, 11, -2, -9};
//+
Plane Surface(20) = {19};
//+
Line Loop(21) = {15, -6, -16, -2};
//+
Plane Surface(22) = {21};
//+
Line Loop(23) = {13, 14, -12, 6};
//+
Plane Surface(24) = {23};
//+
Line Loop(25) = {3, 4, 5, -16};
//+
Plane Surface(26) = {25};



ny1=Floor((0.8*ny)/2.8)+1;
ny=ny-ny1*2;

Transfinite Line {8, 15, 16, 4, 9, 11, 13, 12} = ny Using Progression 1;
Transfinite Line {1, 7, 10, 14,2,6, 3, 5} = nx Using Progression 1;
Transfinite Line {9, 11, 13, 12} = ny1 Using Progression 1;


Recombine Surface {18, 20, 22, 24, 26};
Transfinite Surface {18};
Transfinite Surface {20};
Transfinite Surface {22};
Transfinite Surface {24};
Transfinite Surface {26};

Physical Line(27) = {8};
Physical Line(28) = {4};
Physical Line(29) = {1,9,10,11, 3, 5,12,14,13, 7};
Physical Surface(30) = {18, 20, 22, 24, 26};

