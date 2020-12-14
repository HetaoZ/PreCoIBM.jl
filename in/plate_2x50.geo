lc = 0.00025 ;
width = 0.002 ;
height = 0.05 ;

Point(1) = {0,0,0,lc};
Point(2) = {width,0,0,lc};
Point(3) = {width,height,0,lc};
Point(4) = {0,height,0,lc};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Line Loop(9) = {5,6,7,8};
Plane Surface(10) = {9};
Recombine Surface(10);

Physical Line(101) = {7};
Physical Line(102) = {5};
Physical Line(103) = {8};
Physical Line(104) = {6};

Physical Surface("boundary") = {10};