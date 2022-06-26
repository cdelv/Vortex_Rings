//+
SetFactory("OpenCASCADE");

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
Height = 2;
Radius = 1;
Eps = 0.001;

Nb = 5+1;
Rb = 1;
Nc1 = 5+1;
Rc1 = 1.00;
Nz = 3;

//+
Point(1) = {0, Radius, Radius, 1.0};
//+
Point(2) = {0, 2*Radius, Radius, 1.0};
//+
Point(3) = {0, Radius, 2*Radius, 1.0};
//+
Point(4) = {0, 0, Radius, 1.0};
//+
Point(5) = {0, Radius, 0, 1.0};
//+
Point(6) = {0, Eps + Radius,Radius, 1.0};
Point(7) = {0, Radius, Eps + Radius, 1.0};
Point(8) = {0, -Eps + Radius, Radius, 1.0};
Point(9) = {0, Radius, -Eps +  Radius, 1.0};

Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 2};
//+
Circle(5) = {9, 1, 6};
//+
Circle(6) = {6, 1, 7};
//+
Circle(7) = {7, 1, 8};
//+
Circle(8) = {8, 1, 9};
//+
Line(9) = {4, 8};
//+
Line(10) = {5, 9};
//+
Line(11) = {2, 6};
//+
Line(12) = {3, 7};
//+
Transfinite Curve {2, 7} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {1, 6} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {4, 5} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {3, 8} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {9, 12, 11, 10} = Nb Using Progression Rb;
//+
Curve Loop(1) = {2, 9, -7, -12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 12, -6, -11};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {4, 11, -5, -10};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, 10, -8, -9};
//+
Plane Surface(4) = {4};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Recombine Surface {1, 2, 3, 4};
Extrude {Height, 0, 0} {
  Surface{4}; Surface{1}; Surface{2}; Surface{3};
  Layers{Nz};
  Recombine;
}
//+
Physical Surface("Inlet", 1) = {2, 1, 3, 4};
//+
Physical Surface("Side", 2) = {10, 5, 18, 14};
//+
Physical Surface("Outlet", 3) = {17, 13, 9, 20};
//+
Physical Volume("Sea", 4) = {2, 3, 4, 1};
