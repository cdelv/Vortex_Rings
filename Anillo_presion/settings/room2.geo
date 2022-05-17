// Gmsh project created on Tue Apr 26 20:39:16 2022
SetFactory("OpenCASCADE");
//+
//Dimensiones de la caja grande
Alto1 = 1;
Ancho1 = 1;
Largo1 = 3;

//Dimensiones del tubo
Radio2 = 0.25*Ancho1/2;
Largo2 = 0.2;
//Proporción del diametro contra el radio del tubo
R = 0.5;

Box(1) = {0, 0, 0, Largo1, Ancho1, Alto1};
Cylinder(2) = {-Largo2, Ancho1/2, Alto1/2, Largo2, 0, 0, Radio2, 2*Pi};
Circle(25) = {0, Ancho1/2, Alto1/2, R*Radio2, 0, 2*Pi};
//+
Rotate {{0, 1, 0}, {0, 0, Alto1/2}, Pi/2} {
  Curve{25}; 
}
//+
Curve Loop(13) = {25};
//+
Plane Surface(13) = {13};//+
BooleanDifference{ Surface{1};}{ Surface{13}; Delete; }//+
Physical Surface("Back", 1) = {9};
//+
Physical Surface("Side_small", 2) = {7};
//+
Physical Surface("Side", 3) = {4, 6, 3, 5};
//+
Physical Surface("Border_wall", 4) = {10};
//+
Physical Surface("Front", 5) = {2};
//+
Physical Volume("Space", 6) = {2, 1};
