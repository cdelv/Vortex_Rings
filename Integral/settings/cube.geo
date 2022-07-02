// Gmsh project created on Tue Apr 26 20:39:16 2022
SetFactory("OpenCASCADE");
//+
Alto = 1;
Ancho = 1;
Largo = 2;
Box(1) = {0, 0, 0, Largo, Ancho, Alto};
//+
Physical Surface("left", 1) = {1};
//+
Physical Surface("side", 2) = {6, 4, 3, 5};
//+
Physical Surface("right", 3) = {2};
//+
Physical Volume("environment", 4) = {1};
