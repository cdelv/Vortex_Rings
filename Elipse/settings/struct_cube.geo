// Gmsh project created on Tue Apr 26 20:39:16 2022
SetFactory("OpenCASCADE");
//+
Alto = 1;
Nz = 3;
Ancho = 1;
Ny = 5; Ry = 1.0;
Largo = 2;//+
Nx = 7; Rx = 1.0;
Rectangle(1) = {0, 0, 0, Largo, Ancho, 0};//+
Transfinite Curve {3, 1} = Nx Using Progression Rx;
//+
Transfinite Curve {4, 2} = Ny Using Progression Ry;
Transfinite Surface {1} = {4, 3, 2, 1};
Recombine Surface{1};
//+
Extrude {0, 0, Alto} {
  Surface{1};
  Layers{Nz};
  Recombine;
}
//+
Physical Surface("Inlet", 1) = {5};
//+
Physical Surface("Side", 2) = {6, 4, 2, 1};
//+
Physical Surface("Outlet", 3) = {3};
//+
Physical Volume("Room", 4) = {1};
