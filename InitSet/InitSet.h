#ifndef _InitSet_H
#define _InitSet_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"
#include "header.h"

using namespace std;
using namespace mfem;

void SaveState(ParMesh *pmesh,ParGridFunction *Velocity);
void InitState(ParMesh *pmesh,ParGridFunction *u);
void InitMesh();

#endif
