#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
using namespace std;

int main(void) {
  vector3D a,b,c;
  a.cargue(1,2,3);
  b.cargue(-2,-1,-3);
  c = a+b;
  b = b*3;
  b = 3*b;
  c = b^a;
  c.show();
  
  return 0;
}
