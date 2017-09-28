#include <iostream>
#include <cmath>
#include "Random64.h"

using namespace std;

int main() {
  Crandom ran64(1);
  int i;
  double mu=3,sigma=2;

  for (i=0;i<1e4;i++)
    cout << sigma*sqrt(-2*log(ran64.r()))*cos(2*M_PI*ran64.r())+mu << endl;

  return 0;
}
