#include <iostream>
#include <cmath>
#include "Random64.h"

using namespace std;

const int N = 1000;
class Agente {
private:
  double dinero;
public:
  friend class Mercado;
  void Inicie(double dinero0) {dinero = dinero0;};
  double GetDinero(void) {return dinero;};
};

class Mercado {
private:
  
public:
  void HagaTransaccionEntre( Agente& ag1, Agente& ag2, Crandom& ran64 );
};

void Mercado::HagaTransaccionEntre( Agente& ag1, Agente& ag2, Crandom& ran64 ) {
  double gan1, gan2, res;
  // chitrapati economic model
  res = 2*ran64.r();
  gan1 = res-1;
  gan2 = (2-res)-1;
  if (ag1.dinero+gan1>0 && ag2.dinero+gan2>0) {
    ag1.dinero+=gan1;
    ag2.dinero+=gan2;
  }
}

int main() {
  Mercado Corabastos;
  Agente agentes[N];
  Crandom ran64(1);

  int i,j,t,Ntransacciones=N*10000;
  for (i=0;i<N;i++)
    agentes[i].Inicie(10);
  // Inicie los agentes con 10miles de pesos
  for (t=0;t<Ntransacciones;t++) {
    i = (int) N*ran64.r();
    j = (int) (N-1)*ran64.r(); if (j>=i) j++;
    // escojo dos agentes al azar e interactuan
    Corabastos.HagaTransaccionEntre(agentes[i],agentes[j],ran64);
  }
  // imprimo el dinero de todo
  for (i=0;i<N;i++)
    cout << agentes[i].GetDinero() << endl;
  
  return 0;
}
