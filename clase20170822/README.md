Métodos de integración

Error sistematico vs probabilistico?

* Primer orden no centrada: Euler

  x(t+dt) := x(t) + dt*dx/dt|t + O(dt^3)

* Segundo orden centrada: Leap-Frog

* Tercer orden: Verlet

  x(t+dt) := 2*x(t) - x(t-dt) + dt^2*F(t)/m
  x(-dt)  := x(0) - dt*v(0) + 1/2*dt^2*a(0)

* Segundo orden: Algoritmos simplecticos Velocity-Verlet(VV) y Position-Verlet(PV)

* Omelyan Optimized-(Velocity/Postion)-Verlet 3rd Order (OVV or OPV)

* Forest-Ruth (FR) 199a 4th Order

* VEFRL 2002 4th Order
  
* Algoritmos predictor-corrector