Lattice Gases
  * HPP-GAS: J. Hardy, Y. Pomeau, and O. de Pazzis. Molecular Dynamics of a classical lattice gas: Transport properties and time correlation functions. 1976
     Esta limitación no reproduce las ecuaciones de Navier-Stokes
  * FHP-GAS: U. Frisch, B. Hasslacher and Y. Pomeau. Lattice-gas automata for the Navier-Stokes equation. 1986.
     Se requere muchas, muchas celdas. Dio origen a los  Cellar Automata Machine: CAM6 y CAM8, programadas en FORTH. Ademas se tiene el problema de representar las probabilidades de la malla hexagonal, como una combinatoria de posibilidades.

Lattice-Boltzmann
   * La ecuacion de transporte de Boltzmann
      f:= f(t,x,v)
      la derivada total sería (indicial en i)
      Df/Dt = df/dt+df/dxi*dxi/dt+df/dvi*dvi/dt
            = df/dt+v.grad_x(f)+(F/m).grad_v(f)
	    = Omega(vi,vj)
   * What is a Lattice-Boltzmann BGK? (1988 McNamara, 1989 Higuera-Jimenez, 1992 Benzi)