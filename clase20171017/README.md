The Lattice Boltzmann Method
Springer

Ecuación de evolución:

 Ni(x+Vi*dt, t+dt) - Ni(x,dt) = (1-rho)*[Nj(x,t)-Ni(x,t)]

Expansion de Chapman-Enskog
 1. Serie de Taylor
    dx/dt = dNi/dx*dt + dNi/dt*dt + dt^2*d^2Ni/dt^2 + ...
          =   dt*[d()/dt+Vi*d()/dxi]*Ni
	    + 1/2*dt^2[d^2()/dt^2 + 2*Vi*d^2()/(dt*dx) + Vi^2*d^2()/dx^2]*Ni

    => se define los operadores diferenciales
     Ni(x+Vi*dt, t+dt) = (dt*[d()/dt+Vi*d()/dx]+1/2*dt^2*[d()/dt+Vi*d()/dx]^2)*Ni
     		       = (dt*[d()/dt+Vi.Nabla()]+1/2*dt^2*[d()/dt+Vi.Nabla()]^2)*Ni
		       
 2. Expansion Perturbativa para aproximar los operadores diferenciales
    (1) d()/dx = epsilon*d()/dx1
    (2) d()/dt = epsilon^2*d()/dt2
    (3) Ni = Ni0+Ni1+Ni2  # igualar por ordenes
    	(dt*[d()/dt+Vi*d()/dx]+1/2*dt^2*[d()/dt+Vi*d()/dx]^2)*(Ni0+Ni1+Ni2) = (1-rho)*[(Nj0+Nj1+Nj2)-(Ni0+Ni1+Ni2)]

	=> zero orden, sabiendo que: Ni0+Nj0 = rho
	   Ni0 = rho/2
	=> primer orden, sabiendo que: Ni1+Nj1 = 0
	   Ni1 = dt*Vi/(2*(rho-1))*d()/dx1*Ni0
	=> segundo orden, sabiendo que: Ni2+Nj2 = 0
	   Ni2 =  [dt*d()/dt2*Ni0 + dt^2*Vi^2*1/(2*(rho-1))*d^2()/dx1^2*Ni0 + 1/2*dt^2*Vi^2*d^2()/dx1^2*Ni0]/(2*(rho-1))

	   Observe que, al sumar sobre los i: sabiendo que Vi = Ci*lambda/dt y dt=tau
	   d rho/dt - lambda^2/tau*(rho/(2*1-rho))*d^2 rho/dx^2 = 0
	   
 3. Cantidades Macroscopicas