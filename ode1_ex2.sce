clear, close
exec ode1.sci;
clc

// parametry
EI = 250e9;            // sztywność na zginanie [Nmm^2]
P = 1000;              // wartość siły skupionej [N]
L = 2000;              // długość wspornika [mm]

// warunki początkowe
y = [0;0];  // w(0))=0, fi(0)=0
h = 250;    // krok
x = 0:h:L;  // x <0,L>

// równanie
function dydx = pochodna(x,y,P,EI,L)
  
  M = P*(L-x);  // równanie momentu (siła na końcu belki wspornikowej)
  
  dydx(1,1) = y(2);
  dydx(2,1) = -M/EI;
endfunction

y1 = euler1(y,x,pochodna); // do wyboru: euler2, midpoint, rk2, rk4
yr = ode(y,0,x,pochodna);  // rozw. referencyje

plot(x',y1(1,:)','b');
plot(x',yr(1,:)','r');
