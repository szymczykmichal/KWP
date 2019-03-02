clear, close
exec ode1.sci;
clc

//warunki początkowe
L = 3;
y = 1;      // y(0))=1
h = 0.5;    // krok
x = 0:h:L;  // x <0,L>

// równanie
function dydx = pochodna(x,y)
  dydx(1,1) = 1 - x*y;
endfunction

// rozwiązanie
y1 = euler1(y,x,pochodna);
y2 = euler2(y,x,pochodna);
y3 = midpoint(y,x,pochodna);
y4 = rk2(y,x,pochodna);
y5 = rk4(y,x,pochodna);

yr = ode(y,x(1),x,pochodna);  // rozw. referencyje

plot(x',y1','b');
plot(x',y2','r');
plot(x',y3','m');
plot(x',y4','c');
plot(x',y5','g');
plot(x',yr','k:');

