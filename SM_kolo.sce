clear, close
exec ode1.sci;
clc

                                |         |
    |   q1   |                  |    q2   |
    ---------------------------------------
   /\    L/4        L/2            L/4    /\

// parametry
L = 4000;              // długość belki
E = 210000 //N/mm^2 - szytwnosc na zginanie
I = 1940e4 //mm^4
EI = E*I
_q = 10000;    // obciazenie liniowo rozlozne, kN/m
q1=_q/1000
q2=20;

// warunki początkowe
h = 10;    // krok
x = 0:h:L;  // x <0,L>

function dydx = pochodna(x, y, q1, q2, EI, L)
    RA = q1*(7/32)*L+q2*(L/32)
// równanie momentu 
  if x <= L/4 then
      M = -RA*x + q1*x*x/2
  elseif x > L/4 && x <= (3*L)/4 then
      M = -RA*x + q1*(L/4)*(x-L/8)
  else
      M = -RA*x + q1*(L/4)*(x-L/8) + q2*(x-(3*L)/4)*(x-(3*L)/4)*0.5
  end

  dydx(1,1) = y(2);
  dydx(2,1) = -M/EI;
endfunction

function M = Momentum(x, q1, q2, L)
    RA = q1*(7/32)*L+q2*(L/32)
// równanie momentu 
  if x <= L/4 then
      M = -RA*x + q1*x*x/2
  elseif x > L/4 && x <= (3*L)/4 then
      M = -RA*x + q1*(L/4)*(x-L/8)
  else
      M = -RA*x + q1*(L/4)*(x-L/8) + q2*(x-(3*L)/4)*(x-(3*L)/4)*0.5
  end
endfunction

for i=1:length(x)
    M(i)=Momentum(x(i), q1, q2, L);
end

// pierwszy strzał
fi1 = -1e-3;
//fi1 = 1;
y = rk4([0; fi1], x, pochodna); // do wyboru: euler1, euler2, midpoint, rk2

// błąd 
e1 = y(1,$) - 0;

// drugi strzał
fi2 = -1e-2;
//fi2 = 2;
y = rk4([0; fi2], x, pochodna); // do wyboru: euler1, euler2, midpoint, rk2

// błąd
e2 = y(1,$) - 0;

// skorygowany kąt obrotu
fi = fi2 - e2*(fi1 - fi2)/(e1 - e2)

//disp(y(1,$))
//disp(y(2,$))
disp(fi*180/%pi)

// rozwiązanie
y = rk4([0; fi], x, pochodna); // do wyboru: euler1, euler2, midpoint, rk2

plot(x,y(1,:))
//disp(min(y(1,:)))
//plot(x',M)
