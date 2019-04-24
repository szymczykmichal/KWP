clear, close
exec ode1.sci;
clc

    |   q1           |               q2   |
    ---------------------------------------
   /\     E1                    E2        /\
    <-------L/2-------><------L/2--------->

// parametry
L = 4000;              // długość belki
E = 210000 //N/mm^2 - szytwnosc na zginanie
E2 = 270000 //N/mm^2 - szytwnosc na zginanie
I = 1940e4 //mm^4
EI = E*I
E2I = E2*I
q2=20; // obciazenie liniowo rozlozne, kN/m
q1=10;

// warunki początkowe
h = 10;    // krok
x = 0:h:L;  // x <0,L>

function dydx = pochodna(x, y, q1, q2, EI, E2I, L)
    RA = q1*(3/8)*L+q2*(L/8)
// równanie momentu 
  if x <= L/2 then
      M = -RA*x + q1*x*x/2
  else
      M = -RA*x + q1*(L/2)*(x-L/4) + q2*(x-L/2)*(x-L/2)*0.5
  end

  dydx(1,1) = y(2);
  if x <= L/2 then
      dydx(2,1) = -M/EI;
  else
      dydx(2,1) = -M/E2I;
  end
endfunction

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

plot(x,y(2,:)*1000,'r-','LineWidth',3);
plot(x,y(1,:),'b-','LineWidth',3);
//xlabel("x [mm]]");
xlabel('$x \quad [\text{mm}]$','fontsize',4)
ylabel('$w(x), \phi(x)$','fontsize',4);
title("Porownanie");
legend(['$\phi(x) \quad [\text{rad}]$';'$w(x) \quad [mm]$']);
