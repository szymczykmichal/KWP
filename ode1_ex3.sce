clear, close
exec ode1.sci;
clc

// parametry
EI = 250e9;            // sztywność na zginanie
P = 5000;              // wartość siły skupionej
L = 2000;              // długość belki

// warunki początkowe
h = 100;    // krok
x = 0:h:L;  // x <0,L>

function dydx = pochodna(x,y,P,EI,L)
  
  // równanie momentu (siła w środku rozpiętości belki wolnopodpartej)
  if x<=L/2 then
    M = -P/2*x;
  else
    M = -P/2*x - P*(L/2 - x);
  end
  
  dydx(1,1) = y(2);
  dydx(2,1) = -M/EI;
endfunction

// pierwszy strzał
fi1 = -1e-3;
y = rk4([0; fi1], x, pochodna); // do wyboru: euler1, euler2, midpoint, rk2

// błąd 
e1 = y(1,$) - 0;

// drugi strzał
fi2 = -1e-2;
y = rk4([0; fi2], x, pochodna); // do wyboru: euler1, euler2, midpoint, rk2

// błąd
e2 = y(1,$) - 0;

// skorygowany kąt obrotu
fi = fi2 - e2*(fi1 - fi2)/(e1 - e2)

// rozwiązanie
y = rk4([0; fi], x, pochodna); // do wyboru: euler1, euler2, midpoint, rk2

plot(x,y(1,:))
