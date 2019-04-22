clear, close
exec ode1.sci;
clc

//             EI1                EI2   
//    |            q1   |                
//    ---------------------------------------
//   /\         L/2                  L/2    /\

// parametry
L = 4000;              // długość belki
E = 210000 //N/mm^2 - szytwnosc na zginanie
I = 1940e4 //mm^4
EI = E*I
_q = 10000;    // obciazenie liniowo rozlozne, kN/m
q=_q/1000


// warunki początkowe
h = 10;    // krok
x = 0:h:L;  // x <0,L>

function M = Momentum(x, q, L)
    RA = q*(3/8)*L
// równanie momentu 
  if x <= L/2 then
      M = -RA*x + q*x*x/2
  else
      M = -RA*x + q*(L/2)*(x-L/4)
  end
endfunction

function dydx = pochodna(x, y, q, EI, L)
//    RA = q*(3/8)*L
//// równanie momentu 
//  if x <= L/2 then
//      M = -RA*x + q*x*x/2
//  else
//      M = -RA*x + q*(L/2)*(x-L/4)
//  end
  M = Momentum(x,q,L);

  dydx(1,1) = y(2);
  dydx(2,1) = -M/EI;
endfunction

//for i=1:length(x)
//    M(i)=Momentum(x(i), q1, q2, L);
//end

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
