clc; clear;close;

exec('./ode1.sci');


//EI = 250e9;  // szywnosc na zginanie, Nmm^2
E = 210000 //N/mm^2
I = 1940e4 //mm^4
EI = E*I
_q = 10000;    // obciazenie rozlozne [kN/m], 
q = _q/1000;  // obciazenie rozlozne [kN/m], 
L = 4000;    // dlugosc belki, mm


y=[0;0]
x = 0:20:L;

function dydx = f(x,y,q,L,EI)
    // moment
    if x<=0.5*L
        M = q*0.5*L*(0.75*L-x);
    else
        M = q*0.5*L*(0.75*L-x) + q*(x-0.5*L)*(x-0.5*L)*0.5;
    end
    // uklad rownan rozniczkowych    
    dydx(1,1) = y(2); //
    dydx(2,1) = -M/EI;
endfunction

y=rk4(y,x,f);

disp(y(1,$))
disp(y(2,$))
disp(y(2,$)*180/%pi)

plot(x,y(2,:)*1000,'r-','LineWidth',3);
plot(x,y(1,:),'b-','LineWidth',3);
//xlabel("x [mm]]");
xlabel('$x \quad [\text{mm}]$','fontsize',4)
ylabel('$w(x), \phi(x)$','fontsize',4);
title("Porownanie");
legend(['$\phi(x) \quad [\text{rad}]$';'$w(x) \quad [mm]$']);

