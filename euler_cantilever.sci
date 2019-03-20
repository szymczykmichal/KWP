clc; clear;close;

exec('./ode1.sci');


//EI = 250e9;  // szywnosc na zginanie, Nmm^2
E = 210000 //N/mm^2
I = 1940e4 //mm^4
EI = E*I
P = 10000;    // sila skupiona, N
L = 4000;    // dlugosc belki, mm
xp = 1.0;

y=[0;0]
x = 0:10:L;

function dydx = f(x,y,P,L,EI,xp)
    // moment
    if x<=L*xp
        M = P*(L*xp - x);
    else
        M = 0;
    end
    // uklad rownan rozniczkowych    
    dydx(1,1) = y(2); //
    dydx(2,1) = -M/EI;
endfunction

y=euler1(y,x,f);

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

