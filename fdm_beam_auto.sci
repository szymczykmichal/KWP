clc; clear;close;

exec('./ode1.sci');


//EI = 250e9;  // szywnosc na zginanie, Nmm^2
E = 210000 //N/mm^2
I = 1940e4 //mm^4
EI = E*I
_q = 10000;    // sila skupiona, N
q=_q/1000
L = 4000;    // dlugosc belki, mm


A = [1 0 0 0 0;
     1 -2 1 0 0;
     0 1 -2 1 0;
     0 0 1 -2 1;
     0 0 0 0 1]

b = [0;
    (3*q*L^4)/(512*E*I);
    (q*L^4)/(256*E*I);
    (3*q*L^4)/(512*E*I);
    0;]

w = linsolve(A,-b)
plot(0:L/4:L,w)
//plot(x,y(2,:)*1000,'r-','LineWidth',3);
//plot(x,y(1,:),'b-','LineWidth',3);
////xlabel("x [mm]]");
//xlabel('$x \quad [\text{mm}]$','fontsize',4)
//ylabel('$w(x), \phi(x)$','fontsize',4);
//title("Porownanie");
//legend(['$\phi(x) \quad [\text{rad}]$';'$w(x) \quad [mm]$']);
//
//
