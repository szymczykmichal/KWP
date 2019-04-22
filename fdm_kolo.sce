clc; clear;close;

exec('./ode1.sci');


//EI = 250e9;  // szywnosc na zginanie, Nmm^2
E = 210000 //N/mm^2
I = 1940e4 //mm^4
EI = E*I
_q = 10000;    // sila skupiona, N
q=_q/1000
L = 4000;    // dlugosc belki, mm

n=5;
A = zeros(n,n)

for i=1:n
    if i==1 || i==n then
        A(i,i) = 1
    else
        A(i,i-1) = 1
        A(i,i) = -2
        A(i,i+1) = 1
    end        
end
