clc; clear;close;

exec('./ode1.sci');


//EI = 250e9;  // szywnosc na zginanie, Nmm^2
E = 210000 //N/mm^2
I = 1940e4 //mm^4
EI = E*I
_q = 10000;    // sila skupiona, N
q=_q/1000
L = 4000;    // dlugosc belki, mm

n=10;
A = zeros(n,n);
for i =1:n
    if i == 1 || i == n then
        A(i,i) = 1;
    else
        A(i,i-1) = 1;
        A(i,i)  = -2;
        A(i,i+1) = 1  ;
     end              
end

for i = 1:n-1
    b(i+1,:) = (1/(2*E*I)) * (L/(n-1))^2 * q * L*(i/(n-1)) * (L - L*(i/(n-1)));
end

w = linsolve(A,-b)
plot(0:L/(n-1):L,w,'b-','LineWidth',3);

disp(min(w))
//plot(x,y(2,:)*1000,'r-','LineWidth',3);
//plot(x,y(1,:),'b-','LineWidth',3);
////xlabel("x [mm]]");
//xlabel('$x \quad [\text{mm}]$','fontsize',4)
//ylabel('$w(x), \phi(x)$','fontsize',4);
//title("Porownanie");
//legend(['$\phi(x) \quad [\text{rad}]$';'$w(x) \quad [mm]$']);
//
//
