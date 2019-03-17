clc; clear;close;

EI = 250e9;  // szywnosc na zginanie, Nmm^2
P = 1000;    // sila skupiona, N
L = 2000;    // dlugosc belki, mm
xp = 0.5;

y = [0;0]
x = 0:10:L


// dla i=1
i=1;
printf('----- KROK %i -----\n',i)
h = x(i+1) - x(i);

//k1 = f(x(i), y(:,i));

printf("Wartosc x(i)= %i \n",x(i));
printf("Wartosc y(:,i)= ");
disp(y(:,i));

//wchodzimy do funkcji dydx
//if
M = P*(L*xp - x(i));
printf("Wartosc M= %i \n",M);

//k1 = f(x(i), y(:,i));
dydx(1,1) = y(2);
printf("Wartosc dydx(1,1)= %i \n", dydx(1,1));

dydx(2,1) = -M/EI
printf("Wartosc dydx(2,1)= %i \n", dydx(2,1));

k1 = dydx
printf("Wartosc k1 = \n");
disp(k1);

//y(:,i+1) = y(:,i) + h*k1;
y(:,i+1) = y(:,i) + h*k1;
printf("Wartosc y(:,i+1)= \n");
disp(y(:,i+1));


////////////////////////////////////////////
// dla i=2
i=2;
printf('\n\n----- KROK %i -----\n',i)
h = x(i+1) - x(i);

//k1 = f(x(i), y(:,i));

printf("Wartosc x(i)= %i \n",x(i));
printf("Wartosc y(:,i)= ");
disp(y(:,i));

//wchodzimy do funkcji dydx
//if
M = P*(L*xp - x(i));
printf("Wartosc M= %i \n",M);

//k1 = f(x(i), y(:,i));
dydx(1,1) = y(2,i);
printf("Wartosc dydx(1,1)= %i \n", dydx(1,1));

dydx(2,1) = -M/EI
printf("Wartosc dydx(2,1)= %i \n", dydx(2,1));

k1 = dydx
printf("Wartosc k1 = \n");
disp(k1);

//y(:,i+1) = y(:,i) + h*k1;
y(:,i+1) = y(:,i) + h*k1;
printf("Wartosc y(:,i+1)= \n");
disp(y(:,i+1));


////////////////////////////////////////////
// dla i=3
i=3;
printf('\n\n----- KROK %i -----\n',i)
h = x(i+1) - x(i);

//k1 = f(x(i), y(:,i));

printf("Wartosc x(i)= %i \n",x(i));
printf("Wartosc y(:,i)= ");
disp(y(:,i));

//wchodzimy do funkcji dydx
//if
M = P*(L*xp - x(i));
printf("Wartosc M= %i \n",M);

//k1 = f(x(i), y(:,i));
dydx(1,1) = y(2,i);
printf("Wartosc dydx(1,1)= %i \n", dydx(1,1));

dydx(2,1) = -M/EI
printf("Wartosc dydx(2,1)= %i \n", dydx(2,1));

k1 = dydx
printf("Wartosc k1 = \n");
disp(k1);

//y(:,i+1) = y(:,i) + h*k1;
y(:,i+1) = y(:,i) + h*k1;
printf("Wartosc y(:,i+1)= \n");
disp(y(:,i+1));

