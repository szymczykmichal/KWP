clc; clear;close;

exec('./ode1.sci');
x = 0:0.5:4;

deff('[t]=f(x,y)','t = -2*x^3+12*x^2-20*x+8.5');
deff('[t]=f0(x)', 't = -0.5*x^4 + 4*x^3 -10*x^2 + 8.5*x + 1')

//function t=f(x, y)
//    t = -2*x^3+12*x^2-20*x+8.5
//endfunction

y(1)=1;
y0(1)= f0(x(1));

y=euler1(y,x,f);

for i=1:length(x)-1
    y0(:,i+1) = f0(x(i+1));
end

plot(x,y,'ro-','LineWidth',2);
plot(x,y0,'bs:','LineWidth',2);
xlabel("x axis");
ylabel("y axis");
title("Porownanie");
legend(["Wartość dokładna";"Wartośc przybliżona"]);
//set(gca(), "data bounds", );
