// kolekcja funkcji do rozwiązywania ODE dowolnego rzedu, metodami eulera i rk 

function y = euler1(y,x,f)
  for i=1:length(x)-1
    h = x(i+1) - x(i);
    k1 = f(x(i), y(:,i));
    y(:,i+1) = y(:,i) + h*k1;
  end
endfunction

function y = euler2(y,x,f)
  for i=1:length(x)-1
    h = x(i+1) - x(i);
    k1 = f(x(i), y(:,i));
    k2 = f(x(i)+h, y(:,i)+h*k1);
    y(:,i+1) = y(:,i) + h*k2;
  end
endfunction

function y = midpoint(y,x,f)
  for i=1:length(x)-1
    h = x(i+1) - x(i);
    k1 = f(x(i), y(:,i));
    k2 = f(x(i)+h/2, y(:,i)+h/2*k1);
    y(:,i+1) = y(:,i) + h*k2;
  end
endfunction

function y = rk2(y,x,f)
  for i=1:length(x)-1
    h = x(i+1) - x(i);
    k1 = f(x(i), y(:,i));
    k2 = f(x(i)+h, y(:,i)+h*k1);
    y(:,i+1) = y(:,i) + h*(k1+k2)/2;
  end    
endfunction

function y = rk4(y,x,f)
  for i = 1:length(x)-1
    h = x(i+1) - x(i);
    k1 = f(x(i), y(:,i));
    k2 = f(x(i)+h/2, y(:,i)+h/2*k1);
    k3 = f(x(i)+h/2, y(:,i)+h/2*k2);
    k4 = f(x(i)+h, y(:,i)+h*k3);
    k = (k1 + 2*k2 + 2*k3 + k4)/6;
    y(:,i+1) = y(:,i) + h*k;
  end
endfunction
