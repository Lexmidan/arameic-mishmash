function eu=eulers(x0,y0,x1,h)
% Eulerova metoda pro reseni ODR
global n x y

N = round((x1 - x0)/h);
x = zeros(N,1);
y= zeros(N,2);

x(1) = x0;
y(1, :) = y0;
n = 1;
while x(n) <= x1
  y(n+1, :) = y(n, :) + h*f(x(n),y(n, :));
  x(n+1) = x(n) + h;
  n = n+1;
%  h = min(h, 2/x(n));
end;
plot(x,y(:,1));

% pocitani chyby metody:
error=max(abs(y-exp(-x)))   % y-reseni metody; exp(-x) - analyticke reseni

function ff = f(x,y)
ff = [y(2),-0.2*y(2)-4.01*y(1)];


