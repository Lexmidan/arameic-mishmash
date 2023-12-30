function pp=poison
global dx N x y u un f r B eps C
% dx je prostorovy krok
% N je pocet bunek
% x(j) je prostorova souradnice j-teho bodu
% u(j) je reseni v bode x(j) na iteraci n
% un(j) je reseni v bode x(j) na iteraci n+1
% resime na intervalu (-B,B)x(-B,B)
% resime rovnici  u_xx + u_yy = f
B = 1;
N = 80;
eps = 1.e-6;
C = 2;
solve;

function so=solve
% resi PDR u_t = a u_x
global dx N x y u un f r B eps
init;
n=1;
while max(max(abs(r))) > eps
  step;
%  pcolor(un);
%  colorbar;
%  axis([-B B 0 1]);
%  drawnow;
  rmax = max(max(abs(r)));
  n=n+1;
end
for i=1:N
  for j=1:N
    r(i,j) = u(i,j) - ubc(x(i),y(j));
  end
end
n=n
format long
maxerr=max(max(abs(r)))
l1err = sum(sum(abs(r)))*dx^2
format short
surf(x,y,u)
colorbar

function ii=init
% inicializace, okrajove podminky
global dx N x y u un f r B 

dx = 2*B/(N-1);
x = zeros(N,1);
for j=1:N
  x(j) = -B + (j-1)*dx;
end
y = x;
u = zeros(N,N);
un = zeros(N,N);
f = zeros(N,N);
r = zeros(N,N);
r(2,2)=1;
for j=1:N
  u(1,j) = ubc(x(1),y(j));
  u(N,j) = ubc(x(N),y(j));
  u(j,1) = ubc(x(j),y(1));
  u(j,N) = ubc(x(j),y(N));
  for i=1:N
    f(i,j) = ff(x(i),y(j));
  end
end
un = u;

function lf=step
global dx N x u un f r C
omega = 2/(1 + C*dx);
for i=2:N-1
  for j=2:N-1
% jedna Jacobiho iterace
%    un(i,j) = (u(i,j+1) + u(i,j-1) + u(i+1,j) + u(i-1,j) - ...
%     f(i,j)*dx^2)/4;  % Jacobi
    un(i,j) = (u(i,j+1) + un(i,j-1) + u(i+1,j) + un(i-1,j) - ...
     f(i,j)*dx^2)/4;  % Gauss-Seidel
%    un(i,j) = u(i,j) + omega* ...
%              ((u(i,j+1) + un(i,j-1) + u(i+1,j) + un(i-1,j) - f(i,j)*dx^2)/4 ...
%               - u(i,j));  % SOR
  end
end
u=un;
for i=2:N-1
  for j=2:N-1
    r(i,j) = u(i,j+1) + u(i,j-1) + u(i+1,j) + u(i-1,j) -4*u(i,j) - f(i,j)*dx^2;
  end
end

function ub=ubc(x,y)
%  ub = x^2 + y^2;
  ub = x^4 + y^4;
%  ub = sin(x^2 + y^2);
  
function fb=ff(x,y) %prava strana Poiss rce tvaru uxx+uyy=f(x,y)
%  fb = 4;
   fb = 12*(x^2+y^2);
%    fb = 4*(cos(x^2 + y^2) - sin(x^2 + y^2)*(x^2 + y^2));

  % 0.012962207194092
  % 0.024681879179625
  
