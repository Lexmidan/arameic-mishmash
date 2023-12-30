function pp=poison_cg
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
eps = 1.e-7;
C = 2;
solve;

function so=solve                                                                      x                                                                                               
global dx N x y u un f r B eps p q
init;
n=1;
r=zeros(N,N);
q=r;
p=r;
% r = b - A x
for i=2:N-1
  for j=2:N-1
    r(i,j) = dx^2*f(i,j) - (u(i,j+1) + u(i,j-1) + u(i+1,j) + u(i-1,j) - 4*u(i,j));
  end
end
% q = A r; p=q
for i=2:N-1
  for j=2:N-1
    q(i,j) = (r(i,j+1) + r(i,j-1) + r(i+1,j) + r(i-1,j) - 4*r(i,j));
    p(i,j) = r(i,j);
  end
end

while max(max(abs(r))) > eps
  cgstep;
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

function lf=cgstep
global dx N x u un f r C p q
% alpha_k = |r^k|^2/(p^k,q^k)
rk = sum(sum(r .* r));
pq = sum(sum(p .* q));
alpha = rk/pq;
for i=2:N-1
  for j=2:N-1
      % x^k+1 = x^k + alpha_k p^k
      u(i,j) = u(i,j) + alpha*p(i,j);
      % r^k+1 = r^k - alpha_k q^k
      r(i,j) = r(i,j) - alpha*q(i,j);
  end
end
% beta_k = |r^k+1|^2/|r^k|^2
beta = sum(sum(r .* r))/rk;
for i=2:N-1
  for j=2:N-1
      % p^k+1 = r^k+1 + beta_k p^k
      p(i,j) = r(i,j) + beta*p(i,j);
      % q^k+1 = A r^k+1 + beta_k q^k
      q(i,j) = (r(i,j+1) + r(i,j-1) + r(i+1,j) + r(i-1,j) ...
                - 4*r(i,j)) + beta*q(i,j);
  end
end

function ub=ubc(x,y)
%  ub = x^2 + y^2;
  ub = x^4 + y^4;
%  ub = sin(x^2 + y^2);
  
function fb=ff(x,y)
%  fb = 4;
   fb = 12*(x^2+y^2);
%    fb = 4*(cos(x^2 + y^2) - sin(x^2 + y^2)*(x^2 + y^2));

  % 0.012962207194092
  % 0.024681879179625
  
