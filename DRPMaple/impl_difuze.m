%---------------------------------------------------------------------
function de=impl
global T t dt dx J x u un B C a dc sip sim theta b
% T je konecny cas
% dt,dx jsou casovy a prostorovy krok
% J je pocet bunek
% x(j) je prostorova souradnice j-teho bodu
% u(j) je reseni v bode x(j) na casove vrstve n
% un(j) je reseni v bode x(j) na casove vrstve n+1
% resime na intervalu (-B,B)
% C je CFL cislo
% resime rovnici u_t + a u_x = 0
B = 20;
J = 200;
T = 15;
C = 2;
a = 1;
b = 0.2;


solve;
%---------------------------------------------------------------------
function so=solve
% resi PDR u_t = a u_x
global T t dt dx J x u un B C a M b 

t=0;
init;
M = zeros(J,J);
M(1,1) = 1;
M(1,2) = -1;
%M(1,J-1) = -1;
for j=2:J-1
  M(j,j-1) = -a *dt/(2*dx) - b*dt/(dx^2);
  M(j,j)   = 1 + 2*b*dt/dx^2;
  M(j,j+1) = a *dt/(2*dx)  - b*dt/(dx^2);
end
M(J,J-1) = -1;
M(J,J) = 1;
%M(J,2) = -1;
M = inv(M);
while t < T
  imstep;
  plot(x,un);
%  axis([-B B 0 1]);
  drawnow;
  %pause
  u = un;
  t = t + dt;
end
%---------------------------------------------------------------------
function ii=init
% inicializace, pocatecni podminky
global T dt dx J x u un B C a b

dx = 2*B/(J);
x = zeros(J,1);
for j=1:J
  x(j) = -B + (j)*dx;
end
dt = C*dx/abs(a);
u = zeros(J,1);
un = zeros(J,1);

for j=1:J
    if x(j) < 0
      u(j) = 1;
%      u(j) = 1;
    end
end
%---------------------------------------------------------------------
function im=imstep
% jeden casovy krok implicitni schema
global T dt dx J x u un B C a M b 

uu = u;
uu(1) = 0;
uu(J) = 0;
% for j=2:J-1
%     uu(j) = u(j) + b*dt/(2*dx^2)*(u(j+1)-2*u(j)+u(j-1));
% end
un = M*uu;
%---------------------------------------------------------------------


