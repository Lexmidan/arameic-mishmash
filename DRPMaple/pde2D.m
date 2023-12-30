function pp=pde
global T t dt dx J x u un B C a b dc sip sim theta
% T je konecny cas
% dt,dx jsou casovy a prostorovy krok
% J je pocet bunek
% x(j) je prostorova souradnice j-teho bodu
% u(j) je reseni v bode x(j) na casove vrstve n
% un(j) je reseni v bode x(j) na casove vrstve n+1
% resime na intervalu (-B,B)
% C je CFL cislo
% resime rovnici u_t + a u_x + b u_y = 0
B = 10;
J = 100;
T = 20;
C = 0.5;
a = -1;
b = 2;
solve;

function so=solve
% resi PDR u_t = a u_x
global T t dt dx J x u un B C a

t=0;
init;
du = u;
while t < T
  lfstep;
  %  contour(un',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]);
  contour(un');
  colorbar;
%  axis([-B B 0 1]);
  drawnow;
  u = un;
  t = t + dt;
end

function ii=init
% inicializace, pocatecni podminky
global T dt dx J x u un B C a b

dx = 2*B/J;
x = zeros(J,1);
for j=1:J
  x(j) = -B + j*dx;
end
dt = C*dx/max(abs(a),abs(b));
u = zeros(J,J);
un = zeros(J,J);
for j=1:J
  for i=1:J
    if abs(x(j)) < pi & abs(x(i)) < pi
        %      u(i,j) = (1 + cos(x(j)) + cos(x(i)))/3;
      u(i,j) = 1;
    end
  end
end

function lf=lfstep
% jeden casovy krok Lax-Friendrichsovo schema
global T dt dx J x u un B C a b

for i=2:J-1
  for j=2:J-1
    un(i,j) = (u(i,j+1) + u(i,j-1) + u(i+1,j) + u(i-1,j))/4 ...
	      - a*dt*(u(i+1,j) - u(i-1,j))/(2*dx) ...
	      - b*dt*(u(i,j+1) - u(i,j-1))/(2*dx);  
  end
end
% periodicke okrajove podminky
for i=1:J
  un(i,1) = un(i,J-1);
  un(i,J) = un(i,2);
end
for j=1:J
  un(1,j) = un(J-1,j);
  un(J,j) = un(2,j);
end
