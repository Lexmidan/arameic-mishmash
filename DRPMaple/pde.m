%---------------------------------------------------------------------
function de=pde
global T t dt dx J x u un B C a b dc sip sim theta
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
T = 820;
C = 2; 
a = -1; %nestabilni pro 1 u dopredneho, -1 u zpetneho

b = 5;
solve;
%---------------------------------------------------------------------
function so=solve
% resi PDR u_t = a u_x
global T t dt dx J x u un B C a ermax

t=0;
init;
du = u;
while t < T
  vedeni_tepla
%   dopredu;
%   dozadu;
% central;
%LaxFriedrich;
%Laxwendroff;
  plot(x,un,'x-');
%  axis([-B B 0 1]);
  drawnow; % pro urychleni vypoctu kreslit az kazdy desaty krok
%     pause; %aby bylo videt jak nestabilita vznika
  u = un;
  t = t + dt;
end
u_e = exp(-((x-a*t)./2).^2);
% 
% ue = zeros(J,1);
% for j = 1:J
%     ue(j) = exp(-((x(j)-a*t)/2)^2);
% end
err = sum(abs(u-ue))*dx


%---------------------------------------------------------------------
function ii=init
% inicializace, pocatecni podminky
global T dt dx J x u un B C a b

dx = 2*B/J;
x = zeros(J,1);
for j=1:J
  x(j) = -B + j*dx;
end
dt = C*dx^2/(2*b);
u = zeros(J,1);
un = zeros(J,1);
for j=1:J
  if x(j)<0
      u(j) = 1;
  end
end


%---------------------------------------------------------------------
function dd=dozadu
% jeden casovy krok  schema
global T dt dx J x u un B C a

for j=2:J-1
  un(j) = u(j) - a*dt*(u(j) - u(j-1))/dx;
end
un(1) = un(2);
un(J) = un(J-1);

function dd=dopredu
% jeden casovy krok  schema
global T dt dx J x u un B C a

for j=2:J-1
  un(j) = u(j) - a*dt*(u(j+1) - u(j))/dx;
end
un(1) = un(J-1);
un(J) = un(2);

function dd=central
% jeden casovy krok  schema
global T dt dx J x u un B C a

for j=2:J-1
  un(j) = u(j) - a*dt*(u(j+1) - u(j))/(2*dx);
end

function dd=LaxFriedrich
% jeden casovy krok  schema
global T dt dx J x u un B C a

for j=2:J-1
  un(j) = (u(j+1) + u(j-1))/2 - a*dt*(u(j+1) - u(j-1))/(2*dx);
end
% okrajove podminky
un(1) = un(2);
un(J) = un(J-1);

function dd=Laxwendroff
% jeden casovy krok  schema
global T dt dx J x u un B C a

for j=2:J-1
  un(j) = u(j) - a*dt*(u(j+1) - u(j-1))/(2*dx) + (a*dt)^2/2*(u(j+1) - 2*u(j) + u(j-1))/(dx^2);
  if u(j)<0
      u(j)=1
  end
end

function dd=vedeni_tepla
% jeden casovy krok  schema
global T dt dx J x u un B C b

for j=2:J-1
  un(j) = u(j) + b*dt*(u(j+1) - 2*u(j) + u(j-1)) / dx^2;
end
un(1) = un(2);
un(J) = un(J-1);


% % okrajove podminky
% un(1) = un(2);
% un(J) = un(J-1);


%Periodicke podminky:
% okrajove podminky
un(1) = un(J-1);
un(J) = un(2);
