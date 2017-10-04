%-------------------------------------------------------------------------%
%                           INICIALIZACION
% Dimensiones, cargas, potenciales, etc.
%-------------------------------------------------------------------------%

% Dimensiones
m_mm = 1/1000;
h = 10*m_mm;
a = 5*h;
b = a;
w = 8.25*m_mm;
t = 0.1*m_mm;
s = w;
e0 = 8.854e-12;
e =10*e0;

alto = a;
ancho = 2*b;

% Potenciales, Corriente
Vd = 5;
I = 0.001;

% Cada cuadro será de 0.1 mm x 0.1 mm

resolucion = 0.1*m_mm;

Nx = floor(ancho/resolucion) + 1 ;    % Numero de grillas x
Ny = floor(alto/resolucion) + 1;      % Numero de grillas y
mpx = ceil(Nx/2);              % Punto medio de x
mpy = ceil(Ny/2);              % Punto medio de y


Ni = 1000;      % Numero de iteraciones para Poisson

V = zeros(Nx,Ny);  % Matriz de potencial (voltaje)
A = zeros(Nx,Ny);
J = zeros(Nx,Ny);

T = 0;      % Pared superior
B = 0;      % Pared Inferior
L = 0;      % Pared Izquierda
R = 0;      % Pared Derecha


%-------------------------------------------------------------------------%
% Inicializando condiciones de borde
%-------------------------------------------------------------------------%

V(1,:) = L;
V(Nx,:) = R;
V(:,1) = B;
V(:,Ny) = T;

% A(1,:) = L;
% A(Nx,:) = R;
% A(:,1) = B;
% A(:,Ny) = T;

%-------------------------------------------------------------------------%
% Inicializando potenciales en las esquinas
%-------------------------------------------------------------------------%

V(1,1) = 0.5*(V(1,2)+V(2,1));
V(Nx,1) = 0.5*(V(Nx-1,1)+V(Nx,2));
V(1,Ny) = 0.5*(V(1,Ny-1)+V(2,Ny));
V(Nx,Ny) = 0.5*(V(Nx,Ny-1)+V(Nx-1,Ny));

%-------------------------------------------------------------------------%

largo_placa = 2*w/resolucion;
largo_mitad = floor(largo_placa/2);

posicion_x = floor((s/2 + w)/resolucion);
posicion_y =-floor((3*h/2)/resolucion);

pp1 = mpy+posicion_y;
pp2 = mpy+posicion_y;


for z = 1:Ni        % Numero de iteraciones
    for i = 2:Nx-1
        for j = 2:Ny-1
            
            
            V(mpx-largo_mitad+posicion_x:mpx+largo_mitad+posicion_x,pp1)=-Vd;
            V(mpx-largo_mitad-posicion_x:mpx+largo_mitad-posicion_x,pp2)=+Vd;
            
            if j == pp1 && ((i < mpx-largo_mitad-posicion_x) ...
                    || (i > mpx-posicion_x+largo_mitad ...
                    && i< mpx +posicion_x-largo_mitad) ...
                    || (i>mpx+posicion_x+largo_mitad))
                   

                    V(i,j) = ((V(i+1,j)/4)+(V(i-1,j)/4)+(e0/(2*(e0+e)))*V(i,j+1)+(e/(2*(e0+e)))*V(i,j-1));
            else
                V(i,j)=0.25*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1));
            end
            
        end
    end
end


%Transponer, para orientacion x-y apropiada

V = V';

[Ex,Ey] = gradient(V);
Ex = -Ex;
Ey = -Ey;


% Magnitud del campo electrico

E = sqrt(Ey.^2+Ex.^2);
[Dx,Dy] = gradient(V);
Dx = -Dx;
Dy = -Dy;
for j=1:Ny-1
    if j<ceil((h+1)/resolucion)
        Dx(j,:) = Ex(j,:);
        Dy(j,:) = Ey(j,:)/e;
    end
end

D = sqrt(Dx.^2+Dy.^2);
    


% Potencial magnético y campo magnético

for z = 1:Ni        % Numero de iteraciones
    for i = 2:Nx-1
        for j = 2:Ny-1
            if j == pp1 && i > mpx-largo_mitad-posicion_x && i < mpx+largo_mitad-posicion_x
                J(i,j) = -I/(t*w); 
                A(i,j) = 0.25*(A(i+1,j)+A(i-1,j)+A(i,j+1)+A(i,j-1)-J(i,j)*h*h);
            else
                A(i,j)=0.25*(A(i+1,j)+A(i-1,j)+A(i,j+1)+A(i,j-1));
            end
            
        end
    end
end


A = A';

[dAx,dAy] = gradient(A);
Bx = dAy;
By = -dAx;
B = sqrt(Bx.^2+By.^2);

x = ((1:Nx))*resolucion;
y = ((1:Ny))*resolucion;

paso = 5;

% Contorno para el potencial electrico
figure(1)
clf
colormap('jet')
contour_range_V = -5:0.1:5;
contour(x,y,V,contour_range_V,'linewidth',0.5);
refline([0 h+1*resolucion]);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('eje x [m]','fontsize',14);
ylabel('eje y [m]','fontsize',14);
title('Distribucion de potencial electrico, V(x,y) in volts','fontsize',14);
h1=gca;
set(h1,'fontsize',14);
fh1 = figure(1); 
set(fh1, 'color', 'white')

% Contorno para el campo electrico
figure(2)
clf
contour_range_E = 0:0.05:20;
contour(x,y,E,contour_range_E,'linewidth',0.5);
axis([min(x) max(x) min(y) max(y)]);
refline([0 h+1*resolucion]);
colorbar('location','eastoutside','fontsize',14);
xlabel('eje x [m]','fontsize',14);
ylabel('eje y [m]','fontsize',14);
title('Distribucion de campo Eléctrico, E (x,y)  [V/m]','fontsize',14);
h2=gca;
set(h2,'fontsize',14);
fh2 = figure(2); 
set(fh2, 'color', 'white')

% Lineas de Campo electrico
figure(3)
clf
hold on, quiver(x(1:paso:end,1:paso:end),y(1:paso:end,1:paso:end)...
    ,Ex(1:paso:end,1:paso:end),Ey(1:paso:end,1:paso:end),2)
refline([0 h+1*resolucion]);
title('Lineas de campo Eléctrico, E (x,y) in [V/m]','fontsize',14);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('eje x [m]','fontsize',14);
ylabel('eje y [m]','fontsize',14);
h3=gca;
set(h3,'fontsize',14);
fh3 = figure(3); 
set(fh3, 'color', 'white')

% Lineas de Campo magnético
figure(4)
clf
hold on
qui = quiver(x(1:paso:end,1:paso:end),y(1:paso:end,1:paso:end)...
    ,Bx(1:paso:end,1:paso:end),By(1:paso:end,1:paso:end),2);
refline([0 h+1*resolucion]);
title('Lineas de campo Magnético','fontsize',14);
axis([min(x) max(x) min(y) max(y)]);
xlabel('eje x [m]','fontsize',14);
ylabel('eje y [m]','fontsize',14);
h3=gca;
set(h3,'fontsize',14);
fh4 = figure(4); 
set(fh3, 'color', 'white')


