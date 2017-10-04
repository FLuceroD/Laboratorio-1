clc;
clear all;
close all;

display('The given dimensions are-');
display('h=15,lx=120,ly=60,dx=dy=5 ; everything is in mm ');
display('the permittivities are-');
display('E1=E0=8.854e-12,E2= 3.5*E0');
display('Initial step sizes are');
display('dx=dy=5 ');

lx = 120;  % in mm
ly = 60; % in mm
w = 30; %lx/4; % in mm
h = ly/4; % in mm
d = 5;  % dx=dy  in mm
wd = w / d;
hd = h / d;
ldx = lx / d;
ldy = ly / d;
Eo = 1; %8.854e-12;
E1 = Eo;
E2 = 1000*Eo;
p = zeros(ldy,ldx);
p(1,:)=0;p(ldy,:)=0;p(:,1)=0;p(:,ldx)=0;
p((ldy - hd),(ldx/2-wd/2):((ldx/2)+(wd/2))) = 1;

for N = 1:300
for i= 2:ldy-1
    for j=2:ldx-1
       
        if (i == (ldy - hd))
            if (j<((ldx/2)-(wd/2))+1) || (j>((ldx/2)+(wd/2)))
                p(i,j)= ((p(i,j+1)/4)+(p(i,j-1)/4)+(E1/(2*(E1+E2)))*p(i+1,j)+(E2/(2*(E1+E2)))*p(i-1,j));
            end
            
        else
            p(i,j)= (p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1))/4;
        end
    end
    
end
end
p0 = zeros(ldy,ldx);
p0(1,:)=0;p0(ldy,:)=0;p0(:,1)=0;p0(:,ldx)=0;
p0((ldy - hd),((ldx/2)-(wd/2)):((ldx/2)+(wd/2))) = 1;

for N = 1:300
for i= 2:ldy-1
    for j=2:ldx-1
       if  ~(  (j>=((ldx/2)-(wd/2)+1)  &&  j<=((ldx/2)+(wd/2))) && (i == (ldy - hd))   )
              p0(i,j)= (p0(i+1,j)+p0(i-1,j)+p0(i,j+1)+p0(i,j-1))/4; 
       end
    end
    
end
end
subplot(2,1,1);
imagesc(p0);
title('Potential in the solution domain for the cases where the bottom medium is free space','fontsize',12);
subplot(2,1,2);
imagesc(p);

p0 = Eo*p0;

sm1 = sum(p0((ldy-hd-2),((ldx-wd)/2)-1:((ldx+wd)/2)+1)) - sum(p0((ldy-hd-1),((ldx-wd)/2)-1:((ldx+wd)/2)+1));
sm2 = sum(p0((ldy-hd+2),((ldx-wd)/2)-1:((ldx+wd)/2)+1)) - sum(p0((ldy-hd+1),((ldx-wd)/2)-1:((ldx+wd)/2)+1));
sm3 = sum(p0((ldy-hd-1):(ldy-hd+1),((ldx-wd)/2)-1)) - sum(p0((ldy-hd-1):(ldy-hd+1),((ldx-wd)/2)-0));
sm4 = sum(p0((ldy-hd-1):(ldy-hd+1),((ldx+wd)/2)+2)) - sum(p0((ldy-hd-1):(ldy-hd+1),((ldx+wd)/2)+1));
k = sm1 + sm2 + sm3 + sm4;

sq1 = E1*(sum(p((ldy-hd-2),((ldx-wd)/2)-1:((ldx+wd)/2)+1)) - sum(p((ldy-hd-1),((ldx-wd)/2)-1:((ldx+wd)/2)+1)));
sq2 = E2*(sum(p((ldy-hd+2),((ldx-wd)/2)-1:((ldx+wd)/2)+1)) - sum(p((ldy-hd+1),((ldx-wd)/2)-1:((ldx+wd)/2)+1)));

sq3 =              E1*(p((ldy-hd-1),((ldx-wd)/2)-2) - p((ldy-hd-1),((ldx-wd)/2)-1));
sq4 = ((E1/2)+(E2/2))*(p((ldy-hd+0),((ldx+wd)/2)-2) - p((ldy-hd+0),((ldx-wd)/2)-1));
sq5 =              E2*(p((ldy-hd+1),((ldx-wd)/2)-2) - p((ldy-hd+1),((ldx-wd)/2)-1));

sq6 =              E2*(p((ldy-hd+1),((ldx+wd)/2)+2) - p((ldy-hd+1),((ldx+wd)/2)+1));
sq7 = ((E1/2)+(E2/2))*(p((ldy-hd+0),((ldx+wd)/2)+2) - p((ldy-hd-0),((ldx+wd)/2)+1));
sq8 =              E1*(p((ldy-hd-1),((ldx+wd)/2)+2) - p((ldy-hd-1),((ldx+wd)/2)+1));

sq = sq1 + sq2 + sq3 + sq4 + sq5 + sq6 + sq7 + sq8;

c = sq; c0 = k;
L = (1/(3e+8*sqrt(c0)))^2;
Zo = 1 / (3e+8*sqrt(c*c0));
Vp = 3e+8*sqrt(c0/c);