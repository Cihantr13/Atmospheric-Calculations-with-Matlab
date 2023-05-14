function [T,P,a,rho]=atmos2(h)
gama = 1.4;
mol = 0.0289644; %molar mass of air
R = 8.31446261815324; 
rad = 6371 * 1000;
g = 9.80665;
    
IGF = 9.780327 * (1+0.0053024*sin(45.5/57.2957)^2 - 0.0000058*2*sin(45.5/57.2957)^2);
FAC = -3.086 * 10^(-6) * h;
g = IGF + FAC;

    
if 0<=h<=11000;
P0 = 101325;
T0 = 288.15;
L0 = 0.0065;


H = (h*rad) / (h+rad);

T = T0-L0*(H-0);

P = P0 * (1  - ((L0 * (H-0))/T0)) ^ ((g*mol)/(R*L0));

a = sqrt((gama*R*T)/mol);

rho = (P*mol)/(R*T);

elseif h>11000 && h<=20000;
P0 = 22632.064;
T0 = 216.65;
L0 = 0;


H = (h*rad) / (h+rad);

T = T0-L0*(H-11000);


P = P0*exp(-mol*g*(H-11000))/(R*T0);

a = sqrt((gama*R*T)/mol);

rho = (P*mol)/(R*T);


elseif h>20000 && h<=32000
P0 = 5474.88867;
T0 = 216.65;
L0 = -0.001; 

    
H = (h*rad) / (h+rad);

T = T0-L0*(H-0);

P = P0 * (1  - ((L0 * (H-0))/T0)) ^ ((g*mol)/(R*L0));

a = sqrt((gama*R*T)/mol);

rho = (P*mol)/(R*T);


end
fprintf("T=")
disp(T);
fprintf("P=")
disp(P);
fprintf("a=")
disp(a);
fprintf("rho=")
disp(rho);

end