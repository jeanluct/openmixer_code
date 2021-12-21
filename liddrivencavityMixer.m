function [ y ] = liddrivencavityMixer(T, x1, uwx,uwy,tf, U1,U2,a,b)

f = @(k,y)(2*pi*y/a)*cosh(k*pi*b/a).*sinh(k*pi*y/a)-(2*pi*b/a)*sinh(k*pi*b/a)*cosh(k*pi*y/a);
df = @(k,y) cosh(k*pi*b/a)*((2*pi/a)*sinh(k*pi*y/a)+(2*pi*y/a).*cosh(k*pi*y/a)*(k*pi/a))-(2*pi*b/a)*sinh(k*pi*b/a)*sinh(k*pi*y/a)*(k*pi/a);

C1=(a^2/(2*pi^2*b))/((a/(2*pi*b))*sinh(2*pi*b/a)+1);
C2=(a^2/(4*pi^2*b))/((a/(4*pi*b))*sinh(4*pi*b/a)+1);

dx1 = (x1(:,1)>=0 & x1(:,1)<=a).*(((mod(T,tf)<tf/2)*(-1)+(mod(T,tf)>=tf/2)*(1))*(U1*C1*df(1,x1(:,2)).*sin(pi*x1(:,1)/a))+U2*C2*df(2,x1(:,2)).*sin(2*pi*x1(:,1)/a))+uwx ;
dy1 = (x1(:,1)>=0 & x1(:,1)<=a).*(((mod(T,tf)<tf/2)*(-1)+(mod(T,tf)>=tf/2)*(1))*(-U1*C1*f(1,x1(:,2)).*cos(pi*x1(:,1)/a)*(pi/a))-U2*C2*f(2,x1(:,2)).*cos(2*pi*x1(:,1)/a)*2*pi/a)+uwy;


y=[dx1, dy1];

end