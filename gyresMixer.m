
function [ y ] = gyresMixer(T, x1, A,eps, omega,uwx,uwy)

g = @(t,z) eps*sin(omega*t)*z.^2+(1-2*eps*sin(omega*t))*z;
dg = @(t,z) 2*z.*eps*sin(omega*t)+(1-2*eps*sin(omega*t));


dx1 = (x1(:,1)>=0 & x1(:,1)<=2).*(-pi*A*sin(pi*g(T,x1(:,1))).*cos(pi*x1(:,2)))+uwx;
dy1 = (x1(:,1)>=0 & x1(:,1)<=2).*(pi*A*cos(pi*g(T,x1(:,1))).*sin(pi*x1(:,2)).*dg(T,x1(:,1)))+uwy;

y=[dx1, dy1];
end
