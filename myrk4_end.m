function [y,rt]=rk4(ode, tstart, tfinal, h, iv, A, eps, omega,uwx,uwy) %rk4(ode, tstart, tfinal, h, iv)   %

T=tstart;
dt=T;
j=1;
x=iv;
dy=x;
m=ceil((tfinal-tstart)/h+1);
yarray=zeros(m, size(x,1),size(x,2));
tvector =zeros(1, m);
tvector(1)=T;
yarray(1, :,:)=x;

taus=NaN(size(x,1),1);
tein=NaN(size(x,1),1);



for j=2:m
    
    
 tein(x(:,1)>=0 & isnan(tein))=T;
 taus(x(:,1)>2 & isnan(taus))=T;

    
    
k1=feval(ode, T, x,A,eps, omega,uwx,uwy);
k2=feval(ode, T+0.5*h, x+0.5*h.*k1,A,eps, omega,uwx,uwy); 
k3=feval(ode, T+0.5*h, x+0.5*h.*k2,A,eps, omega,uwx,uwy); 
k4=feval(ode, T+h, x+h.*k3,A,eps, omega,uwx,uwy);
x_temp=x+1/6*h*(k1+2.*k2+2.*k3+k4);


%x(:,1)=max(min(x_temp(:,1),4),0);
%x(:,2)=max(min(x_temp(:,2),1),0);

%x(:,1)=mod(x_temp(:,1),1);
%x(:,2)=mod(x_temp(:,2),1);

%x(:,1)=mod(x_temp(:,1),pi);



x = x_temp;

%x(:,1)=mod(x_temp(:,1),2.5);

T=T+h;
dt=T;
dy=x;
tvector(j)=T;
yarray(j,:,:)=x;





end;

rt=taus-tein;

%t = T;
%y = x;

%t=tvector;
%yarray(:,:,1)=mod(yarray(:,:,1),2.5 );


y=squeeze(yarray(end,:,:));
%y=yarray;


