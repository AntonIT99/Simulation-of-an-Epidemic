  clear all;

T = 200;
R0 = 3;
Ti = 10;
N  = 10^6;
beta= R0/Ti;


gamma = 1/Ti;

X = [1 - 1/N,1/N,0];
f = @(X) [- beta*X(1)*X(2), beta*X(1)*X(2) - gamma*X(2), gamma*X(2)];



h = 0.1;
for i=1:T/h,
   %% Euler 
   %X = X + h*f(X);
   %% Rk2 
   k = X + h*f(X);
   X = X + h/2*(f(X) + f(k));
   
   %%%%%%%ù
   t(i) = i*h;
   I(i) = X(1); S(i) = X(2); R(i) = X(3);
   
end

clf
plot(t,I);
hold on;
plot(t,S,'r');
plot(t,R,'g');

%%%%%%%%%%%%%%%%%%%% quarantaine %%%%%%%%%%%%ù
% 
% 
% T =100;
% R0 = 3;
% tau =6;
% Ti = 12;
% N  = 10^8;
% beta= R0/(Ti-tau);
% 
% 
% gamma = 1/(Ti-tau);
% 
% X = [1 - 1/N,1/N,0,0];
% f = @(X,IS_tau) [- beta*X(1)*X(2),beta*(X(1)*X(2)- IS_tau), beta*IS_tau-gamma*X(3),gamma*X(3)];
% 
% 
% h = 0.01; i_tau = round(tau/h);
% 
% I(1) = X(1); S(1) = X(2); S_q(1) = X(3);
% for i=1:T/h,
%    %% Euler 
%    %X = X + h*f(X);
%    %% Rk2 
%    IStau = (i>i_tau)*I(max(i-i_tau,1))*S(max(i-i_tau,1));
%    IStaup = (i>i_tau)*I(max(i-i_tau+1,1))*S(max(i-i_tau+1,1));
%    k = X + h*f(X, IStau);
%    X = X + h/2*(f(X, IStau) + f(k, IStaup));
%    
%    %%%%%%%ù
%    t(i) = i*h;
%    I(i) = X(1); S(i) = X(2); S_q(i) = X(3); R(i) = X(4); 
%    
% end
% 
% clf
% plot(t,I);
% hold on;
% plot(t,S,'r');
% plot(t,S_q,'m');
% plot(t,R,'g');



% %%%%%%%%%%%%%%%%% confinement  %%%%%%%%ù
% 
% clear all;
% 
% T = 300;
% R0 = 3;
% Ti = 10;
% N  = 10^8;
% beta= R0/Ti;
% 
% gamma = 1/Ti;
% 
% X = [1 - 1/N,1/N,0];
% f = @(X) [- beta*X(1)*X(2), beta*X(1)*X(2) - gamma*X(2), gamma*X(2)];
% 
% Bool_conf= 0;
% 
% 
% h = 0.001;
% for i=1:T/h,
%     
%     if (X(2)>0.1 && Bool_conf==0)
%        alpha = 0.4;
%         beta= (1-alpha)*R0/Ti;
%        f = @(X) [- beta*X(1)*X(2), beta*X(1)*X(2) - gamma*X(2), gamma*X(2)];
%         Bool_conf=1;
%     elseif (X(2)<0.05 && Bool_conf==1)
%        beta= R0/Ti;
%        f = @(X) [- beta*X(1)*X(2), beta*X(1)*X(2) - gamma*X(2), gamma*X(2)];
%         Bool_conf=0;
%     end
%     
%    %% Euler 
%    %X = X + h*f(X);
%    %% Rk2 
%    k = X + h*f(X);
%    X = X + h/2*(f(X) + f(k));
%    
%    %%%%%%%ù
%    t(i) = i*h;
%    I(i) = X(1); S(i) = X(2); R(i) = X(3);
%    
% end
% 
% clf
% plot(t,I);
% hold on;
% plot(t,S,'r');
% plot(t,R,'g');



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Modele avec diffusion %%%%%%%%%%%%%
% T = 300;
% R0 = 3;
% Ti = 10;
% N  = 10^8;
% beta= R0/Ti;
% gamma = 1/Ti;
% dt = 0.01
% 
% 
% Nx = 100;
% posx = linspace(0,1,Nx)
% dx = 1/Nx;
% 
% %%%%% Coef diffusion %%%%%
% Da = 0.00001;
% A = (2*diag(ones(1,Nx)) - diag(ones(1,Nx-1),1) -diag(ones(1,Nx-1),-1));
% A(1,1) = 1; A(Nx,Nx) = 1;
% L = eye(Nx) + dt*Da*A/dx^2;
% 
% 
% X = [ones(size(posx));zeros(size(posx));zeros(size(posx))];
% X(1,50) = 1- 10/N;
% X(2,50) = 10/N;
% 
% 
% f = @(X) [- beta*X(1,:).*X(2,:); beta*X(1,:).*X(2,:) - gamma*X(2,:); gamma*X(2,:)];
% 
% for i=1:T/dt,
%    %% Euler 
%    %X = X + h*f(X);
%    %% Rk2 
%    k = X + dt*f(X);
%    X = X + dt/2*(f(X) + f(k));
%    
%    %%%%% diffusion
%    X(1,:) =  (L\(X(1,:)'))';
%    X(2,:) =  (L\(X(2,:)'))';
%    X(3,:) =  (L\(X(3,:)'))';
%    
%    
%    %%%%%%%ù
%    t(i) = i*dt;
%    I(i) = sum(X(1,:))/Nx; S(i) = sum(X(2,:))/Nx; R(i) = sum(X(3,:))/Nx;
%   
%   if mod(i,100)==1,
%       i
%   clf
%   plot(posx, X(1,:));
%   hold on;
%   plot(posx,X(2,:),'r');
%   plot(posx,X(3,:),'g');
%   pause(0.01)
%   end
%    
% end
% 
% 
% 
%  clf
%  Ii = 1:100:T/dt;
%  plot(t(Ii),I(Ii));
%  hold on;
%  plot(t(Ii),S(Ii),'r');
%  plot(t(Ii),R(Ii),'g');
% 
% % 



