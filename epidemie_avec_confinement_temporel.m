clear all;

N = 6.7*10^6; %population
Ri = 3.41; % le r0 de l'épidémie
Ti = 10; %temps de guérison
Mi = 0.052; %taux de mortalité

%conditions initiales pour S, I, R
S0 = (N-1)/N;
I0 = (1)/N;
R0 = (0)/N;
M0 = (0)/N;

X = [S0, I0, R0, M0];

%%% X'(t) = f(X) %%%
alpha = Ri/Ti;
beta = 1/Ti;
delta = Mi/Ti;
f = @(X) [-alpha*X(1)*X(2), alpha*X(1)*X(2) - (beta+delta)*X(2) , beta*X(2), delta*X(2)];
bool_confinement = 0;
debut_conf = 0;

%%%%%
T = 400; %temps
h = 0.01; %%pas de temps

for i=1:T/h,
    
    %55 jours de confinement à partir de 3%
    if (X(2) > 0.03 && bool_confinement == 0) %%déclenchement confinement
        %dimunition du R0%
        Ri = 0.52;
        alpha = Ri/Ti;
        f = @(X) [-alpha*X(1)*X(2), alpha*X(1)*X(2) - (beta+delta)*X(2) , beta*X(2), delta*X(2)];
        bool_confinement = 1;
        debut_conf = i;
    elseif (i > debut_conf + 55/h && bool_confinement == 1) %%arrêt confinement
        %augmentation du R0%
        Ri = 3.41;
        alpha = Ri/Ti;
        f = @(X) [-alpha*X(1)*X(2), alpha*X(1)*X(2) - (beta+delta)*X(2) , beta*X(2), delta*X(2)];
        bool_confinement = 2;
    end
    %%%%%%%%%%
    k = X + h*f(X);
    X = X + h*(f(X)+f(k))/2;
    
    %%%%%%%%%%
    S(i) = X(1);
    I(i) = X(2);
    R(i) = X(3);
    M(i) = X(4);
    t(i) = i*h;
end

figure('Name','Epidemie avec confinement');
clf;
plot(t, N*S, 'b');
hold on;
plot(t, N*I, 'r');
plot (t, N*R, 'g');
plot (t, N*M, 'k');
plot (t, N*(S+I+R+M), 'c');