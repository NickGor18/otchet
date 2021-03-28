function kurs3

clc, clear all, close all

Nmpc = 10
0;
t0 = 0;
T = 10;
N  = 10;
dt = T/N;
Th = t0:dt:T;
km = ones(1,Nmpc);
um = ones(1,Nmpc);
% параметры системы
k0    = 1;
eps   = 0.01;
rho   = 0.05;
alpha = 0.3;
mu    = 0.03;
k_1 =(alpha/(rho+mu))^(1/(1-alpha));
u_1 = (mu * k_1)/(k_1^alpha);
lambda = 1/((k_1^alpha)*(1-u_1));
for l = 1:Nmpc
    km(l)=k_1;
end
for l = 1:Nmpc
    um(l)=u_1;
end

figure; hold on; grid on;

% вызов МРС
[k, c] = mpcsolve(Nmpc, k0, t0, 0.3*ones(N,1));
disp(k);
disp(c);


for itau = 1:Nmpc
     plot((itau-1)*dt+Th,k(:,itau), '--k', 'LineWidth',1);
end
plot((0:Nmpc)*dt,k(1,:), '-r', 'LineWidth',1);
plot((1:Nmpc)*dt,km, '--m', 'LineWidth',1);
%plot((1:Nmpc)*dt,um, '--m', 'LineWidth',1);
%plot((1:Nmpc)*dt,c(1,:), '-g', 'LineWidth',1);


% вычисление следующего состояния
function x = applycontrol(t,k,c)
    [~, x] = ode45(@(t,k) deq(k,c),[t, t+dt], k);
    x = x(end,:);
end

% алгоритм mpc
function [x, u] = mpcsolve(Nmpc, x0, t0, u0)
    t = t0;
    x = zeros(N + 1,Nmpc);
    x(1,1) = x0; 
    u = zeros(N,Nmpc);
    for i = 1:Nmpc
        u(:, i) = fmincon(@(u) costfun(x(1,i),t,u),...
            u0,[],[],[],[],zeros(N,1),ones(N,1)-eps);
        for j = 1:N
            x(j+1,i) = applycontrol(t + j*dt,x(j,i),u(j,i));
        end
        x(1,i+1) = x(2,i);
        t = t + dt;
    end
end


% критерий качества
function cost = costfun(x0,t0,u) 
    cost = 0; 
    x = x0; 
    t = t0; 
    for i = 1:N
        cost = cost + stagecost(t, x, u(i));
        x = applycontrol(t,x,u(i));
        t = t + dt;
    end
    %cost = cost + termcost(t0, x);
end

% стоимость этапа
function r = stagecost(t, k, c)
    r = (log(1 - c) + log(k^alpha))*...
        (exp(-rho * (t+dt))-exp(-rho * t))/rho;
end

% терминальная стоимость
function z = termcost(t, k)
    z = -(exp(-rho * (t+T-1)))*lambda*k;
end
% система
function dk = deq(k,c)
    dk = c*k^(alpha) - mu*k;
end

end