function [k,t_mean,t_med]=new_chronometric(H_itp,theta,mu,alpha,K,Deltar,do_plots,n)

%Input parameters
%=================

%Time
T       = 80;    %final time
dt      = 0.01; %interval steps

%Define and init variables
%==========================
t   = (.1+dt:dt:3);
t_l = -5:.5:-3;
t_l = 10.^t_l;
t   = [10^-6 t_l (.002:.001:.1) t 3+0.1:0.1:T]';
n_t = length(t);
Phi = 0*t;
u_R = 0*t;
u_L = 0*t;
p_R = 0*t;
p_L = 0*t;
P_R = 0*t;
P_L = 0*t;


%dx    = normrnd(mu*dt*ones(n_t,1),dt*ones(n_t,1));
%dx(1) = 0;
%x     = tril(ones(n_t,n_t))*dx;
%g     = normcdf(bsxfun(@rdivide,x,sqrt(t+alpha)));

x_a     = mu*t;
g_a     = normcdf(bsxfun(@rdivide,x_a,sqrt(t+alpha)));

omega   = Deltar*(2*g_a-1)/2;
Omega_R = exp(omega);
Omega_L = exp(-omega);


for i = 2:n_t
    
    g_i = int32(round(10000*g_a(i))+1);
    h      = H_itp(:,g_i); 
    Phi(i) = new_desirability(h,theta,t(i),K,alpha);
    u_R(i) = Omega_R(i)/Phi(i) * 0.5;
    u_L(i) = Omega_L(i)/Phi(i) * 0.5;
    
    I_R    = trapz(t(1:i),u_R(1:i));
    I_L    = trapz(t(1:i),u_L(1:i));
    I_RL   = I_R + I_L;
    
    p_R(i) = u_R(i)*exp(-I_RL);
    p_L(i) = u_L(i)*exp(-I_RL);
    
    P_R(i) = trapz(t(1:i),p_R(1:i));
    P_L(i) = trapz(t(1:i),p_L(1:i));
    
end

t       = t(2:end);
g_a     = g_a(2:end);
Phi     = Phi(2:end);
omega   = omega(2:end);
Omega_R = Omega_R(2:end);
Omega_L = Omega_L(2:end);
u_R     = u_R(2:end);
u_L     = u_L(2:end);
p_R     = p_R(2:end);
p_L     = p_L(2:end);
P_R     = P_R(2:end);
P_L     = P_L(2:end);

%p_naive     = 0.5*exp(-t);
acc         = bsxfun(@rdivide,P_R,P_R+P_L);
k(1)        = acc(end);
t_mean      = trapz(t,t.*p_R);
[~, i_mean] = min(abs(t-t_mean));
k(2)        = acc(i_mean);

[~, i_med] = min(abs(P_R-.5*P_R(end)));
t_1        = t(i_med);
y_1        = P_R(i_med);
t_2        = 0;
y_2        = 0;
if y_1<0.5
    t_2 = t(i_med+1);
    y_2 = P_R(i_med+1);
else
    t_2 = t(i_med-1);
    y_2 = P_R(i_med-1);
end
t_med = (t_2-t_1)/(y_2-y_1)*(0.5*P_R(end)-y_1)+t_1;

if (do_plots && nargin == 8)

    k          = max(p_R');
    [~, x_cut] = min(abs(p_R/max(p_R)-1e-3));
    t_med      = t(x_cut);
    
    figure(1)
    
    global b1;
    b1(n) = axes('position',[0.4+(n-1)*0.19, 0.68, 0.18, 0.26]);
    plot(t,[p_R,p_L]);
    hold on
    title(['\mu = ', num2str(mu)]);
    set(gca,'XTickLabel',[]);
    
    global b2;
    b2(n) = axes('position',[0.4+(n-1)*0.19, 0.4, 0.18, 0.26]);
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(t,[P_R,P_L]);
    ylim([0 1])
    hold on
    ax.ColorOrderIndex = 4;
    plot([t_mean t_mean],ylim,'--');
    xlabel('t')
    xlabh = get(gca,'XLabel');
    set(xlabh,'Position',get(xlabh,'Position') + [0 .15 0])
    
    global c1;
    axes(c1);
    ax = gca;
    ax.ColorOrderIndex = n;
    plot(t,acc);
    hold on
    ax.ColorOrderIndex = 4;
    plot(t,g_a,'-.');
    
    global c2;
    axes(c2);
    plot(t,(acc-0.5)/(acc(end)-0.5));
    hold on
    
    if n==2
        figure(2)
        plot(t,[Omega_R,2*Phi,u_R,...
            K/2*(1+tanh(omega)),...
            K/2*Omega_R./cosh(Deltar/2)])
        xlim([0 80])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
    end
    


end