function new_runchrono(H_itp,theta,Deltar,alpha,K)

%Input parameters
%=================

% Deltar  = 7;
% alpha   = 1;
% C       = 1;
% lambda  = 1;

%=================

%Load H and theta and interpolate

% filename = ['H_Dr_' num2str(Deltar) '.mat'];
% load(filename);
% 
% g_grid = 0:.002:1;
% g_itp  = 0:.0001:1;
% 
% H_itp = interp2(g_grid,theta,H,g_itp,theta,'spline');

%=================

% mu     = -2:.2:3;
% mu     = 10.^mu;
sigma  = 1/sqrt(alpha);
mu     = logspace(log10(.01*sigma),log10(3*sigma),20);
K1     = 0*mu;
K2     = 0*mu;
T_mean = 0*mu;
T_med  = 0*mu;

for i = 1:length(mu)
    [k, t_mean, t_med] = new_chronometric(H_itp,theta,mu(i),alpha,K,Deltar,0);
    K1(i)              = k(1);
    K2(i)              = k(2);
    T_mean(i)          = t_mean;
    T_med(i)           = t_med;
end

[~, easy]    = min(abs(K1-.95*K1(end)));
[~, medium]  = min(abs(K1-(K1(1)+K1(end))/2));
[~, hard]    = min(abs(K1-1.1*K1(1)));

mu_easy   = round(mu(easy),2);
mu_medium = round(mu(medium),2);
mu_hard   = round(mu(hard-1),2);

close all
figure(1)

values = annotation('textbox',[0.1 0.8 0.3 0.1],...
                   'FitHeightToText','off',...
                   'Margin',0,...
                   'FontSize',13,...
                   'FontWeight','bold',...
                   'LineStyle','none');
               
set(values,'String',['\Deltar = ' num2str(Deltar)...
                     '        \alpha = ' num2str(alpha)...
                     '        C + \lambda = ' num2str(K)]);

a11=axes('position',[0.07 0.45 0.28 0.33]);
semilogx(mu,K1)
xlim([mu(1) mu(end)])
ylim([0.5 1]);
hold on
semilogx(mu,K2)
%a11.ColorOrderIndex = 3;
%get(gca,'YLim')
semilogx([mu_hard mu_hard],[0.5 1],'--k');
text(mu_hard,1.05,sprintf(' \\mu = %.2f \n ( hard)',mu_hard),...
    'HorizontalAlignment','center','FontSize',7);
semilogx([mu_medium mu_medium],[0.5 1],'--k');
text(mu_medium,1.05,sprintf('\\mu = %.2f \n ( medium) ',mu_medium),...
    'HorizontalAlignment','center','FontSize',7);
semilogx([mu_easy mu_easy],[0.5 1],'--k');
text(mu_easy,1.05,sprintf('\\mu = %.2f \n ( easy) ',mu_easy),...
    'HorizontalAlignment','center','FontSize',7);
ylabel('%R_{max}');
legend('t_{\infty}','<t>','Location','southeast');
set(gca,'XTickLabel',[]);

a12=axes('position',[0.07 0.10 0.28 0.33]);
semilogx(mu,T_mean,mu,smooth(T_med,'lowess'))
xlim([mu(1) mu(end)])
T_lim = 1.2*max(max(T_mean,T_med));
ylim([0 T_lim]);
hold on
semilogx([mu_hard mu_hard],[0 T_lim],'--k');
semilogx([mu_medium mu_medium],[0 T_lim],'--k');
semilogx([mu_easy mu_easy],[0 T_lim],'--k');
xlabel('\mu');
ylabel('t');
legend('<t>','t_{median}','Location','southeast');

global b1;
global b2;

global c1;
c1 = axes('position',[0.4, 0.07, 0.26, 0.26]);
global c2;
c2 = axes('position',[0.7, 0.07, 0.26, 0.26]);

[k_hard,t_hard,x_hard]       = new_chronometric(H_itp,theta,mu_hard,alpha,K,Deltar,1,1);
[k_medium,t_medium,x_medium] = new_chronometric(H_itp,theta,mu_medium,alpha,K,Deltar,1,2);
[k_easy,t_easy,x_easy]       = new_chronometric(H_itp,theta,mu_easy,alpha,K,Deltar,1,3);

figure(1)

b_lim = 1.2*max([k_hard(1),k_medium(1),k_easy(1)]);
x_lim = ceil(max([x_hard,x_medium,x_easy]));
if x_lim==1
    x_lim = round(max([x_hard,x_medium,x_easy]),1);
end
if x_lim<0.1
    x_lim = 0.1;
end

axes(b1(1))
xlim([0 x_lim])
ylim([0 b_lim])
ax = gca;
ax.ColorOrderIndex = 4;
plot([t_hard t_hard],[0 b_lim],'--');
ylabel('p*')

axes(b1(2))
xlim([0 x_lim])
ylim([0 b_lim])
ax = gca;
ax.ColorOrderIndex = 4;
plot([t_medium t_medium],[0 b_lim],'--');
set(gca,'YTickLabel',[]);

axes(b1(3))
xlim([0 x_lim])
ylim([0 b_lim])
ax = gca;
ax.ColorOrderIndex = 4;
plot([t_easy t_easy],[0 b_lim],'--');
set(gca,'YTickLabel',[]);
legend('p_{R}','p_{L}','<t_{R}>','Location','northeast')

axes(b2(1))
xlim([0 x_lim])
ylabel('P*')

axes(b2(2))
xlim([0 x_lim])
set(gca,'YTickLabel',[]);

axes(b2(3))
xlim([0 x_lim])
set(gca,'YTickLabel',[]);
legend('P_{R}','P_{L}','Location','east')

axes(c1);
xlim([0 x_lim])
ylim([.5 1])
xlabel('t')
ylabel('%R')
xlab2 = get(gca,'XLabel');
set(xlab2,'Position',get(xlab2,'Position') + [0 .05 0])

axes(c2);
xlim([0 x_lim])
ylim([0 1])
xlabel('t')
ylabel('%R (normalised)')
xlab3 = get(gca,'XLabel');
set(xlab3,'Position',get(xlab3,'Position') + [0 .1 0])
plot([0 0],[0 0],'-.');
legend('hard','medium','easy','<g(\mu)>','Location','southeast');

set(gcf,'paperpositionmode','auto')
%get(gcf,'position')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
figurename = ['plots_chronometrics/figure'...
              '_Dr_'     num2str(Deltar)...
              '_alpha_'  num2str(alpha)...
              '_K_'      num2str(K) '.png'];
saveas(gcf,figurename);
reportfig();

figure(2)
uistack(gcf,'top');
reportfig();

