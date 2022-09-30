%% Benefits of antibiotic resistance in commensals and the scope for resistance optimization
% Kristofer Wollein Waldetoft, Sarah Sundius, Rachel Kuske, Sam P. Brown

% Code by: Sarah Sundius
% September 26, 2022

clear; clc;

%% Calculate + plot optimal f given coexistence

% parameters: 
% rp,rc = maximal growth rate
% kp,kc = carrying capacity
% apc, acp = interspecific interaction coeff
% x = maximal pathogen clearance rate
% A = antibiotic exposure (control parameter)
% f = commensal relative susceptibility (control parameter)
% wp,wc,whb = per-capita risk weightings for each population

% parameters
syms rp rc kp kc apc acp x f A wp wc whb

% pathogen density @ coexistence
Pstar = (rc*kp*(rp-x*A) - apc*rp*kc*(rc-x*f.*A))./(rp*rc*(1-apc*acp));

% commensal density @ coexistence
Cstar = (rp*kc*(rc-x*f.*A) - acp*rc*kp*(rp-x*A))./(rp*rc*(1-apc*acp));

% substitute values for growth, carrying cap, interactions, abx
params = [rp rc kp kc apc acp x A];
paramvals = [0.5 0.5 1 1 0.8 0.8 0.1 1];    %can change interaction type
Pstar1 = subs(Pstar,params,paramvals);
Cstar1 = subs(Cstar,params,paramvals);

% define lines for plotting with weights
Pline = wp*Pstar1;
Cline = wc*Cstar1;
HGTline = whb*Pstar1.*Cstar1;

% substitute values for weights
Pline1 = subs(Pline,wp,1);
Cline1 = subs(Cline,wc,0.1);
HGTline1 = subs(HGTline,whb,5);

% define f and Pstar1, Cstar1 vectors for plotting
fvals = linspace(0,2.1,200);
Pline1_plot = subs(Pline1,f,fvals);
Cline1_plot = subs(Cline1,f,fvals);
Cline1_plot(Cline1_plot<0) = 0;         %exclude negative values for bio relevance
HGTline1_plot = subs(HGTline1,f,fvals);
HGTline1_plot(HGTline1_plot<0) = 0;     %exclude negative values for bio relevance

net1 = Pline1_plot + Cline1_plot;
net2 = net1 + HGTline1_plot;


% create figure
figure
T = tiledlayout(2,1);
set(gcf,'unit','inches','position',[2 2 6 6]);
xlabel(T,'Commensal relative susceptibility, f','FontSize',12,'Color','k')
ylabel(T,'Risk of infection','FontSize',12,'Color','k')

ax1 = nexttile;
hold on
box on
plot(fvals,Pline1_plot,'r--','LineWidth',1.5)
plot(fvals,Cline1_plot,'b--','LineWidth',1.5)
plot(fvals,net1,'k-','LineWidth',1.5)
title('A','FontWeight','normal','FontSize',12,'Color','k')
ax1.TitleHorizontalAlignment = 'left';
ax1.YAxis.FontSize = 9;
ax1.LineWidth = 1;
xlim([0 2.1])
ylim([-0.2 1.5])

ax2 = nexttile;
hold on
box on
plot(fvals,Pline1_plot,'r--','LineWidth',1.5)
plot(fvals,Cline1_plot,'b--','LineWidth',1.5)
plot(fvals,HGTline1_plot,'Color',[0.4940, 0.1840, 0.5560],'LineStyle','--','LineWidth',1.5)
plot(fvals,net2,'k-','LineWidth',1.5)
title('B','FontWeight','normal','FontSize',12,'Color','k')
ax2.TitleHorizontalAlignment = 'left';
ax2.XAxis.FontSize = 9;
ax2.YAxis.FontSize = 9;
ax2.LineWidth = 1;
xlim([0 2.1])
ylim([-0.2 2])

leg_ax1 = legend(ax1,{'Pathogen','Commensal','Net'},'FontSize',9,'TextColor','k','Location','NorthEastOutside');
leg_ax2 = legend(ax2,{'Pathogen','Commensal','HGT pathogen','Net'},'FontSize',9,'TextColor','k','Location','NorthEastOutside');

% save figure
%exportgraphics(T,'OptimizationFig.eps','Resolution',600,'ContentType','vector');


