%% Benefits of antibiotic resistance in commensals and the scope for resistance optimization
% Kristofer Wollein Waldetoft, Sarah Sundius, Rachel Kuske, Sam P. Brown

% Code by: Sarah Sundius
% August 18, 2021


%% Heat Map/Contour Plots
% code to create Figure 3 and Figures S1-4

clear; clc;

% parameters (+ low/high variations for SI))
rp = 0.5;   %rp = 0.25;  %rp = 0.75;  
rc = 0.5;   %rc = 0.25;  %rc = 0.75;
kp = 1;     %kp = 0.5;  %kp = 1.5;
kc = 1;     %kc = 0.5;  %kc = 1.5;
x = 0.1;

% resource competition (apc > 0, acp > 0)
apc1 = 0.8;
acp1 = 0.8;

% commensal exploitation of pathogen (apc > 0, acp < 0)
apc2 = 0.8;
acp2 = -0.8;

% pathogen exploitation of commensal (apc < 0, acp > 0)
apc3 = -0.8;
acp3 = 0.8;

% mutualism (apc < 0, acp < 0)
apc4 = -0.8;
acp4 = -0.8;

% control parameters
a = linspace(0.001,1,100);
f = linspace(0.5,2.1,160);
[A,F] = meshgrid(a,f);

% pathogen density @ coexistence eq
Pstar1 = (rc*kp*(rp-x*A) - apc1*rp*kc*(rc-x*F.*A))./(rp*rc*(1-apc1*acp1));
Pstar2 = (rc*kp*(rp-x*A) - apc2*rp*kc*(rc-x*F.*A))./(rp*rc*(1-apc2*acp2));
Pstar3 = (rc*kp*(rp-x*A) - apc3*rp*kc*(rc-x*F.*A))./(rp*rc*(1-apc3*acp3));
Pstar4 = (rc*kp*(rp-x*A) - apc4*rp*kc*(rc-x*F.*A))./(rp*rc*(1-apc4*acp4));

% calculate competitive release (only in apc > 0 cases)
comp_rel1 = ((rc*kp)/(apc1*rp*kc))*ones(1,length(a));
comp_rel2 = ((rc*kp)/(apc2*rp*kc))*ones(1,length(a));

% calculate stability constraints
f1_1 = (rp*rc*kc - acp1*rc*kp*(rp-x*a))./(rp*kc*x*a);           % upper bound - P dominance
f2_1 = (apc1*rp*rc*kc - rc*kp*(rp-x*a))./(apc1*rp*kc*x*a);      % lower bound`- C dominance

f1_2 = (rp*rc*kc - acp2*rc*kp*(rp-x*a))./(rp*kc*x*a);           % upper bound - P dominance
f2_2 = (apc2*rp*rc*kc - rc*kp*(rp-x*a))./(apc2*rp*kc*x*a);      % lower bound - C dominance

f1_3 = (rp*rc*kc - acp3*rc*kp*(rp-x*a))./(rp*kc*x*a);           % upper bound - P dominance
f2_3 = (apc3*rp*rc*kc - rc*kp*(rp-x*a))./(apc3*rp*kc*x*a);      % upper bound - C dominance
if f1_3 < f2_3
    fhigh_3 = f1_3;
else 
    fhigh_3 = f2_3;
end

f1_4 = (rp*rc*kc - acp4*rc*kp*(rp-x*a))./(rp*kc*x*a);           % upper bound - P dominance
f2_4 = (apc4*rp*rc*kc - rc*kp*(rp-x*a))./(apc4*rp*kc*x*a);      % upper bound - C dominance
if f1_4 < f2_4
    fhigh_4 = f1_4;
else
    fhigh_4 = f2_4;
end

% define colorbar limits for A-C
minPstar = min([min(Pstar1,[],'all'), min(Pstar2,[],'all'), min(Pstar3,[],'all')]);
maxPstar = max([max(Pstar1,[],'all'), max(Pstar2,[],'all'), max(Pstar3,[],'all')]);
if minPstar < 0
    minPstar = 0;
end

% define colormap
d = jet(200);
cmin = 50;          % blue
cmax = 150;         % orange
Ponly = d(180,:);   % red 
Conly = d(10,:);    % dark blue/purple

% make contour plots (FIGURE 3, S1-4)
fig = figure;
set(gcf,'unit','inches','position',[0 0 6 6])

% PANEL (A): COMPETITION
ax(1) = subplot(2,2,1);
set(gca,'Unit','Inches','Position',[0.5 3.25 2.5 2.5]);
contourf(a,f,Pstar1,minPstar:0.03:maxPstar,'ShowText','off','LineWidth',1)
axis square
axis([0 1 0.5 2.1])
colormap(ax(1), d(cmin:cmax,:))
hold on
patch([a fliplr(a)], [f2_1 min(ylim)*ones(size(f2_1))],Conly,'LineStyle','none') % below lower bound
patch([a fliplr(a)], [f1_1 max(ylim)*ones(size(f1_1))],Ponly,'LineStyle','none') % above upper bound
hold off
hold on
line(a,comp_rel1,'Color','white','LineStyle','--','LineWidth',2)
c = colorbar;
c.LineWidth = 1;
caxis([minPstar,maxPstar])
xlim([0 1])
ylim([0.5 2.1])
ax = gca;
ax.FontSize = 8;
ax.LineWidth = 1;
xlabel('antibiotic exposure, A','FontSize',9,'Color','k')
ylabel('commensal relative susceptibility, f','FontSize',9,'Color','k')
title('Competition','FontWeight','Normal','FontSize',10,'Color','k')

% PANEL (B): COMMENSAL EXPLOITS PATHOGEN
ax(2) = subplot(2,2,2);
set(gca,'Unit','Inches','Position',[3.5 3.25 2.5 2.5]);
contourf(a,f,Pstar2,minPstar:0.03:maxPstar,'ShowText','off','LineWidth',1)
axis square
axis([0 1 0.5 2.1])
colormap(ax(2), d(cmin:cmax,:))
hold on
patch([a fliplr(a)], [f2_2 min(ylim)*ones(size(f2_2))],Conly,'LineStyle','none') % below lower bound
patch([a fliplr(a)], [f1_2 max(ylim)*ones(size(f1_2))],Ponly,'LineStyle','none') % above upper bound
hold off
hold on
line(a,comp_rel2,'Color','white','LineStyle','--','LineWidth',2)
c = colorbar;
c.LineWidth = 1;
caxis([minPstar,maxPstar])
xlim([0 1])
ylim([0.5 2.1])
ax = gca;
ax.FontSize = 8;
ax.LineWidth = 1;
xlabel('antibiotic exposure, A','FontSize',9,'Color','k')
ylabel('commensal relative susceptibility, f','FontSize',9,'Color','k')
title({'Commensal'; 'exploits pathogen'},'FontWeight','Normal','FontSize',10,'Color','k')

% PANEL (C): PATHOGEN EXPLOITS COMMENSAL
ax(3) = subplot(2,2,3);
set(gca,'Unit','Inches','Position',[0.5 0.25 2.5 2.5]);
contourf(a,f,Pstar3,minPstar:0.03:maxPstar,'ShowText','off','LineWidth',1)
axis square
axis([0 1 0.5 2.1])
colormap(ax(3),d(cmin:cmax,:))
hold on
patch([a fliplr(a)], [fhigh_3 max(ylim)*ones(size(fhigh_3))],Ponly,'LineStyle','none') % above upper bound - will always be bound for P only
hold off
hold on
c = colorbar;
c.LineWidth = 1;
caxis([minPstar,maxPstar])
xlim([0 1])
ylim([0.5 2.1])
ax = gca;
ax.FontSize = 8;
ax.LineWidth = 1;
xlabel('antibiotic exposure, A','FontSize',9,'Color','k')
ylabel('commensal relative susceptibility, f','FontSize',9,'Color','k')
title({'Pathogen'; 'exploits commensal'},'FontWeight','Normal','FontSize',10,'Color','k')

% PANEL (D): MUTUALISM
ax(4) = subplot(2,2,4);
set(gca,'Unit','Inches','Position',[3.5 0.25 2.5 2.5]);
contourf(a,f,Pstar4,min(min(Pstar4)):0.06:max(max(Pstar4)),'ShowText','off','LineWidth',1)
axis square
axis([0 1 0.5 2.1])
colormap(ax(4),d(cmin:cmax,:))
hold on
patch([a fliplr(a)], [fhigh_4 max(ylim)*ones(size(fhigh_4))],Ponly,'LineStyle','none') % above upper bound - will always be bound for P only
hold off
hold on
c = colorbar;
c.LineWidth = 1;
xlim([0 1])
ylim([0.5 2.1])
ax = gca;
ax.FontSize = 8;
ax.LineWidth = 1;
xlabel('antibiotic exposure, A','FontSize',9,'Color','k')
ylabel('commensal relative susceptibility, f','FontSize',9,'Color','k')
title('Mutualism','FontWeight','Normal','FontSize',10,'Color','k')

% save figure
%exportgraphics(fig,'HeatMap.eps','Resolution',600,'ContentType','vector')
%exportgraphics(fig,'HeatMap_rp.eps','Resolution',600,'ContentType','vector')
%exportgraphics(fig,'HeatMap_rc.eps','Resolution',600,'ContentType','vector')
%exportgraphics(fig,'HeatMap_kp.eps','Resolution',600,'ContentType','vector')
%exportgraphics(fig,'HeatMap_kc.eps','Resolution',600,'ContentType','vector')



