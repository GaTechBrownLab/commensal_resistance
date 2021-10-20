%% Benefits of antibiotic resistance in commensals and the scope for resistance optimization
% Kristofer Wollein Waldetoft, Sarah Sundius, Rachel Kuske, Sam P. Brown

% Code by: Sarah Sundius
% August 18, 2021


%% Qualitative Model (Ecological) Outcomes
% code to create Figure S5

clear; clc;

% parameters
rp = 0.5;
rc = 0.5;
kp = 1;
kc = 1;
x = 0.1;

% control parameters
%A = [1 0.5 0 0.5 1];
%f = [0.5 0.5 1 2 2];
A = [1 0 1];        % used for Fig S5
f = [0.5 1 2];      % used for Fig S5

% alpha range 
a = -1.5:0.001:1.5;

% state bounds 
apc_bound = @(fv,Av) (rc*kp*(rp - x*Av))/(rp*kc*(rc - x*fv*Av));
acp_bound = @(fv,Av) (rp*kc*(rc - x*fv*Av))/(rc*kp*(rp - x*Av));

% set colors
d = jet(200);
Ponly = d(180,:);
Conly = d(10,:);
coexist = [1 0.87 0];
bistable = [0.8 0.6 1];

% plot eco outcomes in (apc,acp)-space
fig = figure;
set(gcf,'unit','inches','position',[0 0 5.5 9]);

% PANEL (A): High commensal resistance, high antibiotic 
ax(1) = subplot(3,1,1);
set(gca,'unit','inches','position',[0.75 6.5 2.25 2.25])
apc = apc_bound(f(1),A(1));
acp = acp_bound(f(1),A(1));
hold on
axis square
axis([-1.5 1.5 -1.5 1.5])
box on
patch([min(xlim) apc apc min(xlim)],[acp acp max(ylim) max(ylim)],Ponly,'LineStyle','none')         % P dominance
patch([apc max(xlim) max(xlim) apc],[min(ylim) min(ylim) acp acp],Conly,'LineStyle','none')         % C dominance
patch([apc max(xlim) max(xlim) apc],[acp acp max(ylim) max(ylim)],bistable,'LineStyle','none')    	% bistability
patch([min(xlim) apc apc min(xlim)],[min(ylim) min(ylim) acp acp],coexist,'LineStyle','none')       % coexistence
set(gca,'Layer','top')
plot(zeros(size(a)),a,'k--','LineWidth',1)  % apc = 0
plot(a,zeros(size(a)),'k--','LineWidth',1)  % acp = 0
ax = gca;
ax.FontSize = 8;
ax.LineWidth = 1;
xticks([-1.5 -1 -0.5 0 0.5 1 1.5])
yticks([-1.5 -1 -0.5 0 0.5 1 1.5])
xticklabels({'-1.5','-1','-0.5','0','0.5','1','1.5'})
yticklabels({'-1.5','-1','-0.5','0','0.5','1','1.5'})
xlim([-1.5,1.5])
ylim([-1.5,1.5])
l1 = legend('pathogen dominance','commensal dominance','bistability','coexistence','FontSize',9,'TextColor','k');
set(l1,'unit','inches','position',[3.233 8.0625 1.875 0.6875])
title('f = 0.5, A = 1','FontWeight','Normal','FontSize',10,'Color','k')

% PANEL (B): Equal susceptibility or no antibiotic
% (equal susceptibility will give this result regardless of antibiotic)
ax(2) = subplot(3,1,2);
set(gca,'unit','inches','position',[0.75 3.5 2.25 2.25])
apc = apc_bound(f(2),A(2));
acp = acp_bound(f(2),A(2));
hold on
axis square
axis([-1.5 1.5 -1.5 1.5])
box on
patch([min(xlim) apc apc min(xlim)],[acp acp max(ylim) max(ylim)],Ponly,'LineStyle','none')         % P dominance
patch([apc max(xlim) max(xlim) apc],[min(ylim) min(ylim) acp acp],Conly,'LineStyle','none')         % C dominance
patch([apc max(xlim) max(xlim) apc],[acp acp max(ylim) max(ylim)],bistable,'LineStyle','none')     	% bistability
patch([min(xlim) apc apc min(xlim)],[min(ylim) min(ylim) acp acp],coexist,'LineStyle','none')       % coexistence
set(gca,'Layer','top')
plot(zeros(size(a)),a,'k--','LineWidth',1)  % apc = 0
plot(a,zeros(size(a)),'k--','LineWidth',1)  % acp = 0
ax = gca;
ax.FontSize = 8;
ax.LineWidth = 1;
xticks([-1.5 -1 -0.5 0 0.5 1 1.5])
yticks([-1.5 -1 -0.5 0 0.5 1 1.5])
xticklabels({'-1.5','-1','-0.5','0','0.5','1','1.5'})
yticklabels({'-1.5','-1','-0.5','0','0.5','1','1.5'})
xlim([-1.5,1.5])
ylim([-1.5,1.5])
title('f = 1 or A = 0','FontWeight','Normal','FontSize',10,'Color','k')
ylabel('a_{cp}','FontSize',10,'Color','k')

% PANEL (C): High commensal susceptibility, high antibiotic 
ax(3) = subplot(3,1,3);
set(gca,'unit','inches','position',[0.75 0.5 2.25 2.25])
apc = apc_bound(f(3),A(3));
acp = acp_bound(f(3),A(3));
hold on
axis square
axis([-1.5 1.5 -1.5 1.5])
box on
patch([min(xlim) apc apc min(xlim)],[acp acp max(ylim) max(ylim)],Ponly,'LineStyle','none')         % P dominance
patch([apc max(xlim) max(xlim) apc],[min(ylim) min(ylim) acp acp],Conly,'LineStyle','none')       	% C dominance
patch([apc max(xlim) max(xlim) apc],[acp acp max(ylim) max(ylim)],bistable,'LineStyle','none')    	% bistability
patch([min(xlim) apc apc min(xlim)],[min(ylim) min(ylim) acp acp],coexist,'LineStyle','none')       % coexistence
set(gca,'Layer','top')
plot(zeros(size(a)),a,'k--','LineWidth',1)  % apc = 0
plot(a,zeros(size(a)),'k--','LineWidth',1)  % acp = 0
ax = gca;
ax.FontSize = 8;
ax.LineWidth = 1;
xticks([-1.5 -1 -0.5 0 0.5 1 1.5])
yticks([-1.5 -1 -0.5 0 0.5 1 1.5])
xticklabels({'-1.5','-1','-0.5','0','0.5','1','1.5'})
yticklabels({'-1.5','-1','-0.5','0','0.5','1','1.5'})
xlim([-1.5,1.5])
ylim([-1.5,1.5])
title('f = 2, A = 1','FontWeight','Normal','FontSize',10,'Color','k')
xlabel('a_{pc}','FontSize',10,'Color','k')

% save figure
%exportgraphics(fig,'EcoOutcomes.eps','Resolution',600,'ContentType','vector')

