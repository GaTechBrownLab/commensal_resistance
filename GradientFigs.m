%% Benefits of antibiotic resistance in commensals and the scope for resistance optimization
% Kristofer Wollein Waldetoft, Sarah Sundius, Rachel Kuske, Sam P. Brown

% Code by: Sarah Sundius
% August 18, 2021


%% dP*/df Gradient
% code to create Figure 4

clear; clc;

% parameters
rp = 0.5;
rc = 0.5;
kp = 1;
kc = 1;
x = 0.1;
A = 1;

% calculate gradient for example interaction coefficients 
% [competition; C exploits P; P exploits C; mutualism]
apc_ex = [0.8 0.8 -0.8 -0.8]'; 
acp_ex = [0.8 -0.8 0.8 -0.8]';  
dPdf_ex = (apc_ex.*kc*x*A)./(rc*(1-apc_ex.*acp_ex));

% alpha range 
a = -0.999:0.001:0.999;
[apc,acp] = meshgrid(a,a);

% gradient function (dP*/df)
dPdf = (apc.*kc*x*A)./(rc*(1-apc.*acp));

% set color limits 
clims_dPdf = [-0.5 0.5];    % restricted range used in manuscript
%clims_dPdf = [-1 1];       
%clims_dPdf = [min(min(dPdf)) max(max(dPdf))];

% create colormaps (red: dP*/df > 0, blue: dP*/df < 0)
redmap = [ones(256,1),linspace(1,0,256)',linspace(1,0,256)'];
bluemap = [linspace(1,0,256)',linspace(1,0,256)',ones(256,1)];

% make heat map
fig_dPdf = figure;
set(gcf,'Unit','Inches','Position',[0 0 4 4])

hold on
axis square
dPdfflip = flipud(dPdf);
imagesc(flipud(dPdfflip),'XData',a,'YData',a);
set(gca,'Layer','top')
box on
caxis(clims_dPdf)
ax = gca;
ax.FontSize = 8;
ax.LineWidth = 1;
colormap([flipud(bluemap);redmap])
c = colorbar;
c.LineWidth = 1;
%plot(a,1./a-0.4,'k-','LineWidth',1)    % pos bound of color scale if using [-0.5,0.5]
%plot(a,1./a+0.4,'k-','LineWidth',1)    % neg bound of color scale if using [-0.5,0.5]
%plot(a,1./a+0.2,'k-','LineWidth',1)   	% pos bound of color scale if using [-1,1]
%plot(a,1./a-0.2,'k-','LineWidth',1)    % neg bound of color scale if using [-1,1]
plot(zeros(size(a)),a,'k--','LineWidth',1)
plot(a,zeros(size(a)),'k--','LineWidth',1)
plot(apc_ex,acp_ex,'k.','MarkerSize',10)
dPdf_exl = {'A','B','C','D'};
dPdf_exl = cellstr(dPdf_exl);
text(apc_ex+0.05,acp_ex+0.05,dPdf_exl,'FontSize',9,'Color','k')
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-0.5','0','0.5','1'})
yticklabels({'-1','-0.5','0','0.5','1'})
xlim([-1 1])
ylim([-1 1])
xlabel('a_{pc}','FontSize',10,'Color','k')
ylabel('a_{cp}','FontSize',10,'Color','k')

% save figure
%exportgraphics(fig_dPdf,'dPdf_grad.eps','Resolution',600,'ContentType','vector')



