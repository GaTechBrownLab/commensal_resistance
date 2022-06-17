%% Benefits of antibiotic resistance in commensals and the scope for resistance optimization
% Kristofer Wollein Waldetoft, Sarah Sundius, Rachel Kuske, Sam P. Brown

% Code by: Sarah Sundius
% May 17, 2022

clear; clc;


%% Find equilibria in no antibiotic case to get initial conditions
% for our model

% Model parameters
rp = 0.5;  
rc = 0.5;  
kp = 1;    
kc = 1;    
x = 0.1;
f = 1;
A = 0;

apc = 0.8;  acp = 0.8;      %Competition

% Equilibria (P, C) in absence of antibiotic
%Eq1 = [0, 0];                      %Joint extinction
Eq2 = [kp*(rp - x*A)/rp, 0];        %Pathogen dominance
%Eq3 = [0, kc*(rc - x*f*A)/rc];     %Commensal dominance 
Eq4 = [(rc*kp*(rp - x*A) - apc*rp*kc*(rc - x*f*A))/(rp*rc*(1 - apc*acp)), ...
    (rp*kc*(rc - x*f*A) - acp*rc*kp*(rp - x*A))/(rp*rc*(1 - apc*acp))]; %Coexistence

% Set initial conditions for simulations (1,1), (1,2), and (1,3)
Peq11 = Eq2(1,1);
Peq123 = Eq4(1,1);
Ceq123 = Eq4(1,2);


%% Simulate no antibiotic case to get initial conditions
% for resource explicit model

% Resource Explicit Model - Pathogen Only 
% Model parameters  
rP = 0.75; 
gP = 1.5;
aP = 0.2;
kP = 0.18;
lP = 0.45;
d = 0.3;
R = 0.0107;

D = 0.05;
S_hat = 1;
x = 0.1;
A = 0;

% Simulation parameters
dt = 0.001;
t0 = 0; 
tf = 10000;     %want long time behavior, i.e. steady state
tsave = t0:dt:tf;

Psave = zeros(1,length(tsave));
Ssave = zeros(1,length(tsave));

% Initial conditions - doesn't matter for now, just needs to go to
% coexistence equilibrium
P0 = 0.1;
S0 = 1;

Psave(1,1) = P0;
Ssave(1,1) = S0;

% Run simulation
for t = 1:length(tsave)-1

    Pcurr = Psave(1,t);
    Scurr = Ssave(1,t);
        
    if Scurr < 0; Scurr = 0; end

    % Dynamics 
    dP = Pcurr*((rP*Scurr)/(aP+Scurr) + (kP*R)/(lP+R) - d*R - x*A - D);
    dS = D*(S_hat - Scurr) - (1/gP)*(Scurr/(Scurr + aP))*rP*Pcurr;

    % Update
    Psave(1,t+1) = Pcurr + dP*dt;
    Ssave(1,t+1) = Scurr + dS*dt;
        
end

% Save end values to use as initial conditions for simulation (2,1)
Peq21 = Psave(1,end);
Seq21 = Ssave(1,end);

% Plot - just to check @ equilibrium
%figure

%subplot(2,1,1)
%hold on
%plot(tsave,Psave(1,:),'k-','LineWidth',1.5)
%xlabel('Time, t')
%ylabel('Pathogen Density, P(t)')
%title('Resource Explicit Model, A=0')

%subplot(2,1,2)
%hold on
%plot(tsave,Ssave(1,:),'k-','LineWidth',1.5)
%xlabel('Time, t')
%ylabel('Resource Concentration, R(t)')


% Resource Explicit Model - Pathogen + Commensal
% Model parameters
rP = 0.75;
rC = 0.89; 
gP = 1.5;
gC = 1.53;
aP = 0.2;
aC = 0.22;
kP = 0.18;
lP = 0.45; 
d = 0.3;
R = 0.0107;

D = 0.05;
S_hat = 1;
x = 0.1;
f = 1;
A = 0;

% Simulation parameters
dt = 0.001;
t0 = 0; 
tf = 10000;     %want long time behavior, i.e. steady state
tsave = t0:dt:tf;

Psave = zeros(1,length(tsave));
Csave = zeros(1,length(tsave));
Ssave = zeros(1,length(tsave));

% Initial conditions - doesn't matter for now, just needs to go to
% coexistence equilibrium
P0 = 1;
C0 = 1;
S0 = 1;

Psave(1,1) = P0;
Csave(1,1) = C0;
Ssave(1,1) = S0;

% Run simulation
for t = 1:length(tsave)-1

    Pcurr = Psave(1,t);
    Ccurr = Csave(1,t);
    Scurr = Ssave(1,t);
        
    if Scurr < 0; Scurr = 0; end

    % Dynamics 
    dP = Pcurr*((rP*Scurr)/(aP+Scurr) + (kP*R)/(lP+R) - d*R - x*A - D);
    dC = Ccurr*((rC*Scurr)/(aC+Scurr) - d*R - x*f*A - D);
    dS = D*(S_hat - Scurr) - (1/gP)*(Scurr/(Scurr + aP))*rP*Pcurr - (1/gC)*(Scurr/(Scurr + aC))*rC*Ccurr;

    % Update
    Psave(1,t+1) = Pcurr + dP*dt;
    Csave(1,t+1) = Ccurr + dC*dt;
    Ssave(1,t+1) = Scurr + dS*dt;
        
end

% Save end values to use as initial conditions for simulations (2,2) and
% (2,3)
Peq223 = Psave(1,end);
Ceq223 = Csave(1,end);
Seq223 = Ssave(1,end);

% Plot - just to check @ equilibrium
%figure

%subplot(3,1,1)
%hold on
%plot(tsave,Psave(1,:),'k-','LineWidth',1.5)
%xlabel('Time, t')
%ylabel('Pathogen Density, P(t)')
%title('Resource Explicit Model, A=0')

%subplot(3,1,2)
%hold on
%plot(tsave,Csave(1,:),'k-','LineWidth',1.5)
%xlabel('Time, t')
%ylabel('Commensal Density, C(t)')

%subplot(3,1,3)
%hold on
%plot(tsave,Ssave(1,:),'k-','LineWidth',1.5)
%xlabel('Time, t')
%ylabel('Resource Concentration, S(t)')


%% Simulate no antibiotic case to get initial conditions (i.e. coexistence)
% for spatial extension model 

% Spatially Extended Model - Pathogen Only 
% Set up grid
xy0=0; xyf=100;     %x,y interval (same in ea direction)
t0=0; tf=1000;      %t interval - expand to get long time behavior

Mxy=200;            %# steps in x,y direction (same in ea direction)
MM=Mxy^2;           %total # of "boxes" in grid
N=10000;            %# time steps

dxy=(xyf-xy0)/Mxy;  
dt=(tf-t0)/N;

% Model parameters
rp = 0.5;  
kp = 1;    
x = 0.1;
A = 0;

DP = 0.01;

RP = DP*dt/(dxy^2);

% Initial conditions 
U0 = zeros(1,2*MM);
ZZZP = randsample(MM,MM/10);    %pick where to put pathogen, 10% of locations
U0(1,ZZZP) = 10;                %put density of 10 at those locations for P

UP = zeros(Mxy,Mxy);
UP(1:MM) = U0(1,1:MM);

% Save initial spatial distribution
P0 = UP;
Pinitial = sum(sum(P0));

Psave = zeros(1,N);
Psave(1,1) = Pinitial/MM;

% Set up coefficient matrices for step 1
P1L=diag(ones(Mxy,1)+RP*ones(Mxy,1))-diag((RP/2)*ones(Mxy-1,1),1);
P1L=P1L+diag((-RP/2)*ones(Mxy-1,1),-1);
P1L(1,2) = -RP; %for BC
P1L(Mxy,Mxy-1) = -RP;
P1R=diag(ones(Mxy,1)-RP*ones(Mxy,1))+diag((RP/2)*ones(Mxy-1,1),1);
P1R=P1R+diag((RP/2)*ones(Mxy-1,1),-1);
P1R(1,2) = RP; %for BC
P1R(Mxy,Mxy-1) = RP;

% Set up coefficient matrices for step 2
P2L=diag(ones(Mxy,1)+RP*ones(Mxy,1))-diag((RP/2)*ones(Mxy-1,1),1);
P2L=P2L+diag((-RP/2)*ones(Mxy-1,1),-1);
P2L(1,2) = -RP; %for BC
P2L(Mxy,Mxy-1) = -RP;
P2R=diag(ones(Mxy,1)-RP*ones(Mxy,1))+diag((RP/2)*ones(Mxy-1,1),1);
P2R=P2R+diag((RP/2)*ones(Mxy-1,1),-1);
P2R(1,2) = RP; %for BC
P2R(Mxy,Mxy-1) = RP;

% Run simulation

% First iteration (n=0) of reaction dynamics
fP = rp.*UP.*(1 - UP./kp) - x.*A.*UP;

for j=1:N-1

    UP(UP<0) = 0;

    % Reaction dyn @ t=n+1/2
    fPnew = rp.*UP.*(1 - UP./kp) - x.*A.*UP;
    fPhalf = (fP+fPnew)/2;
    fPstar = (dt/2)*fPhalf;

    % Step 1
    UPstar = P1L\(UP*P1R + fPstar);

    % Step 2
    UPnext = (P2R*UPstar + fPstar)*inv(P2L);

    % Update U,f
    UP = UPnext;
    fP = fPnew;
  
    Psave(1,j+1) = sum(sum(UP))/MM;

end

% Save end values to use as initial conditions for simulation (3,1) - want
% the grid NOT spatial average
Peq31 = UP;

% Plot - just to check @ equilibrium
%tsave = t0:tf/N:tf;
%figure
%hold on
%plot(tsave(1,1:end-1),Psave(1,:),'k-','LineWidth',1.5)
%xlabel('Time, t')
%ylabel('Pathogen Density, P(t)')
%title('Spatial Extension Model, A=0')


% Spatially Extended Model - Pathogen + Commensal 
% Set up grid
xy0=0; xyf=100;     %x,y interval (same in ea direction)
t0=0; tf=1000;      %t interval - expand to get long time behavior

Mxy=200;            %# steps in x,y direction (same in ea direction)
MM=Mxy^2;           %total # of "boxes" in grid
N=10000;            %# time steps

dxy=(xyf-xy0)/Mxy;  
dt=(tf-t0)/N;

% Model parameters
rp = 0.5;  
rc = 0.5;  
kp = 1;    
kc = 1;    
x = 0.1;
f = 1;
A = 0;

apc = 0.8;  acp = 0.8;      %Competition

DP = 0.01;
DC = 0.01;

RP = DP*dt/(dxy^2);
RC = DC*dt/(dxy^2);

% Initial conditions 
U0 = zeros(1,2*MM);
ZZZP = randsample(MM,MM/10);        %pick where to put pathogen, 10% of locations
ZZZC = randsample(MM+1:2*MM,MM/10); %pick where to put commensal, 10% of locations
U0(1,ZZZP) = 1;                     %put density of 1 at those locations for P
U0(1,ZZZC) = 10;                    %put density of 10 at those locations for C

UP = zeros(Mxy,Mxy);
UC = zeros(Mxy,Mxy);
UP(1:MM) = U0(1,1:MM);
UC(1:MM) = U0(1,MM+1:2*MM);

% Save initial spatial distribution
P0 = UP;
C0 = UC;
Pinitial = sum(sum(P0));
Cinitial = sum(sum(C0));

Psave = zeros(1,N);
Csave = zeros(1,N);
Psave(1,1) = Pinitial/MM;
Csave(1,1) = Cinitial/MM;

% Set up coefficient matrices for step 1
P1L=diag(ones(Mxy,1)+RP*ones(Mxy,1))-diag((RP/2)*ones(Mxy-1,1),1);
P1L=P1L+diag((-RP/2)*ones(Mxy-1,1),-1);
P1L(1,2) = -RP; %for BC
P1L(Mxy,Mxy-1) = -RP;
P1R=diag(ones(Mxy,1)-RP*ones(Mxy,1))+diag((RP/2)*ones(Mxy-1,1),1);
P1R=P1R+diag((RP/2)*ones(Mxy-1,1),-1);
P1R(1,2) = RP; %for BC
P1R(Mxy,Mxy-1) = RP;

C1L=diag(ones(Mxy,1)+RC*ones(Mxy,1))-diag((RC/2)*ones(Mxy-1,1),1);
C1L=C1L+diag((-RC/2)*ones(Mxy-1,1),-1);
C1L(1,2) = -RC; %for BC
C1L(Mxy,Mxy-1) = -RC;
C1R=diag(ones(Mxy,1)-RC*ones(Mxy,1))+diag((RC/2)*ones(Mxy-1,1),1);
C1R=C1R+diag((RC/2)*ones(Mxy-1,1),-1);
C1R(1,2) = RC; %for BC
C1R(Mxy,Mxy-1) = RC;

% Set up coefficient matrices for step 2
P2L=diag(ones(Mxy,1)+RP*ones(Mxy,1))-diag((RP/2)*ones(Mxy-1,1),1);
P2L=P2L+diag((-RP/2)*ones(Mxy-1,1),-1);
P2L(1,2) = -RP; %for BC
P2L(Mxy,Mxy-1) = -RP;
P2R=diag(ones(Mxy,1)-RP*ones(Mxy,1))+diag((RP/2)*ones(Mxy-1,1),1);
P2R=P2R+diag((RP/2)*ones(Mxy-1,1),-1);
P2R(1,2) = RP; %for BC
P2R(Mxy,Mxy-1) = RP;

C2L=diag(ones(Mxy,1)+RC*ones(Mxy,1))-diag((RC/2)*ones(Mxy-1,1),1);
C2L=C2L+diag((-RC/2)*ones(Mxy-1,1),-1);
C2L(1,2) = -RC; %for BC
C2L(Mxy,Mxy-1) = -RC;
C2R=diag(ones(Mxy,1)-RC*ones(Mxy,1))+diag((RC/2)*ones(Mxy-1,1),1);
C2R=C2R+diag((RC/2)*ones(Mxy-1,1),-1);
C2R(1,2) = RC; %for BC
C2R(Mxy,Mxy-1) = RC;

% Run simulation

% First iteration (n=0) of reaction dynamics
fP = rp.*UP.*(1 - (UP + apc.*UC)./kp) - x.*A.*UP;
fC = rc.*UC.*(1 - (UC + acp.*UP)./kc) - x.*f.*A.*UC;

for j=1:N-1

    UP(UP<0) = 0;
    UC(UC<0) = 0;

    % Reaction dyn @ t=n+1/2
    fPnew = rp.*UP.*(1 - (UP + apc.*UC)./kp) - x.*A.*UP;
    fCnew = rc.*UC.*(1 - (UC + acp.*UP)./kc) - x.*f.*A.*UC;
    fPhalf = (fP+fPnew)/2;
    fChalf = (fC+fCnew)/2;
    fPstar = (dt/2)*fPhalf;
    fCstar = (dt/2)*fChalf;

    % Step 1
    UPstar = P1L\(UP*P1R + fPstar);
    UCstar = C1L\(UC*C1R + fCstar);

    % Step 2
    UPnext = (P2R*UPstar + fPstar)*inv(P2L);
    UCnext = (C2R*UCstar + fCstar)*inv(C2L);

    % Update U,f
    UP = UPnext;
    UC = UCnext;
    fP = fPnew;
    fC = fCnew;

    Psave(1,j+1) = sum(sum(UP))/MM;
    Csave(1,j+1) = sum(sum(UC))/MM;

end 

% Save end values to use as initial conditions for simulation (3,2) and (3,3) - want
% the grid NOT spatial average
Peq323 = UP;
Ceq323 = UC;

% Plot - just to check @ equilibrium
%tsave = t0:tf/N:tf;
%figure

%subplot(2,1,1)
%hold on
%plot(tsave(1,1:end-1),Psave(1,:),'k-','LineWidth',1.5)
%xlabel('Time, t')
%ylabel('Pathogen Density, P(t)')
%title('Spatial Extension Model, A=0')

%subplot(2,1,2)
%hold on
%plot(tsave(1,1:end-1),Csave(1,:),'k-','LineWidth',1.5)
%xlabel('Time, t')
%ylabel('Commensal Density, C(t)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create Figure 5
% Panel 1: Our model just pathogen (Figure 5A)
% Panel 2: Our model high commensal susceptibility (f=2) (Figure 5B)
% Panel 3: Our model high commensal resistance (f=0.5) (Figure 5C)
% Repeat 1-3 in next row with resource explicit (Figures 5D-F)
% Repeat 1-3 in next row with spatial (Figures 5G-I)

% All time series are for the competition scenarios

% Initialize figure
T = tiledlayout(3,3,'TileSpacing','Compact');
set(gcf,'Unit','Inches','Position',[0,0,6.75,6.75]);
xlabel(T,'Time, t','FontSize',12,'Color','k');
ylabel(T,'Pathogen density, P(t)','FontSize',12,'Color','k')
c = colormap(jet(200));
colors = c(20:10:80,:);

% Simulation (1,1) - Our Model (just pathogen) 
% Model parameters
rp = 0.5;   
kp = 1;    
x = 0.1;
Avals = 0:0.25:1;

% Simulation parameters
dt = 0.001;
ti = 0;
tf = 100;
tsave = ti:dt:tf;

Psave11 = zeros(length(Avals),length(tsave));
%P0 = 0.1;      %rare P IC
P0 = Peq11;     %equilibrium P IC
Psave11(:,1) = P0;

% Run simulation
for a = 1:length(Avals)
    
    % Set antibiotic exposure
    A = Avals(1,a);
    
    % Simulate
    for i = 1:length(tsave)-1

        % Save current values
        Pcurr = Psave11(a,i);

        % Calculate density @ next time point
        dP = rp*Pcurr*(1 - Pcurr/kp) - x*A*Pcurr;

        Pnew = Pcurr + dP*dt;
   
        % Update
        Psave11(a,i+1) = Pnew;
       
    end

end

% Plot 
ax1 = nexttile;
hold on
box on
plot(tsave,Psave11(1,:),'k-','LineWidth',1.5,'DisplayName',['A =' num2str(Avals(1,1))])
for i = 2:length(Avals)
    plot(tsave,Psave11(i,:),'Color',colors(i,:),'LineWidth',1.5)
end
ylim([0.7 1.1])
title('A','FontWeight','normal','FontSize',12,'Color','k')
ax1.TitleHorizontalAlignment = 'left';
ax1.XTick = [0,25,50,75,100];
ax1.XTickLabel = [];
ax1.YTick = [0.7,0.9,1.1];
ax1.YTickLabel = {'0.7','0.9','1.1'};
ax1.YAxis.FontSize = 9;
ax1.LineWidth = 1;


% Simulation (1,2) - Our Model (susceptible commensal, f=2)
% Model parameters
rp = 0.5;   
rc = 0.5;   
kp = 1;     
kc = 1;    
x = 0.1;
f = 2;      %0.5, 1.2, 1.3, 2
Avals = 0:0.25:1;

apc = 0.8;  acp = 0.8;      %Competition

% Simulation parameters
dt = 0.001;
ti = 0;
tf = 100;
tsave = ti:dt:tf;

Psave12 = zeros(length(Avals),length(tsave));
Csave12 = zeros(length(Avals),length(tsave));
%P0 = 0.1;      %rare P IC 
%C0 = 1;        %rare P IC
P0 = Peq123;    %equilibrium P IC (coexistence)
C0 = Ceq123;    %equilibrium C IC (coexistence)
Psave12(:,1) = P0;
Csave12(:,1) = C0;

% Run simulation
for a = 1:length(Avals)
    
    % Set antibiotic exposure
    A = Avals(1,a);
    
    % Simulate
    for i = 1:length(tsave)-1

        % Save current values
        Pcurr = Psave12(a,i);
        Ccurr = Csave12(a,i);

        % Calculate density @ next time point
        dP = rp*Pcurr*(1 - (Pcurr + apc*Ccurr)/kp) - x*A*Pcurr;
        dC = rc*Ccurr*(1 - (Ccurr + acp*Pcurr)/kc) - x*f*A*Ccurr;

        Pnew = Pcurr + dP*dt;
        Cnew = Ccurr + dC*dt;

        % Update
        Psave12(a,i+1) = Pnew;
        Csave12(a,i+1) = Cnew;

    end

end

% Plot 
ax2 = nexttile;
hold on
box on
plot(tsave,Psave12(1,:),'k-','LineWidth',1.5,'DisplayName',['A =' num2str(Avals(1,1))])
for i = 2:length(Avals)
    plot(tsave,Psave12(i,:),'Color',colors(i,:),'LineWidth',1.5)
end
ylim([0.1 0.9])
title('B','FontWeight','normal','FontSize',12,'Color','k')
ax2.TitleHorizontalAlignment = 'left';
ax2.XTick = [0,25,50,75,100];
ax2.XTickLabel = [];
ax2.YTick = [0.1,0.5,0.9];
ax2.YTickLabel = {'0.1','0.5','0.9'};
ax2.YAxis.FontSize = 9;
ax2.LineWidth = 1;


% Simulation (1,3) - Our Model (resistant commensal, f=0.5)
% Model parameters
rp = 0.5;   
rc = 0.5;   
kp = 1;     
kc = 1;     
x = 0.1;
f = 0.5;    %0.5, 1.2, 1.3, 2
Avals = 0:0.25:1;

apc = 0.8;  acp = 0.8;      %Competition

% Simulation parameters
dt = 0.001;
ti = 0;
tf = 100;
tsave = ti:dt:tf;

Psave13 = zeros(length(Avals),length(tsave));
Csave13 = zeros(length(Avals),length(tsave));
%P0 = 0.1;      %rare P IC
%C0 = 1;        %rare P IC
P0 = Peq123;    %equilibrium P IC (coexistence)
C0 = Ceq123;    %equilibrium C IC (coexistence)
Psave13(:,1) = P0;
Csave13(:,1) = C0;

% Run simulation
for a = 1:length(Avals)
    
    % Set antibiotic exposure
    A = Avals(1,a);
    
    % Simulate
    for i = 1:length(tsave)-1

        % Save current values
        Pcurr = Psave13(a,i);
        Ccurr = Csave13(a,i);

        % Calculate density @ next time point
        dP = rp*Pcurr*(1 - (Pcurr + apc*Ccurr)/kp) - x*A*Pcurr;
        dC = rc*Ccurr*(1 - (Ccurr + acp*Pcurr)/kc) - x*f*A*Ccurr;

        Pnew = Pcurr + dP*dt;
        Cnew = Ccurr + dC*dt;

        % Update
        Psave13(a,i+1) = Pnew;
        Csave13(a,i+1) = Cnew;

    end

end

% Plot 
ax3 = nexttile;
hold on
box on
plot(tsave,Psave13(1,:),'k-','LineWidth',1.5,'DisplayName',['A =' num2str(Avals(1,1))])
for i = 2:length(Avals)
    plot(tsave,Psave13(i,:),'Color',colors(i,:),'LineWidth',1.5)
end
ylim([0.1 0.9])
title('C','FontWeight','normal','FontSize',12,'Color','k')
ax3.TitleHorizontalAlignment = 'left';
ax3.XTick = [0,25,50,75,100];
ax3.XTickLabel = [];
ax3.YTick = [0.1,0.5,0.9];
ax3.YTickLabel = {'0.1','0.5','0.9'};
ax3.YAxis.FontSize = 9;
ax3.LineWidth = 1;


% Simulation (2,1) - Resource Explicit (just pathogen)
% Model parameters
rP = 0.75;         
gP = 1.5;                 
aP = 0.2;
kP = 0.18;
lP = 0.45;
d = 0.3;
R = 0.0107;

D = 0.05;
S_hat = 1;

x = 0.1;
Avals = 0:0.25:1;

% Simulation parameters
dt = 0.001;
t0 = 0; 
tf = 100;
tsave = t0:dt:tf;

Psave21 = zeros(length(Avals),length(tsave));
Ssave21 = zeros(length(Avals),length(tsave));

% Initial conditions
P0 = Peq21; %ran to t=10000 and took end point (from above)
S0 = Seq21;

Psave21(:,1) = P0;
Ssave21(:,1) = S0;

% Run simulation
for a = 1:length(Avals)
    
    % Set antibiotic exposure
    A = Avals(1,a);

    for t = 1:length(tsave)-1

        Pcurr = Psave21(a,t);
        Scurr = Ssave21(a,t);
        
        if Scurr < 0; Scurr = 0; end

        % Dynamics 
        dP = Pcurr*((rP*Scurr)/(aP+Scurr) + (kP*R)/(lP+R) - d*R - x*A - D);
        dS = D*(S_hat - Scurr) - (1/gP)*(Scurr/(Scurr + aP))*rP*Pcurr;

        % Update
        Psave21(a,t+1) = Pcurr + dP*dt;
        Ssave21(a,t+1) = Scurr + dS*dt;
        
    end

end

% Plot 
ax4 = nexttile;
hold on
box on
plot(tsave,Psave21(1,:),'k-','LineWidth',1.5,'DisplayName',['A =' num2str(Avals(1,1))])
for i = 2:length(Avals)
    plot(tsave,Psave21(i,:),'Color',colors(i,:),'LineWidth',1.5)
end
ylim([0.2 1.6])
title('D','FontWeight','normal','FontSize',12,'Color','k')
ax4.TitleHorizontalAlignment = 'left';
ax4.XTick = [0,25,50,75,100];
ax4.XTickLabel = [];
ax4.YTick = [0.2,0.9,1.6];
ax4.YTickLabel = {'0.2','0.9','1.6'};
ax4.YAxis.FontSize = 9;
ax4.LineWidth = 1;


% Simulation (2,2) - Resource Explicit (susceptible commensal, f=2)
% Model parameters
rP = 0.75;     
rC = 0.89;      
gP = 1.5;            
gC = 1.53;       
aP = 0.2;
aC = 0.22;
kP = 0.18;
lP = 0.45;
d = 0.3;
R = 0.0107;

D = 0.05;
S_hat = 1;

x = 0.1;
f = 2;            % 0.5, 1.2, 1.3, 2
Avals = 0:0.25:1;

% Simulation parameters
dt = 0.001;
t0 = 0; 
tf = 100;
tsave = t0:dt:tf;

Psave22 = zeros(length(Avals),length(tsave));
Csave22 = zeros(length(Avals),length(tsave));
Ssave22 = zeros(length(Avals),length(tsave));

% Initial conditions
P0 = Peq223;    %ran to t=10000 and took end point (from above)
C0 = Ceq223;
S0 = Seq223;

Psave22(:,1) = P0;
Csave22(:,1) = C0;
Ssave22(:,1) = S0;

% Run simulation
for a = 1:length(Avals)
    
    % Set antibiotic exposure
    A = Avals(1,a);

    for t = 1:length(tsave)-1

        Pcurr = Psave22(a,t);
        Ccurr = Csave22(a,t);
        Scurr = Ssave22(a,t);
        
        if Scurr < 0; Scurr = 0; end

        % Dynamics 
        dP = Pcurr*((rP*Scurr)/(aP+Scurr) + (kP*R)/(lP+R) - d*R - x*A - D);
        dC = Ccurr*((rC*Scurr)/(aC+Scurr) - d*R - x*f*A - D);
        dS = D*(S_hat - Scurr) - (1/gP)*(Scurr/(Scurr + aP))*rP*Pcurr - (1/gC)*(Scurr/(Scurr + aC))*rC*Ccurr;

        % Update
        Psave22(a,t+1) = Pcurr + dP*dt;
        Csave22(a,t+1) = Ccurr + dC*dt;
        Ssave22(a,t+1) = Scurr + dS*dt;
        
    end

end

% Plot 
ax5 = nexttile;
hold on
box on
plot(tsave,Psave22(1,:),'k-','LineWidth',1.5,'DisplayName',['A =' num2str(Avals(1,1))])
for i = 2:length(Avals)
    plot(tsave,Psave22(i,:),'Color',colors(i,:),'LineWidth',1.5)
end
ylim([0 1])
title('E','FontWeight','normal','FontSize',12,'Color','k')
ax5.TitleHorizontalAlignment = 'left';
ax5.XTick = [0,25,50,75,100];
ax5.XTickLabel = [];
ax5.YTick = [0,0.5,1];
ax5.YTickLabel = {'0','0.5','1'};
ax5.YAxis.FontSize = 9;
ax5.LineWidth = 1;


% Simulation (2,3) - Resource Explicit (resistant commensal, f=0.5)
% Model parameters
rP = 0.75;     
rC = 0.89;      
gP = 1.5;            
gC = 1.53;       
aP = 0.2;
aC = 0.22;
kP = 0.18;
lP = 0.45;
d = 0.3;
R = 0.0107;

D = 0.05;
S_hat = 1;

x = 0.1;
f = 0.5;            % 0.5, 1.2, 1.3, 2
Avals = 0:0.25:1;

% Simulation parameters
dt = 0.0001;
t0 = 0; 
tf = 100;
tsave = t0:dt:tf;

Psave23 = zeros(length(Avals),length(tsave));
Csave23 = zeros(length(Avals),length(tsave));
Ssave23 = zeros(length(Avals),length(tsave));

% Initial conditions
P0 = Peq223;    %ran to t=10000 and took end point (from above)
C0 = Ceq223;
S0 = Seq223;

Psave23(:,1) = P0;
Csave23(:,1) = C0;
Ssave23(:,1) = S0;

% Run simulation
for a = 1:length(Avals)
    
    % Set antibiotic exposure
    A = Avals(1,a);

    for t = 1:length(tsave)-1

        Pcurr = Psave23(a,t);
        Ccurr = Csave23(a,t);
        Scurr = Ssave23(a,t);
        
        if Scurr < 0; Scurr = 0; end

        % Dynamics 
        dP = Pcurr*((rP*Scurr)/(aP+Scurr) + (kP*R)/(lP+R) - d*R - x*A - D);
        dC = Ccurr*((rC*Scurr)/(aC+Scurr) - d*R - x*f*A - D);
        dS = D*(S_hat - Scurr) - (1/gP)*(Scurr/(Scurr + aP))*rP*Pcurr - (1/gC)*(Scurr/(Scurr + aC))*rC*Ccurr;

        % Update
        Psave23(a,t+1) = Pcurr + dP*dt;
        Csave23(a,t+1) = Ccurr + dC*dt;
        Ssave23(a,t+1) = Scurr + dS*dt;
        
    end

end

% Plot 
ax6 = nexttile;
hold on
box on
plot(tsave,Psave23(1,:),'k-','LineWidth',1.5,'DisplayName',['A =' num2str(Avals(1,1))])
for i = 2:length(Avals)
    plot(tsave,Psave23(i,:),'Color',colors(i,:),'LineWidth',1.5)
end
ylim([0 1])
title('F','FontWeight','normal','FontSize',12,'Color','k')
ax6.TitleHorizontalAlignment = 'left';
ax6.XTick = [0,25,50,75,100];
ax6.XTickLabel = [];
ax6.YTick = [0,0.5,1];
ax6.YTickLabel = {'0','0.5','1'};
ax6.YAxis.FontSize = 9;
ax6.LineWidth = 1;


% Simulation (3,1) - Our Model Spatial Extension (pathogen only)
% Crank-Nicolson Scheme
% Note: uses dx=dy

% 2D reaction-diffusion w/ 2 species LV reaction dynamics

% Set up grid
xy0=0; xyf=100;     %x,y interval (same in ea direction)
t0=0; tf=100;       %t interval

Mxy=200;            %# steps in x,y direction (same in ea direction)
MM=Mxy^2;           %total # of "boxes" in grid
N=1000;             %# time steps

dxy=(xyf-xy0)/Mxy;  
dt=(tf-t0)/N;

% Model parameters
rp = 0.5;  
kp = 1;     
x = 0.1;
Avals = 0:0.25:1;

DP = 0.01;

RP = DP*dt/(dxy^2);

% Save initial spatial distribution - from above, run simulation to t=1000
% and take spatial grid
UP = Peq31;
P0 = UP;
Pinitial = sum(sum(P0));

Psave31 = zeros(length(Avals),N);
Psave31(:,1) = Pinitial/MM;

% Set up coefficient matrices for step 1
P1L=diag(ones(Mxy,1)+RP*ones(Mxy,1))-diag((RP/2)*ones(Mxy-1,1),1);
P1L=P1L+diag((-RP/2)*ones(Mxy-1,1),-1);
P1L(1,2) = -RP; %for BC
P1L(Mxy,Mxy-1) = -RP;
P1R=diag(ones(Mxy,1)-RP*ones(Mxy,1))+diag((RP/2)*ones(Mxy-1,1),1);
P1R=P1R+diag((RP/2)*ones(Mxy-1,1),-1);
P1R(1,2) = RP; %for BC
P1R(Mxy,Mxy-1) = RP;

% Set up coefficient matrices for step 2
P2L=diag(ones(Mxy,1)+RP*ones(Mxy,1))-diag((RP/2)*ones(Mxy-1,1),1);
P2L=P2L+diag((-RP/2)*ones(Mxy-1,1),-1);
P2L(1,2) = -RP; %for BC
P2L(Mxy,Mxy-1) = -RP;
P2R=diag(ones(Mxy,1)-RP*ones(Mxy,1))+diag((RP/2)*ones(Mxy-1,1),1);
P2R=P2R+diag((RP/2)*ones(Mxy-1,1),-1);
P2R(1,2) = RP; %for BC
P2R(Mxy,Mxy-1) = RP;

% Run simulation
for a = 1:length(Avals)
    
    % Set antibiotic exposure
    Av = Avals(1,a); 
    if Av > 0 
        Agrad = linspace(1,0,200);      %simple gradient of 0.1 across one direction
        A = Av.*meshgrid(Agrad,1:Mxy);
    else
        A = Av;
    end

    % Reset initial condition
    UP = P0;

    % First iteration (n=0) of reaction dynamics
    fP = rp.*UP.*(1 - UP./kp) - x.*A.*UP;

    for j=1:N-1

        UP(UP<0) = 0;

        % Reaction dyn @ t=n+1/2
        fPnew = rp.*UP.*(1 - UP./kp) - x.*A.*UP;
        fPhalf = (fP+fPnew)/2;
        fPstar = (dt/2)*fPhalf;

        % Step 1
        UPstar = P1L\(UP*P1R + fPstar);

        % Step 2
        UPnext = (P2R*UPstar + fPstar)*inv(P2L);

        % Update U,f
        UP = UPnext;
        fP = fPnew;

        Psave31(a,j+1) = sum(sum(UP))/MM;

    end

end 

% Plot 
tsave = t0:tf/N:tf;
ax7 = nexttile;
hold on
box on
plot(tsave(1,1:end-1),Psave31(1,:),'k-','LineWidth',1.5,'DisplayName',['A =' num2str(Avals(1,1))])
for i = 2:length(Avals)
    plot(tsave(1,1:end-1),Psave31(i,:),'Color',colors(i,:),'LineWidth',1.5)
end
ylim([0.7 1.1])
title('G','FontWeight','normal','FontSize',12,'Color','k')
ax7.TitleHorizontalAlignment = 'left';
ax7.XTick = [0,25,50,75,100];
ax7.XTickLabel = {'0','','50','','100'};
xtickangle(0);
ax7.YTick = [0.7,0.9,1.1];
ax7.YTickLabel = {'0.7','0.9','1.1'};
ax7.XAxis.FontSize = 9;
ax7.YAxis.FontSize = 9;
ax7.LineWidth = 1;


% Simulation (3,2) - Our Model Spatial Extension (susceptible commensal, f=2)
% Crank-Nicolson Scheme
% Note: uses dx=dy

% 2D reaction-diffusion w/ 2 species LV reaction dynamics

% Set up grid
xy0=0; xyf=100;     %x,y interval (same in ea direction)
t0=0; tf=100;       %t interval

Mxy=200;            %# steps in x,y direction (same in ea direction)
MM=Mxy^2;           %total # of "boxes" in grid
N=1000;             %# time steps

dxy=(xyf-xy0)/Mxy;  
dt=(tf-t0)/N;

% Model parameters
rp = 0.5;   
rc = 0.5;   
kp = 1;     
kc = 1;     
x = 0.1;
f = 2;      %0.5, 1.2, 1.3, 2
Avals = 0:0.25:1;

apc = 0.8;  acp = 0.8;      %competition

DP = 0.01;
DC = 0.01;

RP = DP*dt/(dxy^2);
RC = DC*dt/(dxy^2);

% Save initial spatial distribution - from above, run simulation to t=1000
% and take spatial grid
UP = Peq323;
UC = Ceq323;
P0 = UP;
C0 = UC;
Pinitial = sum(sum(P0));
Cinitial = sum(sum(C0));

Psave32 = zeros(length(Avals),N);
Csave32 = zeros(length(Avals),N);
Psave32(:,1) = Pinitial/MM;
Csave32(:,1) = Cinitial/MM;

% Set up coefficient matrices for step 1
P1L=diag(ones(Mxy,1)+RP*ones(Mxy,1))-diag((RP/2)*ones(Mxy-1,1),1);
P1L=P1L+diag((-RP/2)*ones(Mxy-1,1),-1);
P1L(1,2) = -RP; %for BC
P1L(Mxy,Mxy-1) = -RP;
P1R=diag(ones(Mxy,1)-RP*ones(Mxy,1))+diag((RP/2)*ones(Mxy-1,1),1);
P1R=P1R+diag((RP/2)*ones(Mxy-1,1),-1);
P1R(1,2) = RP; %for BC
P1R(Mxy,Mxy-1) = RP;

C1L=diag(ones(Mxy,1)+RC*ones(Mxy,1))-diag((RC/2)*ones(Mxy-1,1),1);
C1L=C1L+diag((-RC/2)*ones(Mxy-1,1),-1);
C1L(1,2) = -RC; %for BC
C1L(Mxy,Mxy-1) = -RC;
C1R=diag(ones(Mxy,1)-RC*ones(Mxy,1))+diag((RC/2)*ones(Mxy-1,1),1);
C1R=C1R+diag((RC/2)*ones(Mxy-1,1),-1);
C1R(1,2) = RC; %for BC
C1R(Mxy,Mxy-1) = RC;

% Set up coefficient matrices for step 2
P2L=diag(ones(Mxy,1)+RP*ones(Mxy,1))-diag((RP/2)*ones(Mxy-1,1),1);
P2L=P2L+diag((-RP/2)*ones(Mxy-1,1),-1);
P2L(1,2) = -RP; %for BC
P2L(Mxy,Mxy-1) = -RP;
P2R=diag(ones(Mxy,1)-RP*ones(Mxy,1))+diag((RP/2)*ones(Mxy-1,1),1);
P2R=P2R+diag((RP/2)*ones(Mxy-1,1),-1);
P2R(1,2) = RP; %for BC
P2R(Mxy,Mxy-1) = RP;

C2L=diag(ones(Mxy,1)+RC*ones(Mxy,1))-diag((RC/2)*ones(Mxy-1,1),1);
C2L=C2L+diag((-RC/2)*ones(Mxy-1,1),-1);
C2L(1,2) = -RC; %for BC
C2L(Mxy,Mxy-1) = -RC;
C2R=diag(ones(Mxy,1)-RC*ones(Mxy,1))+diag((RC/2)*ones(Mxy-1,1),1);
C2R=C2R+diag((RC/2)*ones(Mxy-1,1),-1);
C2R(1,2) = RC; %for BC
C2R(Mxy,Mxy-1) = RC;

% Run simulation
for a = 1:length(Avals)
    
    % Set antibiotic exposure
    Av = Avals(1,a); 
    if Av > 0 
        Agrad = linspace(1,0,200);      %simple gradient of 0.1 across one direction
        A = Av.*meshgrid(Agrad,1:Mxy);
    else
        A = Av;
    end

    % Reset initial condition
    UP = P0;
    UC = C0;

    % First iteration (n=0) of reaction dynamics
    fP = rp.*UP.*(1 - (UP + apc.*UC)./kp) - x.*A.*UP;
    fC = rc.*UC.*(1 - (UC + acp.*UP)./kc) - x.*f.*A.*UC;

    for j=1:N-1

        UP(UP<0) = 0;
        UC(UC<0) = 0;

        % Reaction dyn @ t=n+1/2
        fPnew = rp.*UP.*(1 - (UP + apc.*UC)./kp) - x.*A.*UP;
        fCnew = rc.*UC.*(1 - (UC + acp.*UP)./kc) - x.*f.*A.*UC;
        fPhalf = (fP+fPnew)/2;
        fChalf = (fC+fCnew)/2;
        fPstar = (dt/2)*fPhalf;
        fCstar = (dt/2)*fChalf;

        % Step 1
        UPstar = P1L\(UP*P1R + fPstar);
        UCstar = C1L\(UC*C1R + fCstar);

        % Step 2
        UPnext = (P2R*UPstar + fPstar)*inv(P2L);
        UCnext = (C2R*UCstar + fCstar)*inv(C2L);

        % Update U,f
        UP = UPnext;
        UC = UCnext;
        fP = fPnew;
        fC = fCnew;

        Psave32(a,j+1) = sum(sum(UP))/MM;
        Csave32(a,j+1) = sum(sum(UC))/MM;

    end

end 

% Plot 
tsave = t0:tf/N:tf;
ax8 = nexttile;
hold on
box on
plot(tsave(1,1:end-1),Psave32(1,:),'k-','LineWidth',1.5,'DisplayName',['A =' num2str(Avals(1,1))])
for i = 2:length(Avals)
    plot(tsave(1,1:end-1),Psave32(i,:),'Color',colors(i,:),'LineWidth',1.5)
end
ylim([0.2 0.8])
title('H','FontWeight','normal','FontSize',12,'Color','k')
ax8.TitleHorizontalAlignment = 'left';
ax8.XTick = [0,25,50,75,100];
ax8.XTickLabel = {'0','','50','','100'};
ax8.XAxis.FontSize = 9;
xtickangle(0);
ax8.YTick = [0.2,0.5,0.8];
ax8.YTickLabel = {'0.2','0.5','0.8'};
ax8.YAxis.FontSize = 9;
ax8.LineWidth = 1;


% Simulation (3,3) - Our Model Spatial Extension (resistant commensal, f=0.5)
% Crank-Nicolson Scheme
% Note: uses dx=dy

% 2D reaction-diffusion w/ 2 species LV reaction dynamics

% Set up grid
xy0=0; xyf=100;     %x,y interval (same in ea direction)
t0=0; tf=100;       %t interval

Mxy=200;            %# steps in x,y direction (same in ea direction)
MM=Mxy^2;           %total # of "boxes" in grid
N=1000;             %# time steps

dxy=(xyf-xy0)/Mxy;  
dt=(tf-t0)/N;

% Model parameters
rp = 0.5;  
rc = 0.5;  
kp = 1;    
kc = 1;     
x = 0.1;
f = 0.5;      %0.5, 1.2, 1.3, 2
Avals = 0:0.25:1;

apc = 0.8;  acp = 0.8;      %Competition

DP = 0.01;
DC = 0.01;

RP = DP*dt/(dxy^2);
RC = DC*dt/(dxy^2);

% Save initial spatial distribution - from above, run simulation to t=1000
% and take spatial grid
UP = Peq323;
UC = Ceq323;
P0 = UP;
C0 = UC;
Pinitial = sum(sum(P0));
Cinitial = sum(sum(C0));

Psave33 = zeros(length(Avals),N);
Csave33 = zeros(length(Avals),N);
Psave33(:,1) = Pinitial/MM;
Csave33(:,1) = Cinitial/MM;

% Set up coefficient matrices for step 1
P1L=diag(ones(Mxy,1)+RP*ones(Mxy,1))-diag((RP/2)*ones(Mxy-1,1),1);
P1L=P1L+diag((-RP/2)*ones(Mxy-1,1),-1);
P1L(1,2) = -RP; %for BC
P1L(Mxy,Mxy-1) = -RP;
P1R=diag(ones(Mxy,1)-RP*ones(Mxy,1))+diag((RP/2)*ones(Mxy-1,1),1);
P1R=P1R+diag((RP/2)*ones(Mxy-1,1),-1);
P1R(1,2) = RP; %for BC
P1R(Mxy,Mxy-1) = RP;

C1L=diag(ones(Mxy,1)+RC*ones(Mxy,1))-diag((RC/2)*ones(Mxy-1,1),1);
C1L=C1L+diag((-RC/2)*ones(Mxy-1,1),-1);
C1L(1,2) = -RC; %for BC
C1L(Mxy,Mxy-1) = -RC;
C1R=diag(ones(Mxy,1)-RC*ones(Mxy,1))+diag((RC/2)*ones(Mxy-1,1),1);
C1R=C1R+diag((RC/2)*ones(Mxy-1,1),-1);
C1R(1,2) = RC; %for BC
C1R(Mxy,Mxy-1) = RC;

% Set up coefficient matrices for step 2
P2L=diag(ones(Mxy,1)+RP*ones(Mxy,1))-diag((RP/2)*ones(Mxy-1,1),1);
P2L=P2L+diag((-RP/2)*ones(Mxy-1,1),-1);
P2L(1,2) = -RP; %for BC
P2L(Mxy,Mxy-1) = -RP;
P2R=diag(ones(Mxy,1)-RP*ones(Mxy,1))+diag((RP/2)*ones(Mxy-1,1),1);
P2R=P2R+diag((RP/2)*ones(Mxy-1,1),-1);
P2R(1,2) = RP; %for BC
P2R(Mxy,Mxy-1) = RP;

C2L=diag(ones(Mxy,1)+RC*ones(Mxy,1))-diag((RC/2)*ones(Mxy-1,1),1);
C2L=C2L+diag((-RC/2)*ones(Mxy-1,1),-1);
C2L(1,2) = -RC; %for BC
C2L(Mxy,Mxy-1) = -RC;
C2R=diag(ones(Mxy,1)-RC*ones(Mxy,1))+diag((RC/2)*ones(Mxy-1,1),1);
C2R=C2R+diag((RC/2)*ones(Mxy-1,1),-1);
C2R(1,2) = RC; %for BC
C2R(Mxy,Mxy-1) = RC;

% Run simulation
for a = 1:length(Avals)
    
    % Set antibiotic exposure
    Av = Avals(1,a); 
    if Av > 0 
        Agrad = linspace(1,0,200);      %simple gradient of 0.1 across one direction
        A = Av.*meshgrid(Agrad,1:Mxy);
    else
        A = Av;
    end

    % Reset initial condition
    UP = P0;
    UC = C0;

    % First iteration (n=0) of reaction dynamics
    fP = rp.*UP.*(1 - (UP + apc.*UC)./kp) - x.*A.*UP;
    fC = rc.*UC.*(1 - (UC + acp.*UP)./kc) - x.*f.*A.*UC;

    for j=1:N-1

        UP(UP<0) = 0;
        UC(UC<0) = 0;

        % Reaction dyn @ t=n+1/2
        fPnew = rp.*UP.*(1 - (UP + apc.*UC)./kp) - x.*A.*UP;
        fCnew = rc.*UC.*(1 - (UC + acp.*UP)./kc) - x.*f.*A.*UC;
        fPhalf = (fP+fPnew)/2;
        fChalf = (fC+fCnew)/2;
        fPstar = (dt/2)*fPhalf;
        fCstar = (dt/2)*fChalf;

        % Step 1
        UPstar = P1L\(UP*P1R + fPstar);
        UCstar = C1L\(UC*C1R + fCstar);

        % Step 2
        UPnext = (P2R*UPstar + fPstar)*inv(P2L);
        UCnext = (C2R*UCstar + fCstar)*inv(C2L);

        % Update U,f
        UP = UPnext;
        UC = UCnext;
        fP = fPnew;
        fC = fCnew;

        Psave33(a,j+1) = sum(sum(UP))/MM;
        Csave33(a,j+1) = sum(sum(UC))/MM;

    end

end 

% Plot 
tsave = t0:tf/N:tf;
ax9 = nexttile;
hold on
box on
plot(tsave(1,1:end-1),Psave33(1,:),'k-','LineWidth',1.5,'DisplayName',['A =' num2str(Avals(1,1))])
for i = 2:length(Avals)
    plot(tsave(1,1:end-1),Psave33(i,:),'Color',colors(i,:),'LineWidth',1.5)
end
ylim([0.2 0.8])
title('I','FontWeight','normal','FontSize',12,'Color','k');
ax9.TitleHorizontalAlignment = "left";
ax9.XTick = [0,25,50,75,100];
ax9.XTickLabel = {'0','','50','','100'};
ax9.XAxis.FontSize = 9;
xtickangle(0);
ax9.YTick = [0.2,0.5,0.8];
ax9.YTickLabel = {'0.2','0.5','0.8'};
ax9.YAxis.FontSize = 9;
ax9.LineWidth = 1;


% Add a legend 
leg = legend({'A=0','A=0.25','A=0.5','A=0.75','A=1'},'Orientation','Horizontal','FontSize',9,'TextColor','k');
leg.Layout.Tile = 'North';

% Save figure
%exportgraphics(T,'numsims_0531.eps','Resolution',1200,'ContentType','vector')




