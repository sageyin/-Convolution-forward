%% The main function of the new algorithm, the calculation of the vertical component of the gravity field or its tensors
% Editor：Xianzhe Yin 2022/9/05 China University of Geosciences(Beijing)
clear 
close all
clc

%% ====== Defining the observation grid ======
ObsGrid.dn=2;ObsGrid.dw=2;ObsGrid.dz=2;
ObsGrid.Nmin=-100+ObsGrid.dn/2; ObsGrid.Nmax=100-ObsGrid.dn/2; % North-South
ObsGrid.Wmin=-150+ObsGrid.dw/2; ObsGrid.Wmax=150-ObsGrid.dw/2; % West-East
ObsGrid.zmin=0; ObsGrid.zmax=0;
ObsGrid.n=ObsGrid.Nmin:ObsGrid.dn:ObsGrid.Nmax;
ObsGrid.w=ObsGrid.Wmin:ObsGrid.dw:ObsGrid.Wmax;
ObsGrid.z=ObsGrid.zmin:ObsGrid.dz:ObsGrid.zmax; 
[ObsGrid.W,ObsGrid.N,ObsGrid.Z]=meshgrid(ObsGrid.w,ObsGrid.n,ObsGrid.z);

%% ====== Rectangular model construction ======
SouceGrid.dn=ObsGrid.dn;SouceGrid.dw=ObsGrid.dw;SouceGrid.dz=ObsGrid.dz;
SouceGrid.Nmin=-100; SouceGrid.Nmax=100; % North-South
SouceGrid.Wmin=-150; SouceGrid.Wmax=150; % West-East
SouceGrid.zmin=0; SouceGrid.zmax=100;
SouceGrid.n=SouceGrid.Nmin+SouceGrid.dn/2:SouceGrid.dn:SouceGrid.Nmax-SouceGrid.dn/2;
SouceGrid.w=SouceGrid.Wmin+SouceGrid.dw/2:SouceGrid.dw:SouceGrid.Wmax-SouceGrid.dw/2;
SouceGrid.z=SouceGrid.zmin+SouceGrid.dz/2:SouceGrid.dz:SouceGrid.zmax-SouceGrid.dz/2; 
[SouceGrid.W,SouceGrid.N,SouceGrid.Z]=meshgrid(SouceGrid.w,SouceGrid.n,SouceGrid.z);% Dissection of the underground grid
SouceGrid.density=zeros(size(SouceGrid.W));

logp=logical(SouceGrid.W<=80 & SouceGrid.W>=40 & SouceGrid.N<=80 & SouceGrid.N>=40 ...
                      & SouceGrid.Z<=20 & SouceGrid.Z>=10 );    % Rectangle 50*50*20 Center position: (0, 0, 35)
Souce.W=SouceGrid.W(logp);
Souce.N=SouceGrid.N(logp );
Souce.Z=SouceGrid.Z(logp);
Souce.density=1000; % unit:kg/m^3
SouceGrid.density(logp)=Souce.density;

%% ====== Forward modelling of gravity field by our method ======
tic
dr=[SouceGrid.dw,SouceGrid.dn,SouceGrid.dz];
r=[0,0,SouceGrid.dz/2]; % Upward is positive
t=[0,0,1];
t(1)=size(ObsGrid.W,1);t(2)=size(ObsGrid.W,2);
g =GraconvelP(SouceGrid.density,dr,r,t,'gz');
toc

%% ====== Unit Conversions ======
g=g*10^5;  %   m/s^2 converted to mGal
%g=g*10^6; %   m/s^2 converted to g.u
%g=g*1e9;  %   s^(-2) to E

%% ====== load data from spatial domain calculations ======
g0=load('mode01',"gg");  % Data obtained from spatial domain calculations
er=g0.gg-g;              % Calculation of the difference between the two methods

%% ====== Visualization ======
figure()
contourf(ObsGrid.w,ObsGrid.n,g)
colormap('jet');colorbar
axis equal
xlabel('West-East(m)');
ylabel('South-North(m)')
title('Gravity field calculated by our method')

figure()
contourf(ObsGrid.w,ObsGrid.n,er,50,'LineStyle','none');
colormap('jet');colorbar
axis equal
xlabel('West-East(m)');
ylabel('South-North(m)')
title('The difference between the two methods')

figure()
histogram(er)
xlabel('misfit(mGal)');
ylabel('points');
title('Deviation statistics')




