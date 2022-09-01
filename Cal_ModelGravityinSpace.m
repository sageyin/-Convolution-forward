%% The main function of the calculation of the vertical component of the gravity field or its tensors in space domain
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
                      & SouceGrid.Z<=20 & SouceGrid.Z>=10 );    % Rectangle 50*50*10 Center position: (0, 0, 35)
Souce.W=SouceGrid.W(logp);
Souce.N=SouceGrid.N(logp );
Souce.Z=SouceGrid.Z(logp);
Souce.density=1000; % unit:kg/m^3
SouceGrid.density(logp)=Souce.density;

%% ====== Forward modelling of gravity field by our method in space domain ======

tic
Style='gz';
Souce.Num=length(Souce.W);
gg=0;
for n=1:Souce.Num
    g=Cal_tranGraf(ObsGrid.N,ObsGrid.W,ObsGrid.Z,Souce.N(n),Souce.W(n),Souce.Z(n),SouceGrid.dn,SouceGrid.dw,SouceGrid.dz,Souce.density,Style);
    gg=gg+g;
end
toc

%% ====== Unit Conversions ======
gg=gg*10^5;   %  m/s^2 converted to mGal
%gg=gg*10^6;  %  m/s^2 converted to g.u
%gg=gg*1e9;   %  s^(-2) to E

%% ====== Visualization ======
figure()
contourf(ObsGrid.w,ObsGrid.n,gg)
colormap('jet');colorbar
axis equal
xlabel('West-East(m)');
ylabel('South-North(m)')
title('Gravity field calculated in space domain')

%% ====== Data Storage ======
save('mode01',"gg"); 
