winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels 
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight


% Define Constants

dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);

%Simulation parameters 

tSim = 200e-15
f = 100e12; % define the frequency
lambda = c_c/f; % calculate the wavelength

% Define region 
xMax{1} = 20e-6;
nx{1} = 200;
ny{1} = 0.75*nx{1};


Reg.n = 1;

% create u matrix 

mu{1} = ones(nx{1},ny{1})*c_mu_0;

% Define epsilon structure 
epi{1} = ones(nx{1},ny{1})*c_eps_0;

% Add "inclusion" 
% When this line is commented out, no cmopnent of the wave is reflected
% off of the boundry 

epi{1}(125:150,55:95)= c_eps_0*11.3;
epi{1}(27:30,117:119)= c_eps_0*11.3;
epi{1}(35:45,115:120)= c_eps_0*11.3;
epi{1}(75:100,55:95)= c_eps_0*11.3; 
epi{1}(20:30,115:117)= c_eps_0*11.3;
epi{1}(20:30,123:125)= c_eps_0*11.3;
epi{1}(15:20,115:125)= c_eps_0*11.3;




sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

dx = xMax{1}/nx{1};
dt = 0.25*dx/c_c;
nSteps = round(tSim/dt*2);
yMax = ny{1}*dx;
nsteps_lamda = lambda/dx

% Setup plot

movie = 1;
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];

% Define boundry conditions  
bc{1}.NumS = 4;

% Define EM wave source(s)
% Source 1
bc{1}.s(1).xpos = nx{1}/(4) + 10;
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;

% Source 1
bc{2}.s(2).xpos = nx{1}/(4) + 1;
bc{2}.s(2).type = 'ss';
bc{2}.s(2).fct = @PlaneWaveBC;

% Source 1
bc{3}.s(3).xpos = nx{1}/(4) + 1;
bc{3}.s(3).type = 'ss';
bc{3}.s(3).fct = @PlaneWaveBC;

% Source 1
bc{4}.s(4).xpos = nx{1}/(4) + 1;
bc{4}.s(4).type = 'ss';
bc{4}.s(4).fct = @PlaneWaveBC;

% mag = -1/c_eta_0;
mag = 1;
phi = 0;
omega = f*2*pi;
betap = 0;
t0 = 30e-15;
st = -0.05;
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

Plot.y0 = round(y0/dx);

% Specify the type of boundry 
bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg






