% Plots water level computed from ADCIRC harmonic data and water level
% computed from NOAA harmonic data

% REQUIRES adcirc_util: https://github.com/BrianOBlanton/adcirc_util

% YOU MUST RENAME YOUR FORT.51 TO FORT.53, read_adcirc_fort53
% FUNCTION DOES NOT READ FORT.51

%%
clear; clf; clc;

xlsfilename = 'harmonics.xls';
node = 1; % Node/recording station to plot

%% PLOT TIME PARAMETERS
dt = 1; %one hour timestep
days = 14;

dt2sec = 3600; %coverting dt to seconds
plot_totaltime = 24*days; %converting days to hours
x = (dt:dt:plot_totaltime); %x-axis
tsteps = length(x); %number of timesteps

%% FORMATTING MODELED DATA
fort53 = read_adcirc_fort53('fort.53');
modAmp = fort53.AMP(node,:); %model amplitude
modPhase = deg2rad(fort53.PHA(node,:)); %model phase
modSpeed = fort53.FREQ(1,:); %model frequency
numAdcConstituents = length(modAmp); %number of constituents

%% FORMATTING OBSERVED DATA
NOAA_harmonics = xlsread(xlsfilename, 1, '', 'basic');
obsAmp = NOAA_harmonics(:,1); %observed amplitude
obsPhase = deg2rad(NOAA_harmonics(:,2)); %observed phase
obsSpeed = deg2rad(NOAA_harmonics(:,3)) * (1/3600); %observed frequency
numConstituents = length(obsAmp);

%% CALCULATING ETA
eta_o = zeros(1,numConstituents);
eta_m = zeros(1,numAdcConstituents);
eta_NOAA = zeros(1,tsteps);
eta_model = zeros(1,tsteps);
for t = 1:tsteps
    for i = 1:numConstituents
        eta_o(i) = obsAmp(i) * cos(obsSpeed(i) * (t*dt2sec) + obsPhase(i));
    end
    for i = 1:numAdcConstituents
        eta_m(i) = modAmp(i) * cos(modSpeed(i) * (t*dt2sec) + modPhase(i));
    end
    eta_NOAA(t) = sum(eta_o);
    eta_model(t) = sum(eta_m(~isnan(eta_m)));
end

%% PLOTTING RESYNTHESIS
plot(x,eta_NOAA, 'Linewidth', 1,'Color',[0.00,0.00,0.00]);
hold on;
plot(x,eta_model, 'Linewidth', 1.5)
hold off;

xlabel('Time (days)')
ylabel('Water Surface Elevation (m, NAVD88)')
legend("Observed","Modeled")
grid on
