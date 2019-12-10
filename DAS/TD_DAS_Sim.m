% Michelle Sigona 
% 20191206 - Time-Domain Delay and Sum Beamforming (TD-DAS) Algorithm for
% Passive Acoustic Mapping (DAS-PAM) for simulated data.

close all;
clear;
clc;

load('SimData.mat');

%% Set Parameters
% Setup computational grid
x = -9:0.3:9;           % Lateral location of pixels [mm]
z = 10:1:40;            % Range location of pixel [mm]

% Origin location of emission in mm
Cav_x = 0;              % [mm]
Cav_z = 26;             % [mm]

ncycles = 20;           % Number of cycles in simulated data
Fo = 6;                 % Center frequency of simulated data [MHz]

% Set properties of medium
c = 1.5236;             % speed of sound [mm/us]
rho = 1e-9;             % density of water [kg/mm^3]

%% Setup Transducer
nChannels = 128;                % Number of channels that simultaneously and passively acquired acoustic emissions
ElementPitch = 0.2980;          % Spacing between L7-4 array elements in mm
Aperture = (-18.9230:ElementPitch:18.9230)';     % Position of elements in array
zarray = 0;                     % Range location of array (should always be zero with a linear array)
EleHeight = 6;                  % Element height [mm]
EleWidth = 0.25;                % Element width [mm]
Sl = EleHeight*EleWidth;        % Element surface area [mm^2]
S = Sl*nChannels;               % Active area of array [mm^2]

% Receive parameters
T = 30;                         % Nominal duration of the recorded emission [us]
fs = 40;                        % Sampling frequency [MHz]
dt = 1/fs;                      % Sampling period [us]
NIOI = round(fs/Fo*ncycles);    % Number of points in the interval of interest
TIOI = NIOI*dt;                 % Time duration of interval of interest [us]

TotalTimeI = size(RData,1);     % The total time duration of the simulated received signal (in indicies)
Time = (0:TotalTimeI-1)*dt;     % Time vector

%% Calculate Time Delays
% Distance from point source with position (x(i),z(k)) to element w and
% corresponding time of flight. 

d = single(zeros(length(z),length(x),length(Aperture)));        % Create array to store distances
tof = d;                                                        % Create array to store time of flight
costheta = d;                                                   % Create array to store angles

for w = 1:nChannels         % Index over channels
    for i = 1:length(x)     % Index over lateral pixel index
        for k = 1:length(z) % Index over range pixel index
            d(k,i,w) = sqrt((z(k)-zarray)^2+(x(i)-Aperture(w))^2); % distance between each element and each pixel in mm
            tof(k,i,w) = d(k,i,w)./c;         % time-of-flight between each element and each pixel in us
            costheta(k,i,w) = z(k)./d(k,i,w);   % compute angle between each element and each image point, use for apodization matrix
        end
    end
end

%% Shift Data with Appropriate Time Delays and Create Power Spectrum
delayed_channel = zeros(size(RData,1),nChannels);   % Create array to store delayed data
pow_spec = zeros(length(x),length(z));              % Create array to store power spectrum data

for i = 1:length(x)
    for k = 1:length(z)
        for w = 1:nChannels
            % Delay data by interpolating
            delayed_channel(:,w) = Sl*interp1(Time',RData(:,w),Time'+tof(k,i,w),'linear',0);   
            
            % Apply cosine apodization
            cosmat = squeeze(repmat(costheta(k,i,:),[size(RData,1),1])); % Form cosine apodization matrix if applying cosine apodization
            delayed_channel = delayed_channel.*cosmat;                   % Apply cosine apodization
            
        end
        pow_spec(i,k) = sum(abs(sum(delayed_channel,2)).^2);
        fprintf('Done processing: [%i,%i]\n', i,k);
    end
end

%% Show PAM Image
subplot(122);
imagesc(x,z,pow_spec');
title('Simulated Data - Cosine Apodization');
xlabel('Lateral Location [mm]'); 
ylabel('Range Location [mm]'); % Add axis labels
axis('image');