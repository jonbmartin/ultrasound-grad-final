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
fk = linspace(-fs/2,fs/2-1/(TotalTimeI/fs),TotalTimeI); % Discrete frequency vector, equal to k*Delta-f

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

%% Compute FFT (note that this is not normalized)

% Compute FFT of data
SingleFrameOfData = fftshift(fft(RData,[],1),1);     % Fourier domain data (analgous to X[k]. Use FFT shift so that Fourier-domain data goes from -fs/2 to fs/2

%% Shift Data with Appropriate Time Delays and Create Power Spectrum
delayed_channel = zeros(size(RData,1),nChannels);   % Create array to store delayed data
pow_spec = zeros(length(x),length(z));              % Create array to store power spectrum data
DCOffset = sum(abs(Sl*SingleFrameOfData).^2,2);     % DC bias noted by Norton et al. (2006); J Acoust Soc Am 119(5):2840-2847. This variable represents each term of the second summand within equation 23. For plotting on a dB scale, the user may wish to set this to zero (i.e. "zeros(size(RData,1),1);") to eliminate warning associated with plotting log10 of a negative number (see line 218).

for i = 1:length(x)
    for k = 1:length(z)
        % Compute each term of the first summand in equation 23 by applying
        % phase shifts associated with beamforming.
        shifted = Sl*SingleFrameOfData.*(exp(2i*pi*(squeeze(tof(k,i,:)))*fk)).';  % This line takes advantage of vectorization to eliminate a for loop over frequencies (and thus increase the computational speed). 

        % Apply cosine apodization
        cosmat = squeeze(repmat(costheta(k,i,:),[length(fk),1])); % Form cosine apodization matrix if applying cosine apodization
        shifted = shifted.*cosmat;                   % Apply cosine apodization

        pow_spec = sum(abs(sum(shifted,2)).^2);
        spec_nobias(k,i,:) = squeeze(pow_spec - DCOffset);  % remove DC bias

        fprintf('Done processing: [%i,%i]\n', i,k);
    end
end

spec_nobias = (1/(rho*c*fs^2*S)).*spec_nobias;      % Apply normalization factors from equation (23). Transducer calibration factor M[k] is not included as it depends on the specific array used. Also, no windowing is performed and thus no correction is necessary for windowing. 

%% Show PAM Image
% Create a passive cavitation image for the fundamental near 6 MHz
[~,f6MHzI] = min(abs(fk-6)); f6MHz = fk(f6MHzI);  % Find the index to the frequency closest to 6 MHz
PCI_6MHz = sum(spec_nobias(:,:,f6MHzI-1:f6MHzI+1),3);  % Sum over a frequency band about 6 MHz
PCI_6MHznorm = PCI_6MHz./max(PCI_6MHz(:));  % Normalize the image for plotting on a decibel scale
PCI_6MHznorm(PCI_6MHznorm<0) = eps;         % When removing the DC bias, negative values may occur, set them to a small positive value to enable plotting on using a log scale colormap
imagesc(x,z,10*log10(squeeze(PCI_6MHznorm))');
imagesc(x,z,(squeeze(PCI_6MHznorm)));
xlabel('Lateral Location (mm)'); 
ylabel('Range Location (mm)'); % Add axis labels

%% 
[~,Imax] = max(temp(:));        % find location of maximum pixel amplitude in units of indices
[xmaxI,zmaxI] = ind2sub(size(temp),Imax);   % identify indices x and z location
SpectrumAtMaxLocation = squeeze(spec_nobias(xmaxI,zmaxI,:));
plot(fk,(SpectrumAtMaxLocation), 'LineWidth',2); axis tight

%%
test = sum(spec_nobias,3); 
figure;
imagesc(x,z,temp');
title('Simulated Data - Cosine Apodization');
xlabel('Lateral Location [mm]'); 
ylabel('Range Location [mm]'); % Add axis labels
axis('image');
%%
test = sum(spec_nobias,3); 
figure;
imagesc(x,z,temp');
title('Simulated Data - Cosine Apodization');
xlabel('Lateral Location [mm]'); 
ylabel('Range Location [mm]'); % Add axis labels
axis('image');