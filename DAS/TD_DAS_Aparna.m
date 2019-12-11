% Michelle Sigona 
% 20191206 - Time-Domain Delay and Sum Beamforming (TD-DAS) Algorithm for
% Passive Acoustic Mapping (DAS-PAM) for simulated data.

close all;
clear;
clc;

fstruct = dir('aparna_MBs_3MPa_trial1.bin');
fstruct2 = dir('aparna_MBs_3MPa_trial1_params.mat');

fid = fopen([fstruct.name]);
chandat = fread(fid,Inf,'int16');
fclose(fid);
load([fstruct2.name]);

chandat = reshape(chandat,[params.numRcvSamples params.numacq params.numRcvChannels params.numframes]);
rf_data = zeros([params.numRcvSamples 128 params.numacq*params.numframes]);

for m=1:params.numframes
    m
    for k=1:params.numacq
        my_i = (m-1)*params.numacq+k;
        rf_data(:,:,my_i) = squeeze(chandat(:,k,:,m));
    end
end

%% Set Transducer Parameters

% Material Properties; Do not change these parameters
c =  1.5236;            % speed of sound in mm/us 
rho = 1e-9;             % density of water in kg/mm^3

% Transducer Properties; Do not change these parameters
nChannels = 128;        % Number of channels that simultaneously and passively acquired acoustic emissions
ElementPitch = 0.2980;  % Spacing between L7-4 array elements in mm
Aperture = (-18.9230:ElementPitch:18.9230)';     % Position of elements in array
zarray = 0;             % Range location of array (should always be zero with a linear array)
EleHeight = 6;          % Element height in mm
EleWidth = 0.25;        % Element width in mm
Sl = EleHeight*EleWidth; % Element surface area in mm^2
S = Sl*nChannels;       % Active area of array in mm^2

% Received Data Properties; Do not change these parameters
fs = params.fs/1e6;     % Sampling frequency in MHz
dt = 1/fs;              % Sampling period in us
ncycles = 55;           % Number of cycles in insonation pulse data
Fo = params.f0/1e6;     % Center frequency of insonation ultrasound in MHz
NIOI = round(fs/Fo*ncycles); % Number of points in the interval of interest
TIOI = NIOI*dt;         % Time duration of interval of interest in us
T = (size(rf_data,1)-1)*dt;   % Duration of recorded signal in us
Time = (0:dt:T)';       % Vector of time  in us
fk = linspace(-fs/2,fs/2,size(rf_data,1)); % Discrete frequency vector equal to k*Delta-f

% Passive Cavitation Image Properties; These values may be adjusted to
% modify the image being formed
NFrameStart = 1;        % Frame to start analysis
NFramesAnalyze = 10;    % Number of frames to analyze; up to 100 frames of data are available
x = Aperture(34):ElementPitch/1:Aperture(95);	% Lateral location of pixels in mm
z = 50:1:70;            % Range location of pixel in mm

%% Remove DC Bias 
% Convert from int16 to double
rf_data = double(rf_data); 

% Mean Subtract to remove DC bias on each channel
RData_DC = mean(rf_data,1); 
rf_data = rf_data-repmat(RData_DC,[size(rf_data,1) 1 1]);
NFrames = size(rf_data,3);    % Number of frames

%% Calculate Time Delays
% Distance from point source with position (x(i),z(k)) to element w and
% corresponding time of flight. 

d = single(zeros(length(z),length(x),length(Aperture)));
tof = d;    
costheta = d;

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
delayed_channel = zeros(size(rf_data,1),nChannels);
pow_spec = zeros(length(x),length(z));
full_spec = zeros(NFramesAnalyze,length(x),length(z));
full_specNO = zeros(NFramesAnalyze,length(x),length(z));

for a = 1:2
    for iFrm = NFrameStart:NFrameStart+(NFramesAnalyze-1)
        for i = 1:length(x)
            for k = 1:length(z)
                for w = 1:nChannels
                    delayed_channel(:,w) = Sl*interp1(Time',double(rf_data(:,w,iFrm)),Time'+tof(k,i,w),'linear',0);           
                end
                
                if a == 2
                    % Apply cosine apodization
                    cosmat = squeeze(repmat(costheta(k,i,:),[size(rf_data,1),1])); % Form cosine apodization matrix if applying cosine apodization
                    delayed_channel = delayed_channel.*cosmat;                   % Apply cosine apodization
                end
                
                pow_spec(i,k) = sum(abs(sum(delayed_channel,2)).^2);
                fprintf('Done processing frame #%i: [%i,%i]\n', iFrm, i, k);
            end
        end
        if a == 1
            full_specNO(iFrm,:,:) = pow_spec;
        else
            full_spec(iFrm,:,:) = pow_spec;
        end
    end
end

%% Show PAM Image
nocav = squeeze(mean(full_specNO,1));
cav = squeeze(mean(full_spec,1));

% figure;
subplot(121);
imagesc(x,z,nocav');
title('Experimental Data - No Apodization');
xlabel('Lateral Location [mm]'); 
ylabel('Range Location [mm]'); % Add axis labels
axis('image');
colorbar;

% figure;
subplot(122);
imagesc(x,z,cav');
title('Experimental Data - Cosine Apodization');
xlabel('Lateral Location [mm]'); 
ylabel('Range Location [mm]'); % Add axis labels
axis('image');
colorbar;

%% Look at each frame
for i = 9%1:size(full_spec,1)
    subplot(121);
    imagesc(x,z,squeeze(full_specNO(i,:,:))');
    title('Experimental Data - No Apodization');
    xlabel('Lateral Location [mm]'); 
    ylabel('Range Location [mm]'); % Add axis labels
    axis('image');
    colorbar;

    subplot(122);
    imagesc(x,z,squeeze(full_spec(i,:,:))');
    title('Experimental Data - Cosine Apodization');
    xlabel('Lateral Location [mm]'); 
    ylabel('Range Location [mm]'); % Add axis labels
    axis('image');
    colorbar;
    pause(0.3);
end