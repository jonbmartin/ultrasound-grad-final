function PCIBeamforming_Experiment
%% (C) 2016 Kevin J. Haworth, copyright protected 
% 
% The following code can be used with the MATLAB data file
% ExpData.mat to produce a passive cavitation image using emissions
% recorded on a Vantage 256 ultrasound research scanner (Verasonics Inc.,
% Bothell, WA, USA) using an L7-4 linear array (Philips, Bothell, WA, USA).
% The cavitation emissions were recorded from Definity (Lantheus Medical
% Imaging, Inc., N. Billerica, MA, USA) microbubbles pumped through tubing
% and insonified with 6 MHz ultrasound at a pressure amplitude sufficient
% to produce subharmonic emissions (indicative of stable cavitation). The
% pulse duration was approximately 55 cycles. The system calibration
% factor for each element has not been included.

% The script can be run from the command line or the individual sections
% can be run sequentially.  Each section performs a different step in the
% beamforming algorithm.  The only setup lines the user may wish to modify
% are: 
% 1) Those that specify the x and z coordinates of the pixels, the spacing
% and spatial extent of the coordinates can be modified to increase or
% decrease the processing time.
% 2) Which frames should be analyzed. Up to 100 frames are available.
% 3) The flag which dictates whether or not to save individual frames
% of data or just the average of all the frames. Saving individual frames
% significantly increases the size of the saved .mat files.
% 4) The flag which dictates whether or not cosine apodization will be
% applied to reduce the grating lobes that appear when beamforming the
% fundamental frequency (6 MHz).

% Suggested citations: 
% Haworth KJ, Mast TD, Radhakrishnan K, Burgess MT, Kopechek JA, Huang S-L,
% McPherson DD, Holland CK. Passive imaging with pulsed ultrasound
% insonations. J Acoust Soc Am Acoustical Society of America,
% 2012;132:544?553.    
% Haworth KJ, Bader KB, Rich KT, Holland CK, Mast TD. Quantitative
% Frequency-Domain Passive Cavitation Imaging. IEEE UFFC. DOI:
% 10.1109/TUFFC.2016.2620492

% Note that "us" indicates units of microseconds

clear       % Clear all variables, functions, etc. from memory 

%% Load the experimental data

BASEFolder = './';   % set the path to the data file ExpData.mat
DataFile = 'ExpData.mat';   

[YY,MM,DD,hh,mm,ss] = datevec(now); % This information will be used to create a timestamp that will be affixed on the saved beamformed data

load([BASEFolder DataFile]);        % Load Data
RData = double(RData);              % convert data to double (necessary for later processing), originally save as int16 for reduced file size

%% Transducer, Material, and Data Characteristics

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
fs = 20;                % Sampling frequency in MHz
dt = 1/fs;              % Sampling period in us
ncycles = 55;           % Number of cycles in insonation pulse data
Fo = 6;                 % Center frequency of insonation ultrasound in MHz
NIOI = round(fs/Fo*ncycles); % Number of points in the interval of interest
TIOI = NIOI*dt;         % Time duration of interval of interest in us
T = (size(RData,1)-1)*dt;   % Duration of recorded signal in us
Time = (0:dt:T)';       % Vector of time  in us
fk = linspace(-fs/2,fs/2,size(RData,1)); % Discrete frequency vector equal to k*Delta-f

% Passive Cavitation Image Properties; These values may be adjusted to
% modify the image being formed
NFrameStart = 1;        % Frame to start analysis
NFramesAnalyze = 20;    % Number of frames to analyze; up to 100 frames of data are available
x = Aperture(34):ElementPitch/1:Aperture(95);	% Lateral location of pixels in mm
z = 10:1:40;            % Range location of pixel in mm
SaveIndividFrames = 0;  % Flag to indicate whether to save individual frames, 0 == no, 1 == yes
CosineApodization = 0;  % Flag to indicate whether to apply cosine apodization, 0 == no, 1 == yes

%% Pre-process passively received data

% Mean Subtract to remove DC bias on each channel
RData_DC = mean(RData,1); 
RData = RData-repmat(RData_DC,[size(RData,1) 1 1]);
NFrames = size(RData,3);    % Number of frames

% Compute FFT of data
RFFT = fftshift(fft(RData,[],1),1);

%% Compute time of flight information
% Distance from point source with position (x(i),z(k)) to element w and
% corresponding time of flight. 
% Additionally compute apodization matrix based on transmit angle
d = single(zeros(length(z),length(x),length(Aperture)));
tof = d;    costheta = d;
for w = 1:nChannels         % Index over channels
    for i = 1:length(x)     % Index over lateral pixel index
        for k = 1:length(z) % Index over range pixel index
            d(k,i,w) = sqrt((z(k)-zarray)^2+( x(i)-Aperture(w))^2); % distance between each element and each pixel in mm
            tof(k,i,w) = d(k,i,w)./c;         % time-of-flight between each element and each pixel in us            
            if CosineApodization == 1
                costheta(k,i,w) = z(k)./d(k,i,w);   % compute angle between each element and each image point, use for apodization matrix
            end
        end
    end
end

%% Perform PCI beamforming, pixel amplitude is proportional to the Energy Spectral Density
EDSnobias = double(zeros([length(x) length(z) length(fk)]));  % Allocate memory for PCI that averages across all frames

tic     % Use this to track computational time.

display(sprintf('Starting beamforming, Frame 0001 of %04g.', NFramesAnalyze))
for iFrm = NFrameStart:NFrameStart+(NFramesAnalyze-1)   % Process each frame of interest
    starttoc=toc;
    if SaveIndividFrames == 1   
        PCISpec_IndividFrame = single(zeros([length(x) length(z) length(fk)]));  % Allocate memory for individual frames of data if they will be saved
    end
    
    SingleFrameOfData = squeeze(RFFT(:,:,iFrm));    % Define variable correspond to the current frame of data to be processed
    DCOffset = sum(abs(Sl*SingleFrameOfData).^2,2); % DC bias noted by Norton et al. (2006); J Acoust Soc Am 119(5):2840-2847; for plotting on a dB scale, the user may wish to set this to zero (i.e. "zeros(size(RData,1),1);") to eliminate negative values (See also line 177).
    for iX = 1:length(x)        % index through lateral pixel positions
        for iZ = 1:length(z)    % index through range pixel positions
            
            % Compute each term of the first summand in equation 23 by
            % applying phase shifts associated with beamforming.
            shifted = Sl*SingleFrameOfData.*(exp(2i*pi*(squeeze(tof(iZ,iX,:)))*fk)).';  % This line takes advantage of vectorization to eliminate a for loop over frequencies (and thus increase the computational speed). 
            
            if CosineApodization == 1                            
                cosmat = squeeze(repmat(costheta(iZ,iX,:),length(fk),1)); % Form cosine apodization matrix
                shifted = shifted.*cosmat;      % Apply cosine apodization
            end
            
            PseudoESD = abs(sum(shifted,2)).^2; % Compute the pseudo energy spectral density (i.e., the first term in square brackets in equation 23). It is not the true energy spectral density because we have not included the transducer system calibration factor for each element (M[k]). 
                        
            if SaveIndividFrames == 1
                PCISpec_IndividFrame(iX,iZ,:) = PseudoESD - DCOffset;    % Subtract DC bias and write beamformed data to individual PCI frame
            end
            EDSnobias(iX,iZ,:) = squeeze(EDSnobias(iX,iZ,:)) + PseudoESD - DCOffset;  % Subtract DC bias and compute what will be the cummulative energy density spectra for each pixel            
        end
        
        % If statement used to display the processing time and current data
        % being processed
        if iX == 1
            display(sprintf('Lateral Line %04g of %04g in %05.3g seconds',iX,length(x),toc-starttoc))
        else
            display(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bLateral Line %04g of %04g in %05.3g seconds',iX,length(x),toc-starttoc))
        end
        
    end
    display(sprintf('Frame %04g of %04g in %g seconds',iFrm,NFramesAnalyze,toc-starttoc))   % Display processing time for the most recently processed frame.
end

EDSnobias = (1/(rho*c*fs^2*S)).*EDSnobias;      % Apply normalization factors from equation (23). Transducer calibration factor M[k] is not included as it depends on the specific array used. Also, no windowing is performed and thus no correction is necessary for windowing. 
AveESDAllFrames = (EDSnobias./length(NFramesAnalyze));    % Divide by the number of frames analyzed to compute an average

display('Done Processing')

% Uncomment the following 3 lines to save the result
% clear RFFT RData
% save([BASEFolder DataFile(1:end-4) '_PCI_' num2str(YY) num2str(MM) num2str(DD) num2str(hh) num2str(mm) num2str(round(ss)) '.mat'])  % Save data
% display('Done Saving')

%%  Plots

% Plot the spectra for the pixel location with the maximum total energy (summing over all frequencies)
temp = sum(EDSnobias,3);  % compute matrix proportional to total energy at each pixel
[~,Imax] = max(temp(:));        % find location of maximum pixel amplitude in units of indices
[xmaxI,zmaxI] = ind2sub(size(temp),Imax);   % identify indices x and z location
SpectrumAtMaxLocation = squeeze(EDSnobias(xmaxI,zmaxI,:));
SpectrumAtMaxLocation(SpectrumAtMaxLocation <0) = eps;         % When removing the DC bias, negative values may occur, set them to a small positive value to enable plotting on using a log scale colormap
figure(20); plot(fk, 10*log10(SpectrumAtMaxLocation), 'LineWidth',2);
xlim([1,10]); xlabel('Frequency (MHz)')
ylim([10*log10(max(SpectrumAtMaxLocation))-70 10*log10(max(SpectrumAtMaxLocation))]); ylabel('Pseudo Energy Spectral Density (dB)')
set(gca,'fontsize',18)

% Create a passive cavitation image for the fundamental near 6 MHz 
[~,f6MHzI] = min(abs(fk-6)); f6MHz = fk(f6MHzI);  % Find the index to the frequency closest to 6 MHz
PCI_6MHz = sum(EDSnobias(:,:,f6MHzI-1:f6MHzI+1),3);  % Sum over a frequency band about 6 MHz
PCI_6MHznorm = PCI_6MHz./max(PCI_6MHz(:));  % Normalize the image for plotting on a log scale
PCI_6MHznorm(PCI_6MHznorm<0) = eps;         % When removing the DC bias, negative values may occur, set them to a small positive value to enable plotting on using a log scale colormap
figure(21); imagesc(x,z,10*log10(squeeze(PCI_6MHznorm))'); colorbar; colormap(hot); axis equal; axis tight;
caxis([-20 0]);
xlabel('Lateral Location (mm)'); ylabel('Range Location (mm)'); 
set(gca,'fontsize',18)

% Create a passive cavitation image for the subharmonic near 3 MHz 
[~,f3MHzI] = min(abs(fk-3)); f3MHz = fk(f3MHzI);  % Find the index to the frequency closest to 5 MHz
PCI_3MHz = sum(EDSnobias(:,:,f3MHzI-1:f3MHzI+1),3);  % Sum over a frequency band about 5 MHz
PCI_3MHznorm = PCI_3MHz./max(PCI_3MHz(:));  % Normalize the image for plotting on a log scale
PCI_3MHznorm(PCI_3MHznorm<0) = eps;         % When removing the DC bias, negative values may occur, set them to a small positive value to enable plotting on using a log scale colormap
figure(22); imagesc(x,z,10*log10(squeeze(PCI_3MHznorm))'); colorbar; colormap(hot); axis equal; axis tight;
caxis([-20 0]); % plot over a 15 dB dynamic range; selected based on the approximate increase in signal at 3 MHz relative to nearby inharmonic components
xlabel('Lateral Location (mm)'); ylabel('Range Location (mm)'); 
set(gca,'fontsize',18)


