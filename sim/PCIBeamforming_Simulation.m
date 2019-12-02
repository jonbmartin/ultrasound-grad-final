function PCIBeamforming_Simulation
%% (C) 2016 Kevin J. Haworth, copyright protected
%
% The following code can be used with the MATLAB data file
% SimData.mat to produce a passive cavitation image using simulated
% emissions from a point source located at a depth of 26 mm from the face
% of a linear array and at a lateral position centered about the linear
% array (i.e. lateral position of 0 mm). The source is a 6 MHz, 20-cycle
% tone burst with rectangular windowing. Alternately, the user can set the
% variable DataSource equal to 2 and simulated data will be created based
% on a source specified by the user (details below). The system calibration
% factor for each element has not been included. The array parameters
% specified below are based on a Philips L7-4 linear array.
% The script must be run in its entirety, either from the command line or
% the MATLAB GUI due to the inclusion of local functions.  Each section
% performs a different step in the beamforming algorithm.  The most likely
% lines the user may wish to modify are in the first section.  These lines
% set:  
% 1) The 'DataSource' flag which indicates whether the user wants to use
% the provided simulated data or create new simulated data. 
% 2) If the user choses to create new simulated data, the source
% location (variables 'Cav_x' and 'Cav_z'), the pulse duration ('ncycles'),
% and the scattered frequency ('Fo') can be modified. 
% 3) The x and z coordinates of the pixels in the PCI image. The spacing
% and spatial extent of the coordinates can be modified to increase or
% decrease the processing time.
% 4) The flag which dictates whether or not cosine apodization will be
% applied to reduce the grating lobes that appear when beamforming the
% fundamental frequency (6 MHz).


% Citations: 
% Haworth KJ, Mast TD, Radhakrishnan K, Burgess MT, Kopechek JA, Huang S-L,
% McPherson DD, Holland CK. Passive imaging with pulsed ultrasound
% insonations. J Acoust Soc Am Acoustical Society of America,
% 2012;132:544-553.    
% Haworth KJ, Bader KB, Rich KT, Holland CK, Mast TD. Quantitative
% Frequency-Domain Passive Cavitation Imaging. IEEE UFFC. DOI:
% 10.1109/TUFFC.2016.2620492  

% Note that "us" indicates units of microseconds.

clear all   % Clear all variables, functions, etc. from memory

% Specific pixel locations
x = -9:0.3:9;           % Lateral location of pixels in mm
z = 10:1:40;            % Range location of pixel in mm
% Flag to indicate if the user wants to apply cosine apodization
CosineApodization = 0;  % 0 == no apodization, 1 == apply cosine apodization
% Specify source of data and source characteristics 
DataSource = 1;         % Set DataSource equal to 1 to load the .mat file containing simulated data and set DataSource equal to 2 to create simulated received data
NormalizedDataSource = 0;   % This variable is used when DataSource is set to 2. Set NormalizedDataSource to 0 to use simulated received data as created based on element directivity. Set NormalizedDataSource to 1 to normalize the received elements to have equal amplitude, simulating omnidirectional receivers.
Cav_x = 0; Cav_z = 26;  % Origin location of emission in mm
ncycles = 20;           % Number of cycles in simulated data
Fo = 6;                 % Center frequency of simulated data in MHz


%% Define Transducer, Material, and Data Characteristics

% Material Properties; Do not change these parameters if loading data
c = 1.5236;             % speed of sound in mm/us
rho = 1e-9;             % density of water in kg/mm^3

% Transducer Properties; Do not change these parameters if loading data 
nChannels = 128;        % Number of channels that simultaneously and passively acquired acoustic emissions
ElementPitch = 0.2980;  % Spacing between L7-4 array elements in mm
Aperture = (-18.9230:ElementPitch:18.9230)';     % Position of elements in array
zarray = 0;             % Range location of array (should always be zero with a linear array)
EleHeight = 6;          % Element height in mm
EleWidth = 0.25;        % Element width in mm
Sl = EleHeight*EleWidth; % Element surface area in mm^2
S = Sl*nChannels;       % Active area of array in mm^2

% Received Data Properties; Do not change these parameters if loading data 
T = 30;                 % Nominal duration of the recorded emission in us
fs = 40;                % Sampling frequency in MHz
dt = 1/fs;              % Sampling period in us
NIOI = round(fs/Fo*ncycles); % Number of points in the interval of interest
TIOI = NIOI*dt;         % Time duration of interval of interest in us


%% Load or create the simulated data

BASEFolder = './';   % set the path to the data file SimData.mat
[YY,MM,DD,hh,mm,ss] = datevec(now); % This information will be used to create a timestamp that will be affixed on the saved beamformed data

if DataSource ==1   % Load received data from pre-existing .mat file
    DataFile = 'SimData.mat';
    load([BASEFolder DataFile]);        % Load Data
    
    TotalTimeI = size(RData,1);     % The total time duration of the simulated received signal (in indicies)
    Time = (0:TotalTimeI-1)*dt;     % Time vector
    fk = linspace(-fs/2,fs/2-1/(TotalTimeI/fs),TotalTimeI); % Discrete frequency vector, equal to k*Delta-f
    
elseif DataSource == 2  % Create received simulated data
    [~,Cav_xI]=min(abs(x-Cav_x));  Cav_x1 = x(Cav_xI);   % Make sure the emission originates at a location where a pixel will occur (this command produces a cleaner images, but is not a requirement and may be commented out)
    [~,Cav_zI]=min(abs(z-Cav_z));  Cav_z1 = z(Cav_zI);   % Make sure the emission originates at a location where a pixel will occur (this command produces a cleaner images, but is not a requirement and may be commented out)

    DataFile = sprintf('SimData_x%03.0f_z%03.0f.mat',round(Cav_x),round(Cav_z));    % Define file name for saving
    
    TDelay1 = sqrt((Aperture - Cav_x1).^2 + (zarray - Cav_z1).^2)/c; % Compute approximate delays to assist in computing necessary zero padding to prevent wrap-around (note that this is not strictly needed)
    zpad = ceil(max(TDelay1)./dt);  % Zero pad data so that time-delay does not cause data to wrap around
    TotalTimeI = 2^nextpow2(round(T/dt)+2*zpad);  % Determine the total time (in indicies), then use next highest power of 2 to increase FFT computation speed. The speed increase may be nominal and not worth the additional memory required based on the specifications of the computer used.
    Time = (0:TotalTimeI-1)*dt;     % Time vector
    fk = linspace(-fs/2,fs/2-1/(TotalTimeI/fs),TotalTimeI); % Discrete frequency vector, equal to k*Delta-f
    
    display('Computing propagation to create simulated data')
    
    % Compute matrix of element sensitivities (i.e. complex diffraction
    % pattern amplitude and phase) that are used to simulate received
    % waveforms 
    p1=zeros(length(Aperture), length(fk)); 
    tic
    for jj = 1:nChannels
        if jj == 1
            display(sprintf('Element %03d of %03d after %06.2f sec',jj,nChannels,toc))
        else 
            display(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%03d of %03d after %06.2f sec',jj,nChannels,toc))
        end
        for ll = 1:length(fk)
            % Compute the pressure using the local function pelement; the
            % input parameters are defined based on a Philips L7-4 array.
            % The function uses SI units (m, s, Hz, etc.).  Note that the
            % remainder of the code uses mm, us, and MHz.
            p1(jj, ll)=pelement(2*pi*fk(ll)*1e6/(c*1000),EleWidth*1e-3*0.95/2,EleHeight*1e-3/2,1e6,25e-3,(Cav_x1-Aperture(jj))*1e-3,0,(Cav_z1-zarray)*1e-3);
        end
        
    end
    
    % Next we create the waveform at the source
    MaxCavDurationI = round(ncycles*(fs/Fo));       % Duration of cavitation emission in units of indicies
    CavSig = sin(2*pi*Fo*(0:MaxCavDurationI-1)*dt);   % Emission from source with no windowing applied
%     CavSig = CavSig.*hann(length(CavSig))'; % Simulate ring up and ring down using a hann window.  This line is optional and may be uncommented to crudely simulate ring up and ring down
    CavSig = [zeros(1, floor((TotalTimeI-MaxCavDurationI)/2)) CavSig zeros(1, ceil((TotalTimeI-MaxCavDurationI)/2))];
    
    RectData = repmat(CavSig, [length(Aperture) 1])';   % Create a waveform for each receive channel
    fftDataRect = fftshift(fft(RectData,[],1));         % Compute Fourier tranform
    ReceivedfftDataRect = (fftDataRect.*p1.');          % Apply propagation from emission source to transducer elements
    RData = real(ifft(ifftshift(ReceivedfftDataRect,1),[],1));          % Due to rounding errors, the IFFT creates a complex-valued RData.  The imaginary compnents should be small and can be ignored by taking the real part of the waveforms. An alternate option is to use the IFFT 'symmetric' flag.
    if NormalizedDataSource == 1; % To test the energy equivalence between Parseval's theorem and the PCI algorithm for an approximation of omnidirectional receivers, set this flag to 1.
        RData = RData./sqrt(repmat(sum(abs(RData).^2,1),[size(RData,1) 1]));   
    end
else 
    error('The variable DataSource must be set to either 1 or 2.')
end

%% Compute FFT (note that this is not normalized)

% Compute FFT of data
SingleFrameOfData = fftshift(fft(RData,[],1),1);     % Fourier domain data (analgous to X[k]. Use FFT shift so that Fourier-domain data goes from -fs/2 to fs/2

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
%
starttoc=toc;
display('Performing passive cavitation imaging beamforming')
DCOffset = sum(abs(Sl*SingleFrameOfData).^2,2);    % DC bias noted by Norton et al. (2006); J Acoust Soc Am 119(5):2840-2847. This variable represents each term of the second summand within equation 23. For plotting on a dB scale, the user may wish to set this to zero (i.e. "zeros(size(RData,1),1);") to eliminate warning associated with plotting log10 of a negative number (see line 218).
for iX = 1:length(x)        % index through lateral pixel position in PCI
    for iZ = 1:length(z)    % index through range pixel position in PCI
        
        % Compute each term of the first summand in equation 23 by applying
        % phase shifts associated with beamforming.
        shifted = Sl*SingleFrameOfData.*(exp(2i*pi*(squeeze(tof(iZ,iX,:)))*fk)).';  % This line takes advantage of vectorization to eliminate a for loop over frequencies (and thus increase the computational speed). 

        if CosineApodization == 1
            cosmat = squeeze(repmat(costheta(iZ,iX,:),length(fk),1)); % Form cosine apodization matrix if applying cosine apodization
            shifted = shifted.*cosmat;      % Apply cosine apodization
        end
        
        PseudoESD = abs(sum(shifted,2)).^2; % Compute the pseudo energy spectral density (i.e., the first term in square brackets in equation 23). It is not the true energy spectral density because we have not included the transducer system calibration factor for each element (M[k]). 
        
        EDSnobias(iX,iZ,:) = squeeze(PseudoESD - DCOffset);  % remove DC bias
    end
    
    % If statement used to display the processing time and current data
    % being processed in the MATLAB command window.
    if iX == 1
        display(sprintf('Lateral Line %04g of %04g in %05.3g seconds',iX,length(x),toc-starttoc))
    else
        display(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bLateral Line %04g of %04g in %05.3g seconds',iX,length(x),toc-starttoc))
    end
end

EDSnobias = (1/(rho*c*fs^2*S)).*EDSnobias;      % Apply normalization factors from equation (23). Transducer calibration factor M[k] is not included as it depends on the specific array used. Also, no windowing is performed and thus no correction is necessary for windowing. 

display('Done Processing')

% Uncomment the following lines to save the data set using a time stamp.
% save([BASEFolder DataFile(1:end-4) '_PCI_' num2str(YY) num2str(MM) num2str(DD) num2str(hh) num2str(mm) num2str(round(ss)) '.mat'])  % Save data
% display('Done Saving')


%%  Plots

% Plot the spectra for the pixel location with the maximum total energy (summing over all frequencies)
temp = sum(EDSnobias,3);  % compute matrix proportional to total energy at each pixel
[~,Imax] = max(temp(:));        % find location of maximum pixel amplitude in units of indices (linear)
[xmaxI,zmaxI] = ind2sub(size(temp),Imax);   % identify indices associated with location of maxium pixel amplitude
SpectrumAtMaxLocation = squeeze(EDSnobias(xmaxI,zmaxI,:));
SpectrumAtMaxLocation(SpectrumAtMaxLocation <0) = eps;         % When removing the DC bias, negative values may occur, set them to a small positive value to enable plotting on using a decibel scale
figure(10); plot(fk, 10*log10(SpectrumAtMaxLocation), 'LineWidth',2);   % Plot spectrum at location with maximum total energy
xlim([1,10]); xlabel('Frequency (MHz)');        % Adjust abscissa limits and add label
ylabel('Pseudo Energy Density Spectrum (dB)');  % Add ordinate label
set(gca,'fontsize',18)                          % Adjust font size

% Create a passive cavitation image for the fundamental near 6 MHz
[~,f6MHzI] = min(abs(fk-6)); f6MHz = fk(f6MHzI);  % Find the index to the frequency closest to 6 MHz
PCI_6MHz = sum(EDSnobias(:,:,f6MHzI-1:f6MHzI+1),3);  % Sum over a frequency band about 6 MHz
PCI_6MHznorm = PCI_6MHz./max(PCI_6MHz(:));  % Normalize the image for plotting on a decibel scale
PCI_6MHznorm(PCI_6MHznorm<0) = eps;         % When removing the DC bias, negative values may occur, set them to a small positive value to enable plotting on using a log scale colormap
figure(11); imagesc(x,z,10*log10(squeeze(PCI_6MHznorm))'); colorbar; colormap(hot); axis equal; axis tight; % Plot passive cavitation image
caxis([-15 0]);         % Set colorbar dynamic range to 15 dB
xlabel('Lateral Location (mm)'); ylabel('Range Location (mm)'); % Add axis labels
set(gca,'fontsize',18); % Set font size

%% Uncomment the following code to compare the total incident energy as
% calculated using Parseval's Theorem to the estimate derived from the
% passive cavitation imaging algorithm.  


% Compute shifted and PsuedoESD at source location
shifted = Sl*SingleFrameOfData.*(exp(2i*pi*(squeeze(tof(Cav_zI,Cav_xI,:)))*fk)).';  % apply phase delays associated with beamforming; note that the vector fk is equal to k*Delta-f; This line takes advantage of vectorization to improve processing speed. We also multiply by the area of each element relative to the total area of the array.  In this case all the elements are the same size and there are 128 of them.
PseudoESD = abs(sum(shifted,2)).^2;

% Energy calculated using Parseval's theorem 
E_Pars_TD = (Sl/rho/c*sum(sum(abs(RData.^2),2),1) * 1/fs); % use time-domain data
E_Pars_FD = ((Sl/rho/c*(fk(2) - fk(1))/fs^2)*(sum(sum(abs(SingleFrameOfData).^2,2),1))); % use frequency-domain data
% E_Pars_FD/E_Pars_TD     % This ratio should be equal to 1, neglecting round-off error

% Energy estimated using the passive cavitation imaging algorith
E_PCI = (1/(rho*c*fs^2*S))*sum(PseudoESD,1)*(fk(2) - fk(1));    % Estimate without DC offset
E_PCI2 = ((fk(2) - fk(1)))*sum(EDSnobias(xmaxI,zmaxI,:),3);   % Estimate with DC offset

EnergyError = (E_Pars_TD - E_PCI2)/E_Pars_TD*100;     % percent error; If the m-file is used without modification, this error should be 6.4%. If the only modification to the mfile is setting Fo to 2 MHz, then the error should be 0.8%
display(sprintf('PCI Energy Estimation Error (%%): %g', EnergyError))

end

%% Local functions (i.e. subfunctions) needed to simulate received data
function p=pelement(k,a,b,Fx,Fy,x,y,z)
% pelement(k,a,b,Fx,Fy,xvec,yvec,zvec) returns [p(x,y,z)] under the 
% Fresnel approximation for a unit-amplitude rectangular piston with 
% wavenumber k, half-widths a and b, and geometric focal lengths Fx, Fy. 
% Based on expressions from T. D. Mast, "Fresnel approximations for 
% acoustic fields of rectangularly symmetric sources," J. Acoust. Soc.
% Am. 121(6):3311-3322, 2007, using the Fresnel approximation type
% type 'r'. All parameters must be either scalar variables or 
% vectors/matrices of common size.
% (C) 2007, T. Douglas Mast

eps = 1.e-16;
s2i=sqrt(2i);

% add eps to avoid numerical division by zero
x = x + eps;
y = y + eps;
z = z + eps;
r = sqrt(x.^2+y.^2+z.^2);

ktildex = (k .* (1 - r./Fx)) + eps; 
ktildey = (k .* (1 - r./Fy)) + eps;
   
p = 0.25 .* k .* ...
     exp( 1i*(k.*r - k.^2.*(x.^2./ktildex + y.^2./ktildey)./(2*r)) ) ...
       ./ (sqrt(ktildex).*sqrt(ktildey)) ...
        .* ( cerror((k.*x + ktildex.*a)./(s2i*sqrt(ktildex.*r))) ...
          - cerror((k.*x - ktildex.*a)./(s2i*sqrt(ktildex.*r))) ) ...
            .* ( cerror((k.*y + ktildey.*b)./(s2i*sqrt(ktildey.*r))) ...
              - cerror((k.*y - ktildey.*b)./(s2i*sqrt(ktildey.*r))) );
p = conj(squeeze(p));   % conjugate to obtain exp(+i omega t) time dependence

    function y = cerror(x)
        % compute complex error function using rational approximation,
        % order 16. Based on Table 1 in J.A.C. Weideman,
        % "Computation of the Complex Error Function" SIAM J. Num. Anal.
        % 31(5), 1497-1518 (1994)
        
        xsign = sign(real(x));
        x = x.*xsign;
        L = 3.3635856610148580;
        invsqrtpi = 0.564189583547756287;
        a1 = [ 9.93932253556817358e-07  3.98128757509398406e-06 ...
            -5.58423341110826437e-06 -2.73464046244891389e-05 ...
            2.17098679322337925e-05  0.000210710563965327855 ...
            8.70315842845801376e-05 -0.00152765974012200774  ...
            -0.00388101518902271136   0.00368256731709189645  ...
            0.0518224024316115278    0.191241726746694735    ...
            0.469290900903603592     0.886447830205054577    ...
            1.36224082227195864      1.7483958860819615 ];
        Z = (L-x)./(L+x); p1 = polyval(a1,Z);      % Polynomial evaluation.
        w = 2*p1./(L+x).^2 + invsqrtpi./(L+x);    % Evaluate w(z).
        y = (1 - exp(-x.^2).*w).*xsign;
    end

end

