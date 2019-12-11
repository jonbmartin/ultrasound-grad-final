%% It turns out my portion of project requires das
% implement das
% will write my own to make sure I get this 

load('test_data/SimData.mat')

%% Define grid to image along with a few other properties
% Specific pixel locations
xloc = -9:0.3:9;           % Lateral location of pixels in mm
zloc = 10:1:40;            % depth location of pixel in mm

ncycles = 20;           % Number of cycles in simulated data
Fo = 6;                 % Center frequency of simulated data in MHz

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

TotalTimeI = size(RData,1);     % The total time duration of the simulated received signal (in indicies)
Time = (0:TotalTimeI-1)*dt;     % Time vector
fk = linspace(-fs/2,fs/2-1/(TotalTimeI/fs),TotalTimeI); % Discrete frequency vector, equal to k*Delta-f

%% Optimization woop woop!
%Part 1 : Robust Capon estimation of optimal steering vector

abar = ones(size(Aperture)); %nominal steering vector
a0 = ones(size(Aperture)); %initial estimate of steering vector

%construct matrices D and R (both [numel numel]

d = single(zeros(length(zloc),length(xloc),length(Aperture)));
for w = 1:nChannels         % Index over channels
    for i = 1:length(xloc)     % Index over lateral pixel index
        for k = 1:length(zloc) % Index over range pixel index
            d(k,i,w) = sqrt((zloc(k)-zarray)^2+( xloc(i)-Aperture(w))^2); % distance between each element and each pixel in mm
            tof(k,i,w) = d(k,i,w)./c;         % time-of-flight between each element and each pixel in us
        end
    end
end
d = permute(d, [3, 1, 2]);

%Hard coded, not desirable, but for each voxel...
for ii = 1:length(zloc)
    for jj = 1:length(xloc)
        D = diag(d(:,ii,jj));
    end
end

%pre-steer the data
for i=1:length(xloc) % lateral position
%for i=28:33
    for j=1:length(zloc) %axial position
%    for j=15:19
        for e=1:length(Aperture) %element
            shifted(:,e) = Sl*interp1(Time',RData(:,e),Time'+tof(j,i,e));
        end
    end
end
shifted(isnan(shifted))=0;
shifted = shifted'; % just transpose to be consistent with paper notation
R = shifted * shifted';

%% Now perform RCB
eps = 20;

for ii = 1:length(xloc)
    for kk = 1:length(zloc)
        times = round(squeeze(tof(kk,ii,:)) / (1/(fs))) + 1; % convert time to corresponding pixel indices
        time2ind = times - min(times) + 1; %perform shift
        dist_v = squeeze(d(ii,kk,:));
        shifted_data = zeros(128, size(RData,1) + 1);
        for jj = 1:nChannels
            shifted_data(jj,:) = shifted_data(jj,:) + [RData(1+time2ind(jj):end,jj)' zeros(1, time2ind(jj)+1)]; % shift the pixels up by the number calculated in time2ind and zero pad rest
        end
        
        R = double(shifted_data * shifted_data'); % calculate covariance Rs
        [U, L] = eig(R);
        z_sumsq = transpose(U)*abar;
        
        lam_max = (norm(abar)-sqrt(eps))/(L(128,128)*sqrt(eps));
        lam_min = (norm(abar)-sqrt(eps))/(L(1,1)*sqrt(eps));
        lam0 = [lam_min lam_max]; %initial lambda interval
        fun = @(lam)sum(z_sumsq.^2./((1+(lam*L)).^2))-eps; % eqn (24)
        lam = lsqnonlin(fun,(lam_min+lam_max)/2,lam_min,lam_max);
        
        ahat = abar - pinv(I + lam*R)*abar;
        ahat = (ahat*128)/norm(ahat); %normalize the ahat

        dist_n = sum(dist_v)/128; % find the mean distance of the pixel to transducer element
        
        pwr_rcb_img(ii,kk) = (4*pi*(dist_n^2)/(rho*c))/(ahat'*U*L*pinv((lam^-2*eye(128))+(2*lam^-1*L)+(L^2))*transpose(U)*ahat);

        
    end
end

%% Now perform RCB linalg, find lagrange multiplier
R = double(R);
I = eye(nChannels);
C = D * R * D;
eps = 20;

%lagrange multiplier
[U, L] = eig(R);
z_sumsq = transpose(U)*abar;

lam_max = (norm(abar)-sqrt(eps))/(L(128,128)*sqrt(eps));
lam_min = (norm(abar)-sqrt(eps))/(L(1,1)*sqrt(eps));
lam0 = [lam_min lam_max]; %initial lambda interval
fun = @(lam)sum(z_sumsq.^2./((1+(lam*L)).^2))-eps; % eqn (24)
lam = lsqnonlin(fun,(lam_min+lam_max)/2,lam_min,lam_max);

%% solve steering vector, weight vector

ahat = abar - pinv(I + lam*R)*abar;
%normalize the ahat
ahat = (ahat*128)/norm(ahat);
w = (pinv(R)*ahat)/(conj(ahat)'*pinv(R)*ahat);

%% apply your weight vector to the data

for ii = 1:4096
   image(ii,:) = w.*shifted(:,ii); 
end