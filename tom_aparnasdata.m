%% Implement HOS on aparnas data
clear all
fstruct2 = dir('aparna_MBs_3MPa_trial1_params.mat');
fstruct = dir('aparna_MBs_3MPa_trial1.bin');
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
end;

%% 
frame = 9;
RData = double(squeeze(chandat(:,1,:,frame)));
RData = RData(params.t0:end,:);

imagesc(RData)
%%
for i = 1:50
    plot(rf_data(:,100,i))
    title(i)
    pause
end


%% Define grid to image along with a few other properties
% Transducer Properties; Do not change these parameters if loading data 
nChannels = 128;        % Number of channels that simultaneously and passively acquired acoustic emissions
ElementPitch = 0.2980;  % Spacing between L7-4 array elements in mm
Aperture = (-18.9230:ElementPitch:18.9230)';     % Position of elements in array
zarray = 0;             % Range location of array (should always be zero with a linear array)
EleHeight = 6;          % Element height in mm
EleWidth = 0.25;        % Element width in mm
Sl = EleHeight*EleWidth; % Element surface area in mm^2
S = Sl*nChannels;       % Active area of array in mm^2

% Specific pixel locations
xloc = Aperture(34):ElementPitch/1:Aperture(95);	% Lateral location of pixels in mm
zloc = 50:1:75;            % depth location of pixel in mm

%ncycles = 20;           % Number of cycles in simulated data
Fo = 5.208;                 % Center frequency of simulated data in MHz

% Material Properties; Do not change these parameters if loading data
c = 1.5236;             % speed of sound in mm/us
rho = 1e-9;             % density of water in kg/mm^3

% Received Data Properties; Do not change these parameters if loading data 
%T = 30;                 % Nominal duration of the recorded emission in us
fs = 20.832;                % Sampling frequency in MHz
dt = 1/fs;              % Sampling period in us
%NIOI = round(fs/Fo*ncycles); % Number of points in the interval of interest
%TIOI = NIOI*dt;         % Time duration of interval of interest in us

TotalTimeI = size(RData,1);     % The total time duration of the simulated received signal (in indicies)
Time = (0:TotalTimeI-1)*dt;     % Time vector
fk = linspace(-fs/2,fs/2-1/(TotalTimeI/fs),TotalTimeI); % Discrete frequency vector, equal to k*Delta-f

%% compute delay for each pixel
% we are looking for distance from point source to element w
% Distance from point source with position (x(i),z(k)) to element w and
% corresponding time of flight. 
% Additionally compute apodization matrix based on transmit angle
d = single(zeros(length(zloc),length(xloc),length(Aperture)));
tof = d;    costheta = d;
for w = 1:nChannels         % Index over channels
    for i = 1:length(xloc)     % Index over lateral pixel index
        for k = 1:length(zloc) % Index over range pixel index
            d(k,i,w) = sqrt((zloc(k)-zarray)^2+( xloc(i)-Aperture(w))^2); % distance between each element and each pixel in mm
            tof(k,i,w) = d(k,i,w)./c;         % time-of-flight between each element and each pixel in us
            costheta(k,i,w) = zloc(k)./d(k,i,w);   % compute angle between each element and each image point, use for apodization matrix
        end
    end
end

%% apply delays in time domain
% made need to correct for dc bias prior to summation
DCoff = sum(abs((Sl*RData).^2),2);

shifted =zeros(TotalTimeI,length(Aperture)); %temp var, holds delayed channel data
powIm = zeros(length(zloc),length(xloc));
h = waitbar(0,'applying delay');

for i=1:length(xloc) % lateral position
%for i=26:33
    for j=1:length(zloc) %axial position
    %for j=5:20
        for e=1:length(Aperture) %element
            shifted(:,e) = Sl*interp1(Time',RData(:,e),Time'+tof(j,i,e));
        end
        % account for angular sensitivity (cosine apodization)
        cosmat = repmat(squeeze(costheta(j,i,:))',[TotalTimeI,1]); % Form cosine apodization matrix if applying cosine apodization
%         shifted = shifted.*cosmat;      % Apply cosine apodization
% 
%         % find optimal weights using linprog
%         X = shifted';
%         X(isnan(X))=0;
%         X = double(X);
%         N = TotalTimeI;
%         M = 128;
%         sig = .5; %.04 up to .8 have pretty good results % no sln at 1
%         
%         u = optimvar('u');
%         w = optimvar('w',M);
%         r = optimvar('r');
%         tom = optimvar('tom');
%         prob = optimproblem('Objective',u,'ObjectiveSense','min');
%         prob.Constraints.c1 = X'*w >= -u.*ones([N 1]);
%         prob.Constraints.c2 = X'*w <= u.*ones([N 1]);
%         prob.Constraints.c3 = w >= 1 + sig*r*ones([M 1]);
%         prob.Constraints.c4 = -r*ones([M 1]) <= w;
%         prob.Constraints.c5 = w <= r*ones([M 1]);
% 
%         problem = prob2struct(prob);
% 
%         [sol,fval,exitflag,output] = linprog(problem);
%         
%         % apply weights
%         if exitflag==1 %it worked
%             weights = repmat((sol(2:2+128-1)./sol(1))',[N 1]);
%         else
%             weights = ones([N M]);
%         end
%         shifted = shifted.*weights;
        
        % sum
        latsum = abs(sum(shifted,2,'omitnan')).^2 - 0;
        powIm(j,i) = sum(latsum,'omitnan');
    end
    waitbar(i/length(xloc));
end
close(h)

