% m-file editwav15.m
% This script will:
% 1- Calculate the SPL in 1 s intervals for a SeaExplorer recording 
% 2- Remove two noise patterns created by the SeaExplorer mechanisms based
% on the sound of channel 1
% Two samples of the noise need to be stored in a folder called Noise_DIR
% The files to be analized and cleaned need to be stored in the Mission_DIR
% folder. 
% The software uses spectral matching and the Itakura-Saito distance to
% detect the desired noise patterns. A file segment will be cancelled if it
% matches one of the two patterns provided and lasts at least one second.
% The length of the file will remain unaltered. The undesired noise
% segments will be replaced by NaN values.
% A single clean file containing all the files analyzed will be created and 
% saved in the Mission_DIR folder.
% The parameters have been adjusted to remove all "suspicious" segments.
% Script created by José A. Diaz, Sept 2026.

clear;

% Noise pattern and file directories
Cur_DIR = pwd;
Noise_DIR = 'C:\...'; % There must be two noise patterns in this directory. Edit as needed
Mission_DIR = 'C:\...'; % Place the sound files here. Edit as needed

cd(Noise_DIR);
fs_native = 192000; % Assuming that the SeaExplorer sampling rate is 192 KHz
dec = 8; % Decimation factor
fs_processing = 192000/dec; % Downsampled rate (captures up to 12 kHz)
p = 16; % Prediction order

% Load 192kHz noise samples and downsample
% [file,location] = uigetfile('*.wav')
location = Noise_DIR;
file = 'Noise1.wav';
FileName = fullfile(location,file);
[x,Fs] = audioread(FileName,'double');
xdown = decimate(x,dec);

% Preprocess and obtain noise model A
win = hamming(2400);
xdown = xdown(10001:12400);
xdown = xdown - mean(xdown);
xdown = xdown.*win;
xdown = filter([1, -0.95], 1, xdown);
[a1,g1] = lpc(xdown, p);

% Load 192kHz noise samples and downsample
% [file,location] = uigetfile('*.wav')
location = Noise_DIR;
file = 'Noise2.wav';
FileName = fullfile(location,file);
[x,Fs] = audioread(FileName,'double');
xdown = decimate(x,dec);

% Preprocess and obtain noise model
xdown = xdown(20001:22400).*win;
xdown = xdown - mean(xdown);
xdown = xdown.*win;
xdown = filter([1, -0.95], 1, xdown);
[a2,g2] = lpc(xdown,p);

% Pre-calculate target spectral envelopes
nfft = 512;
[hA,f] = freqz(sqrt(g1),a1,nfft,fs_processing);
[hB,f] = freqz(sqrt(g2),a2,nfft,fs_processing);  

TimeWindowMSFD = 0.1; % 100 ms time window

% SeaExplorer parameters
Sh = -153.6;
FScal = 4*sqrt(2); 
AnalogValueFactorSq = FScal.^2./4;

cd (Mission_DIR);
FileList = dir('**/*.wav');
NumFiles = length(FileList)
yt = []; % File containing all data samples
YdBt = []; % SPL file
distToAt = []; % Distance to Noise 1 (Sound A)
distToBt = []; % Distance to Noise 2 (Sound B)
WinNum = 1; % Window number;

for k =1:NumFiles
        Mission_FILE = FileList(k).name
        Full_Mission_FILE = fullfile(FileList(k).folder,Mission_FILE);
        FILE_Info = audioinfo(Full_Mission_FILE);
        FILE_Info_bin = dir(Full_Mission_FILE);
        Fs = FILE_Info.SampleRate;
        SAMPLES = TimeWindowMSFD*Fs; % 1 sample includes all channels

        [x,Fs] = audioread(Full_Mission_FILE,'double');
        s = size(x);
        y = x(:,1); % y contains only the 1st channel
        
        START = 1;
        i = 1;
        YdB = [];
        distToA = [];
        distToB = [];

        while START <= s(1) - SAMPLES + 1
           yi = y(START:(START + SAMPLES -1),1);
           YdB(i) = 10.*log10(mean(yi.*yi)*AnalogValueFactorSq)-Sh;

           % Downsample to 24 kHz and preprocess
           yidown = decimate(yi,dec);
           yidown = yidown - mean(yidown);
           yidown = yidown.*win;
           yidown = filter([1, -0.95], 1, yidown);

           % Compute LPC model 
           [aT,gT] = lpc(yidown,p);

           % Step C: Generate the live spectral envelope
           [hT,f] = freqz(sqrt(gT),aT,nfft,fs_processing);
  
           PA = g1 ./ (abs(hA).^2 + eps); % Power spectrum of Sound A
           PB = g2 ./ (abs(hB).^2 + eps); % Power spectrum of Sound b
           PT = gT ./ (abs(hT).^2 + eps); % power spectrum of current file

           % Itakura-Saito distance to Sound A
           ratio = PT ./ PA;
           is_vector = ratio - log(ratio + eps) - 1;
           distToA(i) = mean(is_vector);

           % Itakura-Saito distance to Sound B
           ratio = PT ./ PB;
           is_vector = ratio - log(ratio + eps) - 1;
           distToB(i) = mean(is_vector);

           START = START + SAMPLES;
           i = i +1;
           WinNum = WinNum + 1;
        end
        yt =[yt; y]; % Total sample file
        YdBt = [YdBt YdB]; % Total SPL vector
        distToAt = [distToAt distToA]; % Distance to A full vector
        distToBt = [distToBt distToB]; % Distance to B full vector
end

% Smooth distance vectors
distToAt = smoothdata(distToAt,"movmean",6);
distToBt = smoothdata(distToBt,"movmean",6);

contA = 0;
contB = 0;
NumDet = 0;
WinLen = [];
WinPos = [];

for i = 1:(WinNum-1)
    if distToAt(i) < 1.5 % Threshold. Edit if needed
         contA = contA + 1; % Number of detections of A
    else
       if contA >= 10 % Sound detected during 1 s. Edit if needed
           disp(['Sound A identified ',num2str(contA/10),' s']);
           NumDet = NumDet +1;
           WinLen(NumDet) = contA; % Lenghth of detection window
           WinPos(NumDet) = i - contA; % Position of window
       end
       contA = 0;
    end

    if distToBt(i) < 1.5 % Threshold. Edit if needed
       contB = contB + 1;
    else
       if contB >= 10 % Sound detected during 1 s. Edit if needed
           disp(['Sound B identified ',num2str(contB/10),' s']);
           NumDet = NumDet +1;
           WinLen(NumDet) = contB; % Lenghth of detection window
           WinPos(NumDet) = i - contB; % Position of window
       end
       contB = 0;
    end
end

% Remove undesired noise segments
for i=1:NumDet
    iini = (WinPos(i)-1).*fs_native./10 + 1;
    ifin = (WinPos(i)-1+WinLen(i)).*fs_native./10;
    yt(iini:ifin) = NaN;
end

% Save full clean file
audiowrite('FullFileClean.wav',yt,Fs);
cd(Cur_DIR);