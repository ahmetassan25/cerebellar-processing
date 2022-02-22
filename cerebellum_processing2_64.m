
% Cerebellum Signal Processing Master File
clear all
close all 
tic
trial =1;


% This m-file is designed to be used as a single program from which several
% different an array of signal processing functions can be called. Each ofR
% the functions is saved as an individual m-file, and m-files are stored in
% an accompanying folder titled 'sp_mfiles'. Prior to calling each m-file   
% there is a section where necessary parameters are defined. In order to
% save processing time, unnecessary functions and parameter definitions can
% be ignored.

% Functions available in this version are called in the following order:
%     0.001) Play the video file frame-by-frame  
%     0.01) Average the DATA  and break into segments
%     0.02) Search for bad channels and remove it
%     0.03) Cerebellum Common Sourced Signal OUT
%     0.04) Pwelch Analysis
%     0.1) takes the differential of the signal
%     0.5) runs a wavelet denoising algorithim on the signal
%     1) Plot unfiltered data2
%     2) Plot freq response (FFT) of 1
%     3) Filter data for population activity
%     4) Plot filter Transfer Function designed in 3
%     5) Plot filtered Population Activity
%     6) Calculate and plot signal envelope (rectified moving avg of population activity)
%     7) Filter the original signal for LFP 
%     7.01)Plotting filtered LFPs in electrode orientation
%     7.02)Cross-Correlation Across Each Evoked Segments for horizontal and vertical electrode orientation 
%     7.03)Auto -Correlation Across Each Evoked Cycles and Each channels
%     7.1) Normalizing the LFP_data
%     7.5) plot the LFP signals
%     8) Plot rec-averaged version of LFP signals
%     9) Filter the original signal for various frequency bands and plot
%     10) Remove the bad channels from the LFP data
%     20) Calculate principle components and plot
%     30) Plot independent components and plot the rec-ave version
%    40) Coherence Function between contacts of R distance, vertically and horizontally
%     50) Average Coherence between contacts of R distance as a function of time
%    60) Cross correlations on the array as a function of time
%   100) Show Movie
%   110) Show Movie and a running marker on a channel of LFP potentials
%   120) Plot 2D correlation of population activity
%   130) Plot 2D correlation of LFP 
%   140) Movie of 2D correlations as a function of time
%   150) spectrogram for the unfiltered data
%   160) running average of the zero crossing time of the LFPs shown on the array
%   170) calculate force
%   180) Evoked Response Analysis
%   180.1) Simple Way to look at EPs
%   182) CF and MF Analysis, run it with 180


%   190) Spatial Organization Amplitude Analysis
%   200) Temporal Cross-Correlation Analysis
%   201) Half-Wave Rectify Evoked Potential Analysis
%   202) Peak Temporal Delay Analysis
%   203) CF and MF Analysis , run it with 201
%   204) Shading Area plots
%   205) Correlation Analysis Mapped in Electrode Configuration
%   207) Removing 60Hz noise by averaging only for quite data(need to work on that)
%   210) Coherence for desired two channels
%   300) shortcut to plotting functions


% Enter an array containing a list of the corresponding numbers for the
% functions you wish to run. For example, if you would like to only plot
% the unfiltered data (1), filter the data (3), and then plot the filtered
% data (5), enter functions_to_run = [1 3 5]. (Note: if previously existing
% unfiltered, filtered, or envelope data is used, enter names of correct
% files in following section and comment/uncomment appropriately in Data 
% Conditioning section below.)

functions_to_run =[0.01 7] ;

% Begin by setting the name of the raw data file to read in and the name of
% the files which the unfiltered data, filtered data, and signal envelope
% will be written to. All data files are stored within subfolders that are
% in a folder titled 'data'. Raw daq files are stored in a folder titled
% 'raw_daq', unfiltered dat files (which have been compensated for
% amplifier gain) are stored in a folder titled 'unfiltered', filtered datp
% files are stored in a folder titled 'filtered', and envelope dat files
% are stored in a folder titled 'envelope'. If reading in data from an
% existing csv dat file, comment out daqread command and enter name of
% existing dat file (without the .data extension) into file_name.

% file_name = ['CC5_18Apr2008_Trial' num2str(trial)];     % user insert appropriate file name
file_name = ['trial' num2str(trial) ];     % user insert appropriate file name
% num2str(trial)
% raw_data = daqread(['data/raw_daq/' file_name '.daq']);
% unfiltered_data_name = ['data/unfiltered/' file_name '.dat'];
% filtered_data_name = ['data/filtered/' file_name '.dat'];
% envelope_name = ['data/envelope/' file_name '.dat'];


dd=load(['/' file_name] );
raw_data=dd.data;
raw_dataA=raw_data;
% raw_dataA = daqread(['/' file_name 'A.daq']);
% raw_dataB = daqread(['/' file_name 'B.daq']);

% raw_data = [raw_dataA ];
% Set the recording parameters according to what was used for signal
% acquisition.
dur=1;
number_channels =1;
amplifier_gain = 100;
sampling_rate = 25000;  % in samples/sec
video_frame_rate = 30;  % in frames/sec
video_dur        = 1;







% Evoked Potential Parameters
puff_arrival = 5610 ;

% mech_arrival = 800;
% puff_arrival = mech_arrival; 
% 
% elec_arrival = 800;
% puff_arrival = elec_arrival;

signal_end = 5790 ;


EPs1 = 4900;
EPs2 = 5070;

EPs3 = EPs1;
EPs4 = EPs2;

marg =50;

t=1/sampling_rate:1/sampling_rate:dur; % Define the time axis
t=t*1000;

% Enter an array of indices separated by spaces that correspond to bad
% channels (open circuit, corrupted, abnormal response, etc...) For
% example, if channels 3, 8, and 13 are bad enter [3 8 13]. If no channels
% are bad enter bad_channels = [0]. (Bad channels are plotted in 
% different color so as to distinguish them in figures.)

bad_channels = [ ];

% Create an array that matches the number of each channel with its
% corresponding spatial index on the recording array. This uses a 6 x 6
% grid to plot 32 channels (the 4 corners of the grid are set at zero. Each
% index corresponds to electrode number and value corresponds to location
% on the grid. 
%   THESE ARE THE CONTACT NUMBERS ACCORDINGS TO MULTICHANNEL SYSTEMS (from
%   the contact side)
% ===============================  this is the ground contact
%                     32      31      30      29
%
%
%           26      25      24      23      27      28
%
%
%           19      18      17      20      21      22
%
%
%           14      15      16      13      12      11
%
%
%            7       8       9      10       6       5
%
%
%                      1       2       3       4
% =============================== this is the other ground contact

%positions = [32 33 34 35 30 29 25 26 27 28 24 23 22 19 20 21 15 14 13 16 17 18 10 9 8 7 11 12 5 4 3 2];

% THESE ARE THE AMPLIFIER CHANNEL NUMBERS ON THE SAME CONFIGURATION AS ABOVE WITH THE 
% PCB design of the electode connector NOT REVERSED.

%
%                  17      1     18       2 
%
%        20      4       21    5         3     19
%
%         7      24       8     23       6      22
%
%        26      9      25     10      27     11
%
%        13    29     12     28      30     14
%
%                16      32    15       31

% THESE ARE THE AMPLIFIER CHANNEL NUMBERS ON THE SAME CONFIGURATION AS ABOVE WITH THE 
% PCB design of the electode connector IF THE PCB IS REVERSED. IN THIS
% ANIMAL THIS IS THE MAP TO USE.

%                1      17       2        18
%
%        4      20    5       21       19      3
%
%       23     8      24     7         22     6
%
%       10     25    9       26       11     27
%
%        29    13    28     12       14     30
%
%                32    16      31       15

%positions=[6 12 5 11 17 4 10 16 3 9 15 2 8 14 1 7];
%positions= [33 35 29 26 28 23 19 21 14 16 18 9 7 12 4 2 32 34 30 25 27 24 22 20 15 13 17 10 8 11 5 3];
positions_curveMEA=[0 0 0 0 3 2 1 5 6 7;21 20 22 19 23 18 24 17 4 18;30 29 28 31 27 0 26 9 25 10; 0 0 0 0 14 13 15 12 16 11];
positions_curveMEA=positions_curveMEA'
positions= [16 14 12 10 9 11 13 15 23 21 19 17 18 20 22 24 8 6 4 2 1 3 5 7 31 29 27 25 26 28 30 32];
positions_31channels= [8 16 6 14 4 12 2 10 1 9 3 11 5 13 7 15 31 23 29 21 27 19 25 17 26 18 28 20 30 22 24];
pos_x = [21 20 22 19 23 18 24 17 ; 5 4 6 3 7 2 8 1 ; 12 13 11 14 10 15 9 16; 28 29 27 30 26 31 25 32];
pos_y = [21 5 12 28 ; 20 4 13 29 ; 22 6 11 27; 19 3 14 30; 23 7 10 26; 18 2 15 31; 24 8 9 25; 17 1 16 32];
% positions_sorted = [21 20 22 19 23 18 24 17;5 4 6 3 7 2 8 1; 12 13 11 14 10 15 9 16;28 29 27 30 26 31 25 32]' ;
positions_sorted = [9 7 11 5 13 3 15 1; 10 8 12 6 14 4 16 2; 24 26 22 28 20 30   18 31; 23 25 21 27 19 29 17 32]';
positions_sorted(32) = [];
% positions_W_31channels2_splots=[9 16 7 14 5 12 6 11 8 13 10 15 20 17 18 19 30 27 40 25 39 23 38 21 29 22 35 24 36 26 37];

positions_W_31channels2_splots= [39    26    37    24    35    22    36    21    38    23    40    25    30    27    28    29    20    17    10    15     9    13     8    11    19    12     5    14   6    16     7]
positions_new_electrode=[ 35 4 36 8 33 5 34 6 31 9 32 10 29 11 30 12 27 13 28 14 25 15 21 19 23 17 26 16 18];

% data as referred to input. Save this into previously specified dat file
% using csvwrite to create comma separated value file. If unfiltered data
% file already exists, comment out first two lines and read in data file
% using csvread. Data is then scaled according to order of magnitude for
% plotting ('uV', 'mV', or 'V').

 %unfiltered_data = raw_data(:,35:68)/amplifier_gain;
 unfiltered_data = raw_data/amplifier_gain;
 signal_length = min(size(unfiltered_data));
 unfiltered_data = unfiltered_data(:,1:number_channels);
 
 for x = 1:length(functions_to_run)
    if functions_to_run(x) == 0.001
        
        frame_Size = 1;
        
        implay(['videos/trial' num2str(trial) '.avi'],frame_Size)
 
    end
 end
 
 
 for x = 1:length(functions_to_run)
    if functions_to_run(x) == 0.03


unfiltered_data_mean = mean(unfiltered_data(:,1:31)');
% 
[coef, score, latent]=princomp(unfiltered_data(:,1:31));
score1=score(:,1);
coef(:,1)=zeros(number_channels,1);
unfiltered_data=score*coef';
% unfiltered_data = unfiltered_data';


% 
 for i=1:31
%      unfiltered_data(:,i) = unfiltered_data(:,i) - unfiltered_data_mean';
          unfiltered_data(:,i) = unfiltered_data(:,i) - score(:,2);

 end
 
 figure;plot(unfiltered_data_mean)
 
 
 
 
    end
 end
 
 
 
for x = 1:length(functions_to_run)
    if functions_to_run(x) == 0.02
        
%             number_channels=31;
        
            bad_channel = find(var(unfiltered_data(:,1:number_channels)) >= max(var(unfiltered_data(:,1:number_channels))));
            
            meanbar = mean(var(unfiltered_data(:,1:number_channels)));
            stdbar= std(var(unfiltered_data(:,1:number_channels)));
            
            outliers = find(abs(var(unfiltered_data(:,1:number_channels))-meanbar) > 2*stdbar)
            
            figure;stem(var(unfiltered_data(:,1:number_channels)))
            unfiltered_data = unfiltered_data(:,1:number_channels);
            unfiltered_data(:,outliers)=[];
             signal_length = min(size(unfiltered_data));
            number_channels = signal_length;
             
    end
    
end

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 0.01
%         unfiltered_data=unfiltered_data(:,1:signal_length);

duration = 10;
LOOP = 10;
Fs = sampling_rate;
DURATION=duration/LOOP;
N=2048;

%% Taking the Power spectrum first for each period then average it.

% for ii=1:LOOP
%     
% [Pww(:,ii), F] = pwelch(detrend(unfiltered_data((ii-1)*DURATION*Fs+1:ii*DURATION*Fs,:)),hanning(N),N/2,N,sampling_rate);
% Pww(:,ii)  = 5*log10(Pww(:,ii).^2);
% end
% 
% averaged_Pww = sum(Pww(:,1:20)',1)/20;
% averaged_Pww = averaged_Pww';

% figure;plot(F,averaged_Pww)       
% Apparently it does not make any difference taking the power spect from
% the averaged DATA

%% Averaging the data

tt=1/Fs:1/Fs:duration;

DATA = zeros(min(size(unfiltered_data)), Fs*DURATION);
% DATA2 = zeros(min(size(unfiltered_data_mean)), Fs*DURATION);

unfiltered_data = unfiltered_data';
% unfiltered_data_mean_averaged = unfiltered_data_mean';

for Y=1:LOOP,
 
    DATA = DATA+unfiltered_data(:,(Y-1)*DURATION*Fs+1:Y*DURATION*Fs);
%     DATA2 = DATA2+unfiltered_data_mean(:,(Y-1)*DURATION*Fs+1:Y*DURATION*Fs);
end

DATA=DATA'/LOOP;
% DATA2=DATA2'/LOOP;

% unfiltered_data = unfiltered_data';

unfiltered_data = DATA;

%% Break each Evoked Signals into segments  


% for i=1:32
% for p=1:LOOP
% DATA_Segments(:,p,i) = unfiltered_data((p-1)*Fs+1:p*Fs,i);
% end
% end



    end
end


for x = 1:length(functions_to_run)
    if functions_to_run(x) == 0.04
N=2048*4;
% hold on;
[Pww_cer, F] = pwelch(unfiltered_data(:,1:31),hanning(N/2),N/4,N/2,sampling_rate); % Cerebellum Power Spectrum
% [Pww_cor, F] = pwelch(unfiltered_data(:,32:62),hanning(N/2),N/4,N/2,sampling_rate); % Cerebellum Power Spectrum

cc = hsv(50);


Pww_cer = 10*log10(Pww_cer);
% Pww_cor = 10*log10(Pww_cor);

hold on; 
figure(1);
plot(F,mean(Pww_cer'),'color','b');axis([0 1000 -180 -90])
whitebg('w')
% hold on;plot(F,mean(Pww_cor'),'color','r');axis([0 1000 -180 -90])


% hold on;plot(F,Pww_cor,'y');grid;legend('Cerebellum Quiet','Cerebrum Quiet')

% figure;plot(F,Pww_cer);hold on
% [Pww_mean, F] = pwelch(detrend(D ATA2),hanning(N),N/2,N,sampling_rate);
% Pww_mean = 5*log10(Pww_mean.^2);
% plot(F,Pww_mean,'c');legend('Cerebellum Mean Out','Cerebellum  Mean DATA')

    end
end

 
clear raw_data
for x = 1:length(functions_to_run)
    if functions_to_run(x) == 0.1

row1=[zeros(size(unfiltered_data(:,14)))'; unfiltered_data(:,14)'; unfiltered_data(:,11)'; unfiltered_data(:,22)'; unfiltered_data(:,19)'; zeros(size(unfiltered_data(:,14)))'];
row2=[unfiltered_data(:,31)'; unfiltered_data(:,30)'; unfiltered_data(:,27)'; unfiltered_data(:,6)'; unfiltered_data(:,3)'; unfiltered_data(:,2)'];
row3=[unfiltered_data(:,15)'; unfiltered_data(:,28)'; unfiltered_data(:,10)'; unfiltered_data(:,23)'; unfiltered_data(:,5)'; unfiltered_data(:,18)'];
row4=[unfiltered_data(:,32)'; unfiltered_data(:,12)'; unfiltered_data(:,25)'; unfiltered_data(:,8)'; unfiltered_data(:,21)'; unfiltered_data(:,1)'];
row5=[unfiltered_data(:,16)'; unfiltered_data(:,29)'; unfiltered_data(:,9)'; unfiltered_data(:,24)'; unfiltered_data(:,4)'; unfiltered_data(:,17)'];
row6=[zeros(size(unfiltered_data(:,14)))'; unfiltered_data(:,13)'; unfiltered_data(:,26)'; unfiltered_data(:,7)'; unfiltered_data(:,20)'; zeros(size(unfiltered_data(:,14)))'];
%rows=[row1; row2; row3; row4; row5; row6];

s=zeros(150150,3);




   signaldif1(:,1)=row1(3,:)'-row1(2,:)';
   signaldif1(:,2)=row1(5,:)'-row1(4,:)';
 
   s=zeros(150150,3);

signaldif=s;

    for j=1:3
        i=j*2;
        signaldif2(:,j)=row2(i,:)'-row2(i-1,:)';
    end
    
    s=zeros(150150,3);

signaldif=s;

    for j=1:3
        i=j*2;
        signaldif3(:,j)=row3(i,:)'-row3(i-1,:)';
    end
    
    s=zeros(150150,3);

signaldif=s;

    for j=1:3
        i=j*2;
        signaldif4(:,j)=row4(i,:)'-row4(i-1,:)';
    end
    
    s=zeros(150150,3);

signaldif=s;

    for j=1:3
        i=j*2;
        signaldif5(:,j)=row5(i,:)'-row5(i-1,:)';
    end
    
    s=zeros(150150,3);

signaldif=s;

  signaldif6(:,1)=row6(3,:)'-row6(2,:)';
   signaldif6(:,2)=row6(5,:)'-row6(4,:)';
   
    signaldif=[signaldif1 signaldif2 signaldif3 signaldif4 signaldif5 signaldif6];
    unfiltered_data=signaldif;
    number_channels = 16;
    end
end
%csvwrite(unfiltered_data_name, unfiltered_data);
% unfiltered_data = csvread(unfiltered_data_name);
% filtered_data = csvread(filtered_data_name);
% envelope_data = csvread(envelope_name);

scale = 'uV';

if strcmp(scale, 'uV')
    scale_value = 1000000;
elseif strcmp(scale, 'mV')
    scale_value = 1000;
elseif strcmp(scale, 'V')
    scale_value = 1;
end

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 0.5

for i=1:16
%i=500;
Signal=unfiltered_data(:,i);
Signal=Signal-mean(Signal);
Signal=detrend(Signal);
%Signal=Design;
%Signal=scores(:,5);
%Signal=hrf(1:72);

fs = 30000;
Ts = 1/fs;
 x = [0:1/fs:length(Signal)*Ts];
 N = length(Signal);
 index = 0:N-1;


%h0 = [1/sqrt(2) 1/sqrt(2)];
% Daubechies - 4
h0 = [0.4830 0.8365 0.2241 -0.1294];
% Daubechies - 10
%h0 = [0.32580343 1.01094572 0.8922014 -0.03957503 -0.26450717 0.0436163 0.0465036 -0.01498699];
%Daubechies - 2
%h0 = [0.3415    0.5915    0.1585   -0.0915];
%h0 = [0.5000    0.5000];
NSteps = log2(N);
coefs = (2^(-NSteps/2)).*Signal;

j = 1;
[LP1 HP1] = waveanalysis(coefs, h0);
% The output is the low and high pass bands of the signal.
%
% Plot the signal bands.
index = 0:2^(NSteps-j)-1;
time = index.*Ts*2^j;
% figure(1)
% subplot(121)
% plot(time,LP1)
% subplot(122)
% plot(time,HP1)

% This is the second pass. Here the LP range is investigated.
j = 2;
[LP2 HP2] = waveanalysis(LP1, h0);
index2 = 0:2^(NSteps-j)-1;
time = index2.*Ts*2^j;
 %Plot the signal bands.
 %figure(2)
 %subplot(121)
 %plot(time,LP2)
 %subplot(122)
 %plot(time,HP2)

% This is the second pass. Here the LP range is investigated.
j = 3;
[LP3 HP3] = waveanalysis(LP2, h0);
index3 = 0:2^(NSteps-j)-1;
time = index3.*Ts*2^j;
 %Plot the signal bands.
 %figure(3)
 %subplot(121)
 %plot(time,LP3)
 %subplot(122)
 %plot(time,HP3)
 

 
 % This is the second pass. Here the LP range is investigated.
j = 4;
[LP4 HP4] = waveanalysis(LP3, h0);
index4 = 0:2^(NSteps-j);
time = index4.*Ts*2^j;
 %Plot the signal bands.
 %figure(4)
 %subplot(121)
 %plot(time,LP4)
 %subplot(122)
 %plot(time,HP4)
 j = 5;
[LP5 HP5] = waveanalysis(LP4, h0);
index4 = 0:2^(NSteps-j);
time = index4.*Ts*2^j;
 
 
 j = 6;
[LP6 HP6] = waveanalysis(LP5, h0);
index4 = 0:2^(NSteps-j);
time = index4.*Ts*2^j;
 
 j = 7;
[LP7 HP7] = waveanalysis(LP6, h0);
index4 = 0:2^(NSteps-j);
time = index4.*Ts*2^j;
 
  j = 8;
[LP8 HP8] = waveanalysis(LP7, h0);
index4 = 0:2^(NSteps-j);
time = index4.*Ts*2^j;
 
  j = 9;
[LP9 HP9] = waveanalysis(LP8, h0);
index4 = 0:2^(NSteps-j);
time = index4.*Ts*2^j;


  j = 10;
[LP10 HP10] = waveanalysis(LP9, h0);
index4 = 0:2^(NSteps-j);
time = index4.*Ts*2^j;


 sigma1=mad(HP1)/0.6745;
 sigma2=mad(HP2)/0.6745;
 sigma3=mad(HP3)/0.6745;
 sigma4=mad(HP4)/0.6745;
 sigma5=mad(HP5)/0.6745;
 sigma6=mad(HP6)/0.6745;
 sigma7=mad(HP7)/0.6745;
 sigma8=mad(HP8)/0.6745;
 sigma9=mad(HP9)/0.6745;
 lHP1=length(HP1);
 lHP2=length(HP2);
 lHP3=length(HP3);
 lHP4=length(HP4);
 lHP5=length(HP5);
 lHP6=length(HP6);
 lHP7=length(HP7);
 lHP8=length(HP8);
 lHP9=length(HP9);
tHP1=sigma1.*sqrt(2.*log(lHP1));
tHP2=sigma2.*sqrt(2.*log(lHP2));
tHP3=sigma3.*sqrt(2.*log(lHP3));
tHP4=sigma4.*sqrt(2.*log(lHP4));
tHP5=sigma5.*sqrt(2.*log(lHP5));
tHP6=sigma6.*sqrt(2.*log(lHP6));
tHP7=sigma7.*sqrt(2.*log(lHP7));
tHP8=sigma8.*sqrt(2.*log(lHP8));
tHP9=sigma9.*sqrt(2.*log(lHP9));
for k=1:lHP1
    sig=HP1(k);
    sng=sign(sig);
    thresh=abs(sig)-tHP1;
    if thresh < 0
        rep=0;
    else
        rep=thresh;
    end
    dHP1(k)=sng*rep;
end

for k=1:lHP2
    sig=HP2(k);
    sng=sign(sig);
    thresh=abs(sig)-tHP2;
    if thresh < 0
        rep=0;
    else
        rep=thresh;
    end
    dHP2(k)=sng*rep;
end

for k=1:lHP3
    sig=HP3(k);
    sng=sign(sig);
    thresh=abs(sig)-tHP3;
    if thresh < 0
        rep=0;
    else
        rep=thresh;
    end
    dHP3(k)=sng*rep;
end

for k=1:lHP4
    sig=HP4(k);
    sng=sign(sig);
    thresh=abs(sig)-tHP4;
    if thresh < 0
        rep=0;
    else
        rep=thresh;
    end
    dHP4(k)=sng*rep;
end

for k=1:lHP5
    sig=HP5(k);
    sng=sign(sig);
    thresh=abs(sig)-tHP5;
    if thresh < 0
        rep=0;
    else
        rep=thresh;
    end
    dHP5(k)=sng*rep;
end

for k=1:lHP6
    sig=HP6(k);
    sng=sign(sig);
    thresh=abs(sig)-tHP6;
    if thresh < 0
        rep=0;
    else
        rep=thresh;
    end
    dHP6(k)=sng*rep;
end

for k=1:lHP7
    sig=HP7(k);
    sng=sign(sig);
    thresh=abs(sig)-tHP7;
    if thresh < 0
        rep=0;
    else
        rep=thresh;
    end
    dHP7(k)=sng*rep;
end

for k=1:lHP8
    sig=HP8(k);
    sng=sign(sig);
    thresh=abs(sig)-tHP8;
    if thresh < 0
        rep=0;
    else
        rep=thresh;
    end
    dHP8(k)=sng*rep;
end

for k=1:lHP9
    sig=HP9(k);
    sng=sign(sig);
    thresh=abs(sig)-tHP9;
    if thresh < 0
        rep=0;
    else
        rep=thresh;
    end
    dHP9(k)=sng*rep;
end


rLP8 = dewaveanalysis(dHP9',LP9, h0);
rLP7 = dewaveanalysis(dHP8',rLP8(4:end), h0);
rLP6 = dewaveanalysis(dHP7',rLP7(5:end), h0);
rLP5 = dewaveanalysis(dHP6',rLP6(4:end), h0);
rLP4 = dewaveanalysis(dHP5',rLP5(5:end), h0);
rLP3 = dewaveanalysis(dHP4',rLP4(5:end), h0);
%L=length(rLP3)-lHP3;

%rLP2 = dewaveanalysis(dHP3',rLP3(L:length(rLP3)-1)*20, h0);
rLP2 = dewaveanalysis(dHP3',rLP3(4:end), h0);
L=length(rLP2)-lHP2;
%rLP1 = dewaveanalysis(dHP2',rLP2(L+1:length(rLP2))*1, h0);
rLP1 = dewaveanalysis(dHP2',rLP2(4:end), h0);
rrcomp = dewaveanalysis(dHP1',rLP1(1:25000), h0); %%(L+1:length(rLP1))
%size(rrcomp)
 processedsignal(:,i)=rrcomp;
end
size(processedsignal)
unfiltered_data=processedsignal;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create arrays to be used for x-data in future plots. These include time,
% frequency, and video frame number.

number_data_points = length(unfiltered_data);
time = (1:number_data_points)/sampling_rate;
freq = (1:number_data_points)*sampling_rate/number_data_points;
video = (1:number_data_points)*video_frame_rate/sampling_rate;

% 1) Plot unfiltered data
%
% This function plots the unfiltered gain-compensated, scaled data versus
% either a time scale, sample scale, or video frame scale. Specify x_data
% according to which scale is appropriate (time, sample, or video). Also,
% set the appropriate minimum and maximum x and y values to plot.


for x = 1:length(functions_to_run)
    if functions_to_run(x) == 1

        x_data = time;
        
        x_min = 0;
        x_max = max(x_data);
        
        y_min = min(unfiltered_data*scale_value);
        y_max = max(unfiltered_data*scale_value);

        x_label = 'Time (s)';
%         x_label = 'Video Frame';
%         x_label = 'Sample Number';
%         y_label = {'Unfiltered Signal'; ['(' scale ')']};
        y_label = ['Unfiltered Signal (' scale ')'];

        run sp_mfiles/plot_unfiltered_data
%    run sp_mfiles\plot_unfiltered_data_quadrant
    end
end


% 2) Plot frequency response (FFT)
%
% This function computes and plots the fast fourier transform (FFT) of the
% unfiltered gain-compensated data versus the frequency.


for x = 1:length(functions_to_run)
    if functions_to_run(x) == 2
number_of_data_points=25000
        unfiltered_fft = abs(fft(detrend(unfiltered_data)));
        
        x_data = freq;
        
        x_min = 0;
        x_max = sampling_rate/2;
        
        y_min = 0;
%         y_max = max(unfiltered_fft);
        for y = 1:number_channels
            y_max(y) = 0.05; %mean(unfiltered_fft(round(number_data_points/3):round(number_data_points/2),y))+25*std(unfiltered_fft(round(number_data_points/3):round(number_data_points/2),y));
        end
        
        x_label = 'Frequency (Hz)';
        y_label = 'FFT Magnitude';
        
        run sp_mfiles\plot_fft
        
    end
end


% 3) Filter data (Bandpass FIR) for population activity
%
% This function constructs a bandpass FIR filter and filters the
% gain-compensated data accordingly. Filtered data is saved into a
% pre-specified dat file. Filter parameters (number of coefficients, low
% cutoff, and high cutoff) are specified prior to calling the function.

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 3
        
        filter_order = 1000;
        low_cut =200;           % freq in Hz
        high_cut = 600;        % freq in Hz

        high_pass_norm = high_cut/(sampling_rate/2);
         low_pass_norm = low_cut/(sampling_rate/2);
                
        b_filt_band = fir1(filter_order, [ low_pass_norm high_pass_norm  ]);
        a_filt_band=1;
        filtered_data = filtfilt(b_filt_band, a_filt_band, unfiltered_data);
              
%        csvwrite(filtered_data_name, filtered_data);
        
    end
end


% 4) Plot filter
%
% This function plots the frequency response of the filter being used. Set
% the appropriate minimum and maximum x and y values to plot.


for x = 1:length(functions_to_run)
    if functions_to_run(x) == 4

        x_data = freq;
        
        x_min = 0;
        x_max = sampling_rate/2;
        
        y_min = 0;
        y_max = 1.1;

        run sp_mfiles\plot_filter
        
    end
end


% 5) Plot filtered data
%
% This function plots the filtered gain-compensated, scaled data versus
% either a time scale, sample scale, or video frame scale. Specify x_data
% according to which scale is appropriate (time, sample, or video). Also,
% set the appropriate minimum and maximum x and y values to plot.


for x = 1:length(functions_to_run)
    if functions_to_run(x) == 5

        x_data = time;
        
        x_min = 0;
        x_max =max(x_data);
        
        y_min = -70*ones(1,68); %min(filtered_data*scale_value);
        y_max = 70*ones(1,68); % max(filtered_data*scale_value);
        
        x_label = 'Time (s)';
%         x_label = 'Video Frame';
%         x_label = 'Sample Number';
%         y_label = {'Filtered Signal'; ['(' scale ')']};
        y_label = ['Filtered Signal for population activity (' scale ')'];

        run sp_mfiles\plot_filtered_data
        
    end
end


% 6) Calculates and plots signal envelope (rectified moving avg of filtered data)
%
% This function calculates the signal envelope by rectifying and then
% performing a moving average calculation on the filtered data. Set the 
% appropriate bin length for calculating the moving average in seconds.


for x = 1:length(functions_to_run)
    if functions_to_run(x) == 6

        bin_length = 0.05;     % in sec
        
      run sp_mfiles\calculate_envelope
        
%         csvwrite(envelope_name, envelope_data);

        x_data = time;
        
        x_min = 0;
        x_max = max(time); % max(x_data);
        
        y_min = 5e-6; % ones(1,32)* min(min(envelope_data*scale_value));
       y_max = ones(1,32)*40;
        %y_max = ones(1,32)*max(max(envelope_data*scale_value));

        x_label = 'Time (s)';
%         x_label = 'Video Frame';
%         x_label = 'Sample Number';
%         y_label = {'Signal Envelope'; ['(' scale ')']};
        y_label = ['Signal Envelope (' scale ')'];

        run sp_mfiles\plot_envelope
    end
end


% 7) Filter data (Bandpass FIR) for LFP 
%
% This function constructs a bandpass FIR filter and filters the
% gain-compensated data accordingly. Filtered data is saved into a
% pre-specified dat file. Filter parameters (number of coefficients, low
% cutoff, and high cutoff) are specified prior to calling the function.

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 7
        
            
            L=length(unfiltered_data);
            T=1000/sampling_rate;
            % design a low - pass
      
            fch=10;
            fcl=10000;
            
            
            fn2=2*fcl/(sampling_rate);
            fn1=2*fch/(sampling_rate);
          
            Wn=[fn1 fn2];
            [b1,a1]=ellip(4,4,30,Wn);
            
            [b,a]=butter(2,fn1,'high');
            
%            LFP_data_2=filtfilt(b1,a1,raw_dataA);
            LFP_data=filtfilt(b,a,unfiltered_data);
           
%             LFP_data=filtfilt(b,a,DATA);
            [b,a]=butter(2,fn2,'low');
          
            LFP_data=filtfilt(b,a,LFP_data);
            
%             LFP_data2=filtfilt(b,a,LFP_data2);


            fc1=5;
            fc2=12;
            fc3=12;
            fc4=30;
            fc5=20;
            fc6=50;
            fc7=80;
            fc8=300;
            
            
              fn1=2*fc1/(sampling_rate);
            fn2=2*fc2/(sampling_rate);
            
             [b,a]=butter(5,fn2,'low');
            LFP_data1=filtfilt(b,a,unfiltered_data);
            [b,a]=butter(5,fn1,'high');
            LFP_data1=filtfilt(b,a,LFP_data1);
            
            
            fn4=2*fc4/(sampling_rate);
            fn3=2*fc3/(sampling_rate);
%             fn2=[fn3 fn4];              
            [b,a]=butter(5,fn4,'low');
            LFP_data2=filtfilt(b,a,unfiltered_data);
            [b,a]=butter(5,fn3,'high');
            LFP_data2=filtfilt(b,a,LFP_data2);
%             
              fn6=2*fc6/(sampling_rate);
            fn5=2*fc5/(sampling_rate);
            [b,a]=butter(5,fn6,'low');
            LFP_data3=filtfilt(b,a,unfiltered_data);
            [b,a]=butter(5,fn5,'high');
            LFP_data3=filtfilt(b,a,LFP_data3);
%        
            


             fn7=2*fc7/(sampling_rate);
            fn8=2*fc8/(sampling_rate);
            [b,a]=butter(5,fn8,'low');
            LFP_data4=filtfilt(b,a,unfiltered_data);
            [b,a]=butter(5,fn7,'high');
            LFP_data4=filtfilt(b,a,LFP_data4);
           
            % plot the filter amplitude and phase
%             [h,w] = freqz(b,a,10000);	% this calculates the transfer function at N different points of digital frequency; h is the transfer function at w frequencies
%             f=w/(2*pi*T);               % convert the digital frequency array to analog frequency in KHz
%             figure;						
%             subplot(211);
%             plot(f,abs(h));			% amplitide
%             grid
%             subplot(212)
%            % plot(f,angle(h));	% phase
%             grid
%             xlabel('frequency (kHz)');

LFP_data_cer =LFP_data(:,1:number_channels);
% LFP_data_cor = LFP_data(:,number_channels+1:2*number_channels);

unfiltered_data_cer = unfiltered_data(:,1:number_channels);
% unfiltered_data_cor = unfiltered_data(:,number_channels:2*number_channels);

% for p = 1:number_channels
% LFP_data_sorted(:,p) = LFP_data_cer(:,positions_sorted(p));
% unfiltered_DATA(:,p) = unfiltered_data(:,positions_sorted(p));
% 
% end
% for p = 1:number_channels
% LFP_data_sorted(:,p+number_channels) = LFP_data_cor(:,positions_sorted(p));
% unfiltered_DATA(:,p+number_channels) = unfiltered_data_cor(:,positions_sorted(p));
% 
% end
%  



    end
end
    
%     [max_indice max_chan] = find(LFP_data(puff_arrival:signal_end,:) == max(max(LFP_data(puff_arrival:signal_end,:))));
%     max_indice = puff_arrival + max_indice
    
    for x = 1:length(functions_to_run)
    if functions_to_run(x) == 7.01
            figure;
        for l = 1:min(size(pos_x))
            subplot(4,1,l)
            plot(t,1000*LFP_data(:,pos_x(l,:)));
            grid on
%             axis([ puff_arrival-5*marg  puff_arrival+10*marg   -5e-5 3e-5 ])
            axis tight
            title(['trial' num2str(trial) ,' ; ','Contra Arm post Puff-sti 55 min'] )
            ylabel('mV','fontsize',12)
            xlabel('ms','fontsize',12,'fontweight','b')
            
            
            
            
            
        end
        figure;
        for l = 1:length(pos_y)
            subplot(1,8,l) 
            plot(t,1000*LFP_data(:,pos_y(l,:)));
            grid on
            axis([ (EPs1-5*marg)*1000/sampling_rate (EPs2+5*marg)*1000/sampling_rate   -5e-2 10e-2 ])
%             axis tight
            title(['trial' num2str(trial) ,' ; ','CA Puff 55min'] )
            ylabel('mV','fontsize',12)
            xlabel('ms','fontsize',12,'fontweight','b')

        end
    end
    end
        
        for x = 1:length(functions_to_run)
    if functions_to_run(x) == 7.02
        %% Correlation analysis for Averaged DATA
        % 1.part Cross - Correlation is processed for Each loop
        figure;
        for l = 1:min(size(pos_x))
        for i=1:LOOP
            
            subplot(4,1,l)
            surface(corrcoef(squeeze(DATA_Segments(:,i,pos_x(l,:)))))
            axis tight
            title('Corr Distribution across evoked stimulations')
            xlabel(['Line' num2str(l) ] )
            ylabel('8 channels')
            pause(1);end
        end
    figure;
         for l = 1:min(size(pos_y))
         for i=1:LOOP
            
            subplot(1,8,l)
            surface(corrcoef(squeeze(DATA_Segments(:,i,pos_y(l,:)))))
            axis tight
            title('Corr Distribution across evoked stimulations')
            xlabel(['Line' num2str(l) ] )
            ylabel('4 channels')
            pause(1);end
         end
    end
        end
         % 2. Part Autocorrelation is implemented for each channel's
         % averages respect to one selected segment
         % Now upgraded for each 20 LOOP
         % Every time one segment is selected and the autocorr is
         % calculated over the each evoked cycles
         % plotting shows each channels and each cycles(# of loops 20 in
         % our case)auto-correlation
         
         
        for x = 1:length(functions_to_run)
    if functions_to_run(x) == 7.03
         
         for k=1:LOOP
        for z=1:32
         for p=1:LOOP
            cc = corrcoef(DATA_Segments(:,k,z),DATA_Segments(:,p,z));
            cc_corr(k,z,p) = cc(2);
         end
        end
         end
        
        figure;
        for z = 1:20
        subplot(5,4,z);
        surface(squeeze(cc_corr(z,:,:)));axis tight;colorbar % surface auto-correlation of each segments respect all individual channels
        end
        
    end
    end
    

% 7.1) Normalizing the data ; 

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 7.1
        
        
        LFP_data_test = LFP_data(5700:6200,1:32);
        
        min_LFP_data = min(min(LFP_data_test));
        max_LFP_data = max(max(LFP_data_test));
        desired_maxmin = [-1 1];
        dt= max(desired_maxmin) - min(desired_maxmin);
        difference = max_LFP_data - min_LFP_data ;
        n = dt / difference;
        norm_LFP_data = n.*LFP_data_test;
        min_norm_LFP_data = min(min(norm_LFP_data));
        
        padding_data = min_norm_LFP_data - min(desired_maxmin);
        
        norm_LFP_data = norm_LFP_data - padding_data;
%         
%         LFP_data = LFP_data2;
%         LFP_data2 = LFP_data2(100:end,1:31);
%       
%         
%         min_LFP_data2 = min(min(LFP_data2));
%         max_LFP_data2 = max(max(LFP_data2));
%         desired_maxmin = [-1 1];
%         dt= max(desired_maxmin) - min(desired_maxmin);
%         difference = max_LFP_data2 - min_LFP_data2 ;
%         n = dt / difference;
%         norm_LFP_data2 = n.*LFP_data2;
%         min_norm_LFP_data2 = min(min(norm_LFP_data2));
%         
%         padding_data2 = min_norm_LFP_data2 - min(desired_maxmin);
%         
%         norm_LFP_data2 = norm_LFP_data2 - padding_data2;
% %         
%         LFP_data3 = LFP_data3(:,1:32);
%         
%         min_LFP_data3 = min(min(LFP_data3));
%         max_LFP_data3 = max(max(LFP_data3));
%         desired_maxmin = [-1 1];
%         dt= max(desired_maxmin) - min(desired_maxmin);
%         difference = max_LFP_data3 - min_LFP_data3 ;
%         n = dt / difference;
%         norm_LFP_data3 = n.*LFP_data3;
%         min_norm_LFP_data3 = min(min(norm_LFP_data3));
%         
%         padding_data3 = min_norm_LFP_data3 - min(desired_maxmin);
%         
%         norm_LFP_data3 = norm_LFP_data3 - padding_data3;
%         
%         LFP_data4 = LFP_data4(:,1:32);
%         
%         min_LFP_data4 = min(min(LFP_data4));
%         max_LFP_data4 = max(max(LFP_data4));
%         desired_maxmin = [-1 1];
%         dt= max(desired_maxmin) - min(desired_maxmin);
%         difference = max_LFP_data4 - min_LFP_data4 ;
%         n = dt / difference;
%         norm_LFP_data4 = n.*LFP_data4;
%         min_norm_LFP_data4 = min(min(norm_LFP_data4));
%         
%         padding_data4 = min_norm_LFP_data4 - min(desired_maxmin);
%         
%         norm_LFP_data4 = norm_LFP_data4 - padding_data4;
    end
end



% 7.5) plot the LFP_data
for x = 1:length(functions_to_run)
    if functions_to_run(x) == 7.5
            DATAA=1e6*LFP_data;
            
                       
            t=1/sampling_rate:1/sampling_rate:length(DATAA)/sampling_rate;
            [M,N]=size(DATAA);
            figure
            title('LFP data')

                for p=1:number_channels,
                     subplot(6,6,positions(p));
                     temp=DATAA(:,p)-mean(DATAA(:,p));
                    plot(t,temp)
                    y_min = min(min(temp));
                    y_max =max(max(temp));
                    axis([0 5 y_min y_max ]);
%                     head=num2str(p,'ch %d');
%                     title(head)
                end
                
%             t=1/sampling_rate:1/sampling_rate:length(LFP_data)/sampling_rate;
%             plotsize=4;
%             [M,N]=size(LFP_data);
%             numfigs=floor(N/plotsize);
% 
%             for q=1:numfigs,
%                 figure
%                 for p=1:plotsize,
%                     subplot(plotsize,1,p)
%                     plot(t,DATAA(:,p+plotsize*(q-1)))
% %                    amp=5*abs(std((DATAA(:,p+plotsize*(q-1)))));
%                     y_min = min(min(DATAA(:,p+plotsize*(q-1))));
%                     y_max = max(max(DATAA(:, p+plotsize*(q-1))));
%                     axis([0 5 y_min y_max]);
%                     temp=(q-1)*plotsize+p;
%                     head=num2str(temp,'%d');
%                     ylabel(head)
%                     title('LFP data')
%                 end


%             end
    end
end
            
% 8) plot the rec-averaged version of LFP_data

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 8

            fc1=10		% corner frequency of the filter (in Hz)
            p=5;
            fn1=2*fc1/(sampling_rate);		% normalize the corner frequence with fsample
            fn=[fn1];               % low pass filter
           [b,a]=butter(5,fn);
           LFP_ave=filtfilt(b,a,abs(LFP_data));
           
            t=1/sampling_rate:1/sampling_rate:length(LFP_ave)/sampling_rate;
            [M,N]=size(LFP_ave);
            figure
            title('rec-averaged LFP')
            y_min = 0;
            y_max = 20e-6;  %max(max(LFP_ave));
                for p=1:number_channels,
                     subplot(6,6,positions(p));
                    plot(t,LFP_ave(:,p),'r')
                    amp=2*abs(mean((LFP_ave(:,p))));
                    axis([0 5 y_min y_max ]);
                    head=num2str(p,'%d');
                    ylabel(head)
                    
                end
    end
end
     
% 9) Filter data (Bandpass FIR) for LFP and plot the filtered signals
for x = 1:length(functions_to_run)
    if functions_to_run(x) == 9

        % remove the bad channels from unfiltered_data first
data=zeros(length(unfiltered_data),32-length(bad_channels));
p=0;

for q=1:32,
    flag=0;
    for r=1:length(bad_channels),
        if q == bad_channels(r),
            flag=1;
        end
    end
    if flag ==0;
        p=p+1;
        data(:,p)=unfiltered_data(:,q);
    end
end
    unfiltered_data=data;
    clear data
    
%  design the first band - pass       
            L=length(unfiltered_data);
            T=1000/sampling_rate;
            fc2=13
            fn2=2*fc2/(sampling_rate);
            fn=[fn2];              
            [b,a]=butter(5,fn);
            LFP_data=filtfilt(b,a,unfiltered_data);
            fc1=2				% corner frequency of the filter (in kHz)
            fn1=2*fc1/(sampling_rate);		% normalize the corner frequence with fsample
            fn=[fn1];               % band pass filter, so it has two corner frequencies
            [b,a]=butter(3,fn,'high');
            LFP_data1=filtfilt(b,a,LFP_data);
        
% plot the LFP_data
            t=1/sampling_rate:1/sampling_rate:length(LFP_data1)/sampling_rate;
            plotsize=4;
            [M,N]=size(LFP_data1);
            numfigs=floor(N/plotsize);
            for q=1:numfigs,
                figure
                for p=1:plotsize,
                    subplot(plotsize,1,p)
                    plot(t,LFP_data1(:,p+plotsize*(q-1)))
                    amp=3*abs(std((LFP_data1(:,p+plotsize*(q-1)))));
                    axis([0 5 -amp amp]);
                    temp=(q-1)*plotsize+p;
                    head=num2str(temp,'%d');
                    ylabel(head)
                    title('LFP1 data')
                end
            end

%  design the second band - pass       
            fc2=30
            fn2=2*fc2/(sampling_rate);
            fn=[fn2];              
            [b,a]=butter(5,fn);
            LFP_data=filtfilt(b,a,unfiltered_data);
            fc1=13				% corner frequency of the filter (in kHz)
            fn1=2*fc1/(sampling_rate);		% normalize the corner frequence with fsample
            fn=[fn1];               % band pass filter, so it has two corner frequencies
            [b,a]=butter(5,fn,'high');
            LFP_data2=filtfilt(b,a,LFP_data);
        
% plot the LFP_data
            t=1/sampling_rate:1/sampling_rate:length(LFP_data2)/sampling_rate;
            plotsize=4;
            [M,N]=size(LFP_data2);
            numfigs=floor(N/plotsize);
            for q=1:numfigs,
                figure
                for p=1:plotsize,
                    subplot(plotsize,1,p)
                    plot(t,LFP_data2(:,p+plotsize*(q-1)))
                    amp=3*abs(std((LFP_data2(:,p+plotsize*(q-1)))));
                    axis([0 5 -amp amp]);
                    temp=(q-1)*plotsize+p;
                    head=num2str(temp,'%d');
                    ylabel(head)
                    title('LFP2 data')
                end
            end
%  design the third band - pass       
            fc2=100
            fn2=2*fc2/(sampling_rate);
            fn=[fn2];              
            [b,a]=butter(5,fn);
            LFP_data=filtfilt(b,a,unfiltered_data);
            fc1=200				% corner frequency of the filter (in kHz)
            fn1=2*fc1/(sampling_rate);		% normalize the corner frequence with fsample
            fn=[fn1];               % band pass filter, so it has two corner frequencies
            [b,a]=butter(5,fn,'high');
            LFP_data3=filtfilt(b,a,LFP_data);
            figure;plot(LFP_data3(:,1:32));
            axis([3600 4200 -1 1])
% plot the LFP_data
            t=1/sampling_rate:1/sampling_rate:length(LFP_data3)/sampling_rate;
            plotsize=4;
            [M,N]=size(LFP_data3);
            numfigs=floor(N/plotsize);
            for q=1:numfigs,
                figure
                for p=1:plotsize,
                    subplot(plotsize,1,p)
                    plot(t,LFP_data3(:,p+plotsize*(q-1)))
                    amp=3*abs(std((LFP_data3(:,p+plotsize*(q-1)))));
                    axis([0 5 -amp amp]);
                    temp=(q-1)*plotsize+p;
                    head=num2str(temp,'%d');
                    ylabel(head)
                    title('LFP3 data')
                end
            end
%  design the fourth band - pass
            fc2=300
            fn2=2*fc2/(sampling_rate);
            fn=[fn2];              
            [b,a]=butter(5,fn);
            LFP_data=filtfilt(b,a,unfiltered_data);
            fc1=90				% corner frequency of the filter (in kHz)
            fn1=2*fc1/(sampling_rate);		% normalize the corner frequence with fsample
            fn=[fn1];               % band pass filter, so it has two corner frequencies
            [b,a]=butter(5,fn,'high');
            LFP_data4=filtfilt(b,a,LFP_data);
% plot the LFP_data
            t=1/sampling_rate:1/sampling_rate:length(LFP_data4)/sampling_rate;
            plotsize=4;
            [M,N]=size(LFP_data4);
            numfigs=floor(N/plotsize);
            for q=1:numfigs,
                figure
                for p=1:plotsize,
                    subplot(plotsize,1,p)
                    plot(t,LFP_data4(:,p+plotsize*(q-1)))
                    amp=3*abs(std((LFP_data4(:,p+plotsize*(q-1)))));
                    axis([0 5 -amp amp]);
                    temp=(q-1)*plotsize+p;
                    head=num2str(temp,'%d');
                    ylabel(head)
                    title('LFP4 data')
                end
            end
            
%  design the fifth band - pass
            fc2=1000
            fn2=2*fc2/(sampling_rate);
            fn=[fn2];              
            [b,a]=butter(5,fn);
            LFP_data=filtfilt(b,a,unfiltered_data);
            fc1=300				% corner frequency of the filter (in kHz)
            fn1=2*fc1/(sampling_rate);		% normalize the corner frequence with fsample
            fn=[fn1];               % band pass filter, so it has two corner frequencies
            [b,a]=butter(5,fn,'high');
            LFP_data5=filtfilt(b,a,LFP_data);
% plot the LFP_data
            t=1/sampling_rate:1/sampling_rate:length(LFP_data4)/sampling_rate;
            plotsize=4;
            [M,N]=size(LFP_data5);
            numfigs=floor(N/plotsize);
            for q=1:numfigs,
                figure
                for p=1:plotsize,
                    subplot(plotsize,1,p)
                    plot(t,LFP_data5(:,p+plotsize*(q-1)))
                    amp=3*abs(std((LFP_data5(:,p+plotsize*(q-1)))));
                    axis([0 5 -amp amp]);
                    temp=(q-1)*plotsize+p;
                    head=num2str(temp,'%d');
                    ylabel(head)
                    title('LFP5 data')
                end
            end
            
    end
end
        
        
        
        
% 10) Remove bad channels
%
% This function plots the envelope of the signal versus either a time
% scale, sample scale, or video frame scale. Specify x_data according to
% which scale is appropriate (time, sample, or video). Also set the
% appropriate minimum and maximum x and y values to plot.


for x = 1:length(functions_to_run)
    if functions_to_run(x) == 10

    data=zeros(length(LFP_data),32-length(bad_channels));
p=0;
% remove the bad channels
for q=1:32,
    flag=0;
    for r=1:length(bad_channels),
        if q == bad_channels(r),
            flag=1;
        end
    end
    if flag ==0;
        p=p+1;
        data(:,p)=LFP_data(:,q);
    end
end
    LFP_data=data;
    clear data

    end
end


% 20) calculate and plot Principal Components of LFP data

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 20
            [coef, score, latent]=princomp(unfiltered_data);
            latent=100*latent/sum(latent);
            figure; pareto(latent);
            
            plotsize=4;
            [M,N]=size(score);
            numfigs=floor(N/plotsize);
            t=1/sampling_rate:1/sampling_rate:length(score)/sampling_rate;
                for q=1:numfigs,
                figure
                    for p=1:plotsize,
                        subplot(plotsize,1,p)
                        plot(t,score(:,p+plotsize*(q-1)))
                        amp=5*abs(std((score(:,p+plotsize*(q-1)))));
                        axis([0 5 -amp amp]);
                        temp=(q-1)*plotsize+p;
                        head=num2str(temp,'%d');
                        ylabel(head)
                        title('Principle Components of population data')
                    end
                end
     end
end

% 30) Calculate and Plot Independent Components
%
% This function calculates the IC and their correlation coefficients between all pairs
% of channels for a specified time and plots them according to
% the position of each electrode within the array using the specified color
% map. (Typical colormaps are 'jet', 'cool', 'gray', and 'bone'.)


for x = 1:length(functions_to_run)
    if functions_to_run(x) == 30

        ica_data=1000*LFP_data';
        [ica_data,A, W] = fastica(ica_data,'numOFIC',24);
        ica_data=ica_data';
            
        plotsize=4;
        [M,N]=size(ica_data);
        numfigs=floor(N/plotsize);
            t=1/sampling_rate:1/sampling_rate:length(ica_data)/sampling_rate;
                for q=1:numfigs,
                figure
                    for p=1:plotsize,
                    subplot(plotsize,1,p)
                    plot(t,ica_data(:,p+plotsize*(q-1)))
                    amp=5*abs(std((score(:,p+plotsize*(q-1)))));
                    axis([0 5 -amp amp]);
                    end
                end
    end
end

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 40

            %temp = unfiltered_data;
%             temp=DATA(:,1:32);
            temp=unfiltered_data(:,1:32);
           
            N=2048;  
            L=32   % window length for coherence calculation
            [M,L]=size(temp);
% for  q=1:L,
%     temp(:,q)=temp(:,q) - mean(temp(:,q));
% end
                  % pos_x = [3  5  5  2  4  5  1  3  2   4   6   3    1   6   4    2   2   4    6   1   3   6   4    2   3    1   5   4   2   5   5   3];
                   %pos_y = [1  1  2  2  2  3  3  3  4   4   4   5    5   5   6    6   1   1    2   2   2   3   3    3   4    4   4   5   5   5   6   6];
    pos_x=[8 6 4 2 1 3 5 7 7 5 3 1 2 4 6 8 8 6 4 2 1 3 5 7 7 5 3 1 2 4 6 8]
    pos_y=[3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 4 4 4 4 4 4 4 4 1 1 1 1 1 1 1 1]
    
% channel no -        1   2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
%                    pos_x = [2 10 4 5 6 5 1 2 3 4 6 5 4 1 2 3 3 2 1 4 5 6 4 3 2 1 5 6 5 4 3 2];
%                    pos_y = [1 20 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6];
                % The (0,0) corner for the (x,y) coordinate system is the
                % bottom left corner
                T=0;
                COHERE=mscohere((temp(:,1)), (temp(:,1)),hanning(N) ,N/2 ,N ,sampling_rate);
                COHERE=zeros(length(COHERE),1);
                DIST=input('Enter the distance value= ')
                DIRECT=input('Enter the 1 for vertical only, enter 2 for horizantal only= ')
                L=32
                for q=1:L-1,
                    q
                    for r=q+1:L,
                        r
                    X1=pos_x(q);
                    Y1=pos_y(q);
                    X2=pos_x(r);
                    Y2=pos_y(r);
                    X=X2-X1;
                    Y=Y2-Y1;
                    R=round(sqrt(X*X+Y*Y));
                    % distance between contacts that are going to be used for this analysis (averaging)
                    if (DIRECT == 1 & X == 0) | (DIRECT == 2 & Y == 0),
                        if R == DIST,                  
                            T=T+1
                        [Cxy F]= mscohere(temp(:,q), temp(:,r), hanning(N) ,N/2 , N ,sampling_rate);
                        figure(3)
                        plot(F, Cxy)
                        COHERE=COHERE+Cxy;
                        end
                    end
                    end
                end
                 COHERE=COHERE/T;
                 [Cxy F]= mscohere(temp(:,q), temp(:,r),hanning(N) ,N/2,N ,sampling_rate);
                 hold on;figure(2)
                 plot(F,COHERE,'g')
                 xlabel('Frequency (kHz')
                 ylabel('Agregate Coherence')
    end
end


% 40) calculate and plot agregate coherence within a distance of R

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 41

            %temp = unfiltered_data;
            temp=unfiltered_data;
            N=4096  

            [M,L]=size(temp);
            L=32   % max limit of channels
% for  q=1:L,
%     temp(:,q)=temp(:,q) - mean(temp(:,q));
% end

NormBand1L =round( fc1*N/sampling_rate);
NormBand1H =round( fc2*N/sampling_rate);
NormBand2L =round( fc3*N/sampling_rate);
NormBand2H =round( fc4*N/sampling_rate);



for i=1:1
%     h(i) = figure(i);

    for q=1:32
        
[Cxx,F] = mscohere(LFP_data(:,i),LFP_data(:,q),hanning(N) ,N/2 ,N ,sampling_rate);
[Cxx2,F2] = mscohere(LFP_data2(:,i),LFP_data2(:,q),hanning(N) ,N/2 ,N ,sampling_rate);
Cxy = mean(Cxx(1:NormBand1H));
% figure(i);
% 
% subplot(4,8,positions(q));
% imagesc(Cxy,[0.1 0.5]); colormap(gray);axis off
% head=num2str(q,'%d');
% title(head)
figure(i)
subplot(4,8,positions(q));
plot(F(1:NormBand1H),Cxx(1:NormBand1H))
xlabel('Frequency (Hz')
ylabel('Filtered Coherence')
 head=num2str(q,'%d');
title(head)
% axis([0 1000 0 0.5])
figure(i+1)
subplot(4,8,positions(q));
plot(F2(1:NormBand2H),Cxx2(1:NormBand2H))
xlabel('Frequency (Hz')
ylabel('Filtered Coherence')
 head=num2str(q,'%d');
title(head)
% axis([0 1000 0 0.5])

    end

end
%                    pos_x = [3  5  5  2  4  5  1  3  2   4   6   3    1   6   4    2   2   4    6   1   3   6   4    2   3    1   5   4   2   5   5   3];    % horizontal from left to right
%                    pos_y = [1  1  2  2  2  3  3  3  4   4   4   5    5   5   6    6   1   1    2   2   2   3   3    3   4    4   4   5   5   5   6   6];    % vertical from bottom up
%    
%                 T=0;
%                 COHERE=mscohere(detrend(temp(:,1)), detrend(temp(:,1)),hanning(N) ,N/2 ,N ,sampling_rate);
%                 COHERE=zeros(length(COHERE),1);
%                 DIST=input('Enter the distance value= ')
%                 
%                                 figure(101); clf(101); title('cohere with below'); subplot(6,6,1); 
%                                 figure(104); clf(104); title('cohere with right'); subplot(661)       
% 
%                 for q=1:L,
%                     q
%                     for r=1:L,
%                         r
%                     X1=pos_x(q);
%                     Y1=pos_y(q);
%                     X2=pos_x(r);
%                     Y2=pos_y(r);
%                     X=X2-X1;
%                     Y=Y2-Y1;
%                     R=sqrt(X*X+Y*Y);
%                     % distance between contacts that are going to be used for this analysis (averaging)
%                         if R == DIST,                  
%                             T=T+1
%                         [Cxy F]= mscohere(temp(:,q), temp(:,r), hanning(N) ,N/2 , N ,sampling_rate);
% 
%                                 if pos_x(q)==pos_x(r) & pos_y(q) > pos_y(r),   figure(101); subplot(6,6,positions(q)); end % below
%                                 if pos_y(q)==pos_y(r) & pos_x(q) < pos_x(r),   figure(104); subplot(6,6,positions(q)); end % on the right
% 
%                         plot(F,Cxy)
%                         COHERE=COHERE+Cxy;
%                         end
%                     end
%                 end
%                  COHERE=COHERE/T;
%                  [Cxy F]= mscohere(temp(:,q), temp(:,r),hanning(N) ,N/2,N ,sampling_rate);
%                  figure(100)
%                  plot(F,COHERE,'b')
%                  xlabel('Frequency (kHz')
%                  ylabel('Agregate Coherence')
    end
end

% 50) calculate and plot agregate coherence within a distance of R

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 50

            temp = unfiltered_data;
            [M,L]=size(temp);
                % pos_x=[6 6 5 5 5 4 4 4 3 3 3 2 2 2 1 1];
    %pos_y=[1.5 2.5 1 2 3 1 2 3 1 2 3 1 2 3 1.5 2.5];
                   pos_x = [3  5  5  2  4  5  1  3  2   4   6   3    1   6   4    2   2   4    6   1   3   6   4    2   3    1   5   4   2   5   5   3];
                   pos_y = [1  1  2  2  2  3  3  3  4   4   4   5    5   5   6    6   1   1    2   2   2   3   3    3   4    4   4   5   5   5   6   6];
                % pos_x = [2 3 4 5 6 5 1 2 3 4 6 5 4 1 2 3 3 2 1 4 5 6 4 3 2 1 5 6 5 4 3 2];    
                % pos_y = [1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6];
                % The (0,0) corner for the (x,y) coordinate system is the bottom left corner
                
                N=128;               % window length for coherence calculation
                D=0.02;                  % running time window in seconds
                D=D*sampling_rate;
                O=0.01;                 % time iteration in seconds
                O=O*sampling_rate;
                DIST=input('Enter the distance value= ')
                COHEREE=[];
                for u=1:floor((M-D)/O),
                    u
                COHERE=mscohere(detrend(temp((u-1)*O+1:(u-1)*O+D,1)), detrend(temp((u-1)*O+1:(u-1)*O+D,1)), [] ,N/2 ,N ,sampling_rate);
                COHERE=zeros(length(COHERE),1);
                T=0;
                
                for q=16,              % this is the contact
                    q
                    for r=32,
                    X1=pos_x(q);
                    Y1=pos_y(q);
                    X2=pos_x(r);
                    Y2=pos_y(r);
                    X=X2-X1;
                    Y=Y2-Y1;
                    R=round(sqrt(X*X+Y*Y));
                    % distance between contacts that are going to be used for this analysis (averaging)
                        if R == DIST,                  
                            T=T+1;
                            Cxy=mscohere(detrend(temp((u-1)*O+1:(u-1)*O+D,q)), detrend(temp((u-1)*O+1:(u-1)*O+D,r)), hanning(N) ,N/2 ,N ,sampling_rate);
                            % Cxy=mscohere(detrend(temp(:,q)), detrend(temp(:,r)), hanning(N) ,N/2 ,N ,sampling_rate);
                            COHERE=COHERE+Cxy;
                        end
                    end
                end
                 COHERE=COHERE/T;
                 COHEREE=[COHEREE COHERE];
                end
                 %[Cxy F]= mscohere(temp(:,1), temp(:,1),[] ,N/2,N ,sampling_rate);
                 figure
                 map=zeros(64,3);
                for graylevel=0:63; map(graylevel+1,:)=1-graylevel/63;end

                colormap(map)
                 surface(COHEREE)
                 xlabel('Time')
                 ylabel('frequency')
                 
    end
end

% 60) agregate correlation between contacts of distance of R

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 60

            temp = filtered_data;
            [M,L]=size(temp);
pos_y= [6 6 5 5 5 4 4 4 3 3 3 2 2 2 1 1 6 6 5 5 5 4 4 4 3 3 3 2 2 2 1 1];
pos_x=[3 5 5 2 4 5 1 3 2 4 6 3 1 6 4 2 2 4 6 1 3 6 4 2 3 1 5 4 2 5 5 3];
                % The (0,0) corner for the (x,y) coordinate system is the bottom left corner
                
                D=0.1;                  % running time window in seconds
                D=D*sampling_rate;
                O=.01;                 % time iteration in seconds
                time=O:O:(M-D)/sampling_rate;
                O=O*sampling_rate;
                DIST=input('Enter the distance value= ')
                CORR=zeros(1,floor((M-D)/O));
                CORRR=zeros((floor(M-D)/O),L);
                figure(103)
                MAXX=zeros(1,number_channels);
               for q=1:number_channels,              % this is the contact
                    q
                   COHERE=zeros(floor((M-D)/O),2);
                    for r=1:number_channels,
                                X1=pos_x(q);
                                Y1=pos_y(q);
                                X2=pos_x(r);
                                Y2=pos_y(r);
                                X=X2-X1;
                                Y=Y2-Y1;
                                R=round(sqrt(X*X+Y*Y));
                                % distance between contacts that are going to be used for this analysis (averaging)

                                 if R <= DIST & R ~= 0,                  
                                     for u=1:floor((M-D)/O),
                                        Cxy=corrcoef(temp((u-1)*O+1:(u-1)*O+D,q), temp((u-1)*O+1:(u-1)*O+D,r));
                                        COHERE(u,2)=Cxy(1,2);
                                     end
                                     [MAX, I]=max(sum(COHERE));
                                     if I==2,
                                         MAXX(q)=r;
                                        COHERE(:,1)=COHERE(:,I);
                                        CORR=COHERE(:,I);
                                    end
                                end
                                

                    end
                    subplot(6,6,positions(q));
                    plot(time,CORR);
                    CORRR(:,q)=CORR;
                    axis([0 5 -1 1]);
%             
%                  xlabel('Time')
%                  ylabel('Ave Corr')
                 
                end

    end
end


% 100) Show movie

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 100
        
        filename=input('Enter the video file name= ','s')
        cd video
        mov=aviread(filename);
        cd ..
             for r=1:length(mov),
                R=[1 r];
                figure(100)
                movie(mov, R)
                r/30
                pause(0.4)
             end
    end
end

% 110) Show movie and a running marker on a channel of LFP potentials

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 110
        
        filename=input('Enter the video file name= ','s')
        cd video
        mov=aviread(filename);
        cd ..
        FigNo=101;
        ChanNo=21;
        O=1;
        DATAA=LFP_ave;             % modify this variable name as desired
       % time=O/sampling_rate:O/sampling_rate:5.005-D/sampling_rate;
        time=1/10000:1/10000:5;
        figure(FigNo)

        hold off
        plot(time, DATAA(:,ChanNo));
        amp2=max(DATAA(:,ChanNo));
        amp1=min(DATAA(:,ChanNo));
        axis([0 5.005 amp1 amp2]);
        figure(100)
            for r=5:150,
            R=[1 r];
            figure(100)
            movie(mov, R)
            r/30

            figure(FigNo)

            hold on
            temp=(r-1)/30;
            plot(temp,DATAA(floor(temp*sampling_rate/O)+1,ChanNo),'+w')
            plot(temp,DATAA(floor(temp*sampling_rate/O)+1,ChanNo),'.b')
            temp=(r)/30;
            plot(temp,DATAA(floor(temp*sampling_rate/O)+1,ChanNo),'+r')
            hold off
            pause(1)
            end

    end
end



%120) % plot the correlations amongst the population activity
for x = 1:length(functions_to_run)
    if functions_to_run(x) == 120
        DATAA=LFP_data;
%         DATAA=envelope_data;
map=zeros(64,3);
for graylevel=1:32; map(graylevel,2)=graylevel/40;end
for graylevel=33:64; map(graylevel,1)=(64-graylevel)/40;end

colormap(map)
M=[1 9 1 9 ];
C=[-1 1];
AZ=0;
EL=90;
view(AZ,EL);
                    
%     pos_x=[6 6 5 5 5 4 4 4 3 3 3 2 2 2 1 1 ];
%     pos_y=[1.5 2.5 1 2 3 1 2 3 1 2 3 1 2 3 1.5 2.5];
%     
r = [2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1,4,4,4,4,4,4,4,4]
c = [8,6,4,2,1,3,5,7,7,5,3,1,2,4,6,8,8,6,4,2,1,3,5,7,7,5,3,1,2,4,6,8]
     
% positions = [23 19 18 22 24 20 17 21 7 3 2 6 8 4 1 5 14 10 11 15 13 9 12 16 30 26 27 31 29 25 28 32];
positions= [16 14 12 10 9 11 13 15 23 21 19 17 18 20 22 24 8 6 4 2 1 3 5 7 31 29 27 25 26 28 30 32];
sites = [5 4 6 3 7 2 8 1 13 12 14 11 15 10 16 9 20 21 19 22 18 23 17 24 28 29 27 30 26 31 25 32]
s= corrcoef(LFP_data(:,1:32));
figure;
xx=zeros(32,5,9);
for k=1:32
for i=1:32
xx(k,r(i),c(i)) = s(k,i);
end
% b = find(positions==k)
subplot(4,8,positions(k))
head=num2str(sites(positions(k)),'%d');
ylabel(head)
h = surface(squeeze(xx(k,:,:))); caxis([-1 1]);
axis([1 9 1 5]);
colormap(map);

% pause;
end

    
% figure
% correlations=zeros(8,4);
% for k=1:32,
%     for m=1:32,
%     COR=corrcoef(DATAA(:,k),DATAA(:,m));
%     correlations(r(m),c(m))=COR(1,2);
%     end
% %     correlations(7,:)=correlations(6,:);
% %     correlations(:,7)=correlations(:,6);
%     subplot(8,4,positions(k));
%     surface(correlations);
%     axis(M)
%     caxis(C)
%     colormap(map)
%     view(AZ,EL);
% end
%     
    end
end

% 130) plot the correlations amongst the local field potentials

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 130
        
map=zeros(64,3);
for graylevel=1:32; map(graylevel,2)=graylevel/40;end
for graylevel=33:64; map(graylevel,1)=(64-graylevel)/40;end

colormap(map)
M=[1 7 1 7 ];
C=[-1 1];
AZ=0;
EL=90;
view(AZ,EL);
%                16      32    15       31
%
%        13    29     12     28      30     14
%
%        26      9     25     10      27     11
%
%         7      24       8     23       6      22
%
%        20       4     21       5       3     19
%
%                  17      1     18       2 
                                   pos_x=[6 6 5 5 5 4 4 4 3 3 3 2 2 2 1 1];
    pos_y=[1.5 2.5 1 2 3 1 2 3 1 2 3 1 2 3 1.5 2.5];
%pos_y=fliplr(pos_y);
DATA= LFP_data;
figure
subplot(6,6,2)
correlations=zeros(7,7);
for k=1:14,
    for m=1:14,
    COR=corrcoef(DATA(:,k),DATA(:,m));
    correlations(round(pos_y(m)),pos_x(m))=COR(1,2);
    end
    correlations(7,:)=correlations(6,:);
    correlations(:,7)=correlations(:,6);
    subplot(6,6,positions(k));
    surface(correlations);
    axis(M)
    caxis(C)
    colormap(map)
    view(AZ,EL);
end

    end
end
% 140)  Movie of 2D correlations as a function of time

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 140
        
map=zeros(64,3);
for graylevel=1:32; map(graylevel,2)=graylevel/40;end
for graylevel=33:64; map(graylevel,1)=(64-graylevel)/40;end

colormap(map)
M=[1 7 1 7 ];
C=[-1 1];
AZ=0;
EL=90;
view(AZ,EL);
%                16      32    15       31
%
%        13    29     12     28      30     14
%
%        26      9     25     10      27     11
%
%         7      24       8     23       6      22
%
%        20       4     21       5       3     19
%
%                  17      1     18       2 
                     pos_x=[6 6 5 5 5 4 4 4 3 3 3 2 2 2 1 1];
    pos_y=[1.5 2.5 1 2 3 1 2 3 1 2 3 1 2 3 1.5 2.5];
%pos_y=fliplr(pos_y);
DATA=envelope_data;
figure(102)
subplot(3,6,2)
correlations=zeros(4,7);
D=0.2;     % time step
D=D*sampling_rate;

for t=D:D:5.005*sampling_rate-D,
    t/sampling_rate
    clf
    colormap(map)
%     temp=t/sampling_rate;
%     head=num2str(temp,'%d');
%     title(head)
    for k=1:1,  %number_channels,
        for m=1:number_channels,
        COR=corrcoef(DATA(t+1:t+D,k),DATA(t+1:t+D,m));
        correlations(round(pos_y(m)),pos_x(m))=COR(1,2);
        end
        correlations(4,:)=correlations(3,:);
        correlations(:,7)=correlations(:,6);
        figure(102)
        %subplot(3,6,positions(k));
        surface(correlations);
        axis(M)
        caxis(C)
        view(AZ,EL);
        pause(1)
    end
end

    end
end

% 150) spectrogram for the unfiltered data

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 150
        figure; 
        for i=1:number_channels
     ChanNo=i;   
     N=1024;
      [y,f,t,p]= spectrogram(unfiltered_data(:,ChanNo),hanning(N),N/2 ,N,sampling_rate);
         subplot(4,8,positions_31channels(i));
                    
   P=(109 +10*log10(abs(p)))/109;
   %P=10*log10(abs(p))
     surf(t,f,P,'EdgeColor','none');   
      axis xy; axis([0 1 50 500]); %caxis([-1e-15 1e-9]); 
      colormap(hsv); view(0,90);
      xlabel('Time');
      ylabel('Frequency (Hz)');
        end
      end
end

% 160) running average of the zero crossing time of the LFPs shown on the
% array
for x = 1:length(functions_to_run)
    if functions_to_run(x) == 160
        
        DATAA=LFP_data;
        D=0.1;      % running window length
        D=D*sampling_rate;
        t=1/sampling_rate:1/sampling_rate:(length(DATAA)-1-D)/sampling_rate;
        figure(102)
        for k=1:number_channels ,
            k
            temp=(abs(sign(DATAA(1:length(DATAA)-1,k))-sign(DATAA(2:length(DATAA),k))));
            temp_ave=zeros(1,length(temp)-D);

            for T=1:length(temp)-D,
                temp_ave(T)=sum(temp(T:T+D));
            end
            subplot(6,6,positions(k));
            plot(t, temp_ave)
            xlabel('time')
        end

    end
end
% 170) calculate force
for x = 1:length(functions_to_run)
    if functions_to_run(x) == 170
        
FORCE=unfiltered_data(:,70:73) - [ unfiltered_data(:,74)  unfiltered_data(:,74)  unfiltered_data(:,74)  unfiltered_data(:,74)];
figure; plot(time, FORCE)

    end
end

%% 180) Evoked Response Analysis
% close all
if(signal_length > 32),
signal_length = signal_length/2;
end


for x = 1:length(functions_to_run)

    if functions_to_run(x) == 180
        figure;
        DATA2=LFP_data(:,1:signal_length);
        
% for j=1:31
% figure;
% plot(LFP_data2(:,j));axis([1 2000 -0.5e-4 2e-4])        % Monitor the data first channel by channel
% pause;
% end
    for z=1:signal_length;
        
        %% 1.method Graphical 
        
%     clf;plot(DATA(:,z));axis([1 2000 -0.5e-4 2e-4])
%     hold on
%     COOR=ginput(2);
%     X1(z,:)=round(COOR(1,1));
%     X2(z,:)=round(COOR(2,1));
%   
%     Y1=DATA(X1(z,:),z);
%     plot(X1(z,:),Y1,'r*')
%     Y2=DATA(X2(z,:),z);
%     plot(X2(z,:),Y2,'r*')
%   
%     grid on
%     
%     
%     AMP1(z)=Y1-Y2;
%     
%     pause;
    %% 2. method Numerical
    
        for k=1:1;
           
    k = k*10000;
 
    

% Detect the Spike
%     Detect_Start1N = (k+100) - 10000;
%     Detect_End1N = (k+400) - 10000;
% 
%     Detect_Start1P  = (k+100) - 10000;
%     Detect_End1P   = (k+400) - 10000;
% 
%     Detect_Start2N = (k+5800) - 10000;
%     Detect_End2N   = (k+5830) - 10000;
%     
%     Detect_Start2P = (k+5900) - 10000;
%     Detect_End2P   = (k+5970) - 10000;
    
    k = k/10000;
   
    Peak1_Detect_High(z,k) = (max(DATA2(EPs1 - marg:EPs1 + marg,z)));
    Peak1_Detect_Low(z,k) = (min(DATA2(EPs1 - marg:EPs1 + marg,z)));
    
    Peak2_Detect_High(z,k) = (max(DATA2(EPs2 - 2*marg:EPs2 + 2*marg,z)));
    Peak2_Detect_Low(z,k) = (min(DATA2(EPs2 - 2*marg:EPs2+ 2*marg,z)));
    
    
    
    Resp2(z,k) = Peak2_Detect_High(z,k) - Peak2_Detect_Low(z,k);

   
   [rh(z),ch(z)]= find(DATA2 == Peak1_Detect_High(z,k));
    [rl(z),cl(z)]= find(DATA2 == Peak1_Detect_Low(z,k));
     
    
    if (rh(z) > rl(z))
            Resp1(z,k) = Peak1_Detect_High(z,k) - Peak1_Detect_Low(z,k);
    else
    
        Resp1(z,k) = Peak1_Detect_Low(z,k) - Peak1_Detect_High(z,k);
    end
        
        
    [rh2(z),ch2(z)]= find(DATA2 == Peak2_Detect_High(z,k));
    [rl2(z),cl2(z)]= find(DATA2 == Peak2_Detect_Low(z,k));
    
    
    plot(DATA2(:,z));hold on
    
    
     plot(rh(z),Peak1_Detect_High(z,k),'r*')
    plot(rl(z),Peak1_Detect_Low(z,k),'k*')
    plot(rh2(z),Peak2_Detect_High(z,k),'r*')
    plot(rl2(z),Peak2_Detect_Low(z,k),'k*')

        end  
        
end
% 
% figure;plot(Resp1);
% figure;plot(Resp2);
% figure;plot(AMP1)
figure;
aligned_Resp1 = zeros(8,4);
for p = 1:31
subplot(4,8,positions(p));
imagesc(Resp1(p)'); colormap(gray);caxis([min(Resp1) max(Resp1)]);
aligned_Resp1(p) = Resp1(positions_sorted(p));
head=num2str(p,'%d');
timing = num2str(rh(p),'%d');
ylabel(head);
xlabel(timing)
end


figure;
aligned_Resp2 = zeros(8,4);
for p = 1:signal_length
subplot(4,8,positions(p));
imagesc(Resp2(p)'); colormap(gray);caxis([min(Resp2) max(Resp2)]);
aligned_Resp2(p) = Resp2(positions_sorted(p));
head=num2str(p,'%d');
ylabel(head);
timing2 = num2str(rh2(p),'%d');
xlabel(timing2)
end

figure;h=bar3(aligned_Resp1'*1e6);colormap gray
figure;h2=bar3(aligned_Resp2'*1e6);colormap gray

%% Set the Color Shading respect to Z-axis

for i = 1:length(h)
    zdata = get(h(i),'ZData');
    set(h(i),'CData',zdata)
    % Add back edge color removed by interpolating shading
    set(h,'EdgeColor','k') 
       zdata2 = get(h2(i),'ZData');
    set(h2(i),'CData',zdata2)
    % Add back edge color removed by interpolating shading
    set(h2,'EdgeColor','k') 
end


    end
end
toc

%% Simple way for EP analysis

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 180.1

        First = EPs1;
        End = EPs2;
        
Min_LFP = min(LFP_data_sorted(First:End,1:number_channels))';
Max_LFP = max(LFP_data_sorted(First:End,1:number_channels))';
Evoked_diff= Max_LFP -Min_LFP;
% figure;bar(Min_LFP)
% figure;bar(Max_LFP)
% figure;bar(Evoked_diff)

aligned_Resp1 = zeros(8,4);
for p = 1:number_channels
    aligned_Resp1(p) = Evoked_diff(p);
end

figure;h=bar3(aligned_Resp1'*1e6);colormap hot

%% Set the Color Shading respect to Z-axis

for i = 1:length(h)
    zdata = get(h(i),'ZData');
    set(h(i),'CData',zdata)
    % Add back edge color removed by interpolating shading
    set(h,'EdgeColor','k') 

end
% 
% First = EPs3;
% End = EPs4;
% 
% 
% Min_LFP2 = min(LFP_data_sorted(First:End,32:62))';
% Max_LFP2 = max(LFP_data_sorted(First:End,32:62))';
% Evoked_diff2= Max_LFP2 -Min_LFP2;
% % figure;bar(Min_LFP)
% % figure;bar(Max_LFP)
% % figure;bar(Evoked_diff)
% 
% aligned_Resp2 = zeros(8,4);
% for p = 1:31
%     aligned_Resp2(p) = Evoked_diff2(p);
% end
% 
% figure;h2=bar3(aligned_Resp2'*1e6);colormap hot
% figure;subplot(1,2,1);contourf(aligned_Resp2);colormap(jet)
% subplot(1,2,2);contourf(aligned_Resp1);colormap(jet)
% 
% 
% %% Set the Color Shading respect to Z-axis
% 
% for i = 1:length(h2)
%     zdata = get(h2(i),'ZData');
%     set(h2(i),'CData',zdata)
%     % Add back edge color removed by interpolating shading
%     set(h2,'EdgeColor','k') 
% 
% end



    end
end




for x = 1:length(functions_to_run)
    if functions_to_run(x) == 182
        
        MF_baseline  = min(min(rl));
        MF_delays = rl - MF_baseline;
        MF_delays = MF_delays/sampling_rate*1000;
        
        CF_baseline  = min(min(rh2));
        CF_delays = rh2 - CF_baseline;
        CF_delays = CF_delays/sampling_rate*1000;
        
            figure;
        for l = 1:min(size(pos_x))
            subplot(4,1,l)
            bar(MF_delays(:,pos_x(l,:)));
            ylabel('ms','fontsize',12,'fontweight','b')
            title('MF Response Delays')

        end
        
        figure;
        for l = 1:length(pos_y)
            subplot(1,8,l) 
            bar(MF_delays(:,pos_y(l,:)));
            grid on
            axis tight
            title('MF Delays')

        end
        
            figure;
        for l = 1:min(size(pos_x))
            subplot(4,1,l)
            bar(CF_delays(:,pos_x(l,:)));
            ylabel('ms','fontsize',12,'fontweight','b')
            title('CF Response Delays')

        end
        
        figure;
        for l = 1:length(pos_y)
            subplot(1,8,l) 
            bar(CF_delays(:,pos_y(l,:)));
            grid on
            axis tight
            title('CF Delays')
        end
        
    end
end

%% 183) Detection of ONSET latencies

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 183
        
        figure;plot(LFP_data_sorted(:,1:number_channels)) % Ensure the response is exist
        desired_ch = 3;
        
        
        for i=1:number_channels
%         figure(i);plot(LFP_data(:,i));
        ymax =max(LFP_data(:,i)); 
        ymin = min(LFP_data(:,i));
        axis([5550 6100 ymin ymax ])
%         X_s = ginput(2);

        
        min_level1 =5700;
        max_level1 =5770;
        LFP_data_detect = LFP_data(min_level1:max_level1,:);
        
       
        min1 = min(LFP_data_detect(:,i))
        max1 = max(LFP_data_detect(:,i))
        
        
        min_level2 =5800;
        max_level2 =6000;
        
        LFP_data_detect = LFP_data(min_level2:max_level2,:);

         min2 = min(LFP_data_detect(:,i))
         max2 = max(LFP_data_detect(:,i))
         
        detect_max1(i) = find(LFP_data(:,i) == max1)
        detect_min1(i) = find(LFP_data(:,i) == min1)
        
         detect_max2(i) = find(LFP_data(:,i) == max2)
        detect_min2(i) = find(LFP_data(:,i) == min2)
        
    
        figure(i);
        plot(LFP_data(:,i));hold on;
        plot(detect_min1(i),LFP_data(detect_min1(i),i),'r*') 
        plot(detect_max1(i),LFP_data(detect_max1(i),i),'r*')
         plot(detect_min2(i),LFP_data(detect_min2(i),i),'r^') 
        plot(detect_max2(i),LFP_data(detect_max2(i),i),'r^') 
        axis([5500 6000 ymin ymax ])
        pause(1);
        
        a = input('If you didnt like it?,enter "1" to do manually');
        if (a ==1)
             X_s = ginput(4);
             detect_min1(i)  = X_s(1)
             detect_max1(i)  = X_s(2)
             
             detect_min2(i)  = X_s(3)
             detect_max2(i)  = X_s(4)
             
        end
        detect_max1(i) = (detect_max1(i) - puff_arrival)/sampling_rate*1000; %% normalize the delay
        detect_min1(i) = (detect_min1(i) - puff_arrival)/sampling_rate*1000;
        
         detect_max2(i) = (detect_max2(i) - puff_arrival)/sampling_rate*1000; %% normalize the delay
        detect_min2(i) = (detect_min2(i) - puff_arrival)/sampling_rate*1000;

%         L_start(i) = X_s(1);
%         L_start(i) = (L_start(i) - puff_arrival)/sampling_rate*1000;
%      
%         L_peak(i)  = X_s(2);
%         L_peak(i) = (L_peak(i) - puff_arrival)/sampling_rate*1000
%         figure(100*i);bar([L_start(i) L_peak(i)])
        close
        end
 
for p = 1:32    
subplot(4,8,positions(p));

bar([detect_min1(p)  detect_min2(p) ]);
axis([1 2 0 14])
head=num2str(p,'%d');
ylabel(head)
end
    end
end




%% 181)-2 Averaging 

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 181

L=zeros(10000,20);

for c=1:32
c
for  i=1:20
i
L(:,i,c)= LFP_data(((i-1)*10000+1):i*10000,c);
end
end

average = zeros(32,10000);

for i=1:32
average(i,:) = sum(L(1:10000,1:20,i)')/20;
end
average = average';

for i=1:32
odd_average(i,:) = sum(L(1:10000,1:2:20,i)')/10;
end

for i=1:32
even_average(i,:) = sum(L(1:10000,2:2:20,i)')/10;
end

Noise = odd_average' - even_average';
S = odd_average' + even_average';

figure;plot(Noise)
figure;plot(S)
figure;plot(average)


for i=1:30
Resp10_1m(:,i) = Resp10_1(:,i) - mean(Resp10_1(:,i));
end
for i=1:30
Resp10_1mn(:,i)=Resp10_1m(:,i)./std(Resp10_1m(:,i));
end

    end
end





for x = 1:length(functions_to_run)
    if functions_to_run(x) == 190

figure;plot(norm_LFP_data)
pause;





figure;
for i=1:1:500
i=(i)*1;

for p = 1:32
    
subplot(4,8,positions(p));
imagesc(norm_LFP_data(i,p)'); colormap(jet);caxis([-1 1])
yhead_start = num2str(i,'%d');
yhead_end = num2str(i+1,'% d');
title([yhead_start])
head=num2str(p,'%d');
ylabel(head)


end
    
    pause(.1);
    
%      if any(norm_LFP_data2(i,:) >= 0.8 | norm_LFP_data2(i,:) <=-0.8)     % freeze the movie when the amplitude level is above threshold
%         
% 
%         pause;
%     
%     end
    
   i= (i)/1;
end

    end
end





for x = 1:length(functions_to_run)
    if functions_to_run(x) == 200
        
        per = sampling_rate/1000*100;    % Define your temporal window in ms
        figure;
        i=per*3;
 for kkk=1:31
%  subplot(4,8,kkk);
% figure(kkk);
% for i=per*2:per:length(LFP_data)-6*per
for ii=1:31

% clf;
% subplot(1,2,1)
% surface(corrcoef(unfiltered_data(i:i+per,1:number_channels)));caxis([-1 1])
% y_start = num2str(i,'%d');
% title(y_start)
% subplot(1,2,2)
% plot(LFP_data(i:i+per,1:number_channels))
% y_start = num2str(i,'%d');
% title(y_start)


[r,lags]=xcorr(LFP_data_sorted(4800:6500,ii),LFP_data_sorted(4800:6500,kkk),'coeff');

% [r,lags]=xcorr(LFP_data_sorted(i:i+per,ii),LFP_data_sorted(i:i+per,kkk),'coeff');
win_gauss = gausswin(length(r));
az = find(r==max(r));
kk = i/per;
zz(kkk,ii) = lags(az)
win_gauss = gausswin(length(r));
yy = r.*win_gauss;      % Gauss windowing
% 
figure(132);clf;subplot(2,1,1);plot(1e6*LFP_data_sorted(5900:6500,[ii kkk]));
y_start = num2str(kkk,'%d');title(y_start)
subplot(2,1,2); plot(lags,r,'color',cc(4*kk,:));
% pause(.5);

% figure(103);hold on

end
%  end
%  pause;
 %         mean(zz)        % average of lagging or leading across time.
%   figure;plot(zz(kkk,:),'color',cc(4*kk,:))

 end
 figure;stem(zz)

%  for kkk=1:31
%  subplot(4,8,kkk);plot(zz(kkk,:));end
    end
end

        
for x = 1:length(functions_to_run)
    if functions_to_run(x) == 201
%        LFP_data = LFP_data2;
for i=EPs1 - marg:1:EPs1 + marg
    for k=1:32
        
   if (LFP_data(i,k)>0)

       LFP_rectify_exc(i,k)=LFP_data(i,k);
   else
       
       LFP_rectify_inh(i,k) = -LFP_data(i,k);
      
        end
    end
end

for i=EPs2 - marg:1:EPs2 + marg
    for k=1:32
        
   if (LFP_data(i,k)>0)

       LFP_rectify_exc2(i,k)=LFP_data(i,k);
   else
       
       LFP_rectify_inh2(i,k) = -LFP_data(i,k);
      
        end
    end
end

figure;plot(LFP_rectify_exc)
figure;plot(LFP_rectify_exc2)
% figure;plot(LFP_rectify_inh)
% LFP_rectify_exc = LFP_rectify_exc(:,[2:16,19:30]); % Discard the unresponsive channels
LFP_rectify_exc =LFP_rectify_inh;
LFP_rectify_exc2 =LFP_rectify_inh2;

figure;for p = 1:number_channels
subplot(4,8,positions(p));

r= randi(32,32,1);
r=r(1);

maxresp = find(LFP_rectify_exc(:,r) == max(max(LFP_rectify_exc(:,r))));

% maxresp2 = find(LFP_rectify_inh(:,r) == max(max(LFP_rectify_inh(:,r))));

imagesc(LFP_rectify_exc(maxresp,p)'); colormap(gray);caxis([min(min(LFP_rectify_exc(maxresp,1:32))) max(max(LFP_rectify_exc(maxresp,1:32)))])

% imagesc(LFP_rectify_inh(maxresp2,p)'); colormap(jet);caxis([min(min(LFP_rectify_inh(maxresp2,[1:32]))) max(max(LFP_rectify_inh(maxresp2,1:32)))])


yhead_start = num2str(i,'%d');
yhead_end = num2str(i+1,'% d');
title([yhead_start])
head=num2str(p,'%d');
ylabel(head)

end

% figure
% for k=5000:16:7000;
% 
% clf
% for p = 1:32
% subplot(4,8,positions(p));
% imagesc(LFP_rectify_exc(k:k+10,p)'); colormap(jet);caxis([0 max(max(LFP_rectify_exc))])
% yhead_start = num2str(i,'%d');
% yhead_end = num2str(i+1,'% d');
% title([yhead_start])
% head=num2str(p,'%d');
% ylabel(head)
% 
% end
% pause(.1)
% end

    end
    

    end 

% for x = 1:length(functions_to_run)
%     if functions_to_run(x) == 202
%         
%         start_point = max_indice - 50;
%         end_point = max_indice + 50;
%         
% for p = 1:number_channels
% maxresp(p) = find(LFP_data2(start_point:end_point,p) == max(max(LFP_data2(start_point:end_point,p))));
% end
% 
% % maxresp= -maxresp;
% 
% 
% maxresp = maxresp - min(min(maxresp));
% 
% figure;
% 
% for p = 1:number_channels
% subplot(4,8,positions(p));
% imagesc((maxresp(p))'); colormap(gray);caxis([min(min(maxresp)) max(max(maxresp))]);
% head=num2str(p,'%d');
% ylabel(head)
% end

%     end
% end

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 203
        
for k = 1:32
[y_delay(k), x_res(k)] = find(LFP_rectify_exc == max(LFP_rectify_exc(:,k)))
[y_delay2(k), x_res2(k)] = find(LFP_rectify_exc2 == max(LFP_rectify_exc2(:,k)))

end



y_delay = y_delay - puff_arrival;
% y_delay = y_delay - min(y_delay);
% y_delay2 = y_delay2 - min(y_delay2);
y_delay2 = y_delay2 - puff_arrival;


y_delay = y_delay/sampling_rate*1e3;
y_delay2 = y_delay2/sampling_rate*1e3;

MF_delays = y_delay;
CF_delays = y_delay2;


  figure;

        for l = 1:min(size(pos_x))
            subplot(4,1,l)
            bar(MF_delays(:,pos_x(l,:)));
            ylabel('ms','fontsize',12,'fontweight','b')
            title('MF Response Delays')

        end
        
         figure;
        for l = 1:min(size(pos_x))
            subplot(4,1,l)
            bar(CF_delays(:,pos_x(l,:)));
            ylabel('ms','fontsize',12,'fontweight','b')
            title('CF Response Delays')

        end

        figure;               

        for l = 1:length(pos_y)
            subplot(1,8,l) 
            bar(MF_delays(:,pos_y(l,:)));
            grid on
            axis tight
            title('MF Delays')

        end

         figure;
        for l = 1:length(pos_y)
            subplot(1,8,l) 
            bar(CF_delays(:,pos_y(l,:)));
            grid on
            axis tight
            title('CF Delays')

        end
        

    end
end

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 204
        
        
%         x=[0:1:1000];
%         y=x.*x;
%         std1=y*.1;
%         y1=y+std1; y2=y-std1;
        
        f = linspace(1,8000,1025);
%         logf = 10*log10(f);
%         logf(1) = 0;
        xshade = [f,fliplr(f)];
        
        coh   = Ver1_Cere_anest;
        
        meanC1 = mean(coh');
        stdC1 = std(coh');
        y1 = (meanC1  + stdC1);
        y2 = (meanC1  - stdC1);
        
        yshade = [y1,fliplr(y2)];
        
        fill(xshade,yshade,[0.65,0.35,0.55],'EdgeColor',[0.45,0.45,0.45]);
        hold on
        plot(f,meanC1,'b','LineWidth',1)
        
    end
end
        

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 205
% unfiltered_data = DATA;
% for p = 1:number_channels
% unfiltered_data_cere_electrode(:,p) = unfiltered_data(:,positions_sorted(p));
% end
% 
% unfiltered_data_cortex = unfiltered_data(:,33:64);
% for p = 1:32
% unfiltered_data_cor_electrode(:,p) = unfiltered_data_cortex(:,positions_sorted(p));
% end


figure;subplot(2,1,1);
imagesc(corrcoef(unfiltered_DATA(:,1:number_channels)));axis tight

[Corr_Cerebellum,p] = corrcoef(unfiltered_DATA(:,1:number_channels)) ;

% Corr_Cerebellum = Corr_Cerebellum(1,:);
subplot(2,1,2)
imagesc(corrcoef(unfiltered_DATA(:,number_channels:number_channels*2)));axis tight
[Corr_Cortex,p] = corrcoef(unfiltered_DATA(:,number_channels:number_channels*2)) ;



Correlation_Unfiltered_data_cere = zeros(31,31);
Correlation_Unfiltered_data_cor = zeros(31,31);

for i=1:number_channels
    
Correlation_Unfiltered_data_cere(i,i:31) = Corr_Cerebellum(i,i:31);
Correlation_Unfiltered_data_cor(i,i:31) = Corr_Cortex(i,i:31);

end
figure;surface(Correlation_Unfiltered_data_cere) ; caxis([0 1]);axis tight
figure;surface(Correlation_Unfiltered_data_cor) ; caxis([0 1]);axis tight


    end
end



for x = 1:length(functions_to_run)
    if functions_to_run(x) == 206

        LFP_data_sorted = LFP_data_sorted(4900:5500,1:32);
        
% X = input('desired channel')      % to process distinct channel manually

figure;surface(corrcoef(LFP_data_sorted)); % Open the Surface plot to double check 


for X = 1:min(size(LFP_data_sorted)) 
    
row1 = [1:8]; row2 = [9:16]; row3 = [17:24]; row4=[25:32]; % Divide into row vectors

X

if (find(row1 ==X))
    
medial = find(row1<X)
lateral = find(row1>X)          %% Look for medial-lateral direction

rost_caud = [X+8 X+16 X+24]     % Assign the rost-caudal position
row = row1;
y_row = 300;
y_row2 = y_row-20;

x_row = mod(X-1,8)*55;
end

if (find(row2 ==X))
    
medial = find(row2<X)
lateral = find(row2>X)          % 2.row

rost_caud = [X-8 X+8 X+16]      % Intersected Rost-caudal contacts
row = row2;
medial = medial+max(row1);      %% NEEDED for plotting
lateral = lateral+max(row1);
y_row = 225;
y_row2 = y_row-20;

x_row = mod(X-1,8)*55;

end

if (find(row3 ==X))
    
medial = find(row3<X)
lateral = find(row3>X)

rost_caud = [X-16 X-8 X+8]
row = row3;                     %3.row
medial = medial+max(row2);
lateral = lateral+max(row2);
y_row = 125;
y_row2 = y_row-20;

x_row = mod(X-1,8)*55;

end

if (find(row4 ==X))
    
medial = find(row4<X);
lateral = find(row4>X);
                                    %4.row
rost_caud = [X-8 X-16 X-24]
row = row4;
medial = medial+max(row3);
lateral = lateral+max(row3);
y_row=50;
y_row2 = y_row-20;

x_row = mod(X-1,8)*55;

end

R=ones(8,1)                         % Autocorrelation Coeff is defined as '1'.

for i=1:length(medial)
    r = corrcoef(LFP_data_sorted(:,X),LFP_data_sorted(:,medial(i)));
    R(i) = r(2,1) ;   
end

for ii=1:length(lateral)
    r = corrcoef(LFP_data_sorted(:,X),LFP_data_sorted(:,lateral(ii)));
    R(length(medial)+1+ii) = r(2,1); 
  
end


for k=1:length(rost_caud)
   
    r = corrcoef(LFP_data_sorted(:,X),LFP_data_sorted(:,rost_caud(k)))
    Rc(k) = r(2,1) ;   
end

electrode=zeros(8,4);
electrode(row) = R;
electrode(rost_caud) = Rc;
figure(100);clf;imagesc(electrode')

Rmean = mean(R);
Rcmean =mean(Rc);
A = num2str(Rmean);
B = num2str(Rcmean);
text('units','pixels','position',[x_row y_row],'fontsize',12,'string',A) 
text('units','pixels','position',[x_row y_row2],'fontsize',12,'string',B) 

pause(1);

end

    end
end

for x = 1:length(functions_to_run)
    if functions_to_run(x) == 207
  unfiltered_data=unfiltered_data(:,1:signal_length);

duration = 20;
LOOP = 20;
Fs = sampling_rate;
DURATION=duration/LOOP;
% use it for quite data
        
% 
% Ns * Ts = Ncycles * Trej
% 
% Ts = 1/sampling_rate;
% Ncycles = # of cycles will be each averaged data
% Tnoise = 1/60 (60Hz elimination)

Trej = 1/60;
Ncycles = 100; % for 16kHz sampling rate
Ts = 1/sampling_rate;
Ns = Ncycles*Trej/Ts;


duration = 1;
LOOP = 20;
Fs = sampling_rate;
DURATION=duration/LOOP;
N=2048;

%% Taking the Power spectrum first for each period then average it.

% for ii=1:LOOP
%     
% [Pww(:,ii), F] = pwelch(detrend(unfiltered_data((ii-1)*DURATION*Fs+1:ii*DURATION*Fs,:)),hanning(N),N/2,N,sampling_rate);
% Pww(:,ii)  = 5*log10(Pww(:,ii).^2);
% end
% 
% averaged_Pww = sum(Pww(:,1:20)',1)/20;
% averaged_Pww = averaged_Pww';

% figure;plot(F,averaged_Pww)       
% Apparently it does not make any difference taking the power spect from
% the averaged DATA

%% Averaging the data

tt=1/Fs:1/Fs:duration;

DATA = zeros(min(size(unfiltered_data)), Fs*DURATION);
unfiltered_data = unfiltered_data';

for Y=1:LOOP,
 
    DATA = DATA+unfiltered_data(:,(Y-1)*DURATION*Fs+1:Y*DURATION*Fs);
end

DATA=DATA'/LOOP;


unfiltered_data = unfiltered_data';

    end
end


 for x = 1:length(functions_to_run)
    if functions_to_run(x) == 210
        
%         unfiltered_DATA = LFP_data4(:,1:31);
      
        unfiltered_DATA(:,32) = mean(unfiltered_DATA(:,1:31)');
        unfiltered_data(:,32) = mean(unfiltered_data(:,1:31)');       
        N=2048;

%         unfiltered_DATA = LFP_data_sorted;
%       I = input('Like to return desired channels (1) or Electrode Organization (0) or Random (any)'); 
          I=0;

      if (I == 1)
          ch1 = input('1.pair channel 1 ');
          ch2 = input('1.pair channel 2 ');
          ch3 = input('2.pair channel 1 ');
          ch4 = input('2.pair channel 2 ');
          ch5 = input('3.pair channel 1 ');
          ch6 = input('3.pair channel 2 ');
          ch7 = input('4.pair channel 1 ');
          ch8 = input('4.pair channel 2 ');
          
      elseif (I == 0 )
          
          RC1 = randi(8);
          RC2 = randi(16);             % Define your preference whether in Rostro-Caudal (RC) or Medio-Lateral (ML) direction
          RC3 = 24;
          ML1  = 3;
          ML2 = 7 ; 
          
          
          ch1 = randi(8);
          ch2 = ch1+RC1;
          
          ch3 = randi([8 16]);
          ch4 = ch3+RC2;
          
          ch5 = randi([16 24]);
          ch6 = ch5+RC1 ;
          
          ch7 = randi([24 31]);
          
         ch8 = ch7 - ceil(mod(ch7,24)/2);

      else    

ch1 = randi(number_channels);
ch2 = randi(number_channels);

ch3 = randi(number_channels);
ch4 = randi(number_channels);

ch5 = randi(number_channels);
ch6 = randi(number_channels);

ch7 = randi(number_channels);
ch8 = randi(number_channels);
      end
      
%%



%        Windowing and Video segments start here
%        Playing video file simultaneously with temporal coherence analysis
%        Video_data max is taken as 5 sec
       
%       video_data = VideoReader(['videos/trial' num2str(trial) '.avi']);
%       all_frames = read(video_data);
%       
%       clear i
%       Dur = length(unfiltered_DATA)/sampling_rate;              % recalculating the duration of data for temp coherence
%           
%       II = input('Would you like the movie file to be played? (1/0)');  % Run the video file/NOt
%      
%  for i=1:1:Dur;                                               % 1 sec window spacing 
%       i;
%        
%      if (II ==1)
%       if (i<=5)                                             % limit the video file with 5 sec 
%       l =  i*video_frame_rate;                                %% Synch frame video file  to data file   
%       for ll = (l-video_frame_rate) +1 : l
%           
%         figure(112);clf;imagesc(all_frames(:,:,:,ll));colormap(gray);
%         t1 = num2str((i-1),'% d'); t2 = num2str(i,'% d');
%         title([t1,'s - ',t2,'s'])
%         pause(0.1)
%            
%       end
%       end   
%  end
% 
% N=1024;                
% unfiltered_DATA = unfiltered_data(((i-1)*sampling_rate+1):i*sampling_rate,:);      

%%

[COHERE,f] =mscohere(unfiltered_DATA(:,ch1), unfiltered_DATA(:,ch2),hanning(N/4) ,N/8 ,N/2 ,sampling_rate);
[COHERE2,f] =mscohere(unfiltered_DATA(:,ch3), unfiltered_DATA(:,ch4),hanning(N/4) ,N/8 ,N/2 ,sampling_rate);
[COHERE3,f] =mscohere(unfiltered_DATA(:,ch5), unfiltered_DATA(:,ch6),hanning(N/4) ,N/8 ,N/2 ,sampling_rate);
[COHERE4,f] =mscohere(unfiltered_DATA(:,ch7), unfiltered_DATA(:,ch8),hanning(N/4) ,N/8 ,N/2 ,sampling_rate);

figure;
hold on; subplot(2,1,1);
plot(f,COHERE);hold on
plot(f,COHERE2,'r');plot(f,COHERE3,'g');plot(f,COHERE4,'k');
axis([ 0 2000 0 1])
COH1 = num2str([ ch1 ch2 ],'% d');
COH2= num2str([ ch3 ch4 ],'% d');
COH3 = num2str([ ch5 ch6 ],'% d');
COH4 = num2str([ ch7 ch8 ],'% d');
legend(COH1,COH2,COH3,COH4)
t1 = num2str((i-1),'% d'); t2 = num2str(i,'% d');
title([t1,'s - ',t2,'s'])



%% LOAD electode image into matlab

a=imread('electrode_contacts_numbered.png','png');
subplot(2,1,2)
image(a);
hold on

% Specifiy the contact coordinates on the image

refx1 = mod(ch1,8);
if(refx1 == 0);
    refx1 = 8;
end
X_ch1 = 140+((refx1-1) * 50);
Y_ch1 = (ceil(ch1/8) * 50);
plot(X_ch1,Y_ch1,'--bs','LineWidth',10) ;       % 1. pairs

refx2 = mod(ch2,8);
if(refx2 == 0);
    refx2 = 8;
end


X_ch2 = 140+((refx2-1) * 50);
Y_ch2 = (ceil(ch2/8) * 50);
plot(X_ch2,Y_ch2,'--bs','LineWidth',10)   ;  

refx3 = mod(ch3,8);
if(refx3 == 0);
    refx3 = 8;
end
X_ch3 = 140+((refx3-1) * 50);
Y_ch3 = (ceil(ch3/8) * 50);
plot(X_ch3,Y_ch3,'--rs','LineWidth',10)     %2.pairs

refx4 = mod(ch4,8);
if(refx4 == 0);
    refx4 = 8;
end

X_ch4 = 140+((refx4-1) * 50);
Y_ch4 = (ceil(ch4/8) * 50);
plot(X_ch4,Y_ch4,'--rs','LineWidth',10)


refx5 = mod(ch5,8);
if(refx5 == 0);
    refx5 = 8;
end
X_ch5 = 140+((refx5-1) * 50)    ;              %3.pairs
Y_ch5 = (ceil(ch5/8) * 50);
plot(X_ch5,Y_ch5,'--gs','LineWidth',10)

refx6 = mod(ch6,8);
if(refx6 == 0);
    refx6 = 8;
end
X_ch6 = 140+((refx6-1) * 50);
Y_ch6 = (ceil(ch6/8) * 50);
plot(X_ch6,Y_ch6,'--gs','LineWidth',10)

refx7 = mod(ch7,8);                         % 4.pairs
if(refx7 == 0);
    refx7 = 8;
end
X_ch7 = 140+((refx7-1) * 50);
Y_ch7 = (ceil(ch7/8) * 50);
plot(X_ch7,Y_ch7,'--ks','LineWidth',10)

refx8 = mod(ch8,8);
if(refx8 == 0);
    refx8 = 8;  
end
X_ch8 = 140+((refx8-1) * 50);
Y_ch8 = (ceil(ch8/8) * 50);
plot(X_ch8,Y_ch8,'--ks','LineWidth',10)



% N=2048;
% [Pww_cer, F] = pwelch(unfiltered_DATA(:,1:number_channels),hanning(N),N/2,N,sampling_rate); % Cerebellum Power Spectrum
% % [Pww_cor, F] = pwelch(unfiltered_data(:,32:62),hanning(N),N/2,N,sampling_rate); % Cerebellum Power Spectrum
% 
% cc = hsv(50);
% 
% 
% Pww_cer = 10*log10(Pww_cer);
% % Pww_cor = 10*log10(Pww_cor);
% 
% figure(11);clf;plot(F,Pww_cer,'color',cc(trial,:));axis([0 1000 -140 -90])
% pause(3)

%         end
    end
 end
 
 % check the rostro-caudal cross-corr spectrum estimation
 % Enter the desired channels
 
 
 
  for x = 1:length(functions_to_run)
    if functions_to_run(x) == 211
 
 for i=1:6

base_ch1(i) = input('enter first one:'); % enter the desired line starting from first contact
    base_ch2(i) = input('enter second one:');
    
Pwwcere(:,i) = cpsd(unfiltered_DATA(:,base_ch1(i)),unfiltered_DATA(:,base_ch2(i)),hanning(N),N/2,N,sampling_rate);

end

b1 = num2str(base_ch1,' %d') 
b2 = num2str(base_ch2,' %d') 
b = [b1; b2]

Pwwcere = 10*log10(abs(Pwwcere));

figure;plot(f,Pwwcere);axis([0 2000 -140 -90]);hold on

ylabel(b)
    end
  end
  
 

  for x = 1:length(functions_to_run)
    if functions_to_run(x) == 300
 

figure;imagesc(corrcoef(LFP_data(:,1:31)))
figure;plot(t,1e6*LFP_data(:,1:31))
% figure;plot(1e6*LFP_data(:,[48 46]))

% figure;imagesc(corrcoef(LFP_data_sorted(:,32:62)))

 
% figure;plot(mean(LFP_data_sorted(:,1:16)'))
% 
% figure;plot(mean(LFP_data_sorted(:,17:32)'))

Ch = randi(31);
% 
% figure;plot(1e6*LFP_data1(:,1:31),'b')
% hold on;plot(1e6*LFP_data2(:,1:31),'k');
% plot(1e6*LFP_data3(:,1:31),'r');
% plot(1e6*LFP_data4(:,1:31),'m--')



% subplot(1,2,1);plot(1e6*LFP_data_sorted(:,32:62))

%% PLot the LFP potentials in Electrode orientation

%  figure;
%  subplot(4,1,1,'Color','w');plot(LFP_data_sorted(:,[(1:8)])); 
%  
%   subplot(4,1,2);plot(LFP_data_sorted(:,[(9:16)]))      % Medio - Lat Direction
%  subplot(4,1,3);plot(LFP_data_sorted(:,[(17:24)]))
%  subplot(4,1,4);plot(LFP_data_sorted(:,[(25:31)]))
% 
%  figure;
%  subplot(1,3,1);plot(LFP_data_sorted(:,[1 9 17 25]));title('Ch [1 9 17 25]')
% axis ([4800 5000 min(min(LFP_data_sorted(:,[1 9 17 25]))) max(max(LFP_data_sorted(:,[1 9 17 25]))) ])
%   subplot(1,3,2);plot(LFP_data_sorted(:,[4 12 20 28]));title('Ch [4 12 20 28]')
%   axis ([4800 5000 min(min(LFP_data_sorted(:,[4 12 20 28]))) max(max(LFP_data_sorted(:,[4 12 20 28]))) ])   % Rost - Cau Direction
%  subplot(1,3,3);plot(LFP_data_sorted(:,[7 15 23 31]));title('Ch [7 15 23 31]')
%  axis ([4800 5000 min(min(LFP_data_sorted(:,[7 15 23 31]))) max(max(LFP_data_sorted(:,[7 15 23 31]))) ])
% 
%  
% 
% 
% figure;spgrambw(mean(unfiltered_data(:,1:31)'),sampling_rate,'pwjat', [],[200 5 300],10,[0 .001 10],[]);
% 
% figure;spgrambw(mean(unfiltered_data(:,32:62)'),sampling_rate,'pjwat', [],[10 1 500],20,[0 .001 1],[]);

% figure;
% for p = 1:31
% subplot(4,8,positions_31channels(p));
% spgrambw(unfiltered_data(:,p),sampling_rate,'pjwat', [],[10 1 1000],20,[0 .001 1],[]);
% end

%%
 
    end
  end
  
  
  
  for x = 1:length(functions_to_run)
    if functions_to_run(x) == 301


        % Moving Window parameters
        
per=20e-3; %size of window in ms scale
per=1/per;
pp = sampling_rate/per; 
       
       % Applying window in selected freq bands.

for i=pp+1:pp:length(LFP_data1)
k = (i-1)/pp;
C1(k,1:31) = (abs(mean(corrcoef(LFP_data1(i:i+pp-1,1:number_channels)))));
end

for i=pp+1:pp:length(LFP_data2)
k = (i-1)/pp;
C2(k,1:31) = (abs(mean(corrcoef(LFP_data2(i:i+pp-1,1:number_channels)))));
end

for i=pp+1:pp:length(LFP_data3)
k = (i-1)/pp;
C3(k,1:31) = (abs(mean(corrcoef(LFP_data3(i:i+pp-1,1:number_channels)))));
end

for i=pp+1:pp:length(LFP_data4)
k = (i-1)/pp;
C4(k,1:number_channels) = (abs(mean(corrcoef(LFP_data4(i:i+pp-1,1:number_channels)))));
end

% for i=pp+1:pp:length(LFP_data_sorted)
% k = (i-1)/pp;
% C5(k,1:number_channels) = (abs(mean(corrcoef(LFP_data_sorted(i:i+pp-1,1:number_channels)))));
% end
% 




% 
% 
% for i=pp+1:pp:length(LFP_data4)
% k = (i-1)/pp;
% C4(k,1:number_channels) = (abs(mean(corrcoef(LFP_data4(i:i+pp-1,1:number_channels)))));
% end
% 
% for i=pp+1:pp:length(LFP_data4)
% k = (i-1)/pp;
% C11(k,1:number_channels) = min(min(corrcoef(LFP_data1(i:i+pp-1,1:number_channels),LFP_data2(i:i+pp-1,1:number_channels))));
% end
% 
% 
% for i=pp+1:pp:length(LFP_data4)
% k = (i-1)/pp;
% C22(k,1:number_channels) = min(min(corrcoef(LFP_data2(i:i+pp-1,1:number_channels),LFP_data3(i:i+pp-1,1:number_channels))));
% end
% 
% for i=pp+1:pp:length(LFP_data4)
% k = (i-1)/pp;
% C33(k,1:number_channels) = min(min(corrcoef(LFP_data3(i:i+pp-1,1:number_channels),LFP_data4(i:i+pp-1,1:number_channels))));
% end
% 
% for i=pp+1:pp:length(LFP_data4)
% k = (i-1)/pp;
% C44(k,1:number_channels) = min(min(corrcoef(LFP_data1(i:i+pp-1,1:number_channels),LFP_data3(i:i+pp-1,1:number_channels))));
% end


%% This will be run with EPs

TT = pp/sampling_rate*t(end);
% 
figure;plot(t(pp:end)/TT,130+(1e6*abs(LFP_data_sorted(pp:end,1:31))));
hold on


%% This will be with Spontaneous trials.
% % 
% TTT = round(t(end)/length(C1));
% tt = t/TTT;
% figure;plot(tt,130+abs(mean(1e6*LFP_data_cer')),'k')
% hold on
%%

C = [C1';C2';C3';C4'];
% CC = [C11' ; C22'; C33'; C44']; 

contourf(C); colorbar ;caxis([0 1])


% TT = pp/sampling_rate*t(end);
% % 
% figure;plot(t(pp:end)/TT,130+(1e6*abs(LFP_data(pp:end,1:31))));
% % figure;plot(tt,130+abs(1e6*LFP_data_cer),'k')
% hold on
% imagesc(CC); colorbar ;caxis([0 0.5]) ; ylabel('cross-band Corr, (5-30_30-50,30-50_70-220,70-220_220-400Hz)')
% 


[ mean(mean(C1)') mean(mean(C2)') mean(mean(C3)') mean(mean(C4)') ]

figure;plot ([ mean(mean(C1)') mean(mean(C2)') mean(mean(C3)') mean(mean(C4)') ])
% hold;plot(2.5, mean(mean(C5)'),'*r')


    end
  end

% for k=1:31
% 
% t=linspace(0,10,length(unfiltered_data));
% sig1=[t' LFP_data1(:,11)];
% sig2=[t' LFP_data1(:,23)];
% LL=LFPsynch(sig1,sig2,1,0.05);
% sig1=[t' LFP_data2(:,11)];
% sig2=[t' LFP_data2(:,23)];
% LFPsynch(sig1,sig2,1,0.05);
% sig1=[t' LFP_data3(:,11)];
% sig2=[t' LFP_data3(:,23)];
% LL=LFPsynch(sig1,sig2,1,0.05);
% % for ii=1:31
% sig1=[t' LFP_data4(:,11)];
% sig2=[t' LFP_data4(:,23)];
% LL=LFPsynch(sig1,sig2,1,0.01);
% % end
% % pause;
% end

% absLFP=abs(LFP_data);
%  for i=1:31
% amps(i)=mean(absLFP(4801:5280,i))/mean(absLFP(1:4800,i));
% end
% out=[median(amps)]
% %%
% hil_lfp4_amps=median(abs(hilbert(LFP_data4(4700:5200,:))));
% 
% t=linspace(0,dur,length(LFP_data));
% figure(66);
% for i=1:31
% %     subplot(4,10,positions_W_31channels2_splots(i));hold on;h=bar(hil_lfp4_amps(i),'r','EdgeColor','m','FaceAlpha',0.1);axis([0.5 1.5 0 1e-5]);axis off
%     
%         subplot(4,10,positions_W_31channels2_splots(i));hold on;plot(t(4800:6800),LFP_data(4800:6800,i)*1e6,'g','LineWidth',2);axis tight;
%         set(gca,'FontSize',14)
% %         axis([290 330 -30 30])
%            axis off
% end
%  subplot(4,10,40);axis on;axis tight
 
 
 figure;plot(t,LFP_data(:,1),'b')
% hold on;plot(t,LFP_data(:,32:48),'r')
 

%to eliminate delays between Dev-1 and Dev-2 (For EPs only)

 for x = 1:length(functions_to_run)
    if functions_to_run(x) == 315

fs=16000;
temp1=mean(LFP_data(:,16:31)');
temp2=mean(LFP_data(:,32:end)');
cc=xcorr(temp1,temp2);
[~,del]=max(cc)

del=del-fs+1;
% for i=1:31
% temp=ones(size(LFP_data(:,32:end)));
% temp(:,i)=temp(:,i)'.*mean(LFP_data(:,31+i)');
% end
temp=LFP_data(:,32:end);
%temp(del:end,:)=LFP_data(1:end-del+1,32:end);
LFP_data(:,32:end)=temp;
 figure;plot(t,LFP_data(:,1:31),'b')
%hold on;plot(t,LFP_data(:,32:48),'r')
    end
 end