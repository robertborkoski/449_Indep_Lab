%......BME 449 Independent Lab............................................%
%......Dr. Sarles.........................................................%
%......By Robert Borkoski.................................................%
%......Last Revised 11/30/20..............................................%

clear all, close all, clc, format long g
global k % keeps track of what EEG recording is being examined
filename = {'S011R10'}; % currently using one EEG recording; code is largely parametrized to accept 100+ recordings at once
[Data,E,PSD] = importEEG(filename); % general import / preliminary signal processing
[peaks] = findpeaks(PSD) % find locations of PSD peaks

% plotting PSD results - by region

figure(1)
subplot(2,2,1)
surf(PSD.Foot_1_t1.PSD_Delta)
subplot(2,2,2)
surf(PSD.Foot_1_t1.PSD_Theta)
subplot(2,2,3)
surf(PSD.Foot_1_t1.PSD_Alpha)
subplot(2,2,4)
surf(PSD.Foot_1_t1.PSD_SMR)
figure(2)
subplot(2,2,1)
surf(PSD.Foot_1_t2.PSD_Delta)
subplot(2,2,2)
surf(PSD.Foot_1_t2.PSD_Theta)
subplot(2,2,3)
surf(PSD.Foot_1_t2.PSD_Alpha)
subplot(2,2,4)
surf(PSD.Foot_1_t2.PSD_SMR)

function [Data,E,PSD] = importEEG(filename) % Imports edf format EEG data into matlab, perform peliminary signal processing (should be separated into two separate functions but it ain't broke so I ain't fixin' it
global k % keep track of which EEG recording is being examined
records = {'t1','t2'}; % look at EEG 
for k = 1:length(filename) % load in all requested files
    filename{k} = strcat(filename{k},'.edf') % look for edf file
    [hdr, record] = edfread(filename{k}); % import data from session 1
    srate = 1 / hdr.samples(k); % establish sample rate of 160 Hz
    [record_t1, record_t2] = intervals(record,filename,srate); % cut out rest periods and re-splice according to task
    r.t1 = record_t1; % create intermediate r data structs for spliced EEG signals (mainly for testing / dev)
    r.t2 = record_t2;
    for m = 1:length(records) % eval both t1 and t2
        x = size(r.(records{m})); % determine how many samples per EEG capture
        subject_lines = ["Hand_","Foot_"]; % looking at imagined hand, foot movement
        subject_name = strcat(subject_lines(mod(k,2)+1),num2str(ceil(k/2)),'_',records{m}); % form data struct paths
        for n = 1:64 % only read first two channels - change to x(2) when looking at all 64 channels
            time = 0+srate:srate:srate*x(2); % create time domain using sample rate
            chn_name = strcat("Chn_",num2str(n)); % current channel name (1-64)
            Data.(subject_name).(chn_name) = [time;r.(records{m})(n,:)]; % create struct Data containing separated EEG recordings by channel
            [Delta_E,Theta_E,Alpha_E,Beta_E,SMR_E] = eeg2rhythm(Data.(subject_name).(chn_name)(2,:)); % apply hamming window and perform time-frequency domain analysis for each channel individually
            E.(subject_name).(chn_name) = ["Delta",Delta_E(1);"Theta",Theta_E(1);"Alpha",Alpha_E(1);"Beta",Beta_E(1);"SMR",SMR_E(1)]; % pull out PSD values into one table (for testing / dev)
            fprintf('.') % progress bar (testing / dev for runtime)
        end
        grid = headset(); % create blank 2D EEG grid
        s = size(grid); % Should be 11x11 based on current square checkered design
        regions = ["Delta","Theta","Alpha","Beta","SMR"]; % to create struct names
        for o = 1:5 % for each region
            grid = headset(); % create blank 2D EEG grid
            for d = 1:11 % fix to size when parametrizing
                for c = 1:11 % fix to size when parametrizing
                    if isnan(grid(d,c)) == 0 % for all non-empty channel names
                        chn_name = strcat("Chn_",num2str(grid(d,c))); % create channel name header
                        grid(d,c) = E.(subject_name).(chn_name)(o,2); % use grid as ref to build PSD grid
                    end
                end
            end
            subname = strcat("PSD_",regions(o)); % create PSD struct name
            PSD.(subject_name).(subname) = fillmissing(grid,'linear'); % fill in NaN elements by linear column interpolation
        end
    end
end
end

function grid = headset()
   grid = xlsread('grid.xlsx') % import grid template xls
end

function [Delta_E,Theta_E,Alpha_E,Beta_E,SMR_E] = eeg2rhythm(in_signal) % time-frequency analysis of EEG signal, calculate PSD
coder.extrinsic('mapminmax')
coder.extrinsic('display')
Fs = 160; % sample rate = 160 Hz
L = length(in_signal);
eeg_samples = in_signal;
ch = 1; % legacy code bc inputting as single-channel signals anyways

% Implement Hamming window and apply FFT with Hamming window
win = hamming(length(eeg_samples)); % apply hamming window
win = win'; % transpose
fft_eeg_samples = fft(eeg_samples.*repmat(win,1,ch)); % perform fft on EEG signal using hamming window

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum,
% P1 based on P2 and the even-valued signal length L.
P2 = abs(fft_eeg_samples./L); % non-negative
P1 = P2(1,:); % grab first row of data
P1(2:end-1,:) = 2*P1(2:end-1,:); % resize to fit

f = Fs*(0:(L/2))/L; % create frequency vector
freq = [0:(160/length(fft_eeg_samples)):160]; % create sep freq vector for plotting (diff size)
freq = freq(1:end-1); % resize to fit
figure(1) % plotting results (testing / dev)
%Scale the PSD?? for each channel to [0,1]
P1 = mapminmax(P1',0,1)'; % scale transpose
bar(freq,P1) % plot results
% Find indices of Theta, Alpha and Beta band
Delta_idx = find(1 <= f & f < 4 );
Theta_idx = find(4 <= f & f < 7);
Alpha_idx = find(7 <= f & f < 13);
Beta_idx = find(13 <= f & f < 25);
SMR_idx = find(12 <= f & f < 15);

% find total energy in each band
Delta_E = sum(P1(Delta_idx).^2);
Theta_E = sum(P1(Theta_idx).^2);
Alpha_E = sum(P1(Alpha_idx).^2);
Beta_E = sum(P1(Beta_idx).^2);
SMR_E = sum(P1(SMR_idx).^2);
end

function [record_t1, record_t2] = intervals(record,filename,srate) % deconstruct and re-splice signals by task
global k % keep track of which EEG recording is being analyzed
filename{k} = filename{k}(1:end-4); % remove .edf extension from filename
filename{k} = strcat(filename{k},'.xlsx'); % look for xls of imagined movement time intervals
intervals = xlsread(filename{k}); % i'm kinda dumb for naming a var same as func name but I'm too lazy to change it.. if you're reading these comments, hello :)
time = intervals(1,1):srate:intervals(1,end); % create time vector using sample rate
startidx_t1 = find(intervals(2,:) == 1); % finding t1 elements
stopidx_t1 = startidx_t1 + 1; % essentially, grabbing the time frames for each t1/t2 interval
startidx_t2 = find(intervals(2,:) == 2); % finding t2 elements
stopidx_t2 = startidx_t2 + 1;
first = 1; % need for building spliced EEG signals
for s = 1:length(startidx_t1) % splice together all segments
    time_startidx_t1 = find(time == intervals(1,startidx_t1(s))); % grab actual EEG data based on t1/t2 intervals
    time_stopidx_t1 = find(time == intervals(1,stopidx_t1(s)));
    if first == 1 % starting the spliced EEG
        record_t1 = record(:,time_startidx_t1:time_stopidx_t1);
        first = 0;
    else % adding onto spliced EEG
        record_t1 = [record_t1, record(:,time_startidx_t1:time_stopidx_t1)];
    end
end
first = 1;
for u = 1:length(startidx_t2); % same as above loop, but for t2 interval
    time_startidx_t2 = find(time == intervals(1,startidx_t2(u)));
    time_stopidx_t2 = find(time == intervals(1,stopidx_t2(u)));
    if first == 1
        record_t2 = record(:,time_startidx_t2:time_stopidx_t2);
        first = 0;
    else
        record_t2 = [record_t2, record(:,time_startidx_t2:time_stopidx_t2)];
    end
end
end

function [peaks] = findpeaks(PSD) % analyze PSD data to find peaks

struct = 'Foot_1_'; % build struct name
regions = {'Delta','Theta','Alpha','Beta','SMR'}; % region names (for looping)
epochs = {'t1','t2'}; % time names (for looping)
grid = headset(); % create blank template of EEG electrode grid
for i = 1:length(regions) % for each freq region
    for j = 1:length(epochs) % for each task
        nan_absent = 0; % need for checking presence of NaN elements (grid interpolation from above screws up at the edges but this corrects for it)
        while nan_absent == 0; % if there's a NaN
            struct_father = strcat(struct,epochs{j}); % create struct name
            struct_child = strcat('PSD_',regions{i}); % create struct name
            maximum = max2d(PSD.(struct_father).(struct_child)); % for whatever reason I needed to hard-code my own max function to get this to work on a 2d array
            max_idx = find(PSD.(struct_father).(struct_child) == maximum); % find max in actual PSD data
            chn = grid(max_idx); % correlate max in data to what channel it is on grid
            if isnan(chn) == 1 % if there's an erroneous maximum because of an NaN in grid
                PSD.(struct_father).(struct_child)(max_idx) = NaN; % correct PSD matrix to remove false max
            else
                nan_absent = 1;
            end
        end
        peaks.(struct_father).(struct_child) = chn; % create new struct containing peak channel locations
            
    end
end
end

function [maximum] = max2d(table) % needed to hard-code my own little max function for a 2D table??? idk
dim = size(table);
for i = 1:dim(2)
    peak(i) = max(table(:,i));
end
maximum = max(peak);
end