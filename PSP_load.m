%% PSP_load
%  Loading traces, running filters, and run DetectPSP to Analyze EPSCs

%  By Simon Sun; adapted from BSS and SFO detectPSPs and SDS hw_load.m
%  Last edited: 16.02.17

%%   Ask for file numbers to be loaded to use for later
    
    clear all
    close all
    clc

    disp('---------- Determining File Path ----------');
    disp(' ');
    prompt_year = 'Enter year of abf files: ';
    year = input(prompt_year);

    prompt_month = 'Enter month of abf files: ';
    month = input(prompt_month);

    prompt_day = 'Enter day of abf files: ';
    day = input(prompt_day);
    if day<10
        day = ['0',num2str(day)];
    else
    end

    prompt_range_s = 'Enter the beginning of range of file numbers to analyze: ';
    range_s = input(prompt_range_s);

    prompt_range_e = 'Enter the end of the range of file numbers to analyze: ';
    range_e = input(prompt_range_e);
   
    date = [num2str(year),num2str(month),num2str(day)];
    
    prompt_exp_num = 'Enter the experiment number (i.e. rwt106): ';
    exp_num = input(prompt_exp_num);
    
    prompt_date = 'Enter the date in numberical form: YYMMDD: ';
    exp_date = input(prompt_date);
    
    prompt_geno = 'Enter the genotype (i.e. WT or TS2): ';
    geno = input(prompt_geno);
    
    prompt_treat = 'Enter the treatment type (i.e. CTL or TTX): ';
    treat = input(prompt_treat);
    
    prompt_cell = 'Enter the cell number: ';
    cell_num = input(prompt_cell);
    
    prompt_sf = 'Enter the sampling frequency (usually 10000Hz): ';
    sf = input(prompt_sf);
    
    prompt_lsweeps = 'Enter the number of seconds of each sweep: ';
    lsweeps = input(prompt_lsweeps);
    
    directory = ['/Users/simondsaid/Dropbox/Data/Tsien Lab/',exp_num,'/',exp_date,'/',geno,'/',treat,'/',cell_num,'/'];
    
%% Run abfload on range of files interested

RAW = [];
    
pp = range_s;
    for i = range_s:range_e
        if i < 10;
            date_f = [date,num2str(0),num2str(0)];
            file = [directory,date_f,num2str(i),'.abf'];
        elseif i < 100 && i >= 10
            date_f = [date,num2str(0)];
            file = [directory,date_f,num2str(i),'.abf'];
        else
            date_f = date;
            file = [directory,date_f,num2str(i),'.abf'];
        end
        
        RAW_temp = abfload(file);clc
        RAW_temp = squeeze(RAW_temp(:,1,:));
        [x y] = size(RAW_temp);
        if y > 10 % Change this is RAW is coming out empty, this references the number of sweeps
        else
            if (pp-range_s)+1 == 1
                RAW(:,pp-range_s+1:(pp-range_s)+y) = RAW_temp;
                pp = pp+1;
            else 
                RAW(:,((y*(pp-range_s)+1):(y*((pp-range_s)+1)))) = RAW_temp;
                pp = pp+1;
            end
        end
        
        clear RAW_temp
    end   

    prompt_plot_q = 'Should I plot (Y=1/N=0)? ';
    plot_q = input(prompt_plot_q);
    if plot_q == 1
        plot(RAW)
    else
    end
    
    disp('If you need to eliminate traces (for example, from the end) re-run me.');
    
%% Low Pass filtering of the data    
    
clc

    disp('---------- Low Pass filtering ----------');
    figure; hold on;
    plot(RAW(:,1));
    [b_1000,a_1000] = butter(2,1000/sf,'low');
    [b_750,a_750] = butter(2,750/sf,'low');
    [b_400,a_400] = butter(2,400/sf,'low');
    filt1000 = filtfilt(b_1000,a_1000,RAW(:,1)); 
    filt750 = filtfilt(b_750,a_750,RAW(:,1)); 
    filt400 = filtfilt(b_400,a_400,RAW(:,1)); 
    plot(filt1000); plot(filt750); plot(filt400);
    legend('Original','1000Hz','750Hz','400Hz');
    
    prompt_filt = 'Enter frequency of butterworth filter: ';
    filt = input(prompt_filt);
    [b,a] = butter(2,filt/sf,'low');
    
    [length,num_traces] = size(RAW);
    filtered_data = zeros(length,num_traces);
    
    i=1;
    h = waitbar(0,['Trace ',num2str(i),'/',num2str(num_traces)]);
    for i = 1:num_traces
        filtered_data(:,i) = filtfilt(b,a,RAW(:,i));
        waitbar(i/num_traces,h,['Trace ',num2str(i),'/',num2str(num_traces)]);
    end
    
    delete(h);
    
    prompt_plot_q = 'Should I plot (Y=1/N=0)? ';
    plot_q = input(prompt_plot_q);
    
    if plot_q == 1
        prompt_plot_num = 'How Many? ';
        plot_num = input(prompt_plot_num);
        figure; hold on;
        plot(filtered_data(:,1:plot_num));
    else
    end
%% Run DetectPSP
clc
disp('---------- Event Detection ----------');
    downSampleFreq = sf;

    prompt_PSP = 'Which direction are the events (Down = 1, Up = 0)? ';
    PSPsDown = input(prompt_PSP);

% set detection parameters
    detectParams.timesRMS = 2;

    ... minimum allowable amplitude for alpha functions (in units of samples)
        detectParams.minAmp = -750;
    ... maximum allowable amplitude
        detectParams.maxAmp = -5;
    ... minimum allowable tau for alpha functions (in units of samples)
        detectParams.minTau = 0; %original 500e-6,50e-6
    ... maximum allowable tau
        detectParams.maxTau = inf; %original 5e-3,worked 50e-3
    ... minimum allowable yOffset for alpha functions (in units of mV)
        detectParams.minYOffset = -2000;
    ... maximum allowable yOffset
        detectParams.maxYOffset = 2000;
    ... minimum allowable decay tau
        detectParams.minDecay = 1e-3; % CHANGE FROM 1E-3 ON 11-05-19 (Scott had set to 3e-3)1
    ... maximum allowable decay tau
        detectParams.maxDecay = inf; % CHANGE FROM 50E-3 ON 11-05-19  (Scott had set to 20e-3)
    ... threshold used to determine if the change of derivative is great
        ... enough to separate out two alpha functions
        detectParams.derThresh = 10;%original 10,5 worked for cell4, derthresh=5 for 3 and 4
    ... second EPSP is not recognized if it is closer than this value to
        ... the previous one (in units of samples)
        detectParams.closestEPSPs = 1E-3;
    % threshold for standard error above which a fit with multiple alphas is
    % attempted
    detectParams.errThresh = 0.08;
    % 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
    detectParams.dataFilterType = 3;
    % 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
    detectParams.derFilterType = 3;
    % length of data filter
    detectParams.dataFilterLength = 5E-3;
    % length of derivative filter
    detectParams.derFilterLength = 2E-3;
    % if set to 1 then debugging figures appear
    detectParams.debugging = 0;
    % index of first data point
    detectParams.dataStart = 1;
    % forces a graphical output even if other outputs are taken
    detectParams.forceDisplay = 0;
    % turns off best fitting of start time and decay Tau when 1
    detectParams.noFit = 0;
    % set the lowest current value (pA) for event detection.
    detectParams.threshVal = -1; 
    
    
    disp('Running detectPSPs on first trace ... ');
    tic
    detectPSPs(filtered_data(:,1), PSPsDown, ...
                    'minAmp', detectParams.minAmp, ...
                    ... minimum allowable amplitude for alpha functions (in units of samples)
                    'maxAmp', detectParams.threshVal, ...
                    ... maximum allowable amplitude
                    'minTau', (detectParams.minTau * downSampleFreq), ...
                    ... minimum allowable tau for alpha functions (in units of samples)
                    'maxTau', (detectParams.maxTau * downSampleFreq), ...
                    ... maximum allowable tau
                    'minYOffset', detectParams.minYOffset, ...
                    ... minimum allowable yOffset for alpha functions (in units of mV)
                    'maxYOffset', detectParams.maxYOffset, ...
                    ... maximum allowable yOffset
                    'minDecay', (detectParams.minDecay * downSampleFreq), ...
                    ... minimum allowable decay tau
                    'maxDecay', (detectParams.maxDecay * downSampleFreq), ...
                    ... maximum allowable decay tau
                    'derThresh', (detectParams.derThresh), ...
                    ... threshold used to determine if the change of derivative is ...
                    ... great enough to separate out two alpha functions
                    'closestEPSPs', (detectParams.closestEPSPs * downSampleFreq), ...
                    ... second EPSP is not recognized if it is closer than this ...
                    ... value to the previous one (in units of samples)
                    'errThresh', (detectParams.errThresh), ...
                    ... threshold for standard error above which a fit with ...
                    ... multiple alphas is attempted
                    'dataFilterType', detectParams.dataFilterType, ...
                    ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                    'derFilterType', detectParams.derFilterType, ...
                    ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                    'dataFilterLength', (detectParams.dataFilterLength * downSampleFreq), ...
                    ... length of data filter
                    'derFilterLength', (detectParams.derFilterLength * downSampleFreq), ...
                    ... length of derivative filter
                    'debugging', detectParams.debugging, ...
                    ... if set to 1 then debugging figures appear
                    'dataStart', detectParams.dataStart, ...
                    ... index of first data point
                    'forceDisplay', detectParams.forceDisplay,...
                    ... forces a graphical output even if other outputs are taken
                    'noFit', detectParams.noFit)
        disp('The amount of time this single trace took to detect events was: ');
        de_time = toc;
        disp([num2str(de_time/60),' minutes']);
        prompt_thres = 'Is this a good threshold (Y=1/N=0): ';
        thres = input(prompt_thres);
        while thres == 0
            close;
            tic
            prompt_newthres = 'Enter a new threshold: ';
            detectParams.derThresh = input(prompt_newthres);
            detectPSPs(filtered_data(:,1), PSPsDown, ...
                    'minAmp', detectParams.minAmp, ...
                    ... minimum allowable amplitude for alpha functions (in units of samples)
                    'maxAmp', detectParams.threshVal, ...
                    ... maximum allowable amplitude
                    'minTau', (detectParams.minTau * downSampleFreq), ...
                    ... minimum allowable tau for alpha functions (in units of samples)
                    'maxTau', (detectParams.maxTau * downSampleFreq), ...
                    ... maximum allowable tau
                    'minYOffset', detectParams.minYOffset, ...
                    ... minimum allowable yOffset for alpha functions (in units of mV)
                    'maxYOffset', detectParams.maxYOffset, ...
                    ... maximum allowable yOffset
                    'minDecay', (detectParams.minDecay * downSampleFreq), ...
                    ... minimum allowable decay tau
                    'maxDecay', (detectParams.maxDecay * downSampleFreq), ...
                    ... maximum allowable decay tau
                    'derThresh', (detectParams.derThresh), ...
                    ... threshold used to determine if the change of derivative is ...
                    ... great enough to separate out two alpha functions
                    'closestEPSPs', (detectParams.closestEPSPs * downSampleFreq), ...
                    ... second EPSP is not recognized if it is closer than this ...
                    ... value to the previous one (in units of samples)
                    'errThresh', (detectParams.errThresh), ...
                    ... threshold for standard error above which a fit with ...
                    ... multiple alphas is attempted
                    'dataFilterType', detectParams.dataFilterType, ...
                    ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                    'derFilterType', detectParams.derFilterType, ...
                    ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                    'dataFilterLength', (detectParams.dataFilterLength * downSampleFreq), ...
                    ... length of data filter
                    'derFilterLength', (detectParams.derFilterLength * downSampleFreq), ...
                    ... length of derivative filter
                    'debugging', detectParams.debugging, ...
                    ... if set to 1 then debugging figures appear
                    'dataStart', detectParams.dataStart, ...
                    ... index of first data point
                    'forceDisplay', detectParams.forceDisplay,...
                    ... forces a graphical output even if other outputs are taken
                    'noFit', detectParams.noFit)
            disp('The amount of time this single trace took to detect events was: ');
            de_time = toc;
            disp([num2str(de_time/60),' minutes']);
            prompt_thres = 'Is this a good threshold (Y=1/N=0): ';
            thres = input(prompt_thres);
        end

disp('The estimated amount of time this analysis will take is: ');
disp([num2str((de_time*num_traces)/60),' minutes']);
disp('Hit space to proceed');
pause();

%     delete(gcp);
%     parpool('local',4);
    
kk=1;
% h = waitbar(0,['Patience, I am analyzing trace ',num2str(kk),' of ',num2str(num_traces)]);
results = {}; 
mycluster = parcluster('local');
mycluster.NumWorkers = 4;
parpool('local',4);
tic
parfor kk = 1:num_traces
%     disp(num2str(kk));
    results{kk} = detectPSPs(filtered_data(:,kk), PSPsDown, ...
                    'minAmp', detectParams.minAmp, ...
                    ... minimum allowable amplitude for alpha functions (in units of samples)
                    'maxAmp', detectParams.threshVal, ...
                    ... maximum allowable amplitude
                    'minTau', (detectParams.minTau * downSampleFreq), ...
                    ... minimum allowable tau for alpha functions (in units of samples)
                    'maxTau', (detectParams.maxTau * downSampleFreq), ...
                    ... maximum allowable tau
                    'minYOffset', detectParams.minYOffset, ...
                    ... minimum allowable yOffset for alpha functions (in units of mV)
                    'maxYOffset', detectParams.maxYOffset, ...
                    ... maximum allowable yOffset
                    'minDecay', (detectParams.minDecay * downSampleFreq), ...
                    ... minimum allowable decay tau
                    'maxDecay', (detectParams.maxDecay * downSampleFreq), ...
                    ... maximum allowable decay tau
                    'derThresh', (detectParams.derThresh), ...
                    ... threshold used to determine if the change of derivative is ...
                    ... great enough to separate out two alpha functions
                    'closestEPSPs', (detectParams.closestEPSPs * downSampleFreq), ...
                    ... second EPSP is not recognized if it is closer than this ...
                    ... value to the previous one (in units of samples)
                    'errThresh', (detectParams.errThresh), ...
                    ... threshold for standard error above which a fit with ...
                    ... multiple alphas is attempted
                    'dataFilterType', detectParams.dataFilterType, ...
                    ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                    'derFilterType', detectParams.derFilterType, ...
                    ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                    'dataFilterLength', (detectParams.dataFilterLength * downSampleFreq), ...
                    ... length of data filter
                    'derFilterLength', (detectParams.derFilterLength * downSampleFreq), ...
                    ... length of derivative filter
                    'debugging', detectParams.debugging, ...
                    ... if set to 1 then debugging figures appear
                    'dataStart', detectParams.dataStart, ...
                    ... index of first data point
                    'forceDisplay', detectParams.forceDisplay,...
                    ... forces a graphical output even if other outputs are taken
                    'noFit', detectParams.noFit); 
%     waitbar(kk/num_traces,h,['Patience, I am analyzing trace ',num2str(kk),' of ',num2str(num_traces)]);
end
% delete(h);
delete(gcp);
disp('The total amount of time this analysis took was: ');
tot_time = toc;
disp([num2str(tot_time/60),' minutes']);  

savename = [geno,'_',exp_date,'_',treat,'_',cell_num];
    clearvarlist = ['clearvarlist';setdiff(who,{'results';'savename'})];
    clear(clearvarlist{:});
save(savename,'results');

disp('Recall that the output of the results is the following, in each cell: ');
disp('   Params(:,1) = amplitude of psp');
disp('   Params(:,2) = rise time');
disp('   Params(:,3) = x-offset of alpha function');
disp('   Params(:,4) = decay time');
disp(' ');
disp('NEXT CELL PLEASE! :D');

%% Amplitude Data Analysis

%% Frequency Data Analysis

close all
clc

num_traces = max(size(results));
jj = 1;
IEIs = {};

for j = 1:num_traces
    size_trace = max(size(results{j}));
    temp_freqs = [];
    if sum(results{j})~=0    
        for k = 1:size_trace
            if k == 1
                temp_freqs(k,1) = results{j}(k,3);
            else
                temp_freqs(k,1) = results{j}(k,3) - results{j}(k-1,3);
            end
        end
        IEIs{jj} = temp_freqs;
        clear temp_freqs
        jj=jj+1;
    else
    end 
end

IEIs = padcat(IEIs{:});
[l_trace num_traces] = size(IEIs);
IEIs_vec = reshape(IEIs,l_trace*num_traces,1);
IEIs_vec(isnan(IEIs_vec)) = [];
