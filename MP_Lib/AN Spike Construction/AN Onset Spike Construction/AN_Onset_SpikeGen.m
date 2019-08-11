% Program: 
%   AN_Onset_SpikeGen(sig_name_mono,spike_option,save_option,bm_option,norm_option,varargin)
%
% This is the main parameter setup and scripting program for obtaining
% simulations of AN spikes.
%
% ----------------------------
%
% Parameters to use:
%
%
%       calling_mode : 'internal' to use the internal parameters stored
%                      within THIS PROGRAM (e.g. when running a spike
%                      encoding on a single file)
%                       'external' when this program is being called BY
%                       ANOTHER SCRIPTING FILE (usually 'AN_0_SigGen_Script.m')
%
%       sig_name_mono: input sound file name (use a single MONO, .wav file)
%       spike_option : 1 for Mik, 2 for Les (spike encoding program)
%       save_option  : 0 for no data saved, 1 for AN data saved to .mat files
%       bm_option    : 0 for no saving of the BM signal (i.e. the filtered
%                      signal), 1 for it to be saved (in a separate file)
%       plot_option  : 1 for immediate plotting, 0 for no plotting
%       norm_option  : 0 does nothing (no alteration to signal level) 
%                      1 normalises the signal to itself (done to try to
%                      make the analysis level independant).
%       sen_levels   : Input the sensitivity levels
%       n_channels   : Input number of channels
%       SigToRecon_Option:  This is the option to re-assign the spikes to
%            make the reconstruction easy. 0 = No Re-assign; 1 = Re-assign
%       varargin(1)  : any string to be appended to saved data filename
%
% ----------------------------
%
% To run the program:
%       1) First place the .WAV file to be studied in an appropriate 
%       directory (include today's date)  
%       2) Navigate to the directory
%       3) Run this program with appropriate arguments, and having chosen
%       the appropriate parameters within the code (filterbank, ZX detection)
%       4) Output files will also be placed in same directory.
%
% This is a MONO version, requiring a MONO input .WAV file (best to use
% 16bit, 44100Hz). 
%       --> If you supply a STERO file, it will use only the LEFT channel
%
% ----------------------------
%
% MOST CRITICAL PARAMETERS in the program itself: 
%                           n_channels (number of filter channels)
%                           f_low (low limit of filterbank)
%                           f_high (low limit of filterbank)
%                           N_erbs (bandwidth of the filter channels, usually set to 1)
%                           minlevel_zc (minimum signal level at which to produce a spike)
%                           sen_levels (number of sensitivity levels at which to generate spikes)
%
% ----------------------------
%
% This program directly calls: 
%           AN_SpikeGen_Mono_MJN.m 
%           or 
%           zxspiketrain2b_mono_MJNedit
%
% This program indirectly calls: 
%           MakeErbCFs.m        : program to prepare the filterbank frequencies
%           gammatone1.m        : gammatone filter bank
%
% ----------------------------

%function AN_1_SigGen_Mono(calling_mode,sig_name_mono,spike_option,save_option,bm_option,plot_option,norm_option,varargin)
% function [AN,bmSig,output_filename] = AN_Onset_SpikeGen(calling_mode,sig_name_mono,spike_option,save_option,bm_option,plot_option,norm_option, sen_levels, n_channels, Spike_Gen_Option, varargin)


function [AN,bmSig,output_filename] = AN_Onset_SpikeGen(sig_name_mono, norm_option, sen_levels, n_channels, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a simplified auditory processing program.
% It exists to help me develop the ideas of the spike encoding and become
% proficient at writing this kind of software.
%
% It calls several other .m files, most importantly:
%
% 
%
%
% Various other customisable parameters, see 'Declare variables' in the
% code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Limit on total signal length in seconds (including start/pre-transient)
%       --> This will brutally chop off all data after this limit
max_sig_length = 10;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Make sure enough input args have been supplied
% switch nargin
%     case {0,1,2,3}
%         error('Error: not enough input arguments')
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Declare variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Declare filterbank variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Low end of the filter bank (typical: 50Hz)
f_low       = 100;
%f_low       = 1;

% High end of the filter bank (typical: 5000-10000Hz)
f_high      = 10000;
%f_high      = 100;

% Number of channels in the filter bank
% n_channels  = 50;

% Q-value of the channels (typically 1)
N_erbs      = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Declare spike encoder variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the fraction of the cycle over which we will examine the signal
period_frac = 0.25;

% The number of sensitivity levels
% sen_levels  = 16;


% This is the minimum signal threshold level for the zero crossings
% was 0.00014 for Gamma filter
min_level_zc = 0.0002;
% min_level_zc = 0.0005;

% Now for the less sensitivity level, the sen_multiplier shold change
% so that the range of threshold levels, does not change.
% It is inspected that the highest threshold level is about 0.0362
% So,
max_threshold = 0.03620386719651;
sen_multiplier = log10(max_threshold/min_level_zc);
sen_multiplier = sen_multiplier/(sen_levels-1);
sen_multiplier = 10^(sen_multiplier);
disp(' ')
disp('The sensitivity multiplier is - ')
disp(sen_multiplier)
% This is the multiplier between the sensitivity levels (use sqrt(2) for
% 3dB differences, or just 2 for 6dB)
% sen_multiplier = sqrt(2);
% sen_multiplier = 2;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main program starts below here: BE VERY CAREFUL EDITING THIS CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reading in and initial processing of the sound file %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we check the input sound file is a WAV file
str_check = sig_name_mono(end-2:end);
switch strcmp(str_check,'wav')
    case 1
        % Read the wav file into a vector. Remember it will be normalised so the
        % peak amplitudes are 1. Also read the sample rate and bit depth.
        [sig,fs,nbits] = wavread(sig_name_mono);
        % The maximum value of the original signal is necessary to
        % normalise the reconstructed Signal
        PeakofSignal = max(abs(sig));
        length_sig     = length(sig);
        disp(' ')
        disp('Welcome to construct the AN spikes. This performs AN spike encoding in MONO. ')
        disp(' ')
        disp(['Reading wave file: ' sig_name_mono])
        disp(' ')
        
        switch norm_option
            case 0
                disp('You have chosen NOT TO NORMALISE your wave file before spike encoding.')
            case 1
                disp('You have chosen TO NORMALISE your wave file before spike encoding.')
        end
        
        disp(' ')
%         switch spike_option
%             case 1
%                 disp(['Spike encoding using the algorithm by: MJN'])
%             case 2
%                 disp(['Spike encoding using the algorithm by: LSS'])
%         end
        disp(' ')
        disp(['Filtering into: ' num2str(n_channels) ' channels'])
        disp(' ')
        disp(['Over frequency range: ' num2str(f_low) ' - ' num2str(f_high) 'Hz'])
        disp(' ')
        disp(['ZX detection at: ' num2str(sen_levels) ' sensitivity levels'])
        disp(' ')
        disp(['With a minimum (starting) sensitivity level of: ' num2str(min_level_zc)])
        disp(' ')
        disp('... ')
        disp('... ')
        disp('... ')
    case 0
        error('The input file is not a .WAV file: Please convert it to .WAV and try again')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Madhu's Code -
%disp(' ')
%disp('XXXXXXXXX')
%disp('XXXXXXXXX')
%disp('XXXXXXXXX')
%disp(sig)
disp('This is the maximum value in the data')
disp(max(sig))
disp('This is the length of the data')
disp(length_sig)
disp('XXXXXXXXX')
disp('XXXXXXXXX')
disp('XXXXXXXXX')
disp('Rate in Hz is - ')
disp(fs)
disp('XXXXXXXXX')
disp('XXXXXXXXX')
disp('XXXXXXXXX')
disp('No. of bits per sample')
disp(nbits)

%plot(sig)
%grid
disp('XXXXXXXXX')
disp('XXXXXXXXX')
disp('XXXXXXXXX')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Ensure the loaded sound file is mono
if (size(sig, 2) ~= 1)
    disp('Input sound file is not mono: left (chan1) used') ;
    sig = sig(:,1) ;
end ;

%wavwrite(sig,44100,16,'madhu.wav')

% Chop the signal up if required
if length_sig > max_sig_length*fs
    sig=sig(1:max_sig_length*fs);

    % Reassign the cut off signal length as the new signal length
    length_sig = max_sig_length*fs;
end


% If specified, normalise the signal to itself
switch norm_option
    case 0
    case 1
        sig = sig./max(sig);
end


% Create a real-world time vector: note that this runs from 1/fs (i.e.
% starts just AFTER t=0s.
%t     = [1/fs:1/fs:(length_sig/fs)]';

% Create a frequency scale
%f     = linspace(0,fs,length(t));

% Create an empty matrix in which to place the filtered signals (i.e. the
% simulated BM output)
bmSig = zeros(n_channels,length_sig) ; % initialise mono  signal

% Create an empty vector to store the instantaneous frequency outputs from
% the filter (one set for each channel, length same as bmSig)
instf = zeros(n_channels,length_sig) ;

% Create an empty vector to store the filter envelope for each frequency
% (one set for each channel, length same as bmSig)
env   = zeros(n_channels,length_sig) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare the filterbank %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we generate the centre frequencies using the user input data above
% ERB: equivalent rectangular bandwidth (effective bandwidth of a
% square-shaped bandpass filter, centred at the centre frequency of the
% band)
cochCFs     = MakeErbCFs(f_low,f_high,n_channels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp('This is the cochCFs')
%disp(cochCFs(49))
%plot(cochCFs)
%wavwrite(sig,'madhu.wav')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call the gammatone filterbank and return the simulated BM signal %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call gammatone1.m, get back the filtered basilar membrane signal 'bmSig' 
% (each channel should be the same length as the input signal 'sig').
% Also, call the program which calculates the filter delay for each channel
h = waitbar(0,'Please wait. The filterbank signal is being created ...');
for gam_i = 1:n_channels
    [bmSig(gam_i,:),env(gam_i,:),instf(gam_i,:)] = gammatone1(sig',fs,cochCFs(gam_i),N_erbs);    
    % Calculate the channel delay in seconds
    delayVector(gam_i) = gammatoneDelay(cochCFs(gam_i),fs,N_erbs) * (1/fs);
    waitbar(gam_i/n_channels)
end
close(h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save('bmsig_1', 'bmSig');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Now we call the spike generating program %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the detection function or not ...
% 1 - Use it
% 0 - Do not use it
use_DetectFunc = 0;

% Here is where we call the spike encoder, with a choice between Les and Mike's versions 
[zxstructure, assigned_Spikes, bmSig, threshold_level,RMS_Values] = Onset_Spike_Generation(bmSig,n_channels,cochCFs,fs,length_sig,period_frac,sen_levels,sen_multiplier,min_level_zc,use_DetectFunc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Saving the data if required %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the final 'AN' structure with all the parameters and analysed data
AN.Spike_Assign_Sen_Level  = zxstructure;
AN.Spike_Assign_Channels = assigned_Spikes;
AN.cochCFs      = cochCFs;
AN.N_erbs       = N_erbs;
AN.soundlength  = (length_sig/fs);
AN.PeakofSignal = PeakofSignal;
AN.fmin         = f_low;
AN.fmax         = f_high;
AN.channels     = n_channels;
AN.iterations   = sen_levels;
AN.multiplier   = sen_multiplier;
AN.minlevel_zc  = min_level_zc; 
AN.thresh_levels= threshold_level;
AN.RMS_Values= RMS_Values;
AN.fs           = fs;
AN.nbits        = nbits;
AN.datalength   = length_sig;
AN.delayVector  = delayVector;
AN.norm_option  = norm_option;
AN.timestamp    = datestr(clock,30);
        
switch use_DetectFunc
    % Use it
    case 1
        output_filename = [sig_name_mono(1:end-4) '_' num2str(sen_levels) '_' num2str(n_channels) '_OnsetByDetectFunc_AN.mat'];
    % Do not use it
    case 0
        output_filename = [sig_name_mono(1:end-4) '_' num2str(sen_levels) '_' num2str(n_channels) '_Onset_AN.mat'];
end

% Save the file and tell the user
save(output_filename, 'AN')
disp(' ')
disp(['Analysis complete, AN data saved to: ' output_filename])
disp(' ')


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
