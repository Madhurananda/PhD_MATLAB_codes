% Program: 
%   AN_SigGen_Mono_MJN(sig_name_mono,spike_option,save_option,bm_option,norm_option,varargin)
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
function [AN,bmSig,output_filename] = AN_1_SigGen_Mono_EditMP(calling_mode,sig_name_mono,spike_option,save_option,bm_option,plot_option,norm_option, sen_levels, n_channels, varargin)
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
max_sig_length = 100;


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

% First we examine the calling_mode to check whether the parameters come
% from in here (run_mode = 0) or from a parameters file (run_mode = 1)
int_mode = strcmp(calling_mode,'internal');
ext_mode = strcmp(calling_mode,'external');
if int_mode == 1 
    run_mode = 0;
elseif ext_mode == 1
    run_mode = 1;    
else
    error('Error: incorrect assignment of calling mode, check your first input parameter (should be either "internal" or "external"')
end

% Now we either load variables from below, or we load them from a
% parameters file
switch run_mode
    
    % Use variables here
    case 0
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Declare filterbank variables %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [sig,fs] = audioread(sig_name_mono);
        
        % Low end of the filter bank (typical: 50Hz)
        f_low       = 100;
        %f_low       = 1;
        
        % High end of the filter bank (typical: 5000-10000Hz)
        f_high      = min(10000, 0.246875*fs);
        % f_high      = 10000;
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
        
        
    % Load an external parameters file, stored in the local directory
    case 1       
        load('ANparams.mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
        
        % The maximum value of the original signal is necessary to
        % normalise the reconstructed Signal
        PeakofSignal = max(abs(sig));
        length_sig     = length(sig);
        disp(' ')
        disp('Welcome to "AN_1_SigGen_Mono.m". This performs AN spike encoding in MONO. ')
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
        switch spike_option
            case 1
                disp(['Spike encoding using the algorithm by: MJN'])
            case 2
                disp(['Spike encoding using the algorithm by: LSS'])
        end
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
% disp('No. of bits per sample')
% disp(nbits)

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
    % The maximum value of the original signal is necessary to
    % normalise the reconstructed Signal
    PeakofSignal = max(abs(sig));
end ;

%wavwrite(sig,44100,16,'madhu.wav')

% Chop the signal up if required
if length_sig > max_sig_length*fs
    sig=sig(1:max_sig_length*fs);

    % Reassign the cut off signal length as the new signal length
    length_sig = max_sig_length*fs;
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('The maximum length of signal is set as 100 sec.')
    disp('The inputted signal is longer than 100 sec.')
    disp('So, the sound has been chopped off after 100 sec.')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
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

%[ind_zc_fin,spikes] = simple_spike_gen_multilev1(bmSig,n_channels,cochCFs,fs,length_sig,period_frac,sen_levels,sen_multiplier,min_level_zc);
%zc_all = simple_spike_gen_multilev3(bmSig,n_channels,cochCFs,fs,length_sig,period_frac,sen_levels,sen_multiplier,min_level_zc);

% Here is where we call the spike encoder, with a choice between Les and Mike's versions 
switch spike_option
    case 1
        % This is the one whic is edited by MP
        % Mike's spike encoder
        [zxstructure, assigned_Spikes, bmSig, threshold_level] = AN_SpikeGen_Mono_MJN_inNewSpikeForm_EditMP(bmSig,n_channels,cochCFs,fs,length_sig,period_frac,sen_levels,sen_multiplier,min_level_zc);
        
    case 2
        % Leslie's spike encoder
        %[zc_all,bmSig,threshold] = zxspiketrain2b_mono_MJNedit(bmSig, min_level_zc, sen_levels, sen_multiplier, length_sig, cochCFs, fs);
        [zxstructure,bmSig,threshold_level] = zxspiketrain2b_mono_MJNedit(bmSig, min_level_zc, sen_levels, sen_multiplier, length_sig, cochCFs, fs);
    otherwise
        error('Error: please choose a spike encoder')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save('bmsig_2', 'bmSig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Saving the data if required %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the final 'AN' structure with all the parameters and
        % analysed data
        AN.Spike_Assign_Sen_Level  = zxstructure;
        switch spike_option
            case 1
                AN.Spike_Assign_Channels = assigned_Spikes;
            case 2 
                % Do nothing
        end
        
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
        AN.fs           = fs;
%         AN.nbits        = nbits;
        AN.datalength   = length_sig;
        AN.delayVector  = delayVector;
        AN.norm_option  = norm_option;
%         AN.Spike_Gen_Option  = Spike_Gen_Option;
        %AN.t            = t;
        %AN.f            = f;
        AN.timestamp    = datestr(clock,30);
        % Make a note of which ZX encoder was used
        switch spike_option
            case 1                
                AN.encoder = 'MJN';
            case 2                
                AN.encoder = 'LSS';
        end


switch save_option
    case 0   
        % This controls the console output display
        switch run_mode
            case 0
                disp('Analysis complete, but no data file saved.')  
                disp(' ')
            case 1                        
                disp('Analysis complete, data saving handled by calling script.')  
                disp(' ')
        end
    case 1                                
        % Set up the output filename and save the data
        switch nargin
            case 9
                output_filename = [sig_name_mono(1:end-4) '_' num2str(sen_levels) '_' num2str(n_channels) '_AN.mat'];
            case 10
%                 output_filename = [sig_name_mono(1:end-4) '_' char(varargin(1)) '_' num2str(sen_levels) '_' num2str(n_channels) '_AN.mat'];
                output_filename = [sig_name_mono(1:end-4) '_' num2str(sen_levels) '_' num2str(n_channels) '_AN.mat'];
            otherwise
                output_filename = [sig_name_mono(1:end-4) '_' char(varargin(1)) '_' num2str(sen_levels) '_' num2str(n_channels) '_AN.mat'];
                disp('You have specified too many arguments, only the first one has been used to create the output filename.')
        end
        
        % Save the file and tell the user
        save(output_filename, 'AN')
        disp(' ')
        disp(['Analysis complete, AN data saved to: ' output_filename])
        disp(' ')
        
        % Now decide whether or not to save the bmSig in a separate file
        
        % First set up the BM structure
        BM.bmSig        = bmSig;
        BM.delayVector  = delayVector;
        BM.cochCFs      = cochCFs;
        BM.N_erbs       = N_erbs;
        BM.fs           = fs;
        BM.fmin         = f_low;
        BM.fmax         = f_high;
        BM.channels     = n_channels;
        
        switch bm_option
            case 0
                disp('The filtered (BM) signal was not saved')
                disp(' ')
            case 1                                
                % Now save it accordingly
                switch nargin
                    case 7
                        output_filename_BM = [sig_name_mono(1:end-4) '_bmSig.mat'];                        
                    case 8
                        output_filename_BM = [sig_name_mono(1:end-4) '_' char(varargin(1)) '_bmSig.mat'];                                                
                    otherwise
                        output_filename_BM = [sig_name_mono(1:end-4) '_' char(varargin(1)) '_bmSig.mat'];                        
                end                
                save(output_filename_BM,'BM')
                disp(['The filtered (BM) signal was saved into the file: ' output_filename_BM])
                disp(' ')
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=length_sig/fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIAL PLOTTING OF THE DATA, JUST FOR SHOW %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The size of the spike marks in the plots
mark_size   = 3;

% The sensitivity levels to display spikes for in the immediate plotting
disp1 = sen_levels;
disp2 = round(sen_levels/2);
disp3 = 1;
switch plot_option
    
    case 0
        
    case 1
        switch spike_option
            case 1
                figure(11)
                subplot(3,1,1)
                hold on
                plot(zxstructure(disp1).list(:,2),zxstructure(disp1).list(:,1),'.k','MarkerSize', mark_size)
                title(['Spike trains (MJN); Sensitivity level ' num2str(disp1) ' (RMS voltage ' num2str(threshold_level(disp1)) ')'])
                ylabel('Filter channel')
                xlabel('Time (s)')
                ylim([-1 n_channels])
                xlim([0 t(end)])
                
                subplot(3,1,2)
                hold on
                plot(zxstructure(disp2).list(:,2),zxstructure(disp2).list(:,1),'.k','MarkerSize', mark_size)
                title(['Spike trains; Sensitivity level ' num2str(disp2) ' (RMS voltage ' num2str(threshold_level(disp2)) ')'])
                ylabel('Filter channel')
                xlabel('Time (s)')
                ylim([-1 n_channels])
                xlim([0 t(end)])
                
                subplot(3,1,3)
                hold on
                plot(zxstructure(disp3).list(:,2),zxstructure(disp3).list(:,1),'.k','MarkerSize', mark_size)
                title(['Spike trains; Sensitivity level ' num2str(disp3) ' (RMS voltage ' num2str(threshold_level(disp3)) ')'])
                ylabel('Filter channel')
                xlabel('Time (s)')
                ylim([-1 n_channels])
                xlim([0 t(end)])
                
            case 2
                figure(11)
                subplot(3,1,1)
                hold on
                plot(zxstructure(disp1).list(:,2),zxstructure(disp1).list(:,1),'.k','MarkerSize', mark_size)
                title(['Spike trains (LSS); Sensitivity level ' num2str(disp1) ' (RMS voltage ' num2str(threshold_level(disp1)) ')'])
                ylabel('Filter channel')
                xlabel('Time (s)')
                ylim([-1 n_channels])
                xlim([0 t(end)])
                
                subplot(3,1,2)
                hold on
                plot(zxstructure(disp2).list(:,2),zxstructure(disp2).list(:,1),'.k','MarkerSize', mark_size)
                title(['Spike trains; Sensitivity level ' num2str(disp2) ' (RMS voltage ' num2str(threshold_level(disp2)) ')'])
                ylabel('Filter channel')
                xlabel('Time (s)')
                ylim([-1 n_channels])
                xlim([0 t(end)])
                
                subplot(3,1,3)
                hold on
                plot(zxstructure(disp3).list(:,2),zxstructure(disp3).list(:,1),'.k','MarkerSize', mark_size)
                title(['Spike trains; Sensitivity level ' num2str(disp3) ' (RMS voltage ' num2str(threshold_level(disp3)) ')'])
                ylabel('Filter channel')
                xlabel('Time (s)')
                ylim([-1 n_channels])
                xlim([0 t(end)])
                
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






