function [zc_struct,assigned_Spikes,bmSig,threshold_level,RMS_Values] = Onset_Spike_Generation(bmSig,n_channels,cochCFs,fs,length_sig,period_frac,sen_levels,sen_multiplier,min_level_zc,use_DetectFunc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Onset generation has been in two way. One by using the 'Onset Detection
% Function' [use_DetectFunc = 1] , and another is by just simply considering the very first
% occurrence of a continuous train of spikes [use_DetectFunc = 0].  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program has been Edited to Assign the Generated Spikes so that the
% reconstruction work has less latency.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program performs a simplified spike encoding of a simulated BM
% signal (usually this signal is created with a gammatone filterbank).
%
% The gammatone bank outputs multiple filter channels (parallel output, 
% 'n_channels'), each of which is analysed by this program to produce a 
% train of spikes.
%   -->     bmSig(channel,data) where 'channel' is the channel number and 
%           'data' is the data from that channel's filter 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First we work out what the number of samples for the 1/4 cycle (in which 
% to check the sensitivity levels and decide whether to produce a spike) 
% will be for each filter channel

% First create an empty vector to store the sample values
cyc_samples = zeros(n_channels,1);

% Now loop over the number of filterbank channels
for i_cyc = 1:n_channels
    % There should be 'n_channels' number of sample values, one the 1/4
    % period for each filter channel
    %cyc_samples(i_cyc) = floor(period_frac*((1/f(i_cyc))*fs));
    cyc_samples(i_cyc) = floor(period_frac*((1/cochCFs(i_cyc))*fs));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Create an empty vector, then fill in the lowest value threshold
threshold_level = zeros(sen_levels,1);
threshold_level(1) = min_level_zc;

% Loop over the number of sensitivity levels
for i_sen = 2:sen_levels
    % Increment each sensitivity level as the previous one X the
    % multiplier
    threshold_level(i_sen) = threshold_level(i_sen-1)*sen_multiplier;
end
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Detection Function has been applied between the Filter processed
% signal and Onset Localization. 
% The function is the No. 1 & 2 equation in "Trends in Onset Detection" Paper
% [http://delivery.acm.org/10.1145/2020000/2016736/p75-rosao.pdf?ip=139.153.254.154&acc=ACTIVE%20SERVICE&CFID=152750553&CFTOKEN=82060587&__acm__=1354641304_e6d946ec0b84a31662c072062c1a4cd4]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch use_DetectFunc
    
    % Use Detection Function
    case 1
        N = 8;
        E_n = zeros(n_channels, length_sig);

        for ch = 1:n_channels
            for n = 1:length_sig
                e_n = 0;
                for i = (-(N/2)):((N/2)-1)
                    if ((n+i)>=1 && (n+i)<=length_sig)
                        % Eq. 1
        %                 e_n = e_n + (bmSig(1,(n+i)))^2;
                        % Eq. 2
                        e_n = e_n + (bmSig(ch,(n+i)));
                    end
                end

                E_n(ch, n) = e_n;

            end
        end
        bmSig = E_n;
        
    % Do not use Detection Function
    case 0
        % Do nothing %
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This one generates spike in the convenient form for reconstruction.
% This one assigns Spike according to each channel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we prepare the threshold levels to be used for the various
% sensitivities

% Now commence the big loop over each sensitivity level
h1 = waitbar(0,'Assigning Spikes according to Channels ......');

% To set the onset finder on or off
onset_Finder = true;
RMS_Values = zeros(n_channels, length_sig);

% Now we do the spike encoding, once for each sensitivity level
for ch = 1:n_channels

    % Create a cell array in which to store spikes for each channel
    zc_cell = cell(n_channels,1);

    % Create a counter to allow the index to run
    nn=0;
    ind_zc=[];
    % Now run the loop over the whole signal (up to the second to last
    % point to avoid trying to look at a data point that doesn't exist)
    for n = 1:length_sig-1
        % If the current sig level is less than zero and the next sample is
        % greater than zero... it's a zero crossing (count it at the nth
        % sample)
        if (bmSig(ch,n) < 0) && (bmSig(ch,n+1) >= 0)

            % Assign an initial 0 value to the to-be-extracted portion
            % of the signal
            sig_bit=0;

            % Check the zero-crossing isn't too close to the start of
            % the signal (i.e. if there are enough samples over the
            % previous quarter cycle...)
            if (n+1) > cyc_samples(ch)                    
                % Pick out the previous 1/4 cycle (absolute value)...
                %sig_bit = bmSig(ch,n-cyc_samples(ch):n);
                sig_bit = bmSig(ch,n+1-cyc_samples(ch):n+1);
            end

%                 plot(sig_bit)

            % ... and compute the RMS value
            RMS = sqrt(sum(sig_bit.^2)/cyc_samples(ch));
            RMS_Values(ch,n) = RMS;
            
            % Now if the RMS value is greater than any threshold_levels
            % then a spike will be considered to be occurred. 
            spikeLavel = 0;
            for i_thresh = 1:sen_levels
                if RMS >= threshold_level(i_thresh)
                    spikeLavel = spikeLavel+1;
                end
            end


            % Now we search through each sensitivity level, starting with the most sensitive
            if spikeLavel>0
                if (onset_Finder)
                    % If RMS value greater than any threshold level, increment the counter 'nn' by 1...
                    nn=nn+1;

    %                     ind_zc(nn,1) = ch;
    %                     ind_zc(nn,2) = (n+1)/fs;
    %                     ind_zc(nn,3) = i_thresh;

    %                     ind_zc(nn,1) = ch;
                    ind_zc(nn,1) = (n+1)/fs;
                    ind_zc(nn,2) = spikeLavel;
                end
                onset_Finder = false;
            else
                onset_Finder = true;
            end

        end   
    end

    assigned_Spikes{ch} = ind_zc;

    % zc_struct is not usefull here.
%         zc_struct(ch).list = ind_zc;


    waitbar(ch/n_channels)
end
close(h1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is one is the original code writen by Michel Newton, edited by me.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now commence the big loop over each sensitivity level
h2 = waitbar(0,'Assigning Spikes according to Sensitivity Levels ...');

% To set the onset finder on or off
onset_Finder = true;

% Now we do the spike encoding, once for each sensitivity level
for i_thresh = 1:sen_levels

    % Create a cell array in which to store spikes for each channel
    zc_cell = cell(n_channels,1);

    for ch = 1:n_channels

        % Create a counter to allow the index to run
        nn=0;
        ind_zc=[];
        % Now run the loop over the whole signal (up to the second to last
        % point to avoid trying to look at a data point that doesn't exist)
        for n = 1:length_sig-1

            % If the current sig level is less than zero and the next sample is
            % greater than zero... it's a zero crossing (count it at the nth
            % sample)
%             if (bmSig(ch,n) < 0) && (bmSig(ch,n+1) >= 0)
            if (bmSig(ch,n) < 0) && (bmSig(ch,n+1) >= 0)

                % Assign an initial 0 value to the to-be-extracted portion
                % of the signal
                sig_bit=0;

                % Check the zero-crossing isn't too close to the start of
                % the signal (i.e. if there are enough samples over the
                % previous quarter cycle...)
                if (n+1) > cyc_samples(ch)                    
                    % Pick out the previous 1/4 cycle (absolute value)...
                    %sig_bit = bmSig(ch,n-cyc_samples(ch):n);
%                     sig_bit = bmSig(ch,n+1-cyc_samples(ch):n+1);
                    sig_bit = bmSig(ch,n+1-cyc_samples(ch):n+1);
                end

                % ... and compute the RMS value
                RMS = sqrt(sum(sig_bit.^2)/cyc_samples(ch));
                RMS_Values(ch,n) = RMS;

                % Now we search through each sensitivity level, starting with the most sensitive
                if RMS >= threshold_level(i_thresh)
                    
                    % and then put the index of the zero crossing into the zc vector
                    %ind_zc.(['level' num2str(i_thresh)])(nn,1) = n;
                    %ind_zc(nn,1,i_thresh) = n;
                    %ind_zc(nn,2,i_thresh) = ch;
                    
                    if (onset_Finder)
                        % If RMS value greater than any threshold level, increment the counter 'nn' by 1...
                        nn=nn+1;
                        
                        ind_zc(nn,1) = ch;
                        ind_zc(nn,2) = (n+1)/fs;
                        ind_zc(nn,3) = i_thresh;
                    end
                    onset_Finder = false;
                else
                    onset_Finder = true;
                end

            end    

        end

        zc_cell{ch} = ind_zc;        
        zc_all1 = cat(1,zc_cell{:});

        % Here we sort the rows according to spike time, but only if spikes
        % exist
        switch isempty(zc_all1)
            case 0
                zc_all = sortrows(zc_all1,2);
            case 1                
                zc_all = zc_all1;
        end

        % Here we put the data from each sensitivity level into a structure
        % field
        zc_struct(i_thresh).list = zc_all;

%         % assigned_Spikes is not usefull here
%         assigned_Spikes{i_thresh} = zc_all;

        % Clear variables that run on each iteration of the sensitivity
        % level
        clear n nn sig_bit lev ind_zc zc_all1 zc_all
    end
    waitbar(i_thresh/sen_levels)
end

close(h2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end
