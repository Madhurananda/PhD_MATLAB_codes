function [zc_struct,bmSig,threshold_level] = AN_SpikeGen_Mono_MJN_EditMP(bmSig,n_channels,cochCFs,fs,length_sig,period_frac,sen_levels,sen_multiplier,min_level_zc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we prepare the threshold levels to be used for the various
% sensitivities

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

%plot(threshold_level)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First create a big struct to store everything
for ii = 1:sen_levels
    zc_struct(ii).list = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now commence the big loop over each sensitivity level
h1 = waitbar(0,'Finding zero crossings...');
% Now we do the spike encoding, once for each sensitivity level
for i_thresh = 1:sen_levels
    
    % Create a cell array in which to store spikes for each channel
    zc_cell = cell(n_channels,1);
    
    for ch = 1:n_channels
        
        % Initialise the indexing vector
        %ind_zc = [];
        %ind_zc = zeros(length_sig,2);
        
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
                
                % ... and compute the RMS value
                RMS = sqrt(sum(sig_bit.^2)/cyc_samples(ch));
                
                % Now we search through each sensitivity level, starting with the most sensitive
                if RMS >= threshold_level(i_thresh)
                    
                    % If RMS value greater than any threshold level, increment the counter 'nn' by 1...
                    nn=nn+1;
                    % and then put the index of the zero crossing into the zc vector
                    %ind_zc.(['level' num2str(i_thresh)])(nn,1) = n;
                    %ind_zc(nn,1,i_thresh) = n;
                    %ind_zc(nn,2,i_thresh) = ch;
                    
                    ind_zc(nn,1) = ch;
                    ind_zc(nn,2) = (n+1)/fs;
                    ind_zc(nn,3) = i_thresh;
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
        
        % Clear variables that run on each iteration of the sensitivity
        % level
        clear n nn sig_bit lev ind_zc zc_all1 zc_all
    end
    waitbar(i_thresh/sen_levels)
end




close(h1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







end
        
        
        