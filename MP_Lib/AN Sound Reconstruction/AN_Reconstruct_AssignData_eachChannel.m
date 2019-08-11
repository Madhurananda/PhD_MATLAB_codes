function [new_zc_cell] = AN_Reconstruct_AssignData_eachChannel(an_data_files)

% This function converts the spike trains from accroding to sensitivity
% level to channel,
% For 'AN' data. 

% The data should be loaded from the file
big.data = load(an_data_files);

% Pick out some useful information
% How many channels?
num_channels    = big(1).data.AN.channels;
% How many sen levels?
num_sen_levels  = big(1).data.AN.iterations;
% What is the length of the original input signal data?
length_sig = big(1).data.AN.datalength;
% What is the threshhold level values?
thresh_levels = big(1).data.AN.thresh_levels;
% What is the delay vector for this sound signal?
delayVector = big(1).data.AN.delayVector;
% What is the frequency of the original input sound signal?
Fs = big(1).data.AN.fs ;
% What is the peak of the original signal?
PeakofSignal = big(1).data.AN.PeakofSignal;
% What are the center frequencies used for each channel by gammatone filterbank?
cochCFs = big(1).data.AN.cochCFs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This one is constructing a huge big cell array to store all the data.
% zc_cell = cell(num_sen_levels,num_channels);
% h = waitbar(0,'Please wait. First cell array is being created ...');
% for i = 1:num_sen_levels
%     for ch = 1:num_channels
%         zc_cell{i,ch}=AN_FindData(an_data_files,i,ch);
%     end
%     waitbar(i /num_sen_levels)
% end
% close(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for taking all the data from that big cell array and assigning them according to the number of channels 
new_zc_cell = cell(1,num_channels);
log = false;
h2 = waitbar(0,'Please wait. Second cell array is being created ...');
for ch = 1:num_channels
    for i = 1:num_sen_levels

        Time_Value = AN_FindData(an_data_files,(num_sen_levels-i+1),ch);
        %zc_cell{i,ch}=AN_FindData(an_data_files,i,ch);
        %Time_Value = zc_cell{(num_sen_levels-i+1),ch};
        Time_Value(Time_Value==0)=[];
        if (isempty(Time_Value))
            log = true;
            Copy_Time_Value = 0;
        else
            if i==1
                % Do nothing.
%                 if (isempty(Time_Value))
%                     Copy_Time_Value = 0;
%                 end
            else
%                 disp(ch)
%                 disp(i)
                Time_Value = setdiff(Time_Value,Copy_Time_Value);
                
            end
            levels = zeros(length(Time_Value),1);
            for j = 1:length(Time_Value)
                levels(j,1) = (num_sen_levels-i+1);
            end
            if (isempty(Time_Value))
                new_Time_Value = [];
            else
                new_Time_Value = cat(2,Time_Value,levels);
            end
            if i==1
                final_Time_Value = new_Time_Value;
            elseif log
                final_Time_Value = new_Time_Value;
                log = false;
            else
                final_Time_Value = cat(1,final_Time_Value,new_Time_Value);
            end
            % It is necessary to copy the data inside Time_Value at this
            % instance , so that it can be compared to the next Time_Value.
            if (isempty(final_Time_Value))
                
            else
                Copy_Time_Value = final_Time_Value(:,1);
            end
        end
    
    end
    final_Time_Value = sortrows(final_Time_Value,1);
    new_zc_cell{1,ch} = final_Time_Value;
    waitbar(ch/num_channels)
end
close(h2)




end
