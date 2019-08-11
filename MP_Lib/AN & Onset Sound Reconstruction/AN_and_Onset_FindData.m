function [ new_time_value ] = AN_and_Onset_FindData( Onset_data_files,sen_levels,channels)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% %disp(an_data_files);
% % How many files?
% num_files       = length(an_data_files);
% %disp(num_files);
% 
% % Now loop over the number of files, loading them into structures one by one
% for i=1:num_files
%     big(i).data = load(char(an_data_files(i)));        
% end

big.data = load(Onset_data_files); 

% % Pick out some useful information
% % How many channels?
% num_channels    = big(1).data.AN.channels;
% % How many sen levels?
% num_sen_levels  = big(1).data.AN.iterations;
% length_sig = big(1).data.AN.datalength;

% Take the data from the big structure
AN_data = big(1).data.leftonset_wide;
%disp(AN_data(1))

% Extract the data for given sen_level
sen_AN_data = AN_data(sen_levels);
%disp(sen_AN_data)
if (~isempty(sen_AN_data.list))
    channel = sen_AN_data.list(:,1);
    time_value = sen_AN_data.list(:,2);

    new_time_value = zeros(length(time_value),1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is to extract the right data at right level and channel 
    j=1;
    for i = 1:length(channel)
        if channel(i) == channels
            %disp('bingo');
            new_time_value(j,1) = sen_AN_data.list(i,2);
            j=j+1;
        end
    end

new_time_value = new_time_value(1:j-1) ;
else
    new_time_value = zeros(1,1);
end


end

