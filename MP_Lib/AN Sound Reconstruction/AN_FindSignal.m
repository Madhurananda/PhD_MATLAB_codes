function [Y] = AN_FindSignal(an_data_files,sen_levels,channels)
% AN_Reconstruct is re-constructing the input signal on the basis of it's
% produced spike.
%   Detailed explanation goes here

% How many files?
num_files       = length(an_data_files);

% Now loop over the number of files, loading them into structures one by one
for i=1:num_files
    big(i).data = load(char(an_data_files(i)));        
end


% Pick out some useful information
% How many channels?
num_channels    = big(1).data.AN.channels;
% How many sen levels?
num_sen_levels  = big(1).data.AN.iterations;
length_sig = big(1).data.AN.datalength;

% Take the data from the big structure
AN_data = big(1).data.AN.zxstructure;

% Extract the data for given sen_level
sen_AN_data = AN_data(sen_levels);
channel = sen_AN_data.list(:,1);
time_value = sen_AN_data.list(:,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is to extract the right data at right level and channel 
j=1;
for i = 1:length(channel)
    if channel(i) == channels
        %disp('bingo');
        new_time_value(j) = sen_AN_data.list(i,2);
        j=j+1;
    end
end
new_time_value = new_time_value';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is to plot all the spikes at the given channel and channel.
m=1;
for i = 1:length(new_time_value)
    plot_value(m,1) = new_time_value(i);
    plot_value(m,2) = 0;
    m=m+1;
end
figure(1)
plot(plot_value(:,1),plot_value(:,2),'or')
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is to generate the sine curve throgh all the spikes
% By default it is to be assuemd that the frequency of output file will be 44100 Hz
Fs=44100;
h = waitbar(0,'Please wait untill the sine curves are produced...');
%for i = 1:5
for i = 1:length(new_time_value)-1
    x=new_time_value(i):(1/Fs):new_time_value(i+1);
    z=0:(2*pi)/(length(x)-1):(2*pi);
    y=sin(z);
    if (i==1)
        temp_length = new_time_value(i+1)-new_time_value(i);
        temp_x = (new_time_value(i)-temp_length):(1/Fs):new_time_value(i);
        temp_z=0:(2*pi)/(length(temp_x)-1):(2*pi);
        temp_y = sin(temp_z);
        
        X=horzcat(temp_x,x);
        Y=horzcat(temp_y,y);
    else
        X=horzcat(X,x);
        Y=horzcat(Y,y);
    end
    waitbar(i /(length(new_time_value)-1))
end
close(h)
figure(1)
plot(X,Y)
hold on
%grid on
% This is to keep similarity with the audio file reading output
Y=Y';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is to plot the BMSig (Basilar Membrane Signal)
% Take the data from the big structure
BM_data =  load('DbsA1opn_4.57sec_Mono_Right_bmSig.mat');
BM_data = BM_data.BM.bmSig;
BM_data = BM_data(channels,:);
BM_data_axis = 0:(1/Fs):((length_sig/Fs)-(1/Fs));
figure(1)
plot(BM_data_axis,BM_data,'.y')
% hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Create a cell array in which to store spikes for each 

%for ch = 1:num_channels
%     for i = 1:num_sen_levels
%         ind_channel = AN_data(i).list(:,1);
%         j=1;
%         ind_zc=cell(num_channels,1);
%         %ind_timeValue=0;
%         
%         for ch1 = 1:num_channels
% %             disp(ind_channel(1));
%             if (ind_channel(ch1) == ch1)
%                 ind_timeValue(j) = AN_data(i).list(i,2);
%                 j=j+1;
%                 %disp('bingo');
%                 
%             end
%             
%             
% %             if ch1==1
% %                 ind_zc{ch1} = ind_timeValue;
% %             else
% %                 ind_zc = vertcat(ch1,ind_zc,ind_timeValue);
% %             end
%         end
%         disp(ind_timeValue);
%         ind_zc{ch1} = ind_timeValue;
%         
% %         if i==num_sen_levels
% %             if (ind_channel(i) == ch)
% %                 ind_timeValue(j) = AN_data(i).list(i,2);
% %                 j=j+1;
% %             end
% %             for k=(num_sen_levels+1):num_channels
% %                 if (ind_channel(k) == ch)
% %                     ind_timeValue(j) = AN_data(i).list(i,2);
% %                     j=j+1;
% %                 end
% %             end
% %         end
%     end
    %zc_cell{i} = ind_zc;
%end
    
    
%for i = 1:num_sen_levels
    %for ch = 1:num_channels
        %ind_zc = zeros(2,length(AN_data(i).list(:,1)));
%         ind_zc1 = AN_data(i).list(:,1);
%         ind_zc2 = AN_data(i).list(:,2);
        %zc_cell{ch} = cat(2,ind_zc1,ind_zc2);
    %end
    %zc_cell{i} = cat(2,ind_zc1,ind_zc2);
    %zc_cell{i} = cat(2,zc_cell{:});
    %disp('bingo');
% end 
% Gone in vain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% zc_cell = cell(num_sen_levels,num_channels);
% for i = 1:num_sen_levels
%     ind_ch = AN_data(i).list(:,1);
%     ind_time = AN_data(i).list(:,2);
%     j=1;
%     for ch = 1:num_channels
%         if channel(i) == channels
%             %disp('bingo');
%             new_time_value(j) = ind_time.list(i,2);
%             j=j+1;
%         end
%     end
% end






end
