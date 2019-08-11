function [ spf,zf,ef ] = Thomas_code_Work(sig_name_mono,Nb)
% This function is to make Thomas's code to work by calling a single
% function. 
%   Example call - Thomas_code_Work('testfile.wav', 16);

[y, fs]=audioread(sig_name_mono);
sound_length = length(y)/44100;
t=1/44100:1/44100:sound_length;

%% Everything underneath is Koickal's code

row_vector = zeros(2,length(y));
row_vector(1,:)=t;
row_vector(2,:)=y;

% row_vector = [t;y];

%Nb = 4;
[spf,zf,ef] = spikecoder_fixed(row_vector,Nb); %'zf' is the decoded signal from the spike code

% Find out how many spikes has been generated: 
% spikes = spf(2,:);
% real_spikes = spikes(spikes ~= 0);
real_spikes = numel(find(spf>0.5))+numel(find(spf<-0.5));
disp('');
disp(['The total number of spikes ->   ' num2str(length(real_spikes))]);
disp('');
disp('');
ch_time = 1/fs;
expo = abs(min(floor(log10(ch_time))));
mult = 1;
for j = 1:expo
    mult = mult*10;
end
new_sound_length = sound_length*mult;
bit_length = ceil(log2(new_sound_length));
total_bytes = (length(real_spikes)*bit_length)/(8*1024);
disp(' ')
disp(' ------------------------------------------------ ')
disp(['Bytes required for a ' num2str(sound_length) ' sec sound file - '])
disp([num2str(total_bytes), ' Kilobytes'])
disp(' ------------------------------------------------ ')
disp(' ')
disp(' ')

%% UPTO Here was Koickal's code. 

%disp(zf)
new_Name = [sig_name_mono(1:end-4) '_NEW_Thomas_'];
new_Name = [new_Name num2str(Nb)];
new_Name = [new_Name '.wav'];
%zf_new = zf(2,:)';
zf_new = zf(2,:)/max(abs(zf(2,:)));
zf_new = zf_new*0.9;
%disp(zf)
audiowrite(new_Name, zf_new, fs);
disp('Re-Construction complete for Thomas code.');
disp('');
disp(['The audio file has been saved to ' new_Name]);

end
