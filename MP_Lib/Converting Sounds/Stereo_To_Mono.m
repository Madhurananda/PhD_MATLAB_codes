function Stereo_To_Mono(input_Name,output_Name,channel_used)
%Stereo_To_Mono: Summary of this function goes here
%   This function makes a stereo sound file to mono sound file.
%   First parameter: The name of the file you want to convert to Mono (including '.wav')
%   Second parameter: The new name of the converted Mono file 
%   Third parameter: You need to tell that Left or Right channel to be used, if you wish. Or, the left channel will be used by default. 

% First it is necessary to check that the inputted sound file is wav file or not.
str_check = input_Name(end-2:end);
switch strcmp(str_check,'wav')
    case 1
        disp('The inputed file is "wav" file');
    case 0
        disp('The inputed file is not "wav" file');
end

% Now read the data from the sound file.
[sig, Fs] = audioread(input_Name);




% Now do the extraction job.
switch (size(sig, 2) ~= 1)
    case 1
        if nargin ==3
            int_mode = strcmp(channel_used,'Left');
            ext_mode = strcmp(channel_used,'Right');
            if int_mode == 1 
                sig = sig(:,1) ;
            elseif ext_mode == 1
                sig = sig(:,2) ;    
            else
                 error('Error: incorrect assignment of channel mode, check your input parameter (should be either "Left" or "Right"')
            end
        end
        if nargin <=3
            sig = sig(:,1) ;
        end
        disp('The sound has been saved as Mono');
        
        
        
        
    case 0
        disp('Input sound file is already mono: No action has been taken');
end

% Get the maximum value of the sound - 
max_val_sound = max(abs(sig));
MAX_VAL = 0.9; 
sig = sig*(MAX_VAL/max_val_sound);
disp('The amplitude has been normalised to 0.9. ')

% The extracted signal is being written in the new file.
%output_Name = strcat(output_Name, '.wav');
audiowrite(output_Name, sig, Fs);
disp(['The file has saved as ' output_Name])


end