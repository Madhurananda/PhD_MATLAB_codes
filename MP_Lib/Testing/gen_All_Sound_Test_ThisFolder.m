function gen_All_Sound_Test_ThisFolder( type )
% THERE SHOULD NOT BE ANY OTHER FOLDER INSIDE THIS FOLDER. This function regenerates the sounds for all files in a folder. 
%   This is to reconstruct the sounds by all our techniques: - 
%   'AN' Sound Reconstuction
%   'AN' and 'Onset' Sound Reconstruction
%   Thomas's spike technique. 
% The 'type' stands for : 
%   '1' --> Male (channel - 15 to 40)
%   '2' --> Female (channel - 25 to 40)
%   '3' --> Musical (channel - 25 to 45)

% An example call - gen_All_Sound_Test_ThisFolder(3)


dirlist = dir('.');
for i = 1:length(dirlist)
    dirlist(i);
end

% So, the number of files in this folder will be - 'i-2'
% Now run a loop to consider each file in this folder - 
for j=3:i % It starts with 3 , as 1 and 2 repressents . and ..
    sound = [dirlist(j).name];
    gen_All_Sound_Test( sound, type )
end

end

