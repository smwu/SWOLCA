% par_var saves an object to a given file name within a parallel for loop
% Inputs:
%   file_name: name of vile
%   object: name of object to be saved
% Outputs: saves file in the current directory
function par_save(file_name, object)
    save(file_name, 'object');
end