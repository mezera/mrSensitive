function rootPath=mrSensitiveRootPath()
% Return root of the mrQ directory
% 
%        rootPath = mrqRootPath;
%
% This function MUST reside in the directory at the base of the mrSensitive
% directory structure
%
% Wandell & Mezer Copyright Vistasoft Team, 2013

rootPath = fileparts(which(mfilename));

return
