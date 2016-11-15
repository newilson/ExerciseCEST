function fulldestdir = dicom_sort(dicomDir,destdir)
%
% Input:
%   -dicomDir = full path to dicom directory, with all dicoms from
%   different acquisition series in the same directory. Assumes <dicomDir>
%   contains only dicom files.
%
%   note: in some cases, the '.dcm' extension has been removed.  This
%   function will add '.dcm' back to those files
%
% Output:
%  - creates several directories within the <dicomDir>, based on the
%  protocol name and series found in the dicom headers.  Individual dicoms
%  are moved from the dicomDir into their respective protocol and series
%  directory.
%
% Written by Andrew S Bock Feb 2014
%
% Modified by NW Aug 2016

% cur = cd; % save current directory
% cd(dicomDir);
% check if '.dcm' has been removed

if ~exist('dicomDir','var') || isempty(dicomDir)
%     dicomDir = pwd;
    dicomDir = uigetdir('.','Choose the source directory where DICOM files are saved');
end
if ~exist('destdir','var') || isempty(destdir)
    destdir = uigetdir('.','Choose the destination directory where DICOM files are to be saved');
end
files = listdir(dicomDir,'files');
disp(['Found ' num2str(length(files)) ' files in ' dicomDir])
if ~isempty(files)  % if there are any files in this directory
    dcms = nan(length(files),1);
    disp('Determining which files are dicom files');
    for f = 1:length(files)
        dcms(f) = isdicom(fullfile(dicomDir,files{f}));
    end
    disp([num2str(sum(dcms(:))) ' files are dicoms.'])
    w = waitbar(0,'Moving and sorting dicoms');
    for f = 1:length(files)
        if dcms(f) % If file is a dicom file
            clear tmpdir
            info = dicominfo(fullfile(dicomDir,files{f}));
            if isfield(info,'PatientName')
                fulldestdir = [destdir filesep info.PatientName.FamilyName];
            else
                fulldestdir = destdir;
            end
            if isfield(info,'SeriesDescription')
                tmpdir = fullfile(fulldestdir,sprintf('Series_%03d_%s',info.SeriesNumber,info.SeriesDescription));
                if ~exist(tmpdir,'dir');
                    %check name of new folder for possible forbidden
                    %characters
                    forbidden_chars = ':*?"<>|';
                    [tpath,name,ext] = fileparts(tmpdir);
                    tmpdirname = [name ext];
                    for ichar = 1:length(forbidden_chars)
                        I = strfind(tmpdirname,forbidden_chars(ichar));
                        if ~isempty(I)
                            tmpdirname(I) = [];
                        end
                    end
                    tmpdir=fullfile(tpath,tmpdirname);
                    mkdir(tmpdir);
                    disp(['Sorting dicoms into ' info.SeriesDescription '_series_' num2str(info.SeriesNumber)]);
                end
            end
            if ~strcmp(files{f}(end-3:end),'.dcm')
                % If extension is not '.dcm', add that extension
                copyfile(fullfile(dicomDir,files{f}),...
                    ([fullfile(tmpdir,files{f}) '.dcm']));
            else
                % Sort dicom files according to ProtocolName
                copyfile(fullfile(dicomDir,files{f}),...
                    fullfile(tmpdir,files{f}));
            end
        end
        waitbar(f/length(files),w)
    end
    close(w)
end