% cdsort3.m
% Assumes that all the files from a CD are copied to a specific directory
% and then calls stand_alone_sort3 to load the files in the correct
% subdirectories corresponding to series

fclose all;


% Get directory where the dicom files are stored
srcdir = uigetdir('.','Choose the source directory where DICOM files are saved');
% Get destination directory name
destdir = uigetdir('.', 'Choose the destination directory to save sorted data');

if (( srcdir ~= 0) && (destdir ~= 0) )
% Now copy with series as subdirectories full of images from each series
    outdir = standalone_sort3(srcdir,destdir,1);
end





