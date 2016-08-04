% Exercise CEST - mexfiles creation script

disp('Mex files creation make file for CEST COMMON directory');

CC = mex.getCompilerConfigurations('c');
if ( isempty(CC(1).Name) )
    setuprun = 0;
else
    setuprun = 1;
end
if (setuprun == 1)
    disp('mexB0correctedCESTimglfit');
    mex ('mexB0correctedCESTimglfit.c', 'functions/lfit.c', 'functions/nrutil.c',  'functions/gaussj.c','functions/covsrt.c');
    disp('mexwassr_b0map2d');
    mex('mexwassr_b0map2d.c');
%     disp('mexB0corr_interp_map_zspec.c');
%     mex('mexB0corr_interp_map_zspec.c');
    disp('mexB0corr_interp_map_zspecinteg.c');
    mex('mexB0corr_interp_map_zspecinteg.c');
else
    if ( strfind(archstr,'win64') )
        disp(' You may need to download supported 64bit C++ compiler for your MATLAB version');
    end
    disp (' Run mex -setup and choose a C compiler.');
    disp (' Then reexecute this script');
end