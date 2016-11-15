%% CEST images
function [posimages, negimages, hdr, pars, pathname1] = readCEST(pathname1)

if nargin<1 || isempty(pathname1)
    [hdr,images,dicomhdr,pathname1] = readdicomfiles2d;
else
    [hdr,images,dicomhdr] = readdicomfiles2d(pathname1);
end

nfreq = hdr.reps/2;
nfiles = hdr.reps;
seqname = hdr.SeqName;
swversion = hdr.swversion;
if ( strfind( swversion, 'VD13A' ) )
    if (~strfind(seqname,'cest_flash'))
        warning('not cest_flash sequence')
        disp(seqname)
    end
    wipppmindex = 3;
    prepmodeindex = 1;
    cestz2value = 5;
    cestzvalue = 1e5; % arbitrarily large
    CESTpw = hdr.WIPlong(11);
    CESTpw1 = hdr.WIPlong(13);
    CESTb1 = hdr.WIPlong(14);
    CESTdc = hdr.WIPdbl(1);
    reqreps = hdr.WIPlong(8);
    mocoreps = 1;
elseif ( strfind( swversion, 'VD13D' ) )
    seqname = strsplit(seqname, '\');
    seqname = seqname{end};
    if strcmp(seqname,'prep_moco')
        wipppmindex = 2;
        prepmodeindex = 2;
        cestz2value = 3;
        cestzvalue = 11; % NW
        CESTpw = 0.001*hdr.WIPlong(13);
        CESTpw1 = 0.001*hdr.WIPlong(14);
        CESTb1 = hdr.WIPlong(16);
        CESTdc = hdr.WIPdbl(1);
        reqreps = hdr.WIPlong(9); % NW
        mocoreps = hdr.WIPlong(11); % NW
    elseif strcmp(seqname,'prep_moco_iB0')
        wipppmindex = 3; % iB0
        prepmodeindex = 2;
        cestz2value = 3;
        cestzvalue = 11; % NW
        CESTpw = 0.001*hdr.WIPlong(15); % iB0
        CESTpw1 = 0.001*hdr.WIPlong(16); % iB0
        CESTb1 = hdr.WIPlong(18); % iB0
        CESTdc = hdr.WIPdbl(1);
        reqreps = hdr.WIPlong(10); % NW iB0
        mocoreps = hdr.WIPlong(12); % NW iB0
    else
        wipppmindex = input(' Type in the value for wipppmindex : ');
        prepmodeindex = input(' Type in the value for prepmodeindex : ');
        cestzvalue = input(' Type in CESTz value : ');
        CESTb1 = input(' Type in the value for CEST b1 in Hz : ');
        CESTpw = input(' Type in the value for CEST Pulse duration in ms : ');
        CESTpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
        CESTdc = input(' Type in the value for CEST dutycycle in % : ');
    end
else
    if (strfind(seqname,'prep_moco'))
        wipppmindex = 2;
        prepmodeindex = 2;
        cestz2value = 3;
        if ( hdr.WIPlong(19) >= 170 )
            CESTpw = hdr.WIPlong(10);
            CESTb1 = hdr.WIPlong(13);
            CESTpw1 = hdr.WIPlong(11);
            CESTdc = hdr.WIPdbl(1);
        else
            CESTb1 = input(' Type in the value for CEST b1 in Hz : ');
            CESTpw = input(' Type in the value for CEST Pulse duration in ms : ');
            CESTpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
            CESTdc = input(' Type in the value for CEST dutycycle in % : ');
        end
    elseif (strfind(seqname,'prep_cv'))
        wipppmindex = 2;
        prepmodeindex = 1;
        cestz2value = 8;
        if ( hdr.WIPlong(19) >= 170 )
            CESTpw = hdr.WIPlong(10);
            CESTb1 = hdr.WIPlong(13);
            CESTpw1 = hdr.WIPlong(11);
            CESTdc = hdr.WIPdbl(1);
        else
            CESTb1 = input(' Type in the value for CEST b1 in Hz : ');
            CESTpw = input(' Type in the value for CEST Pulse duration in ms : ');
            CESTpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
            CESTdc = input(' Type in the value for CEST dutycycle in % : ');
        end
    elseif (strfind(seqname,'prep_flash'))
        prepmodeindex = 1;
        cestz2value = 8;
        if ( hdr.WIPlong(19) >= 170 )
            CESTpw = hdr.WIPlong(11);
            CESTb1 = hdr.WIPlong(12);
            CESTpw1 = hdr.WIPlong(15);
            CESTdc = hdr.WIPdbl(1);
            wipppmindex = 3;
        else
            CESTb1 = input(' Type in the value for CEST b1 in Hz : ');
            CESTpw = input(' Type in the value for CEST Pulse duration in ms : ');
            CESTpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
            CESTdc = input(' Type in the value for CEST dutycycle in % : ');
            wipppmindex = 5;
        end
    end
end
if ~exist('reqreps','var')
    warning('Required reps and moco reps were not read in')
    reqreps = 1;
    mocoreps = 1;
end
posimages = images(:,:,1:2:end);
negimages = images(:,:,2:2:end);
ppmbegin = hdr.WIPdbl(wipppmindex+1);  % begin
ppmend   = hdr.WIPdbl(wipppmindex+2);  % end
ppmstep  = hdr.WIPdbl(wipppmindex+3);  % step
if (ppmbegin > ppmend)
    ppmstep = -ppmstep;
end
CESTppmlist = ppmbegin:ppmstep:ppmend;
if hdr.WIPlong(prepmodeindex) == cestzvalue
    DSimage = negimages(:,:,1);
    posimages = posimages(:,:,2:2:end);
    negimages = negimages(:,:,2:2:end);
end

pars.CESTpw = CESTpw;
pars.CESTb1 = CESTb1;
pars.CESTpw1 = CESTpw1;
pars.CESTdc = CESTdc;
pars.reqreps = reqreps;
pars.mocoreps = mocoreps;
pars.CESTppmlist = CESTppmlist;

end