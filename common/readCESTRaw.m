function [posimages, negimages, pars, posB0maps, negB0maps, img, WIPlong, WIPdbl, twix_obj, FileName] = readCESTRaw(FileName,flip)

if nargin<1 || isempty(FileName)
    twix_obj = mapVBVD;
else
    twix_obj = mapVBVD(FileName);
end
if nargin<2, flip = false; end

if iscell(twix_obj)
    twix_obj = twix_obj{2};
end

WIPlong = twix_obj.hdr.MeasYaps.sWipMemBlock.alFree;
WIPdbl = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree;

seqname = twix_obj.hdr.Meas.tSequenceFileName;
seqname = strsplit(seqname, '\');
seqname = seqname{end};

if strcmp(seqname,'prep_moco')
    wipppmindex = 2;
    prepmodeindex = 2;
    cestz2value = 3;
    cestzvalue = 11; % NW
    CESTpw = 0.001*WIPlong{13};
    CESTpw1 = 0.001*WIPlong{14};
    CESTb1 = WIPlong{16};
    CESTdc = WIPdbl{1};
    reqreps = WIPlong{9};
    mocoreps = WIPlong{11};
elseif strcmp(seqname,'prep_moco_iB0')
    wipppmindex = 3; % iB0
    prepmodeindex = 2;
    cestz2value = 3;
    cestzvalue = 11; % NW
    CESTpw = 0.001*WIPlong{15}; % iB0
    CESTpw1 = 0.001*WIPlong{16}; % iB0
    CESTb1 = WIPlong{18}; % iB0
    CESTdc = WIPdbl{1};
    reqreps = WIPlong{10}; % iB0
    mocoreps = WIPlong{12}; % iB0
    if length(WIPlong)<27
        presatsense = 1;
    else
        presatsense = WIPlong{27}; % iB0
    end
    prePElines = WIPlong{7}; % iB0
else
    wipppmindex = input(' Type in the value for wipppmindex : ');
    prepmodeindex = input(' Type in the value for prepmodeindex : ');
    cestzvalue = input(' Type in CESTz value : ');
    CESTb1 = input(' Type in the value for CEST b1 in Hz : ');
    CESTpw = input(' Type in the value for CEST Pulse duration in ms : ');
    CESTpw1 = input(' Type in the value for CEST Pulse width1 in ms : ');
    CESTdc = input(' Type in the value for CEST dutycycle in % : ');
end

twix_obj.image.flagRemoveOS = true;
data = 32767*twix_obj.image('');

fe = size(data,1);
ch = size(data,2);
pe = size(data,3);
con = size(data,8);
reps = size(data,9);

data = squeeze(data);
if con>1
    data = permute(data,[3 1 2 4 5]); % pe-fe-ch-con-reps
else
    data = permute(data,[3 1 2 5 4]); % pe-fe-ch-con-reps
end


img = zeros(pe,fe,con,reps);
for jj=1:reps
    temp = data(:,:,:,:,jj);
    tempimg = zeros(pe,fe,con);
    for ii=[con,1:con-1]%con:-1:1
        if ii==con % Don't get S from presat! Why? || (ii==1 && presatsense==1) % same S for each contrast except the last
            [tempimg(:,:,ii),S] = NWcomplex_coil_combine(squeeze(temp(:,:,:,ii)),[],false);
        else
            if presatsense==1
                [tempimg(:,:,ii)] = NWcomplex_coil_combine(squeeze(temp(:,:,:,ii)),S,true,prePElines,prePElines/4);
            else
                tempimg(:,:,ii) = NWsense_coil_combine(squeeze(temp(:,:,:,ii)),S,true,prePElines,prePElines/4,presatsense,1);
            end
        end
    end
    S = [];
    img(:,:,:,jj) = tempimg;
end

if flip
    img = fliplr(img);
end

img = 2000*img/max(abs(img(:)));

posimages = abs(squeeze(img(:,:,con,1:2:end)));
negimages = abs(squeeze(img(:,:,con,2:2:end)));

ppmbegin = WIPdbl{wipppmindex+1};  % begin
ppmend   = WIPdbl{wipppmindex+2};  % end
ppmstep  = WIPdbl{wipppmindex+3};  % step
if (ppmbegin > ppmend)
    ppmstep = -ppmstep;
end
CESTppmlist = ppmbegin:ppmstep:ppmend;
if WIPlong{prepmodeindex} == cestzvalue
    warning('not yet implemented')
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

if con>1 % B0maps
    pars.presatflip = WIPdbl{3};
    pars.prePElines = WIPlong{7};
    pars.TEms = cell2mat(twix_obj.hdr.MeasYaps.alTE)/1000;
    pars.sf = twix_obj.hdr.Dicom.lFrequency*1e-6;
    width = [pe pe];
    slope = round(width/4);
    filt = NWSiemensRawFilter2d([pe fe],width,slope);
    for ii=1:reps % filter high res image
        ktemp = ifft2c(squeeze(img(:,:,con,ii)));
        ktemp = ktemp.*filt;
        Itemp = fft2c(ktemp);
        img(:,:,con,ii) = Itemp;
    end
    phases = angle(img);
    phases = phases(:,:,1:con-1,:); pars.TEms = pars.TEms(1:con-1);
    B0maps = NWcalcB0gre([],phases,pars)/pars.sf;
    for ii=1:reps
        B0maps(:,:,ii) = anisodiff(B0maps(:,:,ii),20,50,0.03,1);
    end
    posB0maps = B0maps(:,:,1:2:end);
    negB0maps = B0maps(:,:,2:2:end);
else
    posB0maps = [];
    negB0maps = [];
end

end

