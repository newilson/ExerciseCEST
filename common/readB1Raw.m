function [B1map, img, WIPlong, WIPdbl, twix_obj, FileName, hdr] = readB1Raw(FileName,flip)

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

alpha = WIPdbl{2}*pi/180;

twix_obj.image.flagRemoveOS = true;
data = 32767*twix_obj.image{''}; % fe-ch-pe-con-reps
data = permute(data,[3 1 2 4 5]); % pe-fe-ch-con-reps

[pe,fe,ch,con,reps] = size(data);

img = zeros(pe,fe,con,reps);
for jj=1:reps
    for ii=1:con
        img(:,:,ii,jj) = NWcomplex_coil_combine(squeeze(data(:,:,:,ii,jj)),[],false);
    end
end

if flip
    img = fliplr(img);
end

img = 2000*img/max(abs(img(:))); % image must be scaled bc calc_B1map uses div by 1+image

B1map = zeros(pe,fe,reps);
if isequal(con,2)
    for ii=1:reps
        temp = calc_B1map(abs(img(:,:,1,ii)),abs(img(:,:,2,ii)),ones(pe,fe),alpha);
        B1map(:,:,ii) = anisodiff(temp,20,50,0.03,1);
    end
elseif isequal(con,3)
    for ii=1:reps
        temp = calc_B1map_new(abs(img(:,:,1,ii)),abs(img(:,:,2,ii)),abs(img(:,:,3,ii)),ones(pe,fe),alpha);
        B1map(:,:,ii) = anisodiff(temp,20,50,0.03,1);
    end
end