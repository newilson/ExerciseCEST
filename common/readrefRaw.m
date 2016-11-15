function [img,WIPlong,WIPdbl,twix_obj,FileName] = readrefRaw(FileName,flip)

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

twix_obj.image.flagRemoveOS = true;
data = 32767*twix_obj.image{''}; % fe-ch-pe-reps
data = permute(data,[3 1 2 4]); % pe-fe-ch-reps

img = zeros(size(data,1),size(data,2),size(data,4));
for ii=1:size(data,4)
    img(:,:,ii) = NWcomplex_coil_combine(squeeze(data(:,:,:,ii)),[],false);
end
    
img = abs(2000*img/max(abs(img(:))));

for ii=1:size(img,3)
    img(:,:,ii) = anisodiff(squeeze(img(:,:,ii)),20,50,0.03,1);
end


if flip
    img = fliplr(img);
end

end