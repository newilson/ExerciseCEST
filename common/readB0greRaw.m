function [B0maps, img, WIPlong, WIPdbl, twix_obj, FileName, hdr] = readB0greRaw(FileName,flip,recon)

if nargin<1 || isempty(FileName)
    twix_obj = mapVBVD;
else
    twix_obj = mapVBVD(FileName);
end

if nargin<2 || isempty(flip), flip = false; end
if nargin<3 || isempty(recon), recon = 0; end
% recon: 0 = linear fit, 1 = weighted mean of coil phase differences, 2 =
% lin fit before coil combination

if iscell(twix_obj)
    twix_obj = twix_obj{2};
end

WIPlong = twix_obj.hdr.MeasYaps.sWipMemBlock.alFree;
WIPdbl = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree;

twix_obj.image.flagRemoveOS = true;
data = 32767*twix_obj.image{''}; % fe-ch-pe-con-reps
data = permute(data,[3 1 2 4 5]); % pe-fe-ch-con-reps

peorig = size(data,1); fe = size(data,2);
flag_zf = false;
if peorig<fe
%     flag_zf = true;
    data = padarray(data,[(fe-peorig)/2 0]);
end

hdr = [];
hdr.TEms = cell2mat(twix_obj.hdr.MeasYaps.alTE)/1000;
hdr.sf = twix_obj.hdr.Dicom.lFrequency*1e-6;

if recon==1 % only 2 contrasts
    data = data(:,:,:,1:2,:);
end

[pe,fe,ch,con,reps] = size(data);

switch recon
    case 0 % linear fit of complex coil combined images
        img = zeros(pe,fe,con,reps);
        for jj=1:reps
            for ii=1:con
                if ii==1 % same S for each contrast
                    if ~flag_zf
                        [img(:,:,ii,jj),S] = NWcomplex_coil_combine(squeeze(data(:,:,:,ii,jj)),[],false);
                    else
                        [img(:,:,ii,jj),S] = NWcomplex_coil_combine(squeeze(data(:,:,:,ii,jj)),[],true,peorig,peorig/4);
                    end
                else
                    if ~flag_zf
                        img(:,:,ii,jj) = NWcomplex_coil_combine(squeeze(data(:,:,:,ii,jj)),S,false);
                    else
                        img(:,:,ii,jj) = NWcomplex_coil_combine(squeeze(data(:,:,:,ii,jj)),S,true,peorig,peorig/4);
                    end
                end
            end
        end
        
        if flip
            img = fliplr(img);
        end
        
        img = 2000*img/max(abs(img(:)));
        phases = angle(img);
        
        B0maps = NWcalcB0gre([],phases,hdr)/hdr.sf;

    case 1 % weight mean of phase difference by coil
        phases = zeros(pe,fe,reps);
        for ii=1:reps
            kc1 = squeeze(data(:,:,:,1,ii)); kc2 = squeeze(data(:,:,:,2,ii));
            phasediff(:,:,ii) = NWphasediff_coil_combine(kc1,kc2,false);
        end
        diffTE = diff(hdr.TEms(1:2)/1000);
        B0maps = phasediff/(2*pi*diffTE*hdr.sf);
        
        if flip
            B0maps = fliplr(B0maps);
        end
        
        img = [];
        
    case 2 % linear fit coil-by-coil and then weighted mean combination
        data = reshape(data,pe,fe,[]);
        img = fft2c(data);
        img = reshape(img,[pe,fe,ch,con,reps]);
        img = permute(img,[1 2 4 3 5]); % pe-fe-con-ch-reps
        if flip
            img = fliplr(img);
        end
        temp_c = reshape(img,[pe,fe,con,ch*reps]);
        phases = angle(temp_c);
        B0maps_c = NWcalcB0gre([],phases,hdr)/hdr.sf;
        B0maps_c = reshape(B0maps_c,[pe,fe,ch,reps]);
        weights = squeeze(abs(img(:,:,1,:,:))); % pe-fe-ch-reps
        sum_weights = squeeze(sum(weights,3)); % pe-fe-reps
        B0maps = squeeze(sum(B0maps_c.*weights,3))./sum_weights;
end

for ii=1:reps
    B0maps(:,:,ii) = anisodiff(B0maps(:,:,ii),20,50,0.03,1);
end

end