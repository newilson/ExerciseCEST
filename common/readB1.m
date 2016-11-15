function [B1map,hdr,alpha,images,pathname1] = readB1(pathname1)

if nargin<1 || isempty(pathname1)
    [hdr,images,dicomhdr,pathname1] = readdicomfiles2d;
else
    [hdr,images,dicomhdr] = readdicomfiles2d(pathname1);
end

if ~isempty(strfind(hdr.SeqName,'prep_moco')) && isequal(size(images,3),3) % NW check this
    alpha = hdr.WIPdbl(2)*pi/180;
    app.txtDisp.Value = 'Calculating B1 map';
    B1map = calc_B1map_new(images(:,:,1),images(:,:,2),images(:,:,3),ones(size(squeeze(images(:,:,1)))),alpha);
elseif ~isempty(strfind(hdr.swversion,'VD13A')) && isequal(size(images,3),2)
    alpha = hdr.WIPdbl(2)*pi/180;
    %                         alpha = 15.2*pi/180; % CLIPPED RF
    app.txtDisp.Value = 'Calculating B1 map';
    B1map = calc_B1map(images(:,:,1),images(:,:,2),ones(size(squeeze(images(:,:,1)))),alpha);
elseif ~isempty(strfind(lower(pathname1),'flip30')) && isequal(size(images,3),1)
    alpha = pi/6.0;
    pathparts = strsplit(pathname1,'_');
    for jj=1:length(pathparts)
        if strcmp(pathparts{jj},'flip30')
            pathparts{jj} = 'flip60';
        end
    end
    pathname2 = strjoin(pathparts,'_');
    [~,image2] = readdicomfiles2d(pathname2);
    app.txtDisp.Value = 'Calculating B1 map';
    B1map = calc_B1map(images,image2,ones(size(squeeze(images(:,:,1)))),alpha);
end
if 1 % Always filter B1map app.fit.isfilt
    B1map = anisodiff(B1map,20,50,0.03,1);
end

end