function roimask = createroimask2d(image,nroi,summode,roimode,threshold,maxval)
%  roimask = output 2D array (H,W,nroi)
%  image = input 2D array (H,W)
%  nroi = number of roi's needed
%  summode = 0 or 1 (returns the individual masks in a 3D
%                    array or the combined mask in the 2D array)
%  roimode : 1 = 'poly', 2 = 'poly_annulus',3 = 'amplitude',4 = 'maxval',5 = 'circle',
%            6 = 'square',7 = 'freehand'
%  threshold = fraction used for roimode = 3.
%
if ( nargin < 1)
    image = imread('cameraman.tif');
    nroi = 3;
    summode = 1;
    threshold = 0.1;
    maxval = max(image(:));
end
if (nargin < 2)
    nroi = 1;
    summode = 0;
    threshold = 0.1;
    maxval = max(image(:));
end
if (nargin < 3)
    nroi = 1;
    summode = 0;
    threshold = 0.1;
    maxval = max(image(:));
end
if (nargin < 4)
    summode = 0;
    threshold = 0.1;
    maxval = max(image(:));
end
if (nargin < 5)
    str1 = sprintf('Choose ROI selection mode ');
    roimode = menu(str1, 'poly', 'poly_annulus','amplitude','maxval','circle','square','freehand');
    threshold = 0.1;
    maxval = max(image(:));
end
if (nargin < 6)
    threshold = 0.1;
    maxval = max(image(:));
end
if (nargin < 7)
    maxval = max(image(:));
end
    

minval = min(image(:));

if (roimode == 3)
    if (nargin < 6)
        nroi = 1;
        prompt = {'Enter threshold value in %'};
        dlgtitle = 'Threshold';
        num_lines = 1;
        def = {'10.0'};
        thresholdstr = inputdlg(prompt,dlgtitle,num_lines,def);
        [threshold, status] = str2num(thresholdstr{1});
        if ( ~status)
            threshold = 10.0;
        end
    end
elseif (roimode == 4)
    nroi = 1;
end

if (nroi == 1)
    summode = 1;
end

H = size(image,1);
W = size(image,2);

if (summode == 0)
    roimask = zeros(H,W,nroi);
else
    roimask = zeros(H,W);
end
mask = zeros(H,W,nroi);
maskt = zeros(H,W);

roih = figure;

redrawmode = 1;

while (redrawmode == 1)
    if (roimode == 1)     % simple polygon
        figure(roih); hold off; imh = imshow(image,[0.5*minval,1.3*maxval]);axis('square'); title(sprintf('Choose %i ROIs',nroi));  hold on;
        for i = 1: nroi
            mask(:,:,i) = maskt;
            maskt = image*0.0;
            V =  roipoly; hold on;
            maskt(V > 0.0) = 1.0;
            mask(:,:,i) = maskt;
        end
    elseif (roimode == 2) % annulus formed by 2 polygons
        figure(roih); hold off; imh = imshow(image,[0.5*minval,1.3*maxval]);axis('square'); title(sprintf('Choose %i ROIs',nroi));  hold on;
        for i = 1: nroi
            maskt = image*0.0;
            V1 = roipoly;
            V2 = roipoly;
            indor = ( find(V1 | V2) );
            indand = ( find(V1 & V2) );
            maskt(indor) = 1.0;
            maskt(indand) = 0.0;
            mask(:,:,i) = maskt;
        end
    elseif (roimode == 3) % Based on amplitude threshold
        for i = 1: nroi
            maskt = image *0.0;
            maskt(image > (threshold * 0.01 *maxval) ) = 1.0;
            mask(:,:,i) = maskt;
        end
    elseif (roimode == 4) % Based on max amplitude
        for i = 1: nroi
            maskt = image *0.0;
            maskt(image == maxval) = 1.0;
            mask(:,:,i) = maskt;
        end
    elseif (roimode == 5) % ellipse
        figure(roih); hold off; imh = imshow(image,[0.5*minval,1.3*maxval]);axis('square'); title(sprintf('Choose %i ROIs',nroi));  hold on;
        for i = 1: nroi
            maskt = image *0.0;
            if (i == 1)
                V =  imellipse;
                P = wait(V);
                pos = getPosition(V);
            else
                V =  imellipse(gca,[pos(1)+5 pos(2)+5 pos(3) pos(4)]);
                P = wait(V);
            end
            maskt = createMask(V,imh);
            mask(:,:,i) = maskt;
        end
    elseif (roimode == 6) % rectangle
        figure(roih); hold off; imh = imshow(image,[0.5*minval,1.3*maxval]);axis('square'); title(sprintf('Choose %i ROIs',nroi));  hold on;
        for i = 1: nroi
            maskt = image *0.0;
            if (i == 1)
                V =  imrect;
                P = wait(V);
                pos = getPosition(V);
            else
                V =  imrect(gca,[pos(1)+5 pos(2)+5 pos(3) pos(4)]);
                P = wait(V);
            end
            maskt = createMask(V,imh);
            mask(:,:,i) = maskt;
        end
    elseif (roimode == 7) % freehand
        figure(roih); hold off; imh = imshow(image,[0.5*minval,1.3*maxval]);axis('square'); title(sprintf('Choose %i ROIs',nroi));  hold on;
        for i = 1: nroi
            maskt = image *0.0;
            V =  imfreehand;
            P = wait(V);
            maskt = createMask(V,imh);
            mask(:,:,i) = maskt;
        end
    else                  % not defined - return warning
        WarnDlg('Mask set to zero since roimode is not correct');
    end
    if (roimode == 3 || roimode == 4)
        redrawmode = 0;
    else
        doneflag = menu('REDO','Done','Redraw');
        redrawmode =doneflag-1;
    end
end

if (summode == 0)
    roimask = mask;
else
    roimask = squeeze(sum(mask,3));
    roimask(roimask>1.0) = 1.0;
end

close(roih);



