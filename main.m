%%% Znib Project
%%% Fabio Maschi
%%% 3 October 2018
%%% Version 0.2

clear all; close all;

%% Parameters
N_REFERENCES = 4;   % number of reference enzymes per patch
B_WITH_CLOSING = 1; % 0 FALSE 1 TRUE
B_WITH_SAVING = 1;  % 0 FALSE 1 TRUE
%
K_SEGMENTATION = 1.5; % between 1.3 and 2
K_ROUND_TOLERANCE = 0.2;
%
CORE_SIZE = 300; % in pixels
K_HITSIZE = 10; % in pixels

%% Open image
img = double(rgb2gray(imread('Picture10.png')));
[w,h] = size(img);

%% Segmentation
% histogram
hstgrm = zeros(1,255);
for x=1:w
    for y=1:h
        hstgrm(img(x,y)+1) = hstgrm(img(x,y)+1)+1;
    end
end
hstgrm = hstgrm / (w*h);

% % omega(i) = \sum_1^i { hstgrm(k) }
omega_0 = cumsum(hstgrm);
omega_1 = 1 - omega_0;

% % mu(i) = \sum_1^i { hstgrm(k)*i }
mu = cumsum(hstgrm .* (1:255));

mu_0 = mu ./ omega_0;
mu_1 = (mu(255) - mu) ./ (omega_1);

var2b = omega_0.*omega_1.*(mu_1 - mu_0).^K_SEGMENTATION;

[~, K] = max(var2b);
img_segmented = logical(img > K);

% figure;
% range1 = [min(img(:)) max(img(:))];
% range2 = [min(img_bw(:)) max(img_bw(:))];
% subplot(1,3,1); imshow(img, range1);
% subplot(1,3,2); imshow(img_bw, range2);
% subplot(1,3,2); plot(hstgrm);
% xlim([0 255]);
% xlabel('Grayscale');
% ylabel('# pixels');
% title('Grayscale histogram');

% figure;
% imshow(img_bw);
clear K mu mu_0 mu_1 var2b omega_0 omega_1 hstgrm;

%% Closing
if B_WITH_CLOSING == 1
    img_connected = single(~(img_segmented/255));
    %kernel = [0 1 0; 1 1 1; 0 1 0];
    %kernel = [0 1 1 1 0; 1 1 1 1 1; 0 1 1 1 0];
    kernel = [0 0 1 0 0; 0 1 1 1 0; 1 1 1 1 1; 0 1 1 1 0; 0 0 1 0 0];
    img_connected = conv2(img_connected, kernel, 'same');
    img_connected = single(~(img_connected));
    img_connected = conv2(img_connected, kernel, 'same');
else
    img_connected = single(img_segmented/255);
end

%% Filling wholes
CC = bwconncomp(img_connected);
numPixels = cellfun(@numel, CC.PixelIdxList);

for i=1:size(numPixels,2)
    if numPixels(i) <= CORE_SIZE
        img_connected(CC.PixelIdxList{i}) = 0;
    end
end

img_connected = ~img_connected;

if B_WITH_SAVING == 1
    mask = single(cat(3,img/255 + (img_connected),img/255,img/255));
    imwrite(mask, 'test.jpg');
end
clear numPixels;


%% ROUND CELLS
CC = bwconncomp(img_connected);
stats = regionprops(CC, 'Area','Perimeter');
pixels = regionprops(CC, 'PixelIdxList');
area_perim = [stats.Area; stats.Perimeter];
C = (4*pi*area_perim(1,:))./(area_perim(2,:).^2);
cellule_ronde = false(w,h);

for i = 1:CC.NumObjects 
   if (C(i) > (1-K_ROUND_TOLERANCE) && C(i) < (1+K_ROUND_TOLERANCE))
       cellule_ronde(CC.PixelIdxList{i}) = true;
   end
end

% mask = single(cat(3,img/255+img_connected,img/255+cellule_ronde,img/255));
% imshow(mask), hold on
% imwrite(mask, 'fill-off.jpg');

%% Grid detection
st = regionprops(cellule_ronde, 'Centroid');
c = vertcat(st.Centroid);
[~, x] = kmeans(c(:,1), 12);
[~, y] = kmeans(c(:,2), 8);
x = sort(x);
y = sort(y);

% I = imoverlay(img, cellule_ronde, [0.9 0.1 0.1]);
% imshow(I, 'InitialMag',200, 'Border','tight'), hold on
% line(c(:,1), c(:,2), 'LineStyle','none', 'Marker','+', 'Color','b')
% for i = 1:length(x)
%     plot([x(i); x(i)], [0; w], 'g-', 'LineWidth',2);
% end
% for i = 1:length(y)
%     plot([0; h], [y(i); y(i)], 'g-', 'LineWidth',2);
% end
% hold off

%% Selecting gridded cells
img = img/255;
img_stdenz = zeros(w, h);
for i = 1:length(x)
    for j = 1:length(y)
        %mask = false(w,h);
        %mask(ceil(y(j)-K_HITSIZE:y(j)+K_HITSIZE), ceil(x(i)-K_HITSIZE:x(i)+K_HITSIZE)) = true;
        my = ceil(y(j)-K_HITSIZE:y(j)+K_HITSIZE);
        my = max(1,my);
        my = min(w,my);
        mx = ceil(x(i)-K_HITSIZE:x(i)+K_HITSIZE);
        mx = max(1,mx);
        mx = min(h,mx);
        xy = sub2ind([w,h], my, mx);
        for k = 1:CC.NumObjects
            if(any(ismember(xy,pixels(k).PixelIdxList)))
                img_stdenz(CC.PixelIdxList{k}) = 1;
                break;
            end
        end
    end
end
mask = single(cat(3,img/255+img_connected,img/255+img_stdenz,img/255));
imshow(mask), hold on
line(c(:,1), c(:,2), 'LineStyle','none', 'Marker','+', 'Color','b')
for i = 1:length(x)
    plot([x(i); x(i)], [0; w], 'g-', 'LineWidth',2);
end
for i = 1:length(y)
    plot([0; h], [y(i); y(i)], 'g-', 'LineWidth',2);
end
hold off
stop
clear i j my mx;

%% Selecting References
img_refenz = zeros(w, h);
ref_area = 0;
for i = 1:N_REFERENCES
    mask = single(cat(3,img+img_refenz,img+img_stdenz,img+img_stdenz));
    imshow(mask);
    [x, y] = ginput(1);
    xy = sub2ind([w,h], ceil(y), ceil(x));
    for j = 1:CC.NumObjects
        if(ismember(xy,pixels(j).PixelIdxList))
            img_stdenz(CC.PixelIdxList{j}) = 0;
            img_refenz(CC.PixelIdxList{j}) = 1;
            ref_area = ref_area + stats(j).Area;
            break;
        end
    end
end
mask = single(cat(3,img+img_refenz,img+img_stdenz,img+img_stdenz));
imshow(mask);
ref_area = ref_area / N_REFERENCES

%% Computing area
CC = bwconncomp(img_stdenz);
stats = regionprops(CC, 'Area','Perimeter');

for i = 1:CC.NumObjects
    if stats(i).Area < ref_area
        img_stdenz(CC.PixelIdxList{i}) = 0;
    else
        % print stats(i).Area / ref_area * 100 "%"
    end
end

mask = single(cat(3,img+img_refenz,img+img_stdenz,img+img_stdenz));
imshow(mask);

% imoverlay function
% http://www.mathworks.com/matlabcentral/fileexchange/10502-image-overlay