
% MB
paramMB = getParam();
paramMB.use_backgroundness = false;
% MB+
paramMBplus = getParam();
paramMBplus.verbose = true; % to show run time for each step 

% Geodesic
paramGeo = getParam();
paramGeo.use_backgroundness = false;
paramGeo.use_geodesic = true;

I = imread('demo.jpg');

paramMBD.use_backgroudness = true;  
[pMap1, dMap1] = doMBS(I, paramMB);
[pMapG, dMapG] = doMBS(I, paramGeo);
[pMap2] = doMBS(I, paramMBplus); 

% display results
fh = figure(1);
set(fh, 'OuterPosition', [500,600,800,300]);

subplot(1,3,1)
imshow(I);
title('Input')

subplot(1,3,2)
imshow(pMap1)
title('MB')

subplot(1,3,3)
imshow(pMap2)
title('MB+')


fh = figure(2);
set(fh, 'OuterPosition', [500,100,800,500]);

subplot(2,3,1)
imshow(I);
title('Input')

subplot(2,3,2)
imshow(dMap1)
title('MBD Map')

subplot(2,3,3)
imshow(pMap1)
title('MB')

subplot(2,3,5)
imshow(dMapG)
title('Geodesic Distance Map')
subplot(2,3,6)
imshow(pMapG)
title('Geodesic Saliency Map')

