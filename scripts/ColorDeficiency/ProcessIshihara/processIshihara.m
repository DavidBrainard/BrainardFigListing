
grayIshihara = imread('GrayIshihara.jpg');

lowThresh = 50;
highThresh = 220;
index = (grayIshihara < highThresh & grayIshihara > lowThresh);
testIndex = grayIshihara;
testIndex(index) = 128;
figure; imshow(testIndex);

redIshihara = double(grayIshihara)/255;
greenIshihara = double(grayIshihara)/255;
blueIshihara = double(grayIshihara)/255;

redIshihara(index) = 1.05*redIshihara(index);
greenIshihara(index) = 1*greenIshihara(index);
blueIshihara(index) = 0.5*blueIshihara(index);

colorIshihara(:,:,1) = redIshihara;
colorIshihara(:,:,2) = greenIshihara;
colorIshihara(:,:,3) = blueIshihara;
figure; imshow(colorIshihara);
imwrite(colorIshihara,'colorIshihara.jpg','jpg');
