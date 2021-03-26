function [ imgout ] = blendTexture(warped_img1, warped_pmap1, warped_img2, warped_pmap2)

% img1Homo：aligned target image
% pmap1Homo: salient target image
% img2Homo：aligned reference image
% pmap2Homo: salient reference image

%%  pre-process of blendtexture
w1 = imfill(imbinarize(rgb2gray(warped_img1), 0),'holes');
w2 = imfill(imbinarize(rgb2gray(warped_img2), 0),'holes');
A = w1;  B = w2;
C = A & B;  % mask of overlapping region
[ sz1, sz2 ]=size(C);
ind = find(C);  % index of overlapping region
nNodes = size(ind,1);
revindC = zeros(sz1*sz2,1);
revindC(C) = 1:length(ind);

[m,~] = find(C);
minrow = min(m);
maxrow = max(m);
%%  terminalWeights, choose source and sink nodes

% border of overlapping region
BR=(B-[B(:,2:end) false(sz1,1)])>0;
BL=(B-[false(sz1,1) B(:,1:end-1)])>0;
BD=(B-[B(2:end,:);false(1,sz2)])>0;
BU=(B-[false(1,sz2);B(1:end-1,:)])>0;

CR=(C-[C(:,2:end) false(sz1,1)])>0;
CL=(C-[false(sz1,1) C(:,1:end-1)])>0;
CD=(C-[C(2:end,:);false(1,sz2)])>0;
CU=(C-[false(1,sz2);C(1:end-1,:)])>0;

imgseedR=(BR|BL|BD|BU)&(CR|CL|CD|CU);
imgseedB = (CR|CL|CD|CU) & ~imgseedR;

% data term
tw=zeros(nNodes,2);
tw(revindC(imgseedB),2)=inf;
tw(revindC(imgseedR),1)=inf;

mask=~(imgseedB | imgseedR);

imglap = 0.5.*imadd(warped_img1.*cat(3,mask&C,mask&C,mask&C), warped_img2.*cat(3,mask&C,mask&C,mask&C));
freeseed = imgseedR;
freeseed(minrow,:) = 0;
freeseed(maxrow,:) = 0;
free_seed = imgseedR & (~freeseed);
% saliency weights
%======================
tw(revindC(free_seed),1) = 0;
tw(revindC(free_seed),2) = 0;
imgseed = warped_img1.*cat(3,A-C,A-C,A-C) + warped_img2.*cat(3,B-C,B-C,B-C) + imglap+cat(3,freeseed,free_seed,imgseedB);
%======================
% figure,imshow(imgseed);
% title('seed image on two warped images');
terminalWeights=tw;       % data term

%% calculate edgeWeights

CL1=C&[C(:,2:end) false(sz1,1)];
CL2=[false(sz1,1) CL1(:,1:end-1)];
CU1=C&[C(2:end,:);false(1,sz2)];
CU2=[false(1,sz2);CU1(1:end-1,:)];
blurred1 = warped_img1; blurred2 = warped_img2;

%  edgeWeights:  Euclidean or sigmoid or saliency-weighted sigmoid norm
ang_1=blurred1(:,:,1); sat_1=blurred1(:,:,2); val_1=blurred1(:,:,3);
ang_2=blurred2(:,:,1); sat_2=blurred2(:,:,2); val_2=blurred2(:,:,3);
% baseline difference map
imgdif = sqrt( ( (ang_1.*C-ang_2.*C).^2+(sat_1.*C-sat_2.*C).^2+ (val_1.*C-val_2.*C).^2 )./3);

% sigmoid-metric difference map
a_rgb = 0.06; % bin of histogram
beta=4/a_rgb; % beta
gamma=exp(1); % base number
para_alpha = histOstu(imgdif(C), a_rgb);  % parameter:tau
imgdif_sig = 1./(1+power(gamma,beta*(-imgdif+para_alpha))); % difference map with logistic function
imgdif_sig = imgdif_sig.*C;   % difference to compute the smoothness term

% saliency-weighted sigmoid difference map
saliency_map = ((warped_pmap1 + warped_pmap2)./2).*C;  % the saliency of overlapping region
imgdif_sal = imgdif_sig.*(1 + saliency_map); % difference to compute the smoothness term
% saliency-weighted sigmoid method
DL = (imgdif_sal(CL1)+imgdif_sal(CL2))./2;
DU = (imgdif_sal(CU1)+imgdif_sal(CU2))./2;


% smoothness term
edgeWeights=[
    revindC(CL1) revindC(CL2) DL+1e-8 DL+1e-8;
    revindC(CU1) revindC(CU2) DU+1e-8 DU+1e-8
    ];

%%  graph-cut labeling

[~, labels] = graphCutMex(terminalWeights,edgeWeights);

As=A;
Bs=B;
As(ind(labels==1))=false;   % mask of target seam
Bs(ind(labels==0))=false;   % mask of reference seam
imgout = gradient_blend(warped_img1, As, warped_img2);
figure,imshow(imgout);

% 在全尺寸拼接结果中显示缝合线
SE_seam = strel('square', 3);
As_seam_ = imdilate(As, SE_seam);
Bs_seam_ = imdilate(Bs, SE_seam);
Cs_seam_ = As_seam_ & Bs_seam_;
imgseam = imgout.*cat(3,(A|B)-Cs_seam_,(A|B)-Cs_seam_,(A|B)-Cs_seam_) + cat(3,Cs_seam_,zeros(sz1,sz2),zeros(sz1,sz2));

% figure,imshow(imgseam);
% title('final stitching seam');

end

function [ alpha ] = histOstu(edge_D, interval_rgb)
% use Ostu's method to compute the parameter of logistic function
xbins = 0+interval_rgb/2:interval_rgb:1-interval_rgb/2;
num_x = size(xbins,2);
counts = hist(edge_D, xbins);
%figure,hist(edge_D, xbins);
num = sum(counts);
pro_c = counts./num;
ut_c = pro_c.*xbins;
sum_ut = sum(ut_c);
energy_max = 0;
for k=1:num_x
    uk_c = pro_c(1:k).*xbins(1:k);
    sum_uk = sum(uk_c);
    sum_wk = sum(pro_c(1:k));
    
    sigma_c = (sum_uk-sum_ut*sum_wk).^2/(sum_wk*(1-sum_wk));
    
    if sigma_c>energy_max
        energy_max = sigma_c;
        threshold = k;
    end
end

alpha = xbins(threshold)+interval_rgb/2;
end