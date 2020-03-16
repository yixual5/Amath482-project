%% simple case
load('cam1_1.mat')
load('cam2_1.mat')
load('cam3_1.mat')
numFrames = size(vidFrames1_1,4);
py=zeros(1,numFrames);
px=zeros(1,numFrames);
for j = 1:numFrames
    X= vidFrames1_1(:,:,:,j);
    X_dgray = double(rgb2gray(X));%we can handle double data in black and white color domain
    X_dgray(:,1:304)=0;%eliminate the range that does not contain the mass
    X_dgray(:,388:end)=0;
    X_dgray([1:169 425:end],305:387)=0;
    Mx = max(X_dgray(:));
    [yrange,xrange] = find(X_dgray>0.95*Mx);
    py(j)=mean(yrange);
    px(j)=mean(xrange);
end
[~,idx] = max(py(1:30));
py=py(idx:end);
px=px(idx:end);

numFrames2 = size(vidFrames2_1,4);
py2=zeros(1,numFrames2);
px2=zeros(1,numFrames2);
for j = 1:numFrames2
    X2= vidFrames2_1(:,:,:,j);
    X2_dgray = double(rgb2gray(X2));%we can handle double data in black and white color domain
    X2_dgray(:,1:241)=0;%eliminate the range that does not contain the mass
    X2_dgray(:,369:end)=0;
    X2_dgray([1:97 371:end],242:368)=0;
    Mx = max(X2_dgray(:));
    [yrange,xrange] = find(X2_dgray>0.95*Mx);
    py2(j)=mean(yrange);
    px2(j)=mean(xrange);
end
[~,idx] = max(py2(1:30));
py2=py2(idx:end);
px2=px2(idx:end);

numFrames3 = size(vidFrames3_1,4);
py3=zeros(1,numFrames3);
px3=zeros(1,numFrames3);
for j = 1:numFrames3
    X3= vidFrames3_1(:,:,:,j);
    X3_dgray = double(rgb2gray(X3));%we can handle double data in black and white color domain
    X3_dgray(:,1:224)=0;%eliminate the range that does not contain the mass
    X3_dgray(:,507:end)=0;
    X3_dgray([1:212 357:end],225:506)=0;
    Mx = max(X3_dgray(:));
    [yrange,xrange] = find(X3_dgray>0.95*Mx);
    py3(j)=mean(xrange); %reverse the coordinate due to the camera angle
    px3(j)=mean(yrange);
end
[~,idx] = max(py3(1:30));
py3=py3(idx:end);
px3=px3(idx:end);

cut=min([length(py) length(py2) length(py3)]);
py=py(1:cut)-mean(py(1:cut));
px=px(1:cut)-mean(px(1:cut));
py2=py2(1:cut)-mean(py2(1:cut));
px2=px2(1:cut)-mean(px2(1:cut));
py3=py3(1:cut)-mean(py3(1:cut));
px3=px3(1:cut)-mean(px3(1:cut));
X=[px;py;px2;py2;px3;py3];
[u,s,~]=svd(X);  % perform the SVD
lambda=diag(s).^2;  % produce diagonal variances
Y=u'*X;  % produce the principal components projection
figure(1)
plot(lambda/sum(lambda),'ko','Linewidth',2);
ylabel('Energy(percentage)')
xlabel('variance')
title('Energy percentage for POD','Fontsize',16)
set(gca,'Fontsize',16,'Xtick',0:5:25)
figure(2)
plot(Y(1,:)); hold on
plot(Y(2,:))
title('principle component','Fontsize',16)

set(gca,'Fontsize',10)
legend('principal component 1','principal component 2','Location','northwest')

%%
clc;clear all;close all
load('cam1_2.mat')
load('cam2_2.mat')
load('cam3_2.mat')
numFrames = size(vidFrames1_2,4);
py=zeros(1,numFrames);
px=zeros(1,numFrames);
for j = 1:numFrames
    X= vidFrames1_2(:,:,:,j);
    X_dgray = double(rgb2gray(X));%we can handle double data in black and white color domain
    X_dgray(:,1:290)=0;%eliminate the range that does not contain the mass
    X_dgray(:,431:end)=0;
    X_dgray([1:203 427:end],291:430)=0;
    Mx = max(X_dgray(:));
    [yrange,xrange] = find(X_dgray>0.95*Mx);
    py(j)=mean(yrange);
    px(j)=mean(xrange);
end
[~,idx] = max(py(1:40));
py=py(idx:end);
px=px(idx:end);

numFrames2 = size(vidFrames2_2,4);
py2=zeros(1,numFrames2);
px2=zeros(1,numFrames2);
for j = 1:numFrames2
    X2= vidFrames2_2(:,:,:,j);
    X2_dgray = double(rgb2gray(X2));%we can handle double data in black and white color domain
    X2_dgray(:,[1:215 448:end])=0;%eliminate the range that does not contain the mass
    X2_dgray([1:61 417:end],216:447)=0;
    Mx = max(X2_dgray(:));
    [yrange,xrange] = find(X2_dgray>0.95*Mx);
    py2(j)=mean(yrange);
    px2(j)=mean(xrange);
end

[~,idx] = max(py2(1:40));
py2=py2(idx:end);
px2=px2(idx:end);

numFrames3 = size(vidFrames3_2,4);
py3=zeros(1,numFrames3);
px3=zeros(1,numFrames3);
for j = 1:numFrames3
    X3= vidFrames3_2(:,:,:,j);
    X3_dgray = double(rgb2gray(X3));%we can handle double data in black and white color domain
    X3_dgray(:,[1:241 514:end])=0;%eliminate the range that does not contain the mass
    X3_dgray([1:189 345:end],242:513)=0;
    Mx = max(X3_dgray(:));
    [yrange,xrange] = find(X3_dgray>0.95*Mx);
    py3(j)=mean(xrange); %reverse the coordinate due to the camera angle
    px3(j)=mean(yrange);
end
[~,idx] = max(py3(1:40));
py3=py3(idx:end);
px3=px3(idx:end);

cut=min([length(py) length(py2) length(py3)]);
py=py(1:cut)-mean(py(1:cut));
px=px(1:cut)-mean(px(1:cut));
py2=py2(1:cut)-mean(py2(1:cut));
px2=px2(1:cut)-mean(px2(1:cut));
py3=py3(1:cut)-mean(py3(1:cut));
px3=px3(1:cut)-mean(px3(1:cut));
X=[px;py;px2;py2;px3;py3];
[u,s,~]=svd(X);  % perform the SVD
lambda=diag(s).^2;  % produce diagonal variances
Y=u'*X;  % produce the principal components projection
figure(1)
plot(lambda/sum(lambda),'ko','Linewidth',2);
ylabel('Energy(percentage)')
xlabel('variance')
title('Energy percentage for POD','Fontsize',16)
set(gca,'Fontsize',16,'Xtick',0:5:25)
figure(2)
plot(Y(1,:)); hold on
plot(Y(2,:))
plot(Y(3,:))
title('principle component','Fontsize',16)

set(gca,'Fontsize',10)
legend('principal component 1','principal component 2','principal component 3','Location','northwest')
%%
clc;clear all;close all
load('cam1_3.mat')
load('cam2_3.mat')
load('cam3_3.mat')
numFrames = size(vidFrames1_3,4);
py=zeros(1,numFrames);
px=zeros(1,numFrames);
for j = 1:numFrames
    X= vidFrames1_3(:,:,:,j);
    X_dgray = double(rgb2gray(X));%we can handle double data in black and white color domain
    X_dgray(:,[1:282 415:end])=0;%eliminate the range that does not contain the mass
    X_dgray([1:225 422:end],283:414)=0;
    Mx = max(X_dgray(:));
    [yrange,xrange] = find(X_dgray>0.9*Mx);
    py(j)=mean(yrange);
    px(j)=mean(xrange);
    
end

py=py(8:end);%manully allign 
px=px(8:end);

numFrames2 = size(vidFrames2_3,4);
py2=zeros(1,numFrames2);
px2=zeros(1,numFrames2);
for j = 1:numFrames2
    X2= vidFrames2_3(:,:,:,j);
    X2_dgray = double(rgb2gray(X2));%we can handle double data in black and white color domain
    X2_dgray(:,[1:220 446:end])=0;%eliminate the range that does not contain the mass
    X2_dgray([1:134 429:end],221:445)=0;
    Mx = max(X2_dgray(:));
    [yrange,xrange] = find(X2_dgray>0.9*Mx);
    py2(j)=mean(yrange);
    px2(j)=mean(xrange);
    
end

py2=py2(16:end);
px2=px2(16:end);

numFrames3 = size(vidFrames3_3,4);
py3=zeros(1,numFrames3);
px3=zeros(1,numFrames3);
for j = 1:numFrames3
    X3= vidFrames3_3(:,:,:,j);
    X3_dgray = double(rgb2gray(X3));%we can handle double data in black and white color domain
    X3_dgray(:,[1:246 519:end])=0;%eliminate the range that does not contain the mass
    X3_dgray([1:145 361:end],247:518)=0;
    Mx = max(X3_dgray(:));
    [yrange,xrange] = find(X3_dgray>0.9*Mx);
    py3(j)=mean(xrange); %reverse the coordinate due to the camera angle
    px3(j)=mean(yrange);
end

py3=py3(43:end);
px3=px3(43:end);

cut=min([length(py) length(py2) length(py3)]);
py=py(1:cut)-mean(py(1:cut));
px=px(1:cut)-mean(px(1:cut));
py2=py2(1:cut)-mean(py2(1:cut));
px2=px2(1:cut)-mean(px2(1:cut));
py3=py3(1:cut)-mean(py3(1:cut));
px3=px3(1:cut)-mean(px3(1:cut));
X=[px;py;px2;py2;px3;py3];
[u,s,~]=svd(X/sqrt(cut-1));  % perform the SVD
lambda=diag(s).^2;  % produce diagonal variances
Y=u'*X;  % produce the principal components projection
figure(1)
plot(lambda/sum(lambda),'ko','Linewidth',2);
ylabel('Energy(percentage)')
xlabel('variance')
title('Energy percentage for POD','Fontsize',16)
set(gca,'Fontsize',16,'Xtick',0:5:25)
figure(2)
plot(Y(1,:)); hold on
plot(Y(2,:))
plot(Y(3,:))
title('principle component','Fontsize',16)

set(gca,'Fontsize',10)
legend('principal component 1','principal component 2','principal component 3','Location','northwest')


%%
clc;clear all;close all
load('cam1_4.mat')
load('cam2_4.mat')
load('cam3_4.mat')
numFrames = size(vidFrames1_4,4);
py=zeros(1,numFrames);
px=zeros(1,numFrames);
for j = 1:numFrames
    X= vidFrames1_4(:,:,:,j);
    X_dgray = double(rgb2gray(X));%we can handle double data in black and white color domain
    X_dgray(:,[1:304 461:end])=0;%eliminate the range that does not contain the mass
    X_dgray([1:216 420:end],305:460)=0;
    Mx = max(X_dgray(:));
    [yrange,xrange] = find(X_dgray>0.9*Mx);
    py(j)=mean(yrange);
    px(j)=mean(xrange);
    
end

py=py(11:end);%manully allign 
px=px(11:end);

numFrames2 = size(vidFrames2_4,4);
py2=zeros(1,numFrames2);
px2=zeros(1,numFrames2);
for j = 1:numFrames2
    X2= vidFrames2_4(:,:,:,j);
    X2_dgray = double(rgb2gray(X2));%we can handle double data in black and white color domain
    X2_dgray(:,[1:220 418:end])=0;%eliminate the range that does not contain the mass
    X2_dgray([1:70 407:end],221:417)=0;
    Mx = max(X2_dgray(:));
    [yrange,xrange] = find(X2_dgray>0.9*Mx);
    py2(j)=mean(yrange);
    px2(j)=mean(xrange);
    
end
py2=py2(11:end);
px2=px2(11:end);

numFrames3 = size(vidFrames3_4,4);
py3=zeros(1,numFrames3);
px3=zeros(1,numFrames3);
for j = 1:numFrames3
    X3= vidFrames3_4(:,:,:,j);
    X3_dgray = double(rgb2gray(X3));%we can handle double data in black and white color domain
    X3_dgray(:,[1:265 533:end])=0;%eliminate the range that does not contain the mass
    X3_dgray([1:131 276:end],266:532)=0;
    Mx = max(X3_dgray(:));
    [yrange,xrange] = find(X3_dgray>0.9*Mx);
    py3(j)=mean(xrange); %reverse the coordinate due to the camera angle
    px3(j)=mean(yrange);
end

py3=py3(19:end);
px3=px3(19:end);

cut=min([length(py) length(py2) length(py3)]);
py=py(1:cut)-mean(py(1:cut));
px=px(1:cut)-mean(px(1:cut));
py2=py2(1:cut)-mean(py2(1:cut));
px2=px2(1:cut)-mean(px2(1:cut));
py3=py3(1:cut)-mean(py3(1:cut));
px3=px3(1:cut)-mean(px3(1:cut));
X=[px;py;px2;py2;px3;py3];
[u,s,~]=svd(X/sqrt(cut-1));  % perform the SVD
lambda=diag(s).^2;  % produce diagonal variances
Y=u'*X;  % produce the principal components projection
figure(1)
plot(lambda/sum(lambda),'ko','Linewidth',2);
ylabel('Energy(percentage)')
xlabel('variance')
title('Energy percentage for POD','Fontsize',16)
set(gca,'Fontsize',16,'Xtick',0:5:25)
figure(2)
plot(Y(1,:)); hold on
plot(Y(2,:))
plot(Y(3,:))
title('principle component','Fontsize',16)
set(gca,'Fontsize',10)
legend('principal component 1','principal component 2','principal component 3','Location','northwest')

