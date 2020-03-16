clc;close all;clear all
%% test 1
e = spec_matrix('./test1',(3:8),(110:5:130));
g = spec_matrix('./test1',(9:14),(110:5:130));
j = spec_matrix('./test1',(15:20),(110:5:130));
[U,S,V] = svd([e g j],'econ');
lambda=diag(S).^2;
figure()
plot(lambda/sum(lambda),'ko','Linewidth',2);
ylabel('Energy(percentage)')
xlabel('mode')
title('Energy percentage for POD','Fontsize',16)
feature = (1:10);%look at the energy diagram
U = U(:,feature);
ne=size(e,2);
ng=size(g,2);
nj=size(j,2);
music = S*V';
es = music(feature,1:ne);
gs = music(feature,ne+1:ng+ne);
js = music(feature,ng+ne+1:ng+ne+nj);
me = mean(es,2);
mg = mean(gs,2);
mj = mean(js,2);
SW = 0;
for k = 1:ne
    SW = SW+(es(:,k)-me)*(es(:,k)-me)';
end
for k = 1:ng
    SW = SW+(gs(:,k)-mg)*(gs(:,k)-mg)';
end
for k = 1:nj
    SW = SW+(js(:,k)-mj)*(js(:,k)-mj)';
end
mut = mean([es gs js],2);
for i = 1:3
    SB = size(es,2)*(me-mut)*(me-mut)' + size(gs,2)*(mg-mut)*(mg-mut)' + size(js,2)*(mj-mut)*(mj-mut)';
end
[V2,D] = eig(SB,SW);
w=V2(:,1:2);%be careful
ve=w'*es;
vg=w'*gs;
vj=w'*js;
figure()
plot(ve(1,:),ve(2,:),'ro'); hold on;
plot(vg(1,:),vg(2,:),'go');
plot(vj(1,:),vj(2,:),'bo');
ylabel('projection onto w2')
xlabel('projection onto w1')
title('Projection of data','Fontsize',16)
thersholdx = (2803+2925)/2;%look at the graph
thersholdy = (3520+4782)/2;%look at the graph
xline(thersholdx)
yline(thersholdy)

teste = U'*spec_matrix('./test1t',(3:6),(110:5:130));
testg = U'*spec_matrix('./test1t',(7:10),(110:5:130));
testj = U'*spec_matrix('./test1t',(11:14),(110:5:130));
teste=w'*teste;
testg=w'*testg;
testj=w'*testj;

ce=0;
cg=0;
cj=0;
for i = 1 :size(teste,2)
    if(teste(1,i)>thersholdx)
       ce = ce+1; 
    end
    if(testg(2,i)<thersholdy)
        cg=cg+1;
    end
    if(testj(2,i)>thersholdy)
        cj=cj+1;
    end
end
acur = (ce+cg+cj)/(3*size(teste,2));

%% test 2
clc;close all;clear all
e = spec_matrix('./test2',(3:6),(30:5:50));
g = spec_matrix('./test2',(7:10),(30:5:50));
j = spec_matrix('./test2',(11:14),(30:5:50));
[U,S,V] = svd([e g j],'econ');
lambda=diag(S).^2;
figure()
plot(lambda/sum(lambda),'ko','Linewidth',2);
ylabel('Energy(percentage)')
xlabel('mode')
title('Energy percentage for POD','Fontsize',16)
feature = (1:10);%look at the energy diagram
U = U(:,feature);
ne=size(e,2);
ng=size(g,2);
nj=size(j,2);
music = S*V';
es = music(feature,1:ne);
gs = music(feature,ne+1:ng+ne);
js = music(feature,ng+ne+1:ng+ne+nj);
me = mean(es,2);
mg = mean(gs,2);
mj = mean(js,2);
SW = 0;
for k = 1:ne
    SW = SW+(es(:,k)-me)*(es(:,k)-me)';
end
for k = 1:ng
    SW = SW+(gs(:,k)-mg)*(gs(:,k)-mg)';
end
for k = 1:nj
    SW = SW+(js(:,k)-mj)*(js(:,k)-mj)';
end
mut = mean([es gs js],2);
for i = 1:3
    SB = size(es,2)*(me-mut)*(me-mut)' + size(gs,2)*(mg-mut)*(mg-mut)' + size(js,2)*(mj-mut)*(mj-mut)';
end
[V2,D] = eig(SB,SW);
w=V2(:,1:2);
ve=w'*es;
vg=w'*gs;
vj=w'*js;
figure()
plot(ve(1,:),ve(2,:),'ro'); hold on;
plot(vg(1,:),vg(2,:),'go');
plot(vj(1,:),vj(2,:),'bo');
ylabel('projection onto w2')
xlabel('projection onto w1')
title('Projection of data','Fontsize',16)
thersholdxr = (8433+9741)/2;%red
thersholdxb = (4750+4988)/2;%blue
thersholdy =(5704+6513)/2;
xline(thersholdxr)
xline(thersholdxb)
teste = U'*spec_matrix('./test2t',(3:6),(20:5:45));
testg = U'*spec_matrix('./test2t',(7:10),(40:5:65));
testj = U'*spec_matrix('./test2t',(11:14),(30:5:55));
teste=w'*teste;
testg=w'*testg;
testj=w'*testj;
% plot(teste(1,:),teste(2,:),'ro')
% plot(testg(1,:),testg(2,:),'go')
% plot(testj(1,:),testj(2,:),'bo')
% ylabel('y-2')
% xlabel('x-1')
ce=0;
cg=0;
cj=0;
for i = 1 :size(teste,2)
    if(teste(1,i)>thersholdxr)
       ce = ce+1; 
    end
    if(testg(1,i)<thersholdxb||testg(2,i)<thersholdy)
        cg=cg+1;
    end
    if(thersholdxb<testj(1,i)&& testj(1,i)<thersholdxr)
        cj=cj+1;
    end
end
acur = (ce+cg+cj)/(3*size(teste,2));

%% test 3
clc;close all;clear all
e = spec_matrix('./test3',(3:12),(30:5:50));
g = spec_matrix('./test3',(13:22),(30:5:50));
j = spec_matrix('./test3',(23:32),(30:5:50));
[U,S,V] = svd([e g j],'econ');
lambda=diag(S).^2;
figure()
plot(lambda/sum(lambda),'ko','Linewidth',2);
ylabel('Energy(percentage)')
xlabel('mode')
title('Energy percentage for POD','Fontsize',16)
feature = (1:2);%look at the energy diagram
U = U(:,feature);
ne=size(e,2);
ng=size(g,2);
nj=size(j,2);
music = S*V';
es = music(feature,1:ne);
gs = music(feature,ne+1:ng+ne);
js = music(feature,ng+ne+1:ng+ne+nj);
me = mean(es,2);
mg = mean(gs,2);
mj = mean(js,2);
SW = 0;
for k = 1:ne
    SW = SW+(es(:,k)-me)*(es(:,k)-me)';
end
for k = 1:ng
    SW = SW+(gs(:,k)-mg)*(gs(:,k)-mg)';
end
for k = 1:nj
    SW = SW+(js(:,k)-mj)*(js(:,k)-mj)';
end
mut = mean([es gs js],2);
for i = 1:3
    SB = size(es,2)*(me-mut)*(me-mut)' + size(gs,2)*(mg-mut)*(mg-mut)' + size(js,2)*(mj-mut)*(mj-mut)';
end
[V2,D] = eig(SB,SW);
w=V2(:,1:2);
ve=w'*es;
vg=w'*gs;
vj=w'*js;
figure()
plot(ve(1,:),ve(2,:),'ro'); hold on;
plot(vg(1,:),vg(2,:),'go');
plot(vj(1,:),vj(2,:),'bo');
ylabel('projection onto w2')
xlabel('projection onto w1')
title('Projection of data','Fontsize',16)
thersholdxb = (0.1921+0.1954)/2;%blue
thersholdy =(0.1635+0.1528)/2;
yline(thersholdxb)
xline(thersholdy)

teste = U'*spec_matrix('./test3',(3:12),(70:5:90));
testg = U'*spec_matrix('./test3',(13:22),(70:5:90));
testj = U'*spec_matrix('./test3',(23:32),(70:5:90));
teste=w'*teste;
testg=w'*testg;
testj=w'*testj;
% plot(teste(1,:),teste(2,:),'ro')
% plot(testg(1,:),testg(2,:),'go')
% plot(testj(1,:),testj(2,:),'bo')
% ylabel('y-2')
% xlabel('x-1')
ce=0;
cg=0;
cj=0;
for i = 1 :size(teste,2)
    if(teste(1,i)>thersholdxb && teste(2,i)<thersholdy)
       ce = ce+1; 
    end
    if(testg(1,i)<thersholdxb && testg(2,i)<thersholdy)
        cg=cg+1;
    end
    if(testj(2,i)>thersholdy)
        cj=cj+1;
    end
end
acur = (ce+cg+cj)/(3*size(teste,2));

function [A] = spec_matrix(path,numofsong,timeclips)
folder = dir(path);
A=[];
for i = numofsong
    filename=folder(i).name;
    [song, Fs] = audioread(strcat(path,'/',filename));
    song = song(:,1);
    song = resample(song,Fs/2,Fs);
    Fs = Fs/2;
    %p8 = audioplayer(song(Fs*100:Fs*110),Fs);playblocking(p8);
    for j=timeclips
        temp = song(Fs*j:Fs*(j+5));
        temp_sp = abs(spectrogram(temp));
        temp_sp = reshape(temp_sp,length(temp_sp(:)),1);
        A = [A temp_sp];
    end
end
end



