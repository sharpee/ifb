tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script parameters
cellnum   = 1;  % number to assign to model cell for saving
saveall   = 1;  % 1 to save the stimuli and responses
repeats   = 1;  % 1 to generate a repeated stimulus
nat       = 1;  % 0 for white noise, 1 for van hateren
Nsamp     = 200000; % # stimulus/response pairs
NsampR    = 200000;  % # stimulus/response pairs in repeated data
Nrep      = 1;  % # repeats
Nd        = 16;    % # pixels per side of image
PSP       = 0.25;  % spike probability of model cell
noisenum  = .25;   % noise added to projections has stdev of noisenum*stdev of projections

% Choose the model
filter = 3;
% filter 1 is 9 blobs of random size and orientation
% filter 2 is a rotationally invariant curved edge
% filter 3 is a translationally invariant center-surround

% Make this directory if it doesn't exist
save_dir = './Model_cell_data/';
% Parameters defined
%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the model filters
Ndim = Nd^2;
if filter==1
    Nfilter = 9;
    F = zeros(Nd,Nd,Nfilter);
    ind=1;
    for n1 = 1:3
        for n2 = 1:3
            M=15;
            s1 = randi([2,M],1);
            s2 = randi([2,M],1);
            s3 = (2*rand()-1)*2;
            x0 = 4.5 + (n1-1)*4 + .5*(2*rand()-1);
            y0 = 4.5 + (n2-1)*4 + .5*(2*rand()-1);
            for i = 1:Nd
                for j = 1:Nd
                    F(i,j,ind) = exp(-[i-x0,j-y0]*inv([s1,s3;s3,s2])*[i-x0;j-y0] + .1*(2*rand()-1) );
                end
            end
            ind=ind+1;
        end
    end
    figure
    for n=1:Nfilter
        subplot(3,3,n)
        imagesc(F(:,:,n))
        colormap(hot)
        colorbar
        axis square
    end
    figure
    imagesc(sum(F,3));
    F_new = zeros(Nd^2,Nfilter);
    for n=1:Nfilter
        F_new(:,n)=reshape(F(:,:,n),[Nd^2,1]);
        F_new(:,n) = F_new(:,n)./norm(F_new(:,n));
    end
    dots=zeros(Nfilter);
    for i=1:Nfilter
        for j=1:Nfilter
            dots(i,j) = dot(F_new(:,i),F_new(:,j));
        end
    end
    clear F;
    F = F_new;
elseif filter==2
    Nfilter = 8;
    F = zeros(Nd,Nd,Nfilter);
    x0 = 8.5;
    y0 = 8.5;
    ind=1;
    for n1 = 1:Nfilter
        t1 = 0+(n1-1)*pi/4;
        for i = 1:Nd
            for j = 1:Nd
                x1new = i-x0;
                x2new = j-y0;
                x1 = cos(t1)*x1new-sin(t1)*x2new;
                x2 = sin(t1)*x1new+cos(t1)*x2new;
                F(i,j,ind) = exp(-x2^2/15-(x1-x2^2/5)^2/6)*sign(sin((x1-x2^2/5)));
            end
        end
        ind=ind+1;
    end
    
    figure
    for n=1:Nfilter
        subplot(3,3,n)
        imagesc(F(:,:,n))
        colormap(hot)
        colorbar
        axis square
    end
    figure
    imagesc(sum(F,3));
    F_new = zeros(Nd^2,Nfilter);
    for n=1:Nfilter
        F_new(:,n)=reshape(F(:,:,n),[Nd^2,1]);
        F_new(:,n) = F_new(:,n)./norm(F_new(:,n));
    end
    dots=zeros(Nfilter);
    for i=1:Nfilter
        for j=1:Nfilter
            dots(i,j) = dot(F_new(:,i),F_new(:,j));
        end
    end
    clear F;
    F = F_new;
elseif filter==3
    Nfilter = 4;
    F = zeros(Nd,Nd,Nfilter);
    x0list = [6.5,10.5];
    ind=1;
    for n1 = 1:length(x0list)
        x0 = x0list(n1);
        for n2 = 1:length(x0list)
            y0 = x0list(n2);
            for i = 1:Nd
                for j = 1:Nd
                    x1 = i-x0;
                    x2 = j-y0;
                    r = sqrt(x1^2+x2^2);
                    F(i,j,ind) = 1.1*exp(-r^2/6)-exp(-r^2/8);
                end
            end
            ind=ind+1;
        end
    end
    c1=max(max(max(F(:,:,1))),abs(min(min(F(:,:,1)))));
    clims=[-c1,c1];
    
    figure
    for n=1:Nfilter
        subplot(2,3,n)
        imagesc(F(:,:,n),clims)
        colormap(hot)
        colorbar
        axis square
    end
    figure
    imagesc(sum(F,3));
    F_new = zeros(Nd^2,Nfilter);
    for n=1:Nfilter
        F_new(:,n)=reshape(F(:,:,n),[Nd^2,1]);
        F_new(:,n) = F_new(:,n)./norm(F_new(:,n));
    end
    dots=zeros(Nfilter);
    for i=1:Nfilter
        for j=1:Nfilter
            dots(i,j) = dot(F_new(:,i),F_new(:,j));
        end
    end
    %dots
    clear F;
    F = F_new;
end
% Model filteres are created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%
% Load/create the stimulus
clear stimulus;
image_dir = '...?';
if nat==1
    % load natural movie stimuli here
else
    % Create a Gaussian stimulus 
    sigma = repmat(rand(1,Ndim),[Ndim,1]).*eye(Ndim);
    sigma = sigma+ rand(Ndim)/100;
    sigma = (sigma+sigma')/2;
    sigma = (2*rand(Ndim)-1)/10;
    sigma = sigma*sigma';
    R= chol(sigma);
    stimulus = repmat((2*rand(1,Ndim)-1)/1000000,[Nsamp+NsampR,1])  +  randn(Nsamp+NsampR,Ndim)*R;
    stimulus = reshape(stimulus,[1,Ndim*(Nsamp+NsampR)]);
    stimulus = stimulus + abs(min(min(stimulus)));
    stimulus = stimulus./max(max(stimulus));
    stimulus = stimulus.*255;
    stimulus = floor(stimulus);
    Rstimulus = stimulus(Ndim*Nsamp+1:Ndim*(Nsamp+NsampR));
    stimulus = stimulus(1:Ndim*Nsamp);
    fid = fopen([save_dir 'SN_' num2str(cellnum) '.raw'],'w');
    fwrite(fid,stimulus);
    fclose(fid);
end
Nsamples=length(stimulus)/Ndim;
stimulus = reshape(stimulus,[Ndim,Nsamples])';
stimulus = 2*(stimulus-255/2)/255;
disp('Unrepeated stimulus loaded');
% Stimuli is loaded
%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%
% Generate responses
proj = stimulus*F;
[d1,d2] = size(proj);
noiselevel = std(reshape(proj,[1,d1*d2]))*noisenum;
proj = proj + noiselevel*randn(d1, d2);
proj = proj';
temp = proj;
proj = max(proj);
avgspkprob=10; %just to set it to something large
a=mean(proj);
tally=0;
resp = zeros(Nsamp,1);
while abs(avgspkprob-PSP)>.05
    for i=1:Nsamp
        if proj(i) >= a
            resp(i) = 1;
        else
            resp(i) = 0;
        end
    end
    %resp = proj>=a;
    avgspkprob = mean(mean(resp));
    tally=tally+1;
    if (avgspkprob-PSP)>.05
        a=a+.01;
    elseif (avgspkprob-PSP)<.05
        a=a-.01;
    end
    if tally==1000
        disp('error')
        break;
    end
end
Filter_and_noiselevel=cell(1,2);
Filter_and_noiselevel{1}=F;
Filter_and_noiselevel{2} = noiselevel;
disp(['average spike probability: ' num2str(mean(mean(resp)))])

temp(temp<a)=-2*abs(a);
temp(temp>=a)=1;
temp(temp<1)=0;
tempsum=sum(temp);
new = zeros(1,Nsamp);
for i=1:Nsamp
    if tempsum(i)==1
        [Y,I]=max(temp(:,i));
        new(i) = I;
    end
end
new(new==0)=[];
figure
hist(new,Nfilter)

if saveall == 1
    save([save_dir 'Cell_' num2str(cellnum) '_filter_and_noise' int2str(marker) '.mat'],'Filter_and_noiselevel');
    fid = fopen([save_dir 'Cell_' num2str(cellnum) '_resp_noise' int2str(marker) '.isk'],'w');
    fprintf(fid,'%u\n',resp);
    fclose(fid);
end





if repeats==1
    %%%%%%%%%%%%%%%%%%%%
    % Generate repeated responses
    
    proj = stimulus*F;
    [d1,d2] = size(proj);
    noiselevel = std(reshape(proj,[1,d1*d2]))*noisenum;
    proj = repmat(proj,[1,1,Nrep]);
    proj = proj + noiselevel*randn(d1, d2, Nrep);
    proj = max(proj,[],2);
    proj = reshape(proj,Nsamp,Nrep);
    proj = reshape(proj,Nsamp*Nrep,1);
    resp = zeros(numel(proj),1);
    for i=1:numel(proj)
        if(proj(i) >= a)
            resp(i) = 1;
        else
            resp(i) = 0;
        end
    end
    if saveall == 1
        fid = fopen([save_dir 'Cell_' num2str(cellnum) '_resp_noise' int2str(marker) '_R.isk'],'w');
        fprintf(fid,'%u\n',resp);
        fclose(fid);
    end
    stimulus = repmat(stimulus,Nrep,1);
    stimconst = stimulus;
end
close all;
marker
toc










