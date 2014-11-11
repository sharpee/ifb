global Nfilter
global Nthres

addpath(genpath('./Functional_basis_code/'))

modelcell=1;     % 1 if model cell, 0 if real data
cellnum = 1;     % model cell number (ignore if using real data)
Nd = 16;         % # pixels per side of image
Nsamp = 200000;  % length of stimulus
logicalOR = 0;   % 1 for logical OR, 0 for logical AND
Nthres = 4;      % # functional basis vectors, >= # STC basis vectors
Nrep = 1;        % # repeats of the stimuli


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load stimulus, responses and STC basis vectors %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if modelcell==1
    % load model cell filters, stimulus and responses
    load(['Cell_' num2str(cellnum) '_filter_and_noise']);
    F = Filter_and_noiselevel{1};
    [Ndim,Nfilter] = size(F);
    fid = fopen(['Model_cell_data/SN_' num2str(cellnum) '.raw'],'rb');
    stimulus=fread(fid,Nsamp*Ndim,'uint8');
    fclose(fid);
    Nsamples = length(stimulus)/Ndim;
    stimulus = reshape(stimulus,[Ndim,Nsamples])';
    stimulus = 2*(stimulus-255/2)/255;
    stimulus = repmat(stimulus,Nrep,1);
    if Nrep > 1
        fid = fopen(['Model_cell_data/Cell_' num2str(cellnum) '_resp_noise_R.isk'],'r');
    else
        fid = fopen(['Model_cell_data/Cell_' num2str(cellnum) '_resp_noise.isk'],'r');
    end
    resp = textscan(fid,'%u\n');
    fclose(fid);
    resp = resp{1,1};
    resp=resp(1:Nsamp*Nrep);
    % separate the spike-triggered ensemble
    Num=sum(resp);
    newstim = zeros(Num,Ndim);
    indx=1;
    for i=1:Nsamples*Nrep 
        if resp(i)==1
            newstim(indx,:) = stimulus(i,:);
            indx=indx+1;
        end
    end
    Nsamp = Nsamp*Nrep; 
    % perform STC (uncorrelated Gaussian stim)
    Cprior = cov(stimulus);
    Cspike = cov(newstim);
    delC = Cspike-Cprior;
    [evecs,evals]=eig(delC);
    [EV,inds] = sort(abs(diag(evals)));
    inds = flipud(inds);
    % plot eigenvalues and eigenvectors
    figure(1)
    plot(sort(diag(evals)),'.');
    axis square
    axis([-10 Ndim+10 min(diag(evals))*1.3 max(diag(evals))*1.3])
    drawnow
    EVECS = evecs(:,inds(1:Nfilter));
    figure(2)
    indx=Nfilter+1;
    c1=max(max(max(EVECS)),abs(min(min(EVECS))));
    clims=[-c1,c1];
    for i=1:Nfilter
        subplot(3,Nfilter,indx)
        a=EVECS(:,indx-Nfilter);
        B = reshape(a,Nd,Nd)';
        imagesc(B,clims);
        axis square
        colormap(hot)
        indx=indx+1;
    end
    drawnow
    basis = EVECS;
else
    
    % load stimulus, responses and the STC basis vectors
    % variables should be named stimulus, resp and basis
    
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the functional basis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AND and OR have same form, but reversed definitions of spike/no spike
if logicalOR==1
    P = 1-resp;
else
    P = resp;
end
P = P(1:Nsamp);
stimulus = stimulus*basis; % convert stimulus to reduced stimulus
Nf = Nfilter;

% initialize parameters
H = .1*(2*rand(1,Nthres)-1);
A = (2*rand(Nf,Nthres)-1);
A = reshape(A,[1,Nf*Nthres]);
A = [H,A];

% find funcitonal basis
[A,f1] = MaxLikelihood(A, @Obj_func_1, @Grad_obj_func_1, H, 1e-10, stimulus, P);
% put the parameters back into the original stimulus space for viewing
H = A(1:Nthres);
Afinal = reshape(A(Nthres+1:length(A)),[Nfilter,Nthres]);
for i=1:Nthres
    Afinal(:,i) = Afinal(:,i)/norm(Afinal(:,i));
end
B=basis*Afinal;
for i=1:Nthres
    B(:,i) = B(:,i)/norm(B(:,i));
end
% save stuff
if logicalOR==1
    if modelcell==1
        save(['Model_cell_data/Model_cell_' num2str(cellnum) '_Nthres' int2str(Nthres) '_OR_parameters.mat'],'A');
        save(['Model_cell_data/Model_cell_' num2str(cellnum) '_Nthres' int2str(Nthres) '_functional_basis_OR.mat'],'B');
    else
        % save the functional basis parameters
    end
else
    if modelcell==1
        save(['Model_cell_data/Model_cell_' num2str(cellnum) '_Nthres' int2str(Nthres) '_AND_parameters.mat'],'A');
        save(['Model_cell_data/Model_cell_' num2str(cellnum) '_Nthres' int2str(Nthres) '_functional_basis_AND.mat'],'B');
    else
        % save the functional basis parameters
    end
end
% plot the vectors
figure(2)
indx=1;
c1=max(max(max(B)),abs(min(min(B))));
clims=[-c1,c1];
for i=1:Nthres
    subplot(3,Nthres,2*Nthres+indx)
    a=B(:,indx);
    imagesc(reshape(a,Nd,Nd),clims)
    axis square
    colormap(hot)
    indx=indx+1;
end
drawnow


