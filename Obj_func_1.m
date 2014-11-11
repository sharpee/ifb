function f = Obj_func_1(p, F, stimulus, P)

global Nfilter
global Nthres

[Nsamples,Ndim]=size(stimulus);
P = double(P); %new

alpha = p(1:Nthres);
A = reshape(p(Nthres+1:length(p)),[Nfilter,Nthres]);

% calculate Cij's
C = repmat(alpha,[Nsamples,1]) + (stimulus*A);
C = 1./(1+exp(C));
C = prod(C,2);  % C = Q_0
%f = mean((P-C).^2);
D=1-C;  % D = Q_1
C = P.*log(C); % P=1 when y=0, so this picks out logC's for y=0
C(isnan(C))=0;
D = (1-P).*log(D); % 1-P=1 when y=1
D(isnan(D))=0;
f = -mean(C+D);  % negative b/c this is min(-likelihood)



