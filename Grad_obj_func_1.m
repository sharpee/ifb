function df = Grad_obj_func_1(p, F, stimulus, P)

global Nfilter
global Nthres

[Nsamples,Nf]=size(stimulus);
P = double(P);

alpha = p(1:Nthres);
A = reshape(p(Nthres+1:length(p)),[Nfilter,Nthres]);

C = repmat(alpha,[Nsamples,1]) + (stimulus*A);
C1 = 1./(1+exp(C));
D1 = 1-C1;
Q0 = prod(C1,2);  % C2 = Q_0
Q1 = 1-Q0;



d1 = repmat(P,[1,Nthres]).*D1;
d2 = -repmat((1-P).*Q0./(Q1),[1,Nthres]).*D1;
d2(isnan(d2))=0;
df1= mean(d1+d2);

df2 = zeros(1,Nf*Nthres);
for k=1:Nthres
    df2(1+(k-1)*Nf:k*Nf) = mean(repmat(d1(:,k)+d2(:,k),[1,Nf]) .* stimulus);
end
df = [df1,df2];

if sum(isnan(df))>0
    pause;
end

