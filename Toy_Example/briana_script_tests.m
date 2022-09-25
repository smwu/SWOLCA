% Testing code combines wtOFM and sLCA

% Input directory
in_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Toy_Example/";   
% Output directory 
out_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Toy_Example/";  
% Population scenario 2 (confounder), iteration 1
pop_data = importdata(strcat(in_dir, 'simdata_scen2_iter1.mat'));
% Sample scenario 14 (unequal sampling), iteration 1, sample 1
samp_data = importdata(strcat(in_dir, 'simdata_scen14_iter1_samp1.mat'));

wt_samp = samp_data.sample_wt;
food_sim = samp_data.X_data;
[n,p] = size(samp_data.X_data);
y_cc=samp_data.Y_data;

kappa = sum(samp_data.sample_wt) / n;
wt_c = wt_samp/kappa;
k_max=50;
d_max=max(food_sim(:));
d=max(food_sim);

wt_theta0=repmat(wt_c,[1,p]);
%vectorization of data
idz = repmat(1:p,n,1); idz = idz(:);
y_d = food_sim(:); lin_idx = sub2ind([p,d_max],idz,y_d);



%% OFMM params
rng(1, 'twister');
sp_k=50;
a_pi=ones(1,k_max)/sp_k;
pi_h=drchrnd(a_pi,1);

C_i=mnrnd(1,pi_h,n); [r, c]=find(C_i); gc=[r c];
    gc=sortrows(gc,1); C_i=gc(:,2);
n_C_i=sum(C_i);

%global theta0/1
eta=ones(1,d_max);
theta0=zeros(p,k_max,d_max);

for k=1:k_max
    for j=1:p
        dj=d(j);
        theta0(j,k,1:dj)=drchrnd(eta(1:dj),1);
    end
end

state=samp_data.true_Si; S=length(unique(state));
pcov=k_max+S;
pdem=S;
X=zeros(n,pcov);
mu0=normrnd(0,1,[pcov,1]);
sig0=1./gamrnd(5/2,2/5,[pcov,1]); %shape=5/2, scale=5/2
Sig0=diag(sig0);
xi_0=mvnrnd(mu0,Sig0);
xi_iter=xi_0;

%subpopulation design matrix: w_sid
w_sid=zeros(n,S); 

for s=1:S
    w_sid(state==(s),s)=1;
end

ep_kp=zeros(n,k_max); 

for k=1:k_max
   w_ip=zeros(n,k_max);
   w_ip(:,k)=ones(n,1);
   W_temp=[w_sid w_ip];
   phi_temp=normcdf(W_temp*transpose(xi_iter));
        phi_temp(phi_temp==1)=1-1e-10;
        phi_temp(phi_temp==0)=1e-10;
   probit_kp=y_cc.*log(phi_temp)+(1-y_cc).*log(1-phi_temp);
   ep_kp(:,k)=exp(probit_kp);
end


nrun=9;
burn=0;
thin=1;
pi_out = zeros(nrun/thin, k_max);
theta0_out = zeros(nrun/thin, p, k_max, d_max);
ci_out = zeros(nrun/thin, n);
xi_out = zeros(nrun/thin, pcov);
loglik_lca = zeros(nrun/thin, 1);
z_probit = zeros(n, 1);  % temp storage

for iter=1:nrun
    %% -- update pi_h -- %%
     for h=1:k_max
%              n_C_i(h)=sum(C_i==h);
        n_C_i(h)=sum(wt_c(C_i==h));
     end

    a_pih=a_pi+n_C_i;
    pi_h=drchrnd(a_pih,1);

%   %% -- Ci ~multinomial(pi_h) -- %%
    Cp_k=zeros(n,k_max);
    for l=1:k_max
        t0h=reshape(theta0(:,l,:),p,d_max);
        tmpmat0=reshape(t0h(lin_idx),[n,p]);
        Cp_k(:,l)=pi_h(l)*prod(tmpmat0,2).*ep_kp(:,l);
    end

    probCi = bsxfun(@times,Cp_k,1./(sum(Cp_k,2)));
    log_lca=log(sum(Cp_k,2));
        x_ci=mnrnd(1,probCi); [r, c]=find(x_ci); x_gc=[r c];
        x_gc=sortrows(x_gc,1); C_i=x_gc(:,2);

    % - update theta - %
    dmat0=zeros(p,d_max);
    for k=1:k_max
        C_is=repmat(C_i,[1,p]);
%          ph0 = (C_is==k); %subj's in global cluster h
           ph0=(C_is==k).*wt_theta0; %global cluster with weights
        for c = 1:d_max
             dmat0(:,c) = sum((food_sim==c).*ph0)';
        end
        for j=1:p
            a_tn0=eta+dmat0(j,:);
            theta0(j,k,:) = drchrnd(a_tn0,1);
        end
    end

    % -- RESPONSE MODEL PARAMETERS UPDATE -- %
    % Wup=[w_sid x_ci(:,2:end)]; %covariate matrix with state/global
    Wup=[w_sid x_ci];
      %add cluster L to referent group
    pcov=size(Wup,2);
    %create latent z_probit
        %create truncated normal for latent z_probit model
        WXi_now=Wup*transpose(xi_iter(1:pcov));
    %truncation for cases (0,inf)
    z_probit(y_cc==1)=truncnormrnd(1,WXi_now(y_cc==1),1,0,inf);
    %truncation for controls (-inf,0)
    z_probit(y_cc==0)=truncnormrnd(1,WXi_now(y_cc==0),1,-inf,0);

    %control extremes;
    if sum(z_probit==Inf)>0
        z_probit(z_probit==Inf)=norminv(1-1e-10);
    end 
    if sum(z_probit==-Inf)>0
        z_probit(z_probit==-Inf)=norminv(1e-10);
    end

    % Response Xi(B) update
    sig0up=Sig0(1:pcov,1:pcov);
    % xi_0up=xi_0(1:pcov);
    mu0up=mu0(1:pcov);
    W_tilde = diag(wt_c);
    xi_sig_up=inv(sig0up)+(transpose(Wup)*W_tilde *Wup);
    xi_mu_up2=(sig0up\mu0up)+(transpose(Wup)*W_tilde*z_probit);
    % xi_mu_up2=(sig0up\transpose(xi_0up))+(transpose(Wup)*z_probit);
    xi_mu_up=xi_sig_up\xi_mu_up2;
    xi_up=mvnrnd(xi_mu_up,inv(xi_sig_up));
    xi_iter(1:length(xi_up))=xi_up;

    for k=1:k_max
       w_ip=zeros(n,k_max);
       w_ip(:,k)=ones(n,1);
       W_temp=[w_sid w_ip];
       phi_temp=normcdf(W_temp*transpose(xi_iter));
            phi_temp(phi_temp==1)=1-1e-10;
            phi_temp(phi_temp==0)=1e-10;
       lprobit_kp=y_cc.*log(phi_temp)+(1-y_cc).*log(1-phi_temp);
       ep_kp(:,k)=exp(lprobit_kp);
    end

   if mod(iter,thin)==0
        pi_out(iter/thin,:)=pi_h;
        ci_out(iter/thin,:)=C_i;
        theta0_out(iter/thin,:,1:size(theta0,2),:)=theta0;
        xi_out(iter/thin,:)=xi_iter;
        loglik_lca(iter/thin)=sum(log_lca);
   end
    
    % RELABELLING STEP TO ENCOURAGE MIXING %
    if mod(iter,10)==0
        new_order=randperm(k_max);
        newCi=C_i;
        for h=1:k_max
            newCi(C_i==h)=new_order(h);
        end
        C_i=newCi;
        theta0=theta0(:,new_order,:);
        ep_kp=ep_kp(:,new_order);
    end
end


nrun=20;
burn=0;
thin=1;
pi_out = zeros(nrun/thin, k_max);
theta0_out = zeros(nrun/thin, p, k_max, d_max);
ci_out = zeros(nrun/thin, n);
xi_out = zeros(nrun/thin, pcov);
loglik_lca = zeros(nrun/thin, 1);
z_probit = zeros(n, 1);  % temp storage

for iter=1:nrun
    %% -- update pi_h -- %%
     for h=1:k_max
%              n_C_i(h)=sum(C_i==h);
        n_C_i(h)=sum(wt_c(C_i==h));
     end

    a_pih=a_pi+n_C_i;
    pi_h=drchrnd(a_pih,1);

%   %% -- Ci ~multinomial(pi_h) -- %%
    Cp_k=zeros(n,k_max);
    for l=1:k_max
        t0h=reshape(theta0(:,l,:),p,d_max);
        tmpmat0=reshape(t0h(lin_idx),[n,p]);
        Cp_k(:,l)=pi_h(l)*prod(tmpmat0,2).*ep_kp(:,l);
    end

    probCi = bsxfun(@times,Cp_k,1./(sum(Cp_k,2)));
    log_lca=log(sum(Cp_k,2));
        x_ci=mnrnd(1,probCi); [r, c]=find(x_ci); x_gc=[r c];
        x_gc=sortrows(x_gc,1); C_i=x_gc(:,2);

    % - update theta - %
    dmat0=zeros(p,d_max);
    for k=1:k_max
        C_is=repmat(C_i,[1,p]);
%          ph0 = (C_is==k); %subj's in global cluster h
           ph0=(C_is==k).*wt_theta0; %global cluster with weights
        for c = 1:d_max
             dmat0(:,c) = sum((food_sim==c).*ph0)';
        end
        for j=1:p
            a_tn0=eta+dmat0(j,:);
            theta0(j,k,:) = drchrnd(a_tn0,1);
        end
    end

    % -- RESPONSE MODEL PARAMETERS UPDATE -- %
    % Wup=[w_sid x_ci(:,2:end)]; %covariate matrix with state/global
    Wup=[w_sid x_ci];
      %add cluster L to referent group
    pcov=size(Wup,2);
    %create latent z_probit
        %create truncated normal for latent z_probit model
        WXi_now=Wup*transpose(xi_iter(1:pcov));
    %truncation for cases (0,inf)
    z_probit(y_cc==1)=truncnormrnd(1,WXi_now(y_cc==1),1,0,inf);
    %truncation for controls (-inf,0)
    z_probit(y_cc==0)=truncnormrnd(1,WXi_now(y_cc==0),1,-inf,0);

    %control extremes;
    if sum(z_probit==Inf)>0
        z_probit(z_probit==Inf)=norminv(1-1e-10);
    end 
    if sum(z_probit==-Inf)>0
        z_probit(z_probit==-Inf)=norminv(1e-10);
    end

    % Response Xi(B) update
    sig0up=Sig0(1:pcov,1:pcov);
    % xi_0up=xi_0(1:pcov);
    mu0up=mu0(1:pcov);
    W_tilde = diag(wt_c);
    xi_sig_up=inv(sig0up)+(transpose(Wup)*W_tilde *Wup);
    xi_mu_up2=(sig0up\mu0up)+(transpose(Wup)*W_tilde*z_probit);
    % xi_mu_up2=(sig0up\transpose(xi_0up))+(transpose(Wup)*z_probit);
    xi_mu_up=xi_sig_up\xi_mu_up2;
    xi_up=mvnrnd(xi_mu_up,inv(xi_sig_up));
    xi_iter(1:length(xi_up))=xi_up;

    for k=1:k_max
       w_ip=zeros(n,k_max);
       w_ip(:,k)=ones(n,1);
       W_temp=[w_sid w_ip];
       phi_temp=normcdf(W_temp*transpose(xi_iter));
            phi_temp(phi_temp==1)=1-1e-10;
            phi_temp(phi_temp==0)=1e-10;
       lprobit_kp=y_cc.*log(phi_temp)+(1-y_cc).*log(1-phi_temp);
       ep_kp(:,k)=exp(lprobit_kp);
    end

   if mod(iter,thin)==0
        pi_out(iter/thin,:)=pi_h;
        ci_out(iter/thin,:)=C_i;
        theta0_out(iter/thin,:,1:size(theta0,2),:)=theta0;
        xi_out(iter/thin,:)=xi_iter;
        loglik_lca(iter/thin)=sum(log_lca);
   end
    
    % RELABELLING STEP TO ENCOURAGE MIXING %
    if mod(iter,10)==0
        new_order=randperm(k_max);
        newCi=C_i;
        for h=1:k_max
            newCi(C_i==h)=new_order(h);
        end
        C_i=newCi;
        theta0=theta0(:,new_order,:);
        ep_kp=ep_kp(:,new_order);
    end
end

pi_burn=pi_out((burn/thin)+1:end,:);
theta0_burn=theta0_out((burn/thin)+1:end,:,:,:);
ci_burn=ci_out((burn/thin)+1:end,:);
k_med=median(sum(pi_burn>0.05,2));

k_max=k_med;

%% Fixed sampler
sp_k=k_max;
a_pi=ones(1,k_max)/sp_k;
pi_h=drchrnd(a_pi,1);

C_i=mnrnd(1,pi_h,n); [r, c]=find(C_i); gc=[r c];
    gc=sortrows(gc,1); C_i=gc(:,2);
n_C_i=sum(C_i);


%global theta0/1
eta=ones(1,d_max);
theta0=zeros(p,k_max,d_max);

for k=1:k_max
    for j=1:p
        dj=d(j);
        theta0(j,k,1:dj)=drchrnd(eta(1:dj),1);
    end
end

state=samp_data.true_Si; S=length(unique(state));
pcov=k_max+S;
pdem=S;
X=zeros(n,pcov);
mu0=normrnd(0,1,[pcov,1]);
sig0=1./gamrnd(5/2,2/5,[pcov,1]); %shape=5/2, scale=5/2
Sig0=diag(sig0);
xi_0=mvnrnd(mu0,Sig0);
xi_iter=xi_0;

%subpopulation design matrix: w_sid
w_sid=zeros(n,S); 

for s=1:S
    w_sid(state==(s),s)=1;
end

ep_kp=zeros(n,k_max); 

for k=1:k_max
   w_ip=zeros(n,k_max);
   w_ip(:,k)=ones(n,1);
   W_temp=[w_sid w_ip];
   phi_temp=normcdf(W_temp*transpose(xi_iter));
        phi_temp(phi_temp==1)=1-1e-10;
        phi_temp(phi_temp==0)=1e-10;
   probit_kp=y_cc.*log(phi_temp)+(1-y_cc).*log(1-phi_temp);
   ep_kp(:,k)=exp(probit_kp);
end

pi_out = zeros(nrun/thin, k_max);
theta0_out = zeros(nrun/thin, p, k_max, d_max);
ci_out = zeros(nrun/thin, n);
xi_out = zeros(nrun/thin, pcov);
loglik_lca = zeros(nrun/thin, 1);
z_probit = zeros(n, 1);  % temp storage
for iter=1:nrun
    %% -- update pi_h -- %%
     for h=1:k_max
%              n_C_i(h)=sum(C_i==h);
        n_C_i(h)=sum(wt_c(C_i==h));
     end

    a_pih=a_pi+n_C_i;
    pi_h=drchrnd(a_pih,1);

%   %% -- Ci ~multinomial(pi_h) -- %%
    Cp_k=zeros(n,k_max);
    for l=1:k_max
        t0h=reshape(theta0(:,l,:),p,d_max);
        tmpmat0=reshape(t0h(lin_idx),[n,p]);
        Cp_k(:,l)=pi_h(l)*prod(tmpmat0,2).*ep_kp(:,l);
    end

    probCi = bsxfun(@times,Cp_k,1./(sum(Cp_k,2)));
    log_lca=log(sum(Cp_k,2));
        x_ci=mnrnd(1,probCi); [r, c]=find(x_ci); x_gc=[r c];
        x_gc=sortrows(x_gc,1); C_i=x_gc(:,2);

    % - update theta - %
    dmat0=zeros(p,d_max);
    for k=1:k_max
        C_is=repmat(C_i,[1,p]);
%          ph0 = (C_is==k); %subj's in global cluster h
           ph0=(C_is==k).*wt_theta0; %global cluster with weights
        for c = 1:d_max
             dmat0(:,c) = sum((food_sim==c).*ph0)';
        end
        for j=1:p
            a_tn0=eta+dmat0(j,:);
            theta0(j,k,:) = drchrnd(a_tn0,1);
        end
    end

    % -- RESPONSE MODEL PARAMETERS UPDATE -- %
    % Wup=[w_sid x_ci(:,2:end)]; %covariate matrix with state/global
    Wup=[w_sid x_ci];
      %add cluster L to referent group
    pcov=size(Wup,2);
    %create latent z_probit
        %create truncated normal for latent z_probit model
        WXi_now=Wup*transpose(xi_iter(1:pcov));
    %truncation for cases (0,inf)
    z_probit(y_cc==1)=truncnormrnd(1,WXi_now(y_cc==1),1,0,inf);
    %truncation for controls (-inf,0)
    z_probit(y_cc==0)=truncnormrnd(1,WXi_now(y_cc==0),1,-inf,0);

    %control extremes;
    if sum(z_probit==Inf)>0
        z_probit(z_probit==Inf)=norminv(1-1e-10);
    end 
    if sum(z_probit==-Inf)>0
        z_probit(z_probit==-Inf)=norminv(1e-10);
    end

    % Response Xi(B) update
    sig0up=Sig0(1:pcov,1:pcov);
    % xi_0up=xi_0(1:pcov);
    mu0up=mu0(1:pcov);
    W_tilde = diag(wt_c);
    xi_sig_up=inv(sig0up)+(transpose(Wup)*W_tilde *Wup);
    xi_mu_up2=(sig0up\mu0up)+(transpose(Wup)*W_tilde*z_probit);
    % xi_mu_up2=(sig0up\transpose(xi_0up))+(transpose(Wup)*z_probit);
    xi_mu_up=xi_sig_up\xi_mu_up2;
    xi_up=mvnrnd(xi_mu_up,inv(xi_sig_up));
    xi_iter(1:length(xi_up))=xi_up;

    for k=1:k_max
       w_ip=zeros(n,k_max);
       w_ip(:,k)=ones(n,1);
       W_temp=[w_sid w_ip];
       phi_temp=normcdf(W_temp*transpose(xi_iter));
            phi_temp(phi_temp==1)=1-1e-10;
            phi_temp(phi_temp==0)=1e-10;
       lprobit_kp=y_cc.*log(phi_temp)+(1-y_cc).*log(1-phi_temp);
       ep_kp(:,k)=exp(lprobit_kp);
    end

   if mod(iter,thin)==0
        pi_out(iter/thin,:)=pi_h;
        ci_out(iter/thin,:)=C_i;
        theta0_out(iter/thin,:,1:size(theta0,2),:)=theta0;
        xi_out(iter/thin,:)=xi_iter;
        loglik_lca(iter/thin)=sum(log_lca);
   end
    
    % RELABELLING STEP TO ENCOURAGE MIXING %
    if mod(iter,10)==0
        new_order=randperm(k_max);
        newCi=C_i;
        for h=1:k_max
            newCi(C_i==h)=new_order(h);
        end
        C_i=newCi;
        theta0=theta0(:,new_order,:);
        ep_kp=ep_kp(:,new_order);
    end
end

%% Post-process
pi_burn=pi_out((burn/thin)+1:end,:);
theta0_burn=theta0_out((burn/thin)+1:end,:,:,:);
ci_burn=ci_out((burn/thin)+1:end,:);
xi_burn = xi_out((burn/thin)+1:end,:);
loglik_lcaburn=loglik_lca((burn/thin)+1:end);

m=size(pi_burn,1);
k_med=median(sum(pi_burn>0.05,2));
pd=pdist(transpose(ci_burn),'hamming'); %prcnt differ
Zci=linkage(pd,'complete');
% cdiff=squareform(pd); %Dij
% Zci=linkage(cdiff,'complete');

k0=k_med;
clust0 = cluster(Zci,'maxclust',k_med); %choose k0=5;


%Ordered MCMC set
ci_relabel=zeros(m,k0); %ci_relabelin=zeros(m,k_in); 
for l=1:k0
    ci_relabel(:,l) = mode(ci_burn(:,clust0==l),2);
end

pi_order=zeros(m,k0);
theta0_order=zeros(m,p,k0,d_max);
xi_order=zeros(m,k0+s);
for iter=1:m
    iter_order=ci_relabel(iter,:);
    pi_h1=pi_burn(iter,iter_order);
    % Normalize for each iteration before taking the median
    pi_order(iter,:)=pi_h1/sum(pi_h1);  %% normalization
    theta0_order(iter,:,:,:)=theta0_burn(iter,:,iter_order,:);
    s_iterorder=pdem+iter_order;
    xi_order(iter,:)=[xi_burn(iter,1:pdem) xi_burn(iter,s_iterorder)];
end


%% Analysis

theta0_med=reshape(median(theta0_order),[p,k0,d_max]);
% pi_med=median(pi_order);  % median over cols, i.e., median over iters
        % pi_med=pi_med/sum(pi_med);
        %     % median, then normalize
        %     temp = [0.2 0.4 0.3; 0.1 0.5 0.2];
        %     med_temp = median(temp);
        %     med_temp ./ sum(med_temp)
        %     % normalize, then median
        %     norm_temp = temp ./ sum(temp, 2);
        %     median(norm_temp)
        %     % Gives different results!!

[p,k0,d]=size(theta0_med);

%%Modal pattern of Global Clusters %%
[val0,ind0]=max(theta0_med,[],3);
t_ind0=transpose(ind0);
[t0,ia,ic]=unique(t_ind0,'rows');
k0=length(ia);

        pi_norm=pi_order(:,ia) ./ sum(pi_order(:,ia), 2);
        pi_med=median(pi_norm);
        theta0_norm=theta0_order(:,:,ia,:) ./ sum(theta0_order(:,:,ia,:), 4);
        theta0_medi=reshape(median(theta0_norm),[p,k0,d_max]);
        xi_med = median(xi_order(:, [1:S (S+ia')]));

% pi_med=pi_med(ia)/sum(pi_med(ia));
% theta0_med=theta0_med(:,ia,:);
% theta0_medi=theta0_med./sum(theta0_med,3);

ep_kpl=zeros(n,k0);
for k=1:k0
   w_ip=zeros(n,k0);
   w_ip(:,k)=ones(n,1);
   W_temp=[w_sid w_ip];
   phi_temp=normcdf(W_temp*transpose(xi_med));
        phi_temp(phi_temp==0)=1e-10;
        phi_temp(phi_temp==1)=1-1e-10;
   probit_kp=y_cc.*log(phi_temp)+(1-y_cc).*log(1-phi_temp);
   ep_kpl(:,k)=exp(probit_kp);
end

Cp_med=zeros(n,k0);
for k=1:k0
    t0h=reshape(theta0_medi(:,k,:),p,d_max);
    tmpmat0=reshape(t0h(lin_idx),[n,p]);
    Cp_med(:,k)=pi_med(k)*prod(tmpmat0,2).*ep_kpl(:,k);
end
pCi = bsxfun(@times,Cp_med,1./(sum(Cp_med,2)));
[z_val,z_max]=max(pCi,[],2);
loglca_mean=log(sum(Cp_med,2));

Wmed=[w_sid dummyvar(z_max)]; %covariate matrix with state/global
%add cluster 1 to referent group
WXi_med=Wmed*transpose(xi_med);     
phiWXimean=normcdf(WXi_med);
% - pull extremes off the bounds - %
phiWXimean(phiWXimean==0)=1e-10;
phiWXimean(phiWXimean==1)=1-1e-10;
py_pred=phiWXimean;
y_mse=immse(py_pred,samp_data.true_Phi);

loglike_mean=y_cc.*log(phiWXimean) + (1-y_cc).*log(1-phiWXimean);
%DIC6 based only on probit model
dic_dem=-6*median(loglik_lcaburn)+4*sum(loglca_mean);

% pi_med(1)
% theta0_medi(1,1,1)
% theta0_medi(4,2,2)
% xi_med(1)
% xi_med(4)
% k0
% ep_kpl(1,2)
% ep_kpl(200,1)
% sum(loglca_mean)
% z_max(1:5)
% z_val(1)
% z_val(200)
% phiWXimean(1)
% phiWXimean(200)
% sum(loglca_mean)




wt_N=sum(wt_samp);
posthoc_pi=[sum(sample_weight(z_max==1))/wt_N sum(sample_weight(z_max==2))/wt_N sum(sample_weight(z_max==3))/wt_N];



