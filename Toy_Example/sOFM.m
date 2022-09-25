
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supervised OFM
%% Programmer: BJKS     
%% Data: Simulated dataset        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

rng(1, 'twister');

load(['sim_sRPCdataB498.mat'])

%remove NaN rows
wt_samp = ones(length(true_y),1);
food_sim = sampledata;
y_cc=true_y;


%normalization constant

[n,p]=size(food_sim);
N=sum(wt_samp);
kappa = sum(wt_samp) / n;
wt_c=wt_samp/kappa;
    k_max=3;
    d_max=max(food_sim(:));
    d=max(food_sim);

wt_theta0=repmat(wt_c,[1,p]);
    %vectorization of data
     idz = repmat(1:p,n,1); idz = idz(:);
     y_d = food_sim(:); lin_idx = sub2ind([p,d_max],idz,y_d);

    %% SET UP PRIORS %%

        %pi_h for all classes
        sp_k=50;   %%%%CHANGE????????
a_pi=ones(1,k_max)/sp_k;
pi_h=drchrnd(a_pi,1);

        %phi - cluster index

        rr = unifrnd(0,1,[n,1]);
       pisums=[0 cumsum(pi_h)];
       C_i=zeros(n,1);
x_ci=zeros(n,k_max);
        for l = 1:k_max
            ind = rr>pisums(l) & rr<=pisums(l+1);
            C_i(ind==1) = l;
            x_ci(:,l)=ind;
        end 
n_C_i=sum(x_ci);


          %global theta0/1
     eta=ones(1,d_max);
    theta0=zeros(p,k_max,d_max);

    for k=1:k_max
        for j=1:p
            dj=d(j);
            theta0(j,k,1:dj)=drchrnd(eta(1:dj),1);
         end
    end

    
state=subpop_samp; S=length(unique(state));
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

    %% ------------ %%
    %% data storage %%
    %% ------------ %%
    nrun=250; burn=150;  thin=5;
    pi_out=zeros(nrun/thin,k_max);
   
    theta0_out=zeros(nrun/thin,p,k_max,d_max);
    ci_out=zeros(nrun/thin,n);
    
    xi_out = zeros(nrun/thin, pcov);
    loglik_lca = zeros(nrun/thin, 1);
    z_probit = zeros(n, 1);  % temp storage

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% POSTERIOR COMPUTATION %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    As=zeros(n,p);
    p_ij=zeros(n,p);
tic
for iter=1:nrun

   
        %% -- update pi_h -- %%
         for h=1:k_max
%              n_C_i(h)=sum(C_i==h);
            n_C_i(h)=sum(wt_c(C_i==h));
         end

        a_pih=a_pi+n_C_i;
        pi_h=drchrnd(a_pih,1);

      %% -- phi ~multinomial(pi_h) -- %%

  Cp_k=zeros(n,k_max);
for k=1:k_max
    t0h=reshape(theta0(:,k,:),p,d_max);
    tmpmat0=reshape(t0h(lin_idx),[n,p]);
    Cp_k(:,k)=pi_h(k)*prod(tmpmat0,2).*ep_kp(:,k);  %%% Should this be k instead of l?
end
probCi = bsxfun(@times,Cp_k,1./(sum(Cp_k,2)));
    log_lca=log(sum(Cp_k,2));
    x_ci=mnrnd(1,probCi); [r, c]=find(x_ci); x_gc=[r c];
    x_gc=sortrows(x_gc,1); C_i=x_gc(:,2);

%store global cluster index for postprocess/relabelling
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
                dj=d(j);
                a_tn0=eta(1:dj)+dmat0(j,1:dj);
                theta0(j,k,1:dj) = drchrnd(a_tn0,1);
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


%% RELABELLING STEP TO ENCOURAGE MIXING %%
    if mod(iter,10)==0
        new_order=randperm(k_max);
        newC_i=C_i;
        
        for k=1:k_max
            newC_i(C_i==k)=new_order(k);
        end
        
        C_i=newC_i;
        theta0=theta0(:,new_order,:);
        
        ep_kp=ep_kp(:,new_order);
    end




end
eltime=toc;
    pi_burn=pi_out((burn/thin)+1:end,:);
    theta0_burn=theta0_out((burn/thin)+1:end,:,:,:);
    ci_burn=ci_out((burn/thin)+1:end,:);
    
    xi_burn = xi_out((burn/thin)+1:end,:);
    loglik_lcaburn=loglik_lca((burn/thin)+1:end);

  
        figure; %check mixing of pi parameter
        plot(pi_burn)
%         saveas(gcf,'wtnhaneslow_pis.png')



% save('wtOFM_nhaneslow_MCMCout','pi_burn','ci_burn','theta0_burn','eltime','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESS: PASAPILIOPOULIS SWITCH %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=size(pi_burn,1); m_perm=size(theta0_burn,1); m_thin=m/m_perm;
k_med=median(sum(pi_burn>0.05,2));
pd=pdist(transpose(ci_burn),'hamming'); %prcnt differ
Zci=linkage(pd,'complete');
% cdiff=squareform(pd); %Dij
% Zci=linkage(cdiff,'complete');

figure; %save dendrogram of hierarchical clustering
dendrogram(Zci);
% saveas(gcf,'wtnhaneslow_dendrogram.png')

k0=k_med;
clust0 = cluster(Zci,'maxclust',k_med); %choose k0=5;


%Ordered MCMC set
ci_relabel=zeros(m,k0); %ci_relabelin=zeros(m,k_in); 
for l=1:k0
    ci_relabel(:,l) = mode(ci_burn(:,clust0==l),2);
end

pi_order=zeros(m,k0);
theta0_order=zeros(m,p,k0,d_max);

xi_order=zeros(m,k0+S);
for iter=1:m
    iter_order=ci_relabel(iter,:);
    pi_h1=pi_burn(iter,iter_order);
    pi_order(iter,:)=pi_h1/sum(pi_h1);
    theta0_order(iter,:,:,:)=theta0_burn(iter,:,iter_order,:);
    
    s_iterorder=pdem+iter_order;
    xi_order(iter,:)=[xi_burn(iter,1:pdem) xi_burn(iter,s_iterorder)];    
    
end

figure;
plot(pi_order)

theta0_med=reshape(median(theta0_order),[p,k0,d_max]);
pi_med=median(pi_order);
% pi_med=pi_med/sum(pi_med);

[p,k0,d]=size(theta0_med);


%%Modal pattern of Global Clusters %%
[val0,ind0]=max(theta0_med,[],3);
t_ind0=transpose(ind0);
[t0,ia,ic]=unique(t_ind0,'rows');
k0=length(ia);

pi_med=pi_med(ia)/sum(pi_med(ia));
theta0_med=theta0_med(:,ia,:);
theta0_medi=theta0_med./sum(theta0_med,3);
xi_med = median(xi_order(:, [1:S (S+ia')]));


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
y_mse=immse(py_pred,phi_WXtrue);
 
wt_N=sum(wt_samp);
posthoc_pi=[sum(wt_samp(z_max==1))/wt_N sum(wt_samp(z_max==2))/wt_N sum(wt_samp(z_max==3))/wt_N];

%  save(strcat('wOFM_simResults'),'pi_med','theta0_medi','sp_k','ind0','val0','eltime','z_max','Zci','-v7.3');
test_true=zeros(3,3);
for zi=1:3
    for ti=1:3
        test_true(zi,ti)=sum(z_max==zi & true_ci==ti);
    end
end
heatmap(test_true)

figure;
plot(theta0_order(:,1,1,1))
hold on
plot(theta0_order(:,1,1,2))
hold on
plot(theta0_order(:,1,1,3))
hold off

figure;
plot(theta0_burn(:,1,1,1))
hold on
plot(theta0_burn(:,1,1,2))
hold on 
plot(theta0_burn(:,1,1,3))
hold off

figure;
plot(xi_burn)

figure;
plot(xi_order)