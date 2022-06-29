library(Matrix)
set.seed(sample(1:100000,1))
library(dplyr)
library(MSGLasso)
library(MASS)

load("C:\\Users\\25110\\Desktop\\IVW&TSLS\\PAPER\\Replication Codes\\Replication Codes\\Simulation\\CorrelatedIV-Lemma-S-16\\locus_42.RData")
rownames(expsoure.org) <- 1:nrow(expsoure.org)

IVW_correct<-function(gamma.true, gamma.exp, gamma.gold, gamma.out, se.exp, se.gold, se.out,  etamean = 0.5, pthr = 5e-8, method = "IVW.original",seed = sample(1:100000,1), lambda = 0) {
  
  # denom, optimal, and seed are for carved IVW only
  
  z = qnorm(0.975, 0, 1)
  flip = which(gamma.exp < 0)
  gamma.exp[flip] = -gamma.exp[flip]
  gamma.out[flip] = -gamma.out[flip]
  gamma.true[flip] = -gamma.true[flip]
  
  if (method == "IVW.original"){ #classical IVW
    C_sel = qnorm(pthr/2,lower.tail = FALSE)
    ind_filter = which(abs(gamma.exp / se.exp) >= C_sel)
    numIV_sel = length(ind_filter)
    gamma.exp_sel = gamma.exp[ind_filter]
    gamma.out_sel = gamma.out[ind_filter]
    se.exp_sel = se.exp[ind_filter]
    se.out_sel = se.out[ind_filter]
    
    gamma.exp.pval_sel = 2 * pnorm(-abs(gamma.exp_sel) / se.exp_sel)
    
    beta = sum(gamma.exp_sel*gamma.out_sel *(1/se.out_sel^2)) / sum(gamma.exp_sel^2/se.out_sel^2)
    sd = sqrt(1 / sum(gamma.exp_sel^2/se.out_sel^2))
    p = pnorm(abs(beta/sd), lower.tail = F) * 2
    
  } 
  else if (method == "IVW.gold") { # gold standard IVW
    C_sel = qnorm(pthr/2,lower.tail = FALSE)
    ind_filter = which(abs(gamma.gold / se.gold) >= C_sel)
    numIV_sel = length(ind_filter)
    gamma.exp_sel = gamma.exp[ind_filter]
    gamma.out_sel = gamma.out[ind_filter]
    se.exp_sel = se.exp[ind_filter]
    se.out_sel = se.out[ind_filter]
    gamma.exp.pval_sel = 2 * pnorm(-abs(gamma.exp_sel) / se.exp_sel)
    
    beta = sum(gamma.exp_sel*gamma.out_sel*(1/se.out_sel^2)) / sum(gamma.exp_sel^2/se.out_sel^2)
    sd = sqrt(1 / sum(gamma.exp_sel^2/se.out_sel^2))
    p = pnorm(abs(beta/sd), lower.tail = F) * 2
    
  }   else if (method == "RIVW") { # carved IVW
    
    set.seed(seed)
    W = rnorm(length(gamma.exp), 0, etamean)
    
    # Step 2: Select significant SNPs:
    C_sel = qnorm(pthr/2,lower.tail = FALSE) #* sqrt(1 + eta^2)
    ind_filter = which(abs(gamma.exp/se.exp + W) >= C_sel)
    
    numIV_sel = length(ind_filter)
    gamma.exp_sel = gamma.exp[ind_filter]
    gamma.out_sel = gamma.out[ind_filter]
    se.exp_sel = se.exp[ind_filter]
    se.out_sel = se.out[ind_filter]
    gamma.exp.pval_sel = 2 * pnorm(-abs(gamma.exp_sel) / se.exp_sel)
    gamma.true_sel = gamma.true[ind_filter]
    
    
    # Step 3. Construct the unbiased carved estimator (also the UMVUE)
    alpha1 = (-C_sel - gamma.exp_sel/se.exp_sel) / etamean
    alpha2 = (C_sel - gamma.exp_sel/se.exp_sel) / etamean
    gamma.carve = gamma.exp_sel - (se.exp_sel/etamean) * ( (dnorm(alpha2) - dnorm(alpha1)) / (pnorm(alpha1) + 1 - pnorm(alpha2)) )
    sigma2.carve = (1 - ((alpha2*dnorm(alpha2) - alpha1*dnorm(alpha1)) / (1 - pnorm(alpha2) + pnorm(alpha1) ) - ((dnorm(alpha2) - dnorm(alpha1))/(1 - pnorm(alpha2) + pnorm(alpha1)))^2) / etamean^2 ) * se.exp_sel^2
    gamma.bias = gamma.true_sel - gamma.carve
    gamma.bias.avg = mean(abs(gamma.bias))
    
    beta = sum(gamma.carve*gamma.out_sel*(1/se.out_sel^2)) / sum((gamma.carve^2 - sigma2.carve)/se.out_sel^2)
    
    # estimation based on regression residuals
    RIVW.var = sum( (gamma.out_sel * gamma.carve - beta * (gamma.carve^2 - sigma2.carve) )^2 / se.out_sel^4) / (sum((gamma.carve^2 - sigma2.carve) / se.out_sel^2) )^2
    sd = sqrt(RIVW.var)
    
    p = pnorm(abs(beta/sd), lower.tail = F) * 2
  }  else if (method == "mr.divw") {
    if(lambda == 0) {
      tmp = mr.divw(gamma.exp, gamma.out, se.exp, se.out, diagnostics=FALSE)
    } else {
      # three-sample dIVW
      C_sel = qnorm(pthr/2,lower.tail = FALSE)
      ind_filter = which(abs(gamma.gold / se.gold) >= C_sel)
      numIV_sel = length(ind_filter)
      gamma.exp_sel = gamma.exp[ind_filter]
      gamma.out_sel = gamma.out[ind_filter]
      se.exp_sel = se.exp[ind_filter]
      se.out_sel = se.out[ind_filter]
      tmp = mr.divw(gamma.exp, gamma.out, se.exp, se.out, lambda = 0, diagnostics=FALSE)
    }
    
    beta = tmp$beta.hat
    sd = tmp$beta.se
    ind_filter = tmp$IV
    gamma.exp_sel = gamma.exp[ind_filter]
    se.exp_sel = se.exp[ind_filter]
    gamma.exp.pval_sel = 2 * pnorm(-abs(gamma.exp_sel) / se.exp_sel)
    numIV_sel = length(ind_filter)
    p = pnorm(abs(beta/sd), lower.tail = F) * 2
  }
  
  gamma.true_sel = gamma.true[ind_filter]
  kappa.sel = sum(gamma.true_sel^2 / se.exp_sel^2) / length(gamma.true_sel)
  
  kappa.est = sum(gamma.exp_sel^2 / se.exp_sel^2) / length(gamma.exp_sel) - 1
  
  lower = beta - z * sd
  upper = beta + z * sd
  length = upper - lower
  f=sum(gamma.exp_sel^2 / se.exp_sel^2)/ length(se.exp_sel) - 1
  
  bias = theta - beta
  prop_psmall = round(sum(gamma.exp.pval_sel < 5e-8 & gamma.exp.pval_sel > 5e-10)/numIV_sel * 100, 2)
  
  output = data.frame(beta, sd, p, lower, upper, length, f, numIV_sel, bias, prop_psmall,kappa.sel, kappa.est)
  
  
    return(list(res = output))
  
}



LD_clumping<-function(gamma_summary_x,gamma_summary_y,reference_panel,sd_gamma_x,scale_cluster,nchormosome=20,etamean=0.5,sel.pthr=5e-8,rerand=FALSE)
{ 
  
  snp_index=NULL
  cont=1
  for(i in (1:(nchormosome))){
    left=scale_cluster[i]
    right=scale_cluster[i+1]-1
    reference_panel_i=reference_panel[left:right,left:right]
    gamma_x=gamma_summary_x[left:right]
    gamma_y=gamma_summary_y[left:right]
    sd=sd_gamma_x[left:right]
    len=length(sd)
    W = rnorm(length(gamma_x), 0, etamean)
    
    # Step 2: Select significant SNPs:
    C_sel = qnorm(sel.pthr/2,lower.tail = FALSE) #* sqrt(1 + eta^2)
    
    if(rerand) { # recalculated p value when rerandomized
      tmpZ2 = gamma_x/sd + W
      tmpZ = gamma_x/sd
      
      tmpP2 = pnorm(abs(tmpZ2/sqrt(1 + etamean^2)),lower.tail = FALSE) * 2
      ind_filter = which(abs(tmpZ2) >= C_sel ) # & abs(tmpZ) >= C_sel this should not be used, which will biased the results
    } else {
      # use the original ones
      tmpZ =gamma_x/sd
      ind_filter = which(abs(tmpZ) >= C_sel)
    }
    
    
    
    if(length(ind_filter)==0 ) {
      next
    }
    
    
    # select the SNPs
    sd =sd[ind_filter]
    gamma_x = gamma_x[ind_filter]
    gamma_y = gamma_y[ind_filter]
    p_gamma= pnorm(abs(gamma_x/sd), lower.tail = F) *2
    marked=rep(1,length(gamma_x))
    while(sum(marked)>0)
    { 
      index=which.min(p_gamma)
      snp_index[cont]=scale_cluster[i]+ind_filter[index]-1
      cont=cont+1
      marked[index]=0
      for(k in (1:length(ind_filter)))
      {
        if((k!=index)&(reference_panel_i[ind_filter[k],ind_filter[index]]>0.01))
        {
          marked[k]=0
        }
        
      }
      p_gamma[index]=1
    }
  }
  if(is.null(snp_index))
    print("No valid IV")
  else return(snp_index)
}
LD_sigma_pruning<-function(gamma_summary_x,gamma_summary_y,reference_panel,sd_gamma_x,scale_cluster,nchormosome=20,etamean=0.5,sel.pthr=5e-8,rerand=FALSE)
{
  snp_index=NULL
  cont=1
  for(i in (1:nchormosome)){
    
    left=scale_cluster[i]
    right=scale_cluster[i+1]-1
    reference_panel_i=reference_panel[left:right,left:right]
    gamma_x=gamma_summary_x[left:right]
    gamma_y=gamma_summary_y[left:right]
    sd=sd_gamma_x[left:right]
    len=length(sd)
    W = rnorm(length(gamma_x), 0, etamean)
    
    # Step 2: Select significant SNPs:
    C_sel = qnorm(sel.pthr/2,lower.tail = FALSE) #* sqrt(1 + eta^2)
    
    if(rerand) { # recalculated p value when rerandomized
      tmpZ2 = gamma_x/sd + W
      tmpZ = gamma_x/sd
      
      tmpP2 = pnorm(abs(tmpZ2/sqrt(1 + etamean^2)),lower.tail = FALSE) * 2
      ind_filter = which(abs(tmpZ2) >= C_sel ) # & abs(tmpZ) >= C_sel this should not be used, which will biased the results
    } else {
      # use the original ones
      tmpZ =gamma_x/sd
      ind_filter = which(abs(tmpZ) >= C_sel)
    }
    
    
    
    if(length(ind_filter)==0 ) {
      next
    }
    
    # select the SNPs
    sd =sd[ind_filter]
    gamma_x = gamma_x[ind_filter]
    gamma_y = gamma_y[ind_filter]
    p_gamma= 1/(1+exp(-sd))-0.5
    marked=rep(1,length(gamma_x))
    while(sum(marked)>0)
    { 
      index=which.min(p_gamma)
      snp_index[cont]=scale_cluster[i]+ind_filter[index]-1
      cont=cont+1
      marked[index]=0
      for(k in (1:length(ind_filter)))
      {
        if((k!=index)&(reference_panel_i[ind_filter[k],ind_filter[index]]>0.01))
        {
          marked[k]=0
        }
        
      }
      p_gamma[index]=1
    }
    
  }
  if(is.null(snp_index))
    print("No valid IV")
  return(snp_index)
}
LD_pruning<-function(gamma_summary_x,gamma_summary_y,reference_panel,sd_gamma_x,scale_cluster,nchormosome=20,etamean=0.5,sel.pthr=5e-8,rerand=FALSE)
{
  snp_index=NULL
  cont=1
  EAF=expsoure.org$EAF
  for(i in (1:nchormosome)){
    left=scale_cluster[i]
    right=scale_cluster[i+1]-1
    reference_panel_i=reference_panel[left:right,left:right]
    gamma_x=gamma_summary_x[left:right]
    gamma_y=gamma_summary_y[left:right]
    sd=sd_gamma_x[left:right]
    len=length(sd)
    marked=rep(1,len)
    W = rnorm(length(gamma_x), 0, etamean)
    MAF=EAF[left:right]
    # Step 2: Select significant SNPs:
    C_sel = qnorm(sel.pthr/2,lower.tail = FALSE) #* sqrt(1 + eta^2)
    
    if(rerand) { # recalculated p value when rerandomized
      tmpZ2 = gamma_x/sd + W
      tmpZ = gamma_x/sd
      
      tmpP2 = pnorm(abs(tmpZ2/sqrt(1 + etamean^2)),lower.tail = FALSE) * 2
      ind_filter = which(abs(tmpZ2) >= C_sel ) # & abs(tmpZ) >= C_sel this should not be used, which will biased the results
    } else {
      # use the original ones
      tmpZ =gamma_x/sd
      ind_filter = which(abs(tmpZ) >= C_sel)
      
      
    }
    
    
    
    #if(length(ind_filter)==0 ) {
    #  next
    #}
    
    # select the SNPs
    sd =sd[ind_filter]
    gamma_x = gamma_x[ind_filter]
    gamma_y = gamma_y[ind_filter]
    marked=rep(1,length(gamma_x))
    if(length(ind_filter)<3)
      next
    for (m in (1:10))
    { 
      index=sample(1:length(gamma_x),1)
      if(marked[index]!=0){
        for(k in (1:length(gamma_x)))
        {
          if((k!=index)&(reference_panel_i[ind_filter[k],ind_filter[index]]>0.01))
          {
            if(MAF[ind_filter[k]]<MAF[ind_filter[index]])
              marked[k]=0
            else
              marked[index]=0
            
          }
          
        }}
      else
        next;
    }
    for (k in (1:length(gamma_x))){
      if(marked[k]!=0)
      {snp_index[cont]=scale_cluster[i]+ind_filter[k]-1
      cont=cont+1
      }}
  }
  if(is.null(snp_index))
    print("No valid IV")
  else return(snp_index)
}


Estimation<-function(gamma_true,gamma_summary_X,gamma_summary_Y,gamma_X_sd,gamma_Y_sd,inde.method, reference_panel,scale_cluster,used.method = c("RIVW","IVW","dIVW","RAPS","sRIVW"),etamean=0.5,sel.pthr=1,rerand=FALSE,gamma_summary_Z=rep(0,3963),gamma_Z_sd=rep(0,3693))
{
  if(inde.method=="clumping")
    snp_set=LD_clumping(gamma_summary_Z,gamma_summary_Y,reference_panel,gamma_Z_sd,scale_cluster,nchormosome=19,etamean,sel.pthr,rerand)
  if(inde.method=="Revised_Pruning")
    snp_set=LD_sigma_pruning(gamma_summary_Z,gamma_summary_Y,reference_panel,gamma_Z_sd,scale_cluster,nchormosome=19,etamean,sel.pthr,rerand)
  if(inde.method=="Pruning")
    snp_set=LD_pruning(gamma_summary_Z,gamma_summary_Y,reference_panel,gamma_Z_sd,scale_cluster,nchormosome=19,etamean,sel.pthr,rerand)
  if (!is.null(snp_set))
  {
    if(inde.method=="clumping")
      snp_set=LD_clumping(gamma_summary_Z,gamma_summary_Y,reference_panel,gamma_Z_sd,scale_cluster,nchormosome=19,etamean,sel.pthr,rerand)
    if(inde.method=="Revised_Pruning")
      snp_set=LD_sigma_pruning(gamma_summary_Z,gamma_summary_Y,reference_panel,gamma_Z_sd,scale_cluster,nchormosome=19,etamean,sel.pthr,rerand)
    if(inde.method=="Pruning")
      snp_set=LD_pruning(gamma_summary_Z,gamma_summary_Y,reference_panel,gamma_Z_sd,scale_cluster,nchormosome=19,etamean,sel.pthr,rerand)
    gamma_x=gamma_summary_X[snp_set]
    gamma_y=gamma_summary_Y[snp_set]
    gamma_z=gamma_summary_Z[snp_set]
    gamma_sd_x=gamma_X_sd[snp_set]
    gamma_sd_y=gamma_Y_sd[snp_set]
    gamma_sd_z=gamma_Z_sd[snp_set]
    final.out = list()
    gamma_true=gamma_true[snp_set]
    #RIVW
    #two sample using x information
    #USING Z to select
    RIVW = IVW_correct(gamma_true,gamma_z,gamma_x,gamma_y,gamma_sd_x,gamma_sd_z,gamma_sd_y,etamean = 0.5, pthr = 5e-5, method = "RIVW",seed = 99999, lambda = 0)
  
    final.out[["RIVW"]] = RIVW
    
    # dIVW
    dIVW = IVW_correct(gamma_true,gamma_z,gamma_x,gamma_y,gamma_sd_x,gamma_sd_z,gamma_sd_y,etamean = 0.5, pthr = 5e-8, method = "mr.divw",seed = 99999, lambda = 0)
    
    final.out[["dIVW"]] = dIVW
    
    # IVW (two sample MR results)
    
    IVW.original=IVW_correct(gamma_true,gamma_z,gamma_x,gamma_y,gamma_sd_x,gamma_sd_z,gamma_sd_y,etamean = 0.5, pthr = 5e-8, method = "IVW.original",seed = 99999, lambda = 0)
   
    final.out[["IVW"]] = IVW.original
    
    #threesample
    
    
    # dIVW
    
    dIVW_1 = IVW_correct(gamma_true,gamma_x,gamma_z,gamma_y,gamma_sd_x,gamma_sd_z,gamma_sd_y,etamean = 0.5, pthr = 1 ,method = "mr.divw",seed = 99999, lambda = 5.45)
    
    final.out[["threedIVW"]] = dIVW_1
    
    # IVW (two sample MR results)
    #selection from clumping x but selection bias from gammax
    IVW_1=IVW_correct(gamma_true,gamma_z,gamma_x,gamma_y,gamma_sd_x,gamma_sd_z,gamma_sd_y,etamean = 0.5, pthr = 5e-8, method = "IVW.gold",seed = 99999, lambda = 0)
   
    #selection from both clumping and thresholding x
    IVW_2=IVW_correct(gamma_true,gamma_x,gamma_z,gamma_y,gamma_sd_x,gamma_sd_z,gamma_sd_y,etamean = 0.5, pthr = 5e-8, method = "IVW.gold",seed = 99999, lambda = 0)
    final.out[["threeIVW_part"]] = IVW_1
    final.out[["threeIVW_all"]] = IVW_2
    return(list(res=final.out))
  }
  else
    return(-1)
}

#reference_panel=simBeta_diagonal(G_num=20,seed=sample(1:100000,1),num_row=scale,num_col=scale,0)

#geting GWAS sd directly
require(TwoSampleMR)
require(data.table)
require(dplyr)
require(mr.divw)
iter=1
LD_CL=list()
LD_R_PR=list()
LD_PR=list()
LD_R_PR_RE=list()
LD_PR_RE=list()
LD_CL_RE=list()

repeat_times=2000
jobs=seq(3,9)
setting<-function(job){
  
  if(job==1)
  {
    thre=1
  }
  
  if(job==2)
  {
    thre=1e-1
  }
  if(job==3)
  {
    thre=1
  }
  if(job==4)
  {
    thre=5e-1
  }
  if(job==5)
  {
    thre=5e-2
  }
  if(job==6)
  {
    thre=5e-3
  }
  if(job==7)
  {
    thre=5e-6
  }
  if(job==8)
  {
    thre=5e-7
  }
  if(job==9)
  {
    thre=5e-8
  }
  
  return(thre)
}
res_processing<-function(LD_CL,iter){
  LD_CL_beta_RIVW=NULL
  LD_CL_RIVW_F=NULL
  LD_CL_RIVW_nIV=NULL
  LD_CL_RIVW_sd=NULL

  for ( m in(1:(iter-1)) )
  {
    LD_CL_beta_RIVW[m]=LD_CL[[m]]$RIVW$RIVW$res$beta
    LD_CL_RIVW_F[m]=LD_CL[[m]]$RIVW$RIVW$res$f
    LD_CL_RIVW_nIV[m]=LD_CL[[m]]$RIVW$RIVW$res$numIV_sel
    LD_CL_RIVW_sd[m]=LD_CL[[m]]$RIVW$RIVW$res$sd

  }
  
  LD_CL_beta_IVW=NULL
  LD_CL_IVW_F=NULL
  LD_CL_IVW_nIV=NULL
  LD_CL_IVW_sd=NULL
  LD_CL_beta_IVW_3=NULL
  LD_CL_IVW_F_3=NULL
  LD_CL_IVW_nIV_3=NULL
  LD_CL_IVW_sd_3=NULL
  LD_CL_beta_IVW_31=NULL
  LD_CL_IVW_F_31=NULL
  LD_CL_IVW_nIV_31=NULL
  LD_CL_IVW_sd_31=NULL
  for ( m in(1:(iter-1)) )
  {
    LD_CL_beta_IVW[m]=LD_CL[[m]]$IVW$IVW$res$beta
    LD_CL_IVW_F[m]=LD_CL[[m]]$IVW$IVW$res$sd
    LD_CL_IVW_nIV[m]=LD_CL[[m]]$IVW$IVW$res$numIV_sel
    LD_CL_IVW_sd[m]=LD_CL[[m]]$IVW$IVW$res$sd
    LD_CL_beta_IVW_3[m]=LD_CL[[m]]$threeIVW_part$threeIVW_part$res$beta
    LD_CL_IVW_F_3[m]=LD_CL[[m]]$threeIVW_part$threeIVW_part$res$f
    LD_CL_IVW_nIV_3[m]=LD_CL[[m]]$threeIVW_part$threeIVW_part$res$numIV_sel
    LD_CL_IVW_sd_3[m]=LD_CL[[m]]$threeIVW_part$threeIVW_part$res$sd
    LD_CL_beta_IVW_31[m]=LD_CL[[m]]$threeIVW_all$threeIVW_all$res$beta
    LD_CL_IVW_F_31[m]=LD_CL[[m]]$threeIVW_all$threeIVW_all$res$f
    LD_CL_IVW_nIV_31[m]=LD_CL[[m]]$threeIVW_all$threeIVW_all$res$numIV_sel
    LD_CL_IVW_sd_31[m]=LD_CL[[m]]$threeIVW_all$threeIVW_all$res$sd
  }
  
  LD_CL_beta_dIVW=NULL
  LD_CL_dIVW_F=NULL
  LD_CL_dIVW_nIV=NULL
  LD_CL_dIVW_sd=NULL
  LD_CL_beta_dIVW_3=NULL
  LD_CL_dIVW_F_3=NULL
  LD_CL_dIVW_nIV_3=NULL
  LD_CL_dIVW_sd_3=NULL
  for ( m in(1:(iter-1)) )
  {
    LD_CL_beta_dIVW[m]=LD_CL[[m]]$DIVW$dIVW$res$beta
    LD_CL_dIVW_F[m]=LD_CL[[m]]$DIVW$dIVW$res$f
    LD_CL_dIVW_nIV[m]=LD_CL[[m]]$DIVW$dIVW$res$numIV_sel
    LD_CL_dIVW_sd[m]=LD_CL[[m]]$DIVW$dIVW$res$sd
    LD_CL_beta_dIVW_3[m]=LD_CL[[m]]$threeDIVW$threedIVW$res$beta
    LD_CL_dIVW_F_3[m]=LD_CL[[m]]$threeDIVW$threedIVW$res$f
    LD_CL_dIVW_nIV_3[m]=LD_CL[[m]]$threeDIVW$threedIVW$res$numIV_sel
    LD_CL_dIVW_sd_3[m]=LD_CL[[m]]$threeDIVW$threedIVW$res$sd
  }
  Beta_IVW=mean(na.omit(LD_CL_beta_IVW))
  Beta_IVW_3=mean(na.omit(LD_CL_beta_IVW_3))
  Beta_IVW_31=mean(na.omit(LD_CL_beta_IVW_31))
  sd_IVW=sd(na.omit(LD_CL_beta_IVW))
  sd_IVW_3=sd(na.omit(LD_CL_beta_IVW_3))
  sd_IVW_31=sd(na.omit(LD_CL_beta_IVW_31))
  se_IVW=mean(na.omit(LD_CL_IVW_sd))
  se_IVW_3=mean(na.omit(LD_CL_IVW_sd_3))
  se_IVW_31=mean(na.omit(LD_CL_IVW_sd_31))
  nIV_IVW=mean(na.omit(LD_CL_IVW_nIV))
  F_IVW=mean(na.omit(LD_CL_IVW_F))
  nIV_IVW_3=mean(na.omit(LD_CL_IVW_nIV_3))
  F_IVW_3=mean(na.omit(LD_CL_IVW_F_3))
  nIV_IVW_31=mean(na.omit(LD_CL_IVW_nIV_31))
  F_IVW_31=mean(na.omit(LD_CL_IVW_F_31))
  Beta_dIVW=mean(na.omit(LD_CL_beta_dIVW))
  sd_dIVW=sd(na.omit(LD_CL_beta_dIVW))
  nIV_dIVW=mean(na.omit(LD_CL_dIVW_nIV))
  F_dIVW=mean(na.omit(LD_CL_dIVW_F))
  se_dIVW=mean(na.omit(LD_CL_dIVW_sd))
  Beta_dIVW_3=mean(na.omit(LD_CL_beta_dIVW_3))
  sd_dIVW_3=sd(na.omit(LD_CL_beta_dIVW_3))
  nIV_dIVW_3=mean(na.omit(LD_CL_dIVW_nIV_3))
  F_dIVW_3=mean(na.omit(LD_CL_dIVW_F_3))
  se_dIVW_3=mean(na.omit(LD_CL_dIVW_sd_3))
  
  Beta_RIVW=mean(na.omit(LD_CL_beta_RIVW))
  sd_RIVW=sd(na.omit(LD_CL_beta_RIVW))
  nIV_RIVW=mean(na.omit(LD_CL_RIVW_nIV))
  F_RIVW=mean(na.omit(LD_CL_RIVW_F))
  se_RIVW=mean(na.omit(LD_CL_RIVW_sd))
 
  result_IVW=c(Beta_IVW,sd_IVW,se_IVW,nIV_IVW,F_IVW)
  result_IVW_3=c(Beta_IVW_3,sd_IVW_3,se_IVW_3,nIV_IVW_3,F_IVW_3)
  result_IVW_31=c(Beta_IVW_31,sd_IVW_31,se_IVW_31,nIV_IVW_31,F_IVW_31)
  result_dIVW=c(Beta_dIVW,sd_dIVW,se_dIVW,nIV_dIVW,F_dIVW)
  result_RIVW=c(Beta_RIVW,sd_RIVW,se_RIVW,nIV_RIVW,F_RIVW)
  result_dIVW_3=c(Beta_dIVW_3,sd_dIVW_3,se_dIVW_3,nIV_dIVW_3,F_dIVW_3)
 
  df=rbind(result_IVW,result_IVW_3,result_IVW_31,result_dIVW,result_dIVW_3,result_RIVW)
  df=data.frame(df)
  names(df)=c("BETA","MCSD","SE","NIV","F")
  return(df)
}
#scale=50000
#scale_cluster=as.integer(scale/20)
#reference_panel_cluster=simBeta_diagonal(G_num=20,seed=sample(1:100000,1),num_row=scale,num_col=scale,0)
# read it through gwasvcf
correlation_matrix=matrix(rep(0,3693*3693),3693,3693)
for (j in (1:3693))
{for(k in (1:3693))
{
  correlation_matrix[j,k]=Sigma[j,k]/expsoure.org["SE"]$SE[j]/expsoure.org["SE"]$SE[k]
}
  
  print(j)}
#for (job in jobs){
beta=as.vector(expsoure.org["BETA"])
se=as.vector(expsoure.org["SE"])
indexs=seq(1,3693,200)
indexs[20]=3694
pi=0
#thre=setting(job)
iter=1

while(iter<repeat_times){
  gamma_X_hat=NULL
  gamma_Y_hat=NULL
  gamma_X_sd=NULL
  gamma_Y_sd=NULL
  gamma_Z_hat=NULL
  gamma_Z_sd=NULL
 
  for (i in (1:(length(indexs)-1))){
    left=indexs[i]
    right=indexs[i+1]-1
    covariance_matrix=Sigma[left:right,left:right]
    gamma_x=beta$BETA[left:right]
    #  index=sample((1:length(gamma_x)),length(gamma_x)*pi,replace=FALSE)
    # gamma_x[index]=0
    sd_x=se$SE[left:right]
    gamma_x_hat=mvrnorm(1,gamma_x,covariance_matrix)
    gamma_y_hat=mvrnorm(1,gamma_x,covariance_matrix)
    gamma_z_hat=mvrnorm(1,gamma_x,covariance_matrix)
  
    gamma_X_hat[left:right]=gamma_x_hat
    gamma_Y_hat[left:right]=gamma_y_hat
    gamma_Z_hat[left:right]=gamma_z_hat
 
    gamma_Y_sd[left:right]=sd_x
    gamma_X_sd[left:right]=sd_x
    gamma_Z_sd[left:right]=sd_x

  } 
  
  
  # Conduct analysis
  #RIVW
 
  MR.revclump = Estimation(beta$BETA,gamma_X_hat,gamma_Y_hat,gamma_X_sd,gamma_Y_sd, inde.method="Revised_Pruning",reference_panel=correlation_matrix,scale_cluster=indexs,etamean = 0.5, 
                            sel.pthr = 5e-5,rerand = TRUE,gamma_summary_Z=gamma_Z_hat,gamma_Z_sd=gamma_Z_sd )
  #IVW

 
  MR.clump = Estimation(beta$BETA,gamma_X_hat,gamma_Y_hat,gamma_X_sd,gamma_Y_sd, inde.method="clumping",reference_panel=correlation_matrix,scale_cluster=indexs,etamean = 0.5, 
                         sel.pthr =5e-5,rerand = TRUE,gamma_summary_Z=gamma_Z_hat,gamma_Z_sd=gamma_Z_sd )
 
  #pruning
  #RIVW
  MR.prune =  Estimation(beta$BETA,gamma_X_hat,gamma_Y_hat,gamma_X_sd,gamma_Y_sd, inde.method="Pruning",reference_panel=correlation_matrix,scale_cluster=indexs,etamean = 0.5, 
                         sel.pthr =5e-5,rerand = TRUE,gamma_summary_Z=gamma_Z_hat,gamma_Z_sd=gamma_Z_sd )
                         #revising
 
  #RIVW
  MR.revclump_re = Estimation(beta$BETA,gamma_Y_hat,gamma_X_hat,gamma_Y_sd,gamma_X_sd, inde.method="Revised_Pruning",reference_panel=correlation_matrix,scale_cluster=indexs,etamean = 0.5, 
                               sel.pthr = 5e-5,rerand = TRUE,gamma_summary_Z=gamma_Z_hat,gamma_Z_sd=gamma_Z_sd )
 
  #RIVW
  MR.clump_re = Estimation(beta$BETA,gamma_Y_hat,gamma_X_hat,gamma_Y_sd,gamma_X_sd, inde.method="clumping",reference_panel=correlation_matrix,scale_cluster=indexs,etamean = 0.5, 
                            sel.pthr = 5e-5,rerand = TRUE,gamma_summary_Z=gamma_Z_hat,gamma_Z_sd=gamma_Z_sd )
 
  MR.prune_re =  Estimation(beta$BETA,gamma_Y_hat,gamma_X_hat,gamma_Y_sd,gamma_X_sd, inde.method="Pruning",reference_panel=correlation_matrix,scale_cluster=indexs,etamean = 0.5, 
                             sel.pthr = 5e-5,rerand = TRUE,gamma_summary_Z=gamma_Z_hat,gamma_Z_sd=gamma_Z_sd )
  
  
  
  LD_CL[[iter]]=list(RIVW=MR.clump[["res"]]["RIVW"],IVW=MR.clump[["res"]]["IVW"],DIVW=MR.clump[["res"]]["dIVW"],
                     threeIVW_part=MR.clump[["res"]]["threeIVW_part"],threeIVW_all=MR.clump[["res"]]["threeIVW_all"],threeDIVW=MR.clump[["res"]]["threedIVW"])
  
  LD_PR[[iter]]=list(RIVW=MR.prune[["res"]]["RIVW"],IVW=MR.prune[["res"]]["IVW"],DIVW=MR.prune[["res"]]["dIVW"],
                     threeIVW_part=MR.prune[["res"]]["threeIVW_part"],threeIVW_all=MR.prune[["res"]]["threeIVW_all"],threeDIVW=MR.prune[["res"]]["threedIVW"])
  
  LD_R_PR[[iter]]=list(RIVW=MR.revclump[["res"]]["RIVW"],IVW=MR.revclump[["res"]]["IVW"],DIVW=MR.revclump[["res"]]["dIVW"],
                       threeIVW_part=MR.revclump[["res"]]["threeIVW_part"],threeIVW_all=MR.revclump[["res"]]["threeIVW_all"],threeDIVW=MR.revclump[["res"]]["threedIVW"])
  
  LD_CL_RE[[iter]]=list(RIVW=MR.clump_re[["res"]]["RIVW"],IVW=MR.clump_re[["res"]]["IVW"],DIVW=MR.clump_re[["res"]]["dIVW"],
                        threeIVW_part=MR.clump_re[["res"]]["threeIVW_part"],threeIVW_all=MR.clump_re[["res"]]["threeIVW_all"],threeDIVW=MR.clump_re[["res"]]["threedIVW"])
  
  LD_PR_RE[[iter]]=list(RIVW=MR.prune_re[["res"]]["RIVW"],IVW=MR.prune_re[["res"]]["IVW"],DIVW=MR.prune_re[["res"]]["dIVW"],
                        threeIVW_part=MR.prune_re[["res"]]["threeIVW_part"],threeIVW_all=MR.prune_re[["res"]]["threeIVW_all"],threeDIVW=MR.prune_re[["res"]]["threedIVW"])
  
  LD_R_PR_RE[[iter]]=list(RIVW=MR.revclump_re[["res"]]["RIVW"],IVW=MR.revclump_re[["res"]]["IVW"],DIVW=MR.revclump_re[["res"]]["dIVW"],
                          threeIVW_part=MR.revclump_re[["res"]]["threeIVW_part"],threeIVW_all=MR.revclump_re[["res"]]["threeIVW_all"],threeDIVW=MR.revclump_re[["res"]]["threedIVW"])
  iter=iter+1
  print(iter)
}
file_path=paste0("C:\\Users\\25110\\Desktop\\threrealdata2",".RData")
save.image(file_path)
file_path2=paste0("C:\\Users\\25110\\Desktop\\threrealdata2",".csv")
clumping=res_processing(LD_CL,iter)
pruning=res_processing(LD_PR,iter)
sigma_pruning=res_processing(LD_R_PR,iter)
clumping_re=res_processing(LD_CL_RE,iter)
pruning_re=res_processing(LD_PR_RE,iter)
sigma_pruning_re=res_processing(LD_R_PR_RE,iter)
result=rbind(clumping,clumping_re,pruning,pruning_re,sigma_pruning,sigma_pruning_re)
write.csv(result,file_path2)
#}
