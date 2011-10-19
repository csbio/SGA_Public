function res = adjust_model_fit(vect, vect_std),

%assume vect(1) has the overall mean, and mean + sum(pi) = c for all
%parameters, pi

%adjust mean such that the total number of significantly non-zero
%parameters is minimized

% res = fminsearch(@(x) num_sig_params(x,vect,vect_std),0,1);
 %options = optimset('Display','iter','TolFun',1e-8);
 options = optimset('TolFun',1e-8);
 %init = 5;
 %if(rand > .5) init =-5; end    
 %res = fminsearch(@(x) num_sig_params(x,vect,vect_std),5,options);
 res = fminbnd(@(x) num_sig_params(x,vect,vect_std),-50,50);



function tot = num_sig_params(adj,vect,vect_std),
    vect(1) = vect(1)+adj;
    vect(2:end)=vect(2:end)-adj;
    
    tot1 = sum(vect > 0 & vect - vect_std > 0);
    tot2 = sum(vect < 0 & vect + vect_std < 0);
    
    tot=tot1+tot2;