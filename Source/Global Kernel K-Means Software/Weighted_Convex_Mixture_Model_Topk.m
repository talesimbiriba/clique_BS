function [Top_k_idx,Likelihood]=Weighted_Convex_Mixture_Model_Topk(S,Data_distribution,Convits,Top_k)
%[Top_k_idx,Likelihood]=Weighted_Convex_Mixture_Model_Topk(S,Data_distribution,Convits,Top_k)
%
%This function implements the clustering algorithm proposed by
%D. Lashkari and P. Golland, "Convex Clustering with Exemplar-Based Methods", NIPS 2008
%but appropriately modified to handle weighted data points in feature space as proposed by
%G. Tzortzis and A. Likas, "The Global Kernel k-Means Algorithm for Clustering in Feature Space", to be published in IEEE TNN.
%
%Function Inputs
%===============
%
%S is the matrix containing the pairwise similarities where S=exp(-b*dist(f(x(i)),f(x(j))).
%
%Data_distribution is the empirical distribution of the dataset based on weights.
%
%Convits is the number of consecutive iterations that the Top_k exemplars in the CMM method
%must remain the same in order the CMM method to converge (usually a value of 1000 suffices).
%
%Top_k is the number of returned exemplars. 
%
%Function Outputs
%================
%
%Top_k_idx is the index of the k points with the largest prior probabilities that are selected as exemplars.
%
%Likelihood is the dataset likelihood under the CMM model.
%
%
%Courtesy of G. Tzortzis

%Dataset size.
Data_num=size(S,1);

%Q_thres is a threshold on the components prior probabilities. Each
%such prior that is below this threshold is considered to be zero.
Q_thres=(10^-3)/Data_num; 

%Vector Q contains the prior probabilities for each component and are
%initialized to the value 1/Data_num.
Q=(1/Data_num)*ones(Data_num,1); 

%Store the Top_k exemplars of the last Convits iterations.
Converge=zeros(Top_k,Convits);
Check_conv=zeros(Top_k,1);
Iter=1;

while 1
    
    %Calculate the z values.
    Z=S*Q;
    
    %Calculate the n values.
    Z_tmp=1./Z;
    Z_tmp=Z_tmp';
    N=(Data_distribution.*Z_tmp)*S;
    N=N';
        
    [Top_Q,Top_k_idx]=sort(Q,'descend');
    Top_k_idx=Top_k_idx(1:Top_k);
    
    %Converge when the exemplars stay fixed for Convits iterations.
    Converge(:,mod(Iter-1,Convits)+1)=Top_k_idx;
    if Iter>=Convits
        
        for i=1:Top_k
            Check_conv(i)=sum(Converge(i,:)==Converge(i,1));
        end
        
        Check_conv_sum=sum(Check_conv==Convits);
        
        if Check_conv_sum==Top_k
            break;
        end
    end        
    
    %Calculate the new prior probabilities Q.
    Q=N.*Q;
    
    %Eliminate priors smaller than Q_thres and renormalize.
    Q(Q<Q_thres)=0;
    Q=Q/sum(Q);
        
    Iter=Iter+1;    
end

%Final exemplars.
[Top_Q,Top_k_idx]=sort(Q,'descend');
Top_k_idx=Top_k_idx(1:Top_k);

%Calculate the likelihood.
Likelihood=log(S*Q);
Likelihood=Data_distribution*Likelihood;
   
    
    
    