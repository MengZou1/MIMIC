function [x,selected,accuracy]=MIMIC(S,genes_name,k,lambda,mu,theta,y_0,T)

%[x,selected,accuracy]=MIMIC(S,genes_name,k,lambda,mu,theta,y_0,T)
%
%Input matrix for MIMIC
% S: Expression matrix, i-th row denotes the cell type i and j-th
     %column denotes gene j. i=1,2,...m, j=1,2,...,n
% genes_name: Cell array, one column and j-th row denotes gene j, 
              % the same size with the column of S.
% k: The number of selected genes, can vary from 1 to the number of total
     % genes
%lambda: The balance parameter for differences between using selected genes
        % and entire genes
% mu: The balance parameter for specficity score 
% theta: A parameter of sigmoid function to decide the tendency of 
         %this function, positive value, default 10.
% y_0: The midpoint of sigmid function. positive value, default 3.
% T:  The threshold for the evaluation, positive value, default 1.
%
% Output for MIMIC
% x: The solution of mixed integer programming, the first n elements denote
     % the gene is selected or not, the following m*(m-1)/2 elements denote
     % the difference between using selected and entire samples, the last
     % m*(m-1)/2 elements denote the normailized difference between using
     %selected and entire genes.
% Selected: The cell type specific genes selected by MIMIC
% Accuracy: The model accuracy evaluated by threshlod T, here we set T=1 as
            %defalut
%
% Example for MIMIC
%for Surface markers' (SMs) expression data to identify 4 genes with
% with lambda=10 mu=1 using default parameter
%[x,selected,accuracy]=MIMIC(SMs_expression,SMs,4,10,1,[],[],[]);
%or
%[x,selected,accuracy]=MIMIC(SMs_expression,SMs,,4,10,3,1,10,1);


%% MIMIC in the procedure
%% Initial the parameter
if isempty(theta);
    theta=10;
end
if isempty(T);
    T=1;
end
if isempty(y_0);
    y_0=3;
end
% Calculate the specificity score for each gene in each cell type
JSD_S=JSD(S);
[m,n]=size(S);
%% initilize the inputs for MILP

f=[mu*min(JSD_S)';-ones(m*(m-1)/2,1);lambda*ones(m*(m-1)/2,1)];
Aeq=zeros(m*(m-1)/2+1,n+m*(m-1));
Aineq=zeros(m*(m-1)/2,n+m*(m-1));

%%Aeq for definition of d
q=1;
Z=zeros(m*(m-1)/2,n);
beq=zeros(m*(m-1)/2+1,1);
bineq=zeros(m*(m-1)/2,1);
for i=1:m-1
    for j=i+1:m
        z=distance(S(i,:),S(j,:),theta,y_0);
        Z(q,1:n)=z;
        Aeq(q,1:n)=z;
        Aeq(q,n+q)=-1;
        beq(q,1)=0;
        Aineq(q,n+q)=-1/k;
        Aineq(q,n+m*(m-1)/2+q)=-1;
        bineq(q,1)=-mean(z);
        q=q+1;
    end
end
%% fix the number of genes
Aeq(m*(m-1)/2+1,:)=[ones(1,n) zeros(1,m*(m-1))];

beq(m*(m-1)/2+1,1)=k;

%% mixed integer linear programming
    ctype=strcat(repmat('B',1,n),repmat('C',1,m*(m-1)));
    x=cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],zeros(n+m*(m-1),1),[],ctype);
    % Accuracy
    X=(x(1:n)>0.5);
    selected=genes_name(X)';
    accuracy=sum(x(n+1:n+m*(m-1)/2)>T)/(m*(m-1)/2);   