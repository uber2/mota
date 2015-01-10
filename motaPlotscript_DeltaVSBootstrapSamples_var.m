% DESCRIPTION: This is a pure helper function with the sole purpose
% of reproducing Figure 2 of the paper "S.Hengl, C. Kreutz,
% J.Timmer, T. Maiwald, Data-Based Identifiability Analysis of
% Nonlinear Dynamical Models, Bioinfomatics, 2007". 
%
function MeanDelta=motaPlotscript_DeltaVSBootstrapSamples_var(X)

NOI = 10; % Number of iterations
NoB = [1 2 5 10 20 30 40 50 60]; % number of bootstrap samples
k   = 0; % loop variable
dim=length(X(1,:))-1;
Delta=zeros(NOI,dim+1,length(NoB));

% main loop
for j=NoB
  disp(j)
  k=k+1;
  
  for i=1:NOI     
    % main worker function
    [xsorts,PhiBoot]=motaMexACEofBootstrap(X',j,100);
    
    % testfunction
    Delta(i,:,k)=var(PhiBoot,[],2);
  end  

end

% take the average over all iterations
MeanDelta=zeros(dim+1,length(NoB));
for i=1:length(NoB)
    for j=1:dim+1
        MeanDelta(j,i)=mean(Delta(:,j,i));
    end
end


% plot the results   
figure(2)
set(gca,'FontSize',14)
plot(NoB,MeanDelta','.-','LineWidth',1,'MarkerSize',14);
xlabel('number of bootstrap samples')
ylabel('mean(H_K)')


end
