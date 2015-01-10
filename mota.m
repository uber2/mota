% Usage:
% out = mota(X,IO,threshold1,threshold2,threshold3,sampleSize,numOfBootSamp,MAXNOP,par)
%
% Title:       MOTA (Mean Optimal Tranformations Approach)
% Description: MOTA detects functional relations between parameters. 
% Arguments:   X : (n X p) matrix with n estimates for p parameters
%
%              IO: Output-option
%                  (0) no output (default)
%                  (1) text output in command line. Loop procedure is displayed
%                  (2) calculation of multiple rSquared is displteayed as well
%
%              Threshold 1: (default=0.01) If the H-value of a paramter is smaller than threshold1
%               it is supposed not to contribute at all to a given set of Parameters. The
%               higher Threshold1, the more conservative gets mota.
% 
%              Threshold 2: (default=0.07) If the H-value of a paramter is smaller than threshold
%               2 but greater than threshold 1, it is supposed to contribute
%               to a given set of Parameters. More parameter have to be added until
%               Threshold2 is exceeded
%
%              Threshold 3: (default=0.08) If a new added parameter yields a delta value greater than
%               threshold 3, the testing sequence continues, even if the new MeanOfDelta
%               values is smaller than the previous one. Threshold3 is extremeley close
%               to the analytically determined maximum value of the testfunction.
%               Threshold 3 is not necessary for the functionallity of the algorithm
%
%              sampleSize: (default=N/2) Size of boostrap-samples     
%              
%              numOfBootSamp: (default=35) Number of bootstrap-samples
%
%              MAXNOP: (default=p-1) Maximum number of paramters to
%               be added. This parameter can improve the performance
%               drastically if one is anyway only interested in
%               functional relation of a few (MAXNOP) parameters.
%
%              par: (default=1:p) List of parameters that should be
%               investigated as response variables.
%
% Value:
%              S  : (p X p) matrix with ones and zeros indicating functional relations  
%              r^2: percentage of explained variance 
% 
% Note:        The only input required is the (n X p) matrix X. All other input
%              variables are set to default values.
%
function out=mota(X,IO,threshold1,threshold2,threshold3,sampleSize,numOfBootSamp,MAXNOP,par)

%% defaults
     
  % Output
  if~exist('IO','var')||isempty(IO)
    IO=0;
  end
  
  % Thresholds
  if~exist('threshold1','var')||isempty(threshold1)
    threshold1=0.01;
  end
  
  if~exist('threshold2','var')||isempty(threshold2)
    threshold2=0.07;
  end
  
  if~exist('threshold3','var')||isempty(threshold3)
    threshold3=0.08;
  end
  
  % number of bootstrap samples    
  if~exist('numOfBootSamp','var')||isempty(numOfBootSamp)
    numOfBootSamp=35;
  end
  
  % bootstrap samples size
  if~exist('sampleSize','var')||isempty(sampleSize)
    sampleSize=floor(length(X(:,1))/2);
  end
  
  % MAXNOP
  if~exist('MAXNOP','var')||isempty(MAXNOP)
    MAXNOP=length(X(1,:))-1;
  end
  
  % par
  if~exist('par','var')||isempty(par)
    par=1:length(X(1,:));
  end
    
%% Warnings
  if threshold1>0.09
    warning('threshold1 should not exceed 1/12'); 
  end   
  
  if threshold1>threshold2 || threshold1>threshold3
    error('threshold1 must be smaller than threshold2 and threshold3'); 
  end
  
  if threshold2>0.09
    warning('threshold2 should not exceed 1/12'); 
  end 
  
  if threshold3>0.09
    warning('threshold3 should not exceed 1/12'); 
  end   
  
  if sampleSize>length(X(:,1))
    error('sampleSize must not exceed the number of rows of X') 
  end
  
    
%% Variables  

  % put paramters in output struct
  out.date = date;
  
  out.threshold1 = threshold1;
  out.threshold2 = threshold2;
  out.threshold3 = threshold3;
  
  out.sampleSize = sampleSize;
  out.numOfBootSamp = numOfBootSamp;
  out.MAXNOP = MAXNOP;
  
  out.par = par;


  % backup file
  filename = [date 'mota_backup.mat'];

  % Set history H to zero
  H=0;

  % Number Of parameters
  NOP=length(X(1,:));
  
  % (S)trong functional relation    
  out.S=zeros(length(par),NOP);    
    
  % Delta2 is used to save all current delta values for the new available
  % parameters. 
  Delta2=zeros(1,NOP);
  fprintf('            \n')
  out.rSquared=zeros(2,NOP,length(par));

  out.X=X;

    
%% Main 
for i=par
  
   % counter for number of extension steps     
   m=0;     

   % number of predictors
   dim=0;
   
   % Take i.th paramter
   Xcopy2=X(:,i);
   
   % Which paramters are available
   available=ones(1,NOP);
   available(i)=0;
   
   % Number of available parameters
   NOavailableP=sum(available);
   
   
   % Combine the i.th paramter with ALL other paramters. The
   % maxium number of paramters which can be combined is NOP-1.
   for j=1:MAXNOP
       
            dim=dim+1 ;   
            
            %%  try out all available parameters and determine the one with maximum Delta 
            for k=1:NOavailableP
                
                % add a new one
                    Help=toPn(available,'1');
                    Xcopy=[Xcopy2 X(:,Help(k))];
                    NewNOP=length(Xcopy(1,:));
                
                % test of functional relationship
                    
                    [xsorts,PhiBoot]=motaMexACEofBootstrap(Xcopy',numOfBootSamp,sampleSize);
                                       
                    Delta=var(PhiBoot,[],2);
                
                % The list of available paramaters is filled with the
                % corresponding delta-vallues for each parameter
                    Delta2(Help(k))=Delta(NewNOP);
              
            end

            % Set to zero those parameters which are not discussed at this
            % very moment
            L=toPn(available,'0');
            Delta2(L)=0;
            
            if IO==2
                Delta2
            end

            [MaxOfDelta2,IxOfMaxOfDelta2]=max(Delta2);
            if IO==2
                IxOfMaxOfDelta2
            end
            % Stop if you can not find an additional paramter to establish
            % a functional relation
            
            
            %% "break" if there is no functional relation        
            if MaxOfDelta2<threshold1 && H==0
                
                if IO==1||IO==2
                    disp(sprintf('No strong functional relations of parameter %d !!\n',i))
                end

                out.S(find(par==i),:)=zeros(1,NOP);
                % If there is no functional relation there is at least a 1
                % on the diagonal
                out.S(find(par==i),i)=1; 
                save(filename,'out')
                break
            end
            
            
            %% calculate multiple rSquared (if one parameter has already been added)
            if m>0
                out.rSquared(1,m,find(par==i)) = length(Xcopy2(1,:));
                [psi,phi]=ace(Xcopy2');
                out.rSquared(2,m,find(par==i))=motaR2(phi(1,:)',phi(2:end,:)');
                if IO==2
                    out.rSquared(:,:,find(par==i))
                end
            end
            
          
            %% you've found a functional relation. Update the list of available para.
            % You have found an additonal paramter
            % ->add it to "Xcopy2" ...
            Xcopy2=[Xcopy2 X(:,IxOfMaxOfDelta2)];
            % Counter for number of extension steps
            m=m+1;
                                    
            % ... and update the list of "available" Parameters
            available(IxOfMaxOfDelta2)=0;
            NOavailableP=sum(available);
            
            [xsorts,PhiBoot]=motaMexACEofBootstrap(Xcopy2',numOfBootSamp,sampleSize);
            
            % test function
            Delta=var(PhiBoot,[],2);
                    
            MeanOfDelta=mean(Delta);
            if IO==2
                MeanOfDelta
            end
            
            %% until a first hit, this section has to be run through each time
            if H==0
                if MeanOfDelta>threshold2
                    MeanOfDeltaHistory=MeanOfDelta;
                    availableHistory=available;
                    
                    if j==MAXNOP
                            if IO==1||IO==2
                                disp(sprintf('First hit for parameter %d',i))
                                toPn(availableHistory,'0')
                                disp(sprintf('There are no more parameters to add\n'))
                                disp(sprintf('Strong functional relations of %d: \n',i))
                                toPn(available,'0')
                            end

                            % Create proper entry in output matrix S
                            v=zeros(1,NOP);
                            v(toPn(available,'0'))=1;
                            out.S(find(par==i),:)=v;
                            save(filename,'out')
                            H=0;
                            break 
                    end
                                        
                    if IO==1||IO==2
                        disp(sprintf('First hit for parameter %d',i))
                        toPn(availableHistory,'0')
                    end

                    % set history H to H=1
                    H=1;
                else
                    if j==MAXNOP
                       if IO==1||IO==2
                            disp(sprintf('No strong functional relations of parameter %d !!\n',i))
                       end

                       out.S(find(par==i),:)=zeros(1,NOP);
                       % If there is no functional relation there is at least a 1
                       % on the diagonal
                       out.S(find(par==i),i)=1; 
                       save(filename,'out')
                       break
                    end
                end
                
                
            %%  check if a further improvement is possible   (H=1)          
            else
                    if (MeanOfDelta>MeanOfDeltaHistory||MeanOfDelta>threshold3)
                        
                        if j==MAXNOP
                            
                            if IO==1||IO==2
                                
                                disp(sprintf('Strong functional relations of %d: \n',i))
                                toPn(available,'0')
                            end

                            % Create proper entry in output matrix S
                            v=zeros(1,NOP);
                            v(toPn(available,'0'))=1;
                            out.S(find(par==i),:)=v;
                            save(filename,'out')
                            H=0;
                            break 
                        else
                                            
                            % Keep in mind what you have just done
                            MeanOfDeltaHistory=MeanOfDelta;
                            availableHistory=available;
                            
                            if IO==1||IO==2
                                disp('further improvement')
                                toPn(availableHistory,'0')
                                
                            end

                            H=1; % Set History H to 1
                        end
                    else
                        if IO==1||IO==2
                            disp(sprintf('Strong functional relations of %d: \n',i))
                            toPn(availableHistory,'0')
                        end

                        % Create proper entry in output matrix S
                        v=zeros(1,NOP);
                        v(toPn(availableHistory,'0'))=1;
                        out.S(find(par==i),:)=v;
                        save(filename,'out')
                        H=0;
                        break
                    end
            end
   end
   
   % show progress in command line (thanks to Thomas Maiwald)
   if IO==0
     fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b progress %2.0f%%',find(par==i)/length(par)*100)
   end
   
end
 
fprintf('\n') 
save(filename,'out')

end


%--------------------------------------------------------------------------
% Frequently Used. It converts a vector of ones and zeros to the
% corresponding number of the parameter, e.g. 
% [0 1 0 1 0 0 0 1] -> [2 4 8]
function Pn=toPn(IUL,s)
Pn=(regexp(num2str(IUL),s')+2)/3;
end