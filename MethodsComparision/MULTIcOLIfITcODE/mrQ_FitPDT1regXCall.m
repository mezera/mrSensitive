function [X_valdationErrF, best1, best2,X_valdationErrSN PDfit gEst resnorm exitflag G] =   mrQ_FitPDT1regXCall ...
            (opt,M0_v,BM1,pBasis,Clist,R1basis,g0,Segmask )

% copy write Aviv mezer        
% writen by AM Stanford 2014
% use  X validation to find the best Lamda for T1 regularization
[X_valdationErrF,  X_gEstF]=pdX_valdationLoop_2(opt.lambda,3,M0_v(BM1,Clist), pBasis(BM1,:),R1basis(BM1,:),g0,Segmask(BM1));
        
      
        best1 = find(X_valdationErrF(1,:)==min(X_valdationErrF(1,:))) % sum of abs err
        best2 = find(X_valdationErrF(2,:)==min(X_valdationErrF(2,:)))% RMSE
        %figure;  semilogy(opt.lambda,X_valdationErrF(2,:),'*-'); X_valdationErrF(2,:)./min(X_valdationErrF(2,:))
        %figure;  semilogy(opt.lambda,X_valdationErrF(1,:),'*-'); X_valdationErrF(1,:)./min(X_valdationErrF(1,:))
        % i cases we get no meaningfull Xvalidation with the
        % regularization we should be worry!!!
        % if there is not enght data we many X_valdation Err are
        % almost similar we can cause to take the begest one. this
        % need to be implamented
       
        
         % i cases we get no meaningfull Xvalidation with the
        % regularization we should be worry!!!
        % if there is not enght data we many X_valdation Err are
        % almost similar we can cause to take the begest one. this
        % need to be implamented
        if min(X_valdationErrF(2,:))*2 > X_valdationErrF(2,end)
            X_valdationErrSN=1;
            disp(['no clear  X-validation improvment with regularization'])
        else
            X_valdationErrSN=0;
        end
        
        
        
        
        %  [PDfit,RMSE1]=pdCoilSearch_T1reg( lambda1(best1),M0_v(:,Clist),pBasis,R1basis,X_gEstF(:,:,1,best),[],[],PDsim);
        [PDfit,~,G,gEst, resnorm,exitflag ]=pdCoilSearch_T1reg( opt.lambda(best2),M0_v(BM1,Clist),pBasis(BM1,:),R1basis(BM1,:),X_gEstF(:,:,1,best2),Segmask(BM1));
        