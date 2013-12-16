% =======================================================================
%                          生成VI文件
% =======================================================================

% iCase=409;
% nFile=35;
% nElecNode=66;
% nElecBranch=87;
% nDFIGNode=14;
% nDFIGBranch=25;
% nDFIG=5;

% iCase=413;
% nFile=10;
% nElecNode=26;
% nElecBranch=35;
% nDFIGNode=14;
% nDFIGBranch=25;
% nDFIG=1;

iCase=417;
nFile=21;
nElecNode=3;
nElecBranch=3;
nDFIGNode=14;
nDFIGBranch=25;
nDFIG=5;

% iCase=4171;
% nFile=5;
% nElecNode=3;
% nElecBranch=3;
% nDFIGNode=14;
% nDFIGBranch=25;
% nDFIG=1;
%%

result_PSC = [];
for k = 1:nFile
    if k<9.5
        stg = sprintf('./result/pscadData/case%d_0%d.out',iCase,k);
    else
        stg = sprintf('./result/pscadData/case%d_%d.out',iCase,k);
    end
    result_temp = load(stg);
    t_PSC=result_temp(:,1);
    result_PSC = [result_PSC result_temp(:,2:end)];
end
nT=length(t_PSC);
resultV=zeros(nT,nElecNode+nDFIG*nDFIGNode);
resultI=zeros(nT,nElecBranch+nDFIG*nDFIGBranch);
resultV(:,1:nElecNode)=result_PSC(:,1:nElecNode);
resultI(:,1:nElecBranch)=result_PSC(:,nElecNode+1:nElecNode+nElecBranch);
for i=1:nDFIG
    resultV(:,nElecNode+1+(i-1)*nDFIGNode:nElecNode+i*nDFIGNode)=...
        result_PSC(:,nElecNode+nElecBranch+1+(i-1)*(nDFIGNode+nDFIGBranch):nElecNode+nElecBranch+i*nDFIGNode+(i-1)*nDFIGBranch);
    resultI(:,nElecBranch+1+(i-1)*nDFIGBranch:nElecBranch+i*nDFIGBranch)=...
        result_PSC(:,nElecNode+nElecBranch+1+i*nDFIGNode+(i-1)*nDFIGBranch:nElecNode+nElecBranch+i*(nDFIGNode+nDFIGBranch));
end
resultV=[t_PSC resultV];
resultI=[t_PSC resultI];
stg = sprintf('./result/pscadData/case%d_V.dat',iCase);
save(stg,'-ascii','resultV');
stg = sprintf('./result/pscadData/case%d_I.dat',iCase);
save(stg,'-ascii','resultI');