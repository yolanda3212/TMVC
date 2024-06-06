function [res]= myNMIACC(U,Y,numclass)
stream = RandStream.getGlobalStream;
reset(stream);
U_normalized = U;
MAXiter = 100; 
REPlic = 20; 
for idx = 1:20
    indx = kmeans(U_normalized,numclass,'maxiter',MAXiter,'replicates',REPlic,'emptyaction','singleton');
    Res(idx,:) = ClusteringMeasure(Y, indx);
end
res = [mean(Res)];

