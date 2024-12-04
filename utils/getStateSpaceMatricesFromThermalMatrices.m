function [Amat, Bmat] = getStateSpaceMatricesFromThermalMatrices(CapMat, ResMat)

% Copyright 2022 The MathWorks, Inc.

%%  Create matrixes A B C and D for state space model
[N,~]=size(ResMat);
R_inv = zeros(N);
C_inv = zeros(N,1);
if numel(CapMat)>N % CapMat is a NxN matrix
    CapVec = diag(CapMat);
    CapVec = CapVec(:); % column vec
else
    CapVec = CapMat(:);
end

% Setup model parameters
% Build inverse R,C matrixes
for row=1:1:N
    for col=1:1:N
        R_inv(row,col)=1/ResMat(row,col);
        if row==col
            R_inv(row,col)=0;
        end
    end
end

for i=1:1:N
    C_inv(i,1)=1/CapVec(i,1);
end

% Build Amat matrix

Amat=zeros(N);

% Top right half of Amat matrix
for row=1:1:N
    for col=row+1:1:N
      if C_inv(row,1) < 10^6           
        Amat(row,col)=C_inv(row,1)*R_inv(row,col);
      else
        Amat(row,col)=(10^6)*R_inv(row,col);
      end
    end
end

% Bottom left half of Amat matrix
for row=2:1:N
    for col=1:row-1    
        if C_inv(row,1)< 10^6
            Amat(row,col)=C_inv(row,1)*R_inv(row,col);
        else
            Amat(row,col)=(10^6)*R_inv(row,col);
        end
    end
end

% Diagonal of Amat Matrix
for row=1:1:N 
  if C_inv(row,1)< 10^6
    Amat(row,row)=-(C_inv(row,1)*(sum(R_inv(row,:))));
  else
    Amat(row,row)=(-(10^6)*(sum(R_inv(row,:))));
  end
end

% Build Bmat matrix
Bmat=zeros(N);

for row=1:1:N  
    if C_inv(row,1)< 10^6
        Bmat(row,row)=C_inv(row,1);   

    else 
        Bmat(row,row)= 10^6;
    end   
end

end

