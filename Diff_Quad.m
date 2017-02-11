function [D]=Diff_Quad(N)
%
%  Differential quadrature matrix for N
%  points based on Lobatto grid points.
%  D matrix converts the diff. equation into eigenvalue problem.
%  y'=A*y  --> Dy=Ay
%
% 
%  For any questions and critiques 
%  bicakme@gmail.com
%  Altug Bicak @ 2005

for ii=1:N 
X(ii)=0.5*(1-cos((ii-1)*pi/(N-1)));
end
SagTaraf=1;
Taraf=0;
for ii=1:N
    for jj=1:N
        if ii~=jj
        a1=(1/(X(jj)-X(ii)));
        for k=1:N
            if k~=ii & k~=jj
              SagTaraf=SagTaraf*(X(ii)-X(k))/(X(jj)-X(k));    
           end
        end         
        a(ii,jj)=a1*SagTaraf;
        SagTaraf=1;
        end
        
        if ii==jj
            for k=1:N
                if k~=ii
                   Taraf=Taraf+(1/(X(ii)-X(k))); 
               end
            end
            
            a(ii,jj)=Taraf;
            Taraf=0;
        end
    end
end
D=a;
tic

