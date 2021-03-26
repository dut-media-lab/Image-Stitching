function T = y(Pts,isotropic)

% VGG_CONDITIONER_FROM_PTS - Returns a conditioning matrix for points
%
%	T = vgg_conditioner_from_pts(Pts [,isotropic])
%
% Returns a DxD matrix that normalizes Pts to have mean 0 and stddev sqrt(2)
%
%
%IN:
%	Pts - DxK list of K projective points. Last row is ignored.
%       isotropic - optional; if present then T(1,1)==T(2,2)==...==T(D-1,D-1).
%
%
%OUT:
%	T - DxD conditioning matrix

% Yoni, Thu Feb 14 12:24:48 2002

Dim=size(Pts,1);

Pts=vgg_get_nonhomg(Pts);
Pts=Pts(1:Dim-1,:);

m=mean(Pts,2);
s=std(Pts');
s=s+(s==0);

if nargin==1
  T=[ diag(sqrt(2)./s) -diag(sqrt(2)./s)*m];
else % isotropic; added by TW
  T=[ diag(sqrt(2)./(ones(1,Dim-1)*mean(s))) -diag(sqrt(2)./s)*m];
end
T(Dim,:)=0;
T(Dim,Dim)=1;
