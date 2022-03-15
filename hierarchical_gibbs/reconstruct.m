function electrodes = reconstruct(factors,coeff,latent,means)
%RECONSTRUCT Converts principal components back to a (milli?)volts-based
%EEG topology.
% @ factors: bump magnitudes [nPC x nBumps]
% @ coeff: PCA coeficients [ch x ch] loadings
% @ latent: PCA eiginevalues of the covariance matrix of data [ch x 1]
% @ means: Grand mean [1 x ch]

transform=coeff';

factors=factors';
nscans = size(factors,1);
width = size(factors,2);

electrodes = (factors.*repmat(sqrt(latent(1:width))',nscans,1)) * transform(1:width,:);
electrodes = electrodes + repmat(means,nscans,1); %add means for each electrode
