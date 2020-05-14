using StatGeochem

# Return bin centers, averages, errors
# input: x (bin by), y (bin), y errors, y std
# TODO gailin - xerrors?
function bin(x,y,min,max,oversamplingratio,nbins)
    # [bincenters,averages,errors,varargout]=bin(x,y,min,max,oversamplingratio,nbins,varargin)
    # Return the average values for independent variable y binned as a funciton
    # of independent variable x.

    averages=fill(NaN,nbins)
    errors=fill(NaN,nbins)
    xerrors=fill(NaN,nbins)
    ydevs = fill(NaN,nbins)

    binwidth=(max-min)/nbins
    binedges=range(min, stop=max, length=nbins+1)
    bincenters= range(min + binwidth/2, stop=max - binwidth/2, length= nbins)

    for i=1:nbins
        test = (x .> binedges[i]) .& (x .< binedges[i+1]) .& (.~ map(isnan, y))
        averages[i]=nanmean(y[test])
        # From Brenhin: 
        #multiplying the uncertainty by 2 for two-sigma, and then dividing by the
        #sqrt of the number of samples in each rho bin.
        #If this is resampled data, then you just have to also multiply
        #all the uncertainties by sqrt(N_original_samples) / sqrt(N_resampled_samples)
        #to avoid double-counting.
        errors[i]=nanstd(y[test]) * sqrt(oversamplingratio/sum(test));
        xerrors[i]=bincenters[i] - nanmean(x[test]);
    end
    return bincenters, averages, errors, ydevs
end

"""
   bin 
No oversamplingratio, so return standard error of mean  
"""
function bin(x,y,min,max,nbins)
    # [bincenters,averages,errors,varargout]=bin(x,y,min,max,oversamplingratio,nbins,varargin)
    # Return the average values for independent variable y binned as a funciton
    # of independent variable x.

    averages=fill(NaN,nbins)
    errors=fill(NaN,nbins)
    xerrors=fill(NaN,nbins)

    binwidth=(max-min)/nbins
    binedges=range(min, stop=max, length=nbins+1)
    bincenters= range(min + binwidth/2, stop=max - binwidth/2, length= nbins)

    for i=1:nbins
        test = (x .> binedges[i]) .& (x .< binedges[i+1]) .& (.~ map(isnan, y))
        averages[i]=nanmean(y[test])
        errors[i]=nanstd(y[test])/sqrt(length(y[test]))
        xerrors[i]=bincenters[i] - nanmean(x[test])
    end
    return bincenters, averages, errors
end

export bin
