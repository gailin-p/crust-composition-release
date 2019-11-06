# Return bin centers, averages, errors
# input: x (bin by), y (bin), errors
# TODO gailin - xerrors?
function bin(x,y,min,max,oversamplingratio,nbins)
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
        errors[i]=std(y[test]) * sqrt(oversamplingratio/sum(t));
        xerrors[i]=bincenters[i] - nanmean(x[test]);
    end
    return bincenters, averages, errors
end

export bin
