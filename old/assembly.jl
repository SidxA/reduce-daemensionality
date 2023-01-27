function embed_lag(x,win_l)
    T = length(x)
    P = win_l
    N = floor(Int,(2(T+1)/P))-1
    data = zeros(N-1,P)
    ind = 1
    for t in 1:N
        ind_l = ind + P - 1
        if ind_l < T
            data[t,:] = x[ind:ind_l]
            ind += Int(P/2)
        else
            break
        end
    end
    return data
end

using CSV
using DataFrames

function fluxnet_preprocess(name)
    df = DataFrame(CSV.File(name,delim=",",header=2))
    #remove timestamp
    df = df[:,2:end]
    #remove constant rows
    il=[]
    for i in 1:ncol(df)
        if df[2,i]==df[3,i]
            il = push!(il,false)
        else
            il = push!(il,true)
        end
    end
    return Array(df)[:,Array(Bool.(il))]
end


#the old days of mistrust but its fine!

function toeplitz(X,W)
    C = zeros(W,W)
    for i=1:W,j=1:W
        lag = abs(i-j)
        C[j,i] = 1/(N-lag) * sum([X[t]*X[t+lag] for t=1:N-lag])
    end
return C
end


#normalize RC components AMPLITUDe STILL NEEDS IGENVALUES?
normalize_by_components(array) = hcat([array[:,i]./std(sum(array,dims=2)) for i in 1:size(array)[2]]...)