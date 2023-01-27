using Distributed
@everywhere begin
    include("/net/home/lschulz/scripts/toolbox.jl")

    #return name,B,W,k,Ws,ks
    function recover_para(title)
        title*="_"
        name = title[1:6]
        title = title[8:end]
        ind = findfirst(isequal('_'), title)
        B = parse(Int64,title[2:ind-1])
        title = title[ind+1:end]
        ind = findfirst(isequal('_'), title)
        W = parse(Int64,title[2:ind-1])
        title = title[ind+1:end]
        ind = findfirst(isequal('_'), title)
        k = parse(Int64,title[2:ind-1])
        title = title[ind+1:end]
        ind = findfirst(isequal('_'), title)
        Ws = parse(Int64,title[3:ind-1])
        title = title[ind+1:end]
        ind = findfirst(isequal('_'), title)
        ks = parse(Int64,title[3:ind-1])
        title = title[ind+1:end]
        return name,B,W,k,Ws,ks
    end

    #give an array of fixed size with the first ones beeing the sorted maxmin entries of input array
    function sortarray!(ar::Vector,S)
        L = length(ar)
        if L>=S
            return sort!(ar,rev=true)[1:S]
        else
            return append!(sort!(ar,rev=true),zeros(S-L))
        end
    end
    #read in one in the list, return string with all the information in a fixed k manner
    function readout(name)

        epsi_T = 20
        epsi_A = 0.1
        max_k = 64

        savedir = dir*"sample_ssa2/"
        EOF = npzread(savedir*name)
        signal = npzread(savedir*name*"_signal")

        N = size(signal)[1]
        name,B,W,k,Ws,ks = recover_para(name)
        pairs,periods_F = pairing_full(EOF,epsi_T,epsi_A,false)
        periods_HHT = Float64[0.5 .* N ./ iterated_hilbert(EOF[:,i],64)[1] for i=1:k]

        T = zeros(2*max_k)
        T[1:max_k] = round.(sortarray!(periods_F,max_k),digits=3)
        T[1+max_k:end] = round.(sortarray!(periods_HHT,max_k),digits=3)

        wlist = ""
        for word in [name,N,B,W,k,Ws,ks,T...]
            wlist *=string(word)
            wlist *=",\t"
        end
        wlist *= "\n"
        return wlist
    end
end

#loop that reads the directory and filters out the ones with signal to only give the name
savedir= readdir("/net/home/lschulz/logs/KW_27/sample_ssa2")
for name in savedir[contains.(savedir,"signal")]
    logfile = open(dir*"sample_ssa2_log.txt", "a") 
    write(logfile, readout(name[1:end-7]))
    close(logfile)
end