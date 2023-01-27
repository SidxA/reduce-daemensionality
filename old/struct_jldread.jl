include("/net/home/lschulz/scripts/toolbox.jl")


# this needs a redo include N!!!!



# read a file, calculate PC,RC and count !
function RC_and_protophase(filename,data,indir,outdir)
    filenamesplit = split(filename,"_")
    method = filenamesplit[2]
    i = parse(Int64,filenamesplit[3])
    name = filenamesplit[4]
    W = parse(Int64,filenamesplit[5])
    k = parse(Int64,filenamesplit[6])
    a = parse(Float64,filenamesplit[7])
    e = parse(Float64,filenamesplit[8])
    t = parse(Float64,filenamesplit[9])
    dictfile = load(indir*filename)

    N = 12052
    P = N - W +1

    signal = data[:,i]
    EOFs = dictfile["data"]
    PC = hcat([pc!(signal,EOFs[:,i],P,W) for i in 1:k]...)'
    RC = hcat([reconstructor(PC[i,:],EOFs[:,i],N,W) for i in 1:k]...)
    HilbertT = hcat([hilbert_transform(RC[:,i]) for i=1:k]...)
    protophases = hcat([atan.(imag(HilbertT[:,i]),real(HilbertT[:,i])) for i=1:k]...)

    jldsave(outdir*"RC_Diff_$(i)_$(name)_$(W)_$(k)_$(aet[1])_$(aet[2])_$(aet[3])_$(rand(1000:9999)).jld2",
    rc = RC, protophase = protophases)

end


#read all the files and put them into one file in the KW !

    function recover_keys_ssa(key)#return name
        ind = findfirst(isequal('_'), key)
        name=key[ind-6:ind-1]
        return name
    end

    function recover_keys_diff(key)#return name,a,e,t
        ind = findfirst(isequal('_'), key)
        name=key[ind-6:ind-1]
        key=key[ind+1:end]
        aet = parse(Tuple{Float64},key[find(isequal('_'), key)[1]+1,find(isequal('_'), key)[2]-1])
        a=aet[1]
        e=aet[2]
        t=aet[3]
        return name,a,e,t
    end
    
    function count_T(EOF)#return one-line string with the T,delim = ","

            #give an array of fixed size with the first ones beeing the sorted maxmin entries of input array
        function sortarray!(ar::Vector,S)
            L = length(ar)
            if L>=S
                return sort!(ar,rev=true)[1:S]
            else
                return append!(sort!(ar,rev=true),zeros(S-L))
            end
        end
        
        epsi_T = 20
        epsi_A = 0.1
        max_k = 32
        N = 12052

        pairs,periods_F = pairing_full(EOF,epsi_T,epsi_A,false)
        periods_HHT = Float64[0.5 .* N ./ iterated_hilbert(EOF[:,i],64)[1] for i=1:max_k]

        T = zeros(2*max_k)
        T[1:max_k] = round.(sortarray!(periods_F,max_k),digits=3)
        T[1+max_k:end] = round.(sortarray!(periods_HHT,max_k),digits=3)

        wlist = ""
        for word in T
            wlist *=string(word)
            wlist *=",\t"
        end
        return wlist
    end

    function readfiles_countT(datadir,log_loc)
        filenames = readdir(datadir)
        logfile = open(log_loc, "a") 
        for filename = datadir.*filenames
            filenamesplit = split(filename,"_")
            method = filenamesplit[2]
            name = filenamesplit[3]
            W = parse(Int64,filenamesplit[4])
            k = parse(Int64,filenamesplit[5])
            a = parse(Float64,filenamesplit[6])
            e = parse(Float64,filenamesplit[7])
            t = parse(Float64,filenamesplit[8])
        
            dictfile = load(filename)
            dictkeys = keys(dictfile)
            for dictkey in dictkeys
                data = dictfile[dictkey]
                if method == "ssa"
                    name = recover_keys_ssa(dictkey)
                    a,e,t = 0.0
                elseif method == "diff"
                    name,a,e,t = recover_keys_ssa(dictkey)
                end

                for word in [name,method,W,a,e,t,count_T(data)]
                    write(logfile, string(word))
                    write(logfile, ",\t")
                end
                write(logfile, "\n")
            end
        end
        close(logfile)
    end

    function list_T(log_loc)
        data = Array(CSV.read(log_loc, DataFrame))[:,1:end-1]
        samples = size(data)[1]

        name_f = []
        method_f = []
        W_f = []
        a_f = []
        e_f = []
        t_f = []
        name_h = []
        method_h = []
        W_h = []
        a_h = []
        e_h = []
        t_h = []


        for i in 1:samples
            name = data[i,1]
            method = data[i,2]
            W = data[i,3]
            a = data[i,4]
            e = data[i,5]
            t = data[i,6]
            Tf = filter(!iszero, data[i,7:7+k_max])
            Th = filter(!iszero, data[i,8+k_max:end-1])

            for T in Tf
                name_f = push!(name_f,(name,T))
                method_f = push!(method_f,(method,T))
                W_f = push!(W_f,(W,T))
                a_f = push!(a_f,(a,T))
                e_f = push!(e_f,(e,T))
                t_f = push!(t_f,(t,T))
            end

            for T in Th
                name_h = push!(name_h,(name,T))
                method_h = push!(method_h,(method,T))
                W_h = push!(W_h,(W,T))
                a_h = push!(a_h,(a,T))
                e_h = push!(e_h,(e,T))
                t_h = push!(t_h,(t,T))
            end
        end

        name_f = [name_f[i] for i in 1:length(name_f)]
        method_f = [method_f[i] for i in 1:length(method_f)]
        W_f = [W_f[i] for i in 1:length(W_f)]
        a_f = [a_f[i] for i in 1:length(a_f)]
        e_f = [e_f[i] for i in 1:length(e_f)]
        t_f = [t_f[i] for i in 1:length(t_f)]
        name_h = [name_h[i] for i in 1:length(name_h)]
        method_h = [method_h[i] for i in 1:length(method_h)]
        W_h = [W_h[i] for i in 1:length(W_h)]
        a_h = [a_h[i] for i in 1:length(a_h)]
        e_h = [e_h[i] for i in 1:length(e_h)]
        t_h = [t_h[i] for i in 1:length(t_h)]

        a_h = filter(!iszero, [a_h[i][1] for i=1:length(a_h)])
        e_h = filter(!iszero, [e_h[i][1] for i=1:length(e_h)])
        t_h = filter(!iszero, [t_h[i][1] for i=1:length(t_h)])
        a_f = filter(!iszero, [a_f[i][1] for i=1:length(a_f)])
        e_f = filter(!iszero, [e_f[i][1] for i=1:length(e_f)])
        t_f = filter(!iszero, [t_f[i][1] for i=1:length(t_f)])

        return name_f,method_f,W_f,a_f,e_f,t_f,name_h,method_h,W_h,a_h,e_h,t_h
    end
    #name_f,method_f,W_f,a_f,e_f,t_f,name_h,method_h,W_h,a_h,e_h,t_h

    function plotsave_bothlist(list1,list2,name,title,xlabel)
        p = scatter(list1,c=:black,marker = :x,xlabel=xlabel,ylabel="detected T[days]",title=title)
        p = scatter!(list2,c=:red,marker = :x)
        savefig(plot(p,dpi=800,legend=:none),dir*name)
    end

    function readfiles_filterT(datadir,pp)
        filenames = readdir(datadir)
        modelist=[]
        symbollist=[]
        for filename = filenames
            filenamesplit = split(filename,"_")
            method = filenamesplit[2]
            name = filenamesplit[3]
            W = parse(Int64,filenamesplit[4])
            k = parse(Int64,filenamesplit[5])
            a = parse(Float64,filenamesplit[6])
            e = parse(Float64,filenamesplit[7])
            t = parse(Float64,filenamesplit[8])
            dictfile = load(datadir.*filename)
            EOF = dictfile["data"]

            epsi_T = 20
            epsi_A = 0.1
            N = 12052
            k = 16
    
            pairs,periods_F = pairing_full(EOF,epsi_T,epsi_A,false)
            periods_HHT = Float64[0.5 .* W ./ iterated_hilbert(EOF[:,i],4)[1] for i=1:k]

         #   println("$(method)_$(name)_$(W)_$(k)_$(a)_$(e)_$(t)")

            for (i,period) in enumerate(periods_F)
                if period >= pp
                    symbol = "$(method)_$(name)_$(W)_$(k)_$(a)_$(e)_$(t)_$(period)_CC"
                    println(symbol)
                    modelist=push!(modelist,EOF[:,pairs[i][1]] .+ EOF[:,pairs[i][2]])
                    symbollist=push!(symbollist,symbol)
                end
            end


            for (i,period) in enumerate(periods_HHT)
                if period >= pp
                    symbol = "$(method)_$(name)_$(W)_$(k)_$(a)_$(e)_$(t)_$(period)_HHT"
                    println(symbol)
                    modelist=push!(modelist,EOF[:,i])
                    symbollist=push!(symbollist,symbol)
                end
            end

        end
        return modelist,symbollist
    end

"""
datadir = "/net/scratch/lschulz/testsampling/"
log_loc = dir*"loc"
readfiles_countT(datadir,log_loc)
name_f,method_f,W_f,a_f,e_f,t_f,name_h,method_h,W_h,a_h,e_h,t_h = list_T(log_loc)
plotsave_bothlist(name_f,name_h,"locations","Mode T detection","locations")
plotsave_bothlist(method_f,method_h,"method","Mode T detection","method")
plotsave_bothlist(W_f,W_h,"W","Mode T detection","embedding[days]")
plotsave_bothlist(a_f,a_h,"a","Mode T detection","diffusion sampling importance")
plotsave_bothlist(e_f,e_h,"e","Mode T detection","diffusion kernel strength")
plotsave_bothlist(t_f,t_h,"t","Mode T detection","diffusion iteration number")
"""