#using FFTW
#using EmpiricalModeDecomposition

#phase estimation
#take array of N x k
#first take only ne signal and perform EMD
function perform_EMD(data)
    data = centralizer(data)
    #variables
    n=3
    Sifting_n=10
    S=4
    n_IMF=4
    data_emd = emd(data,EMDSetting(n,Sifting_n,S,n_IMF))
    println("performing sifting such with standard deviations")
    [println(std(data_emd[:,i])) for i in 1:n]
    phase = hcat([atan.(real(hilbert_transform(data_emd[:,i])),imag(hilbert_transform(data_emd[:,i]))) for i in 1:n]...)
    return data_emd,phase
end

function compare_EMD(table)
    name = "RC_phases_"
    for j in 1:10
        x = centralizer(table[:,j])
        p_sig = plot(x,title="component",legend=:none)
        p_sing_phase = plot(phases(x),title="phase of component",legend=:none)
        data_emd,phase = perform_EMD(x)
        p_imf = plot(data_emd,title="different harmonic signals",legend=:none)
        p_phase = plot(phase,title = "phases",legend=:none)
        savefig(plot(p_sig,p_sing_phase,p_imf,p_phase),dir*name*string(j))
    end
    return nothing
end

function individual_plot(list,name)
    l = length(list)
    p=[]
    for i in 1:l
        p=push!(p,plot(list[i]))
    end
    savefig(plot(p...,layout=(l,1)),dir*name)
end

#------------------------------------------------------
#what i need is an algorithm that performs the hilbert trafo and spits out the phase

function phases(signal)
    x = hilbert_transform(signal)
    return atan.(imag(x),real(x))
end


function quickplot()
    name="compact"
    x = centralizer(load_data())
    p_sig = plot(x,title="component",legend=:none)
    p_sing_phase = plot(phases(x),title="phase of component",legend=:none)
    #variables
    n=3
    Sifting_n=10
    S=4
    n_IMF=10
    data_emd = emd(data,EMDSetting(n,Sifting_n,S,n_IMF))
    println("performing sifting such with standard deviations")
    [println(std(data_emd[:,i])) for i in 1:n]
    phase = hcat([atan.(real(hilbert_transform(data_emd[:,i])),imag(hilbert_transform(data_emd[:,i]))) for i in 1:n]...)

    p_imf = plot(data_emd,title="different harmonic signals",legend=:none)
    p_phase = plot(phase,title = "phases",legend=:none)
    savefig(plot(p_sig,p_sing_phase,p_imf,p_phase),dir*name*string(j))
end

function quick_imf()
    x = centralizer(load_data())
    p_sig = plot(x,title="signal",legend=:none)
    p=[]
    #variables
    n=3
    Sifting_n=8
    S=8
    n_IMF=10
    data_emd = emd(data,EMDSetting(n,Sifting_n,S,n_IMF))
    for i in 1:5
        p = push!(p,plot(data_emd[:,i],legend=:none))
    end
    savefig(plot(p_sig,p...,layout=(6,1),dpi=600),dir*"IMF")
end

function quick_RC()
    x = centralizer(load_data())
    p_sig = plot(x,title="signal",legend=:none)
    p=[]
    for i in 1:5
        p = push!(p,plot(rc[:,i],legend=:none))
    end
    savefig(plot(p_sig,p...,layout=(6,1),dpi=600),dir*"RC")
end

function RC_phase_synch()
    z = zeros(10,10)
    for i = 1:10,j = 1:10
        z[i,j] = norm(phases(rc[:,i]).-phases(rc[:,j]))
    end
    return z 
end