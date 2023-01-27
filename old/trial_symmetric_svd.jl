C = 2001
d*s = 144
(T-C*s)/(C-1)=d*s

function symmetric_embedding()
    D = [1,2,3,4,6,8,9,12,16,18,24,36,48]
    S = 144 / D
    C = 2001
    
function sampler(signal::Vector{Float64},d::Int,s::Int,C::Int)
    #samples s
    samples = mean.([signal[Int.(Array(1:s).+i*s)] for i=0:length(signal)/s-1])
    #delays d
    delays = hcat([samples[i*d+1:C+i*d] for i=0:C-1]...)
    return delays
end