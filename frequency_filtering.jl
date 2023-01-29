"""
code execution: extract EOF,lambda,signal from calculated files and jldsave in scratch
"""

W = Int(floor(7*365.25))

for method=["ssa","diff"], vari=1:2, spot=1:20
    read_modes_single_file(method,W,vari,spot)
end




data = load("/net/scratch/lschulz/modes/"*"name")


protophases = rec_protophase(RC,1)
protofreq = hcat([protofrequency(EOF[:,kk],yearsamples,1) for kk in 1:k]...)
