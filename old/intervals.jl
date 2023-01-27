# interval creation from time series X

#X = rand(537)

# structural: sin fit and take the quarter(?) period with axiom DeltaT = 1
# what about the cases where we dont get an integer: rational number might be divided into r/q integers

#phi = phase_variable
#T = period_length
#A = amplitude

#W = floor(T/4)


# smoothing interval finder
# input is the desired window length W

# first take the whole array and create a number of subarrays of the right length and position (loop)

# K = N - W +1

#while the end is not reached
#only admissable are even W !
#create a list of starting,ending points
function embed_half(X,W)
    N = length(X)
    st_pts = 1:Int(W/2):N-W+1
    Y = []
    for i in st_pts
        Y = push!(Y,X[i:i+W-1])
    end
    return Y
end


function multi_dim_split(X,W)

    function reshape_nested(X,W)
        dims = size(X)[1]
        windows = size(X[1])[1]
        Z = zeros(windows,W,dims)
        for i in 1:windows, j in 1:W, k in 1:dims 
            Z[i,j,k] = X[k][i][j]
        end
        return Z 
    end

    dims = size(X)[1]
    Y = []
    for i in 1:dims
        Y = push!(Y,embed_half(X[i,:],W))
    end
    return reshape_nested(Y,W)
end





#take 0*w:1*w
#take 1/2*w:3/2*w


# then create a function that takes one intervall and one desired length and performs sin,bell in order to smoothen the outlines