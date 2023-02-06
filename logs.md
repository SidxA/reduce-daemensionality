# logfile

## 01/28    1h
    created repository structure in cluster
    added authentification
    added chatgpt functionality for commenting

    aims
    quick restructure of the processing functions
    hilbert trafo iteration convergence behavior comparison

## 01/29    1h
    added clearer mode extraction from dimensionality reduction results for the iteration convergence routine
    started cleaning up the toolbox

## 01/30    2h
    mode extraction not working properly
    just use the directly saved data, maybe build some quick struct to work on
    protoperiod not working properly, possibly wrong vari
    iteration scheme itself works

    naive iteration of hilbert trafo does not converge
    there is little amplitude modulation on all modes

## 02/01    2h
    embedding length of 1 year works better
    minor amplitude modulation on harmonics 2,3,...
    starts to incorporate different information then resolution of oscillation
    need some plots to explain this
    a lot of 'trend' modes in NLSA, sometimes split up of oscillatory  into 2 'linear trend'

    maybe new approach to code logging: new analysis script, instead of rework

## 02/02    2h
    strategic paper: the main thematic issue is the disentangling of the non-sinuisoid oscillatory trend
    (that is distributed among different harmonics and is assumed to have persistent amplitude)
    from the rest of the information (such as trend, amplitude modulations)

    BUT: this is possibly not the goal of the work on the paper

    it is confirmed now that pure hilbert iterations do not converge
    there is a tiny tiny portion of modes, mainly first and second harmonic
    that have persistent zero_counts up to 3 iterations, e.g. the amplitude modulation is low enough

    dirty hard coded filter to find these modes works
    this 'kind-of-works', but really shuldnt be published ?
    

## 02/06 1h
    combined plotting on the results of the naive hard filtering
