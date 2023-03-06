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

## 02/08    1h
    since iterated hilbert trafo does not naturally converge as a frequency estimation criterion,
    lets use spectral power by fft  to show the harmonic structure

    the question of ridicolous modes in NLSA: compare W=7a with W=1a for all spectral results

    single fft implemented
    repaired access to the broken metadata file

    need to create large figures for individual spots w. signal, and multiple modes, spectra


## 02/09    2h
    put the fourier spectra for individual modes together into large overview plots for W=7 and W=1

## 02/13    8h
    created the cauchy based sorting algorithm for frequency harmonicity estimation based on threshold
    due to the steepness of the cauchy peaks, the maximum of the residuals couldnt be used for threshold differentiation
    finally, the un-normalized gauss peaks with free width yielded the best spectral fits
    finished some threshold estimation table, mediocre results, keep it at conservative 10%
    finished the sorting and plotting routine to compare SSA with NLSA

    results are pretty miserable: the goal of steady-amplitude trend with individual trends is missed by the selected modes
    there is still so much amplitude modulation in the RC
    are the actual EOF responsible for it?
    jump: they are, lets do a high resolution W=1 re-calculation with half hourly data ATTENTION this would also induce some day night dynamic?
    so maybe check them manually

    the FFt of the NLSA trends have the additional harmonics, this is looking good so far
    NLSA sucks in the case of some signals, the parameter tuning might not have been good enough

##  02/15   8h
    gaussian based harmonicity algorithm:
    unnormalized spectral power
    and free third gaussian parameter
    lead to pretty stable results

## 02/22  2h
    presentation pictures

##  03/02   4h
    figures work: fft peak & seasonality behavior

##  03/04   2h
    figures work: comparison ssa/nlsa and fundamental/harmonic

## 03/05    2h
    bigf figure
    modes figure
    reconstructions figure

## 03/06    1h
    started working on the year diagram for amplitude and phase