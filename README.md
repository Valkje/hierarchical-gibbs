# A hierarchical Gibbs sampler for discovering cognitive stages

This sampler will estimate the parameter values for the hierarchy displayed in the figure below from EEG data.

![parameter hierarchy](Graphical\ models/graphical_model_2.png)

## Requirements

1. MATLAB 2021a (other versions may work just as fine but have simply not been tested)
2. EEGLAB 2021.0 (again, other versions may also work)
3. Preprocessed EEG data that has been fed through a principal component analysis (PCA), contained in the following variables:
    * `all_signal`: One large matrix that holds the chosen number of principal components (PCs) over time, concatenated for all subjects. ([timeSteps * numTrials] x numPCs)
    * `all_x`: The onset indices of every trial.
    * `all_y`: The final indices of every trial.
    * `all_subjects`: An array which provides a subject number for every trial.
    * `all_conds`: A condition number for every trial.

## Setup

Open `main.m`.

### Paths

1. Change the `path` variable to the directory `main.m` is placed in.
2. Replace `genpath([path 'synth_data'])` with the directory your data are placed in.
3. Replace `[path 'eeglab2021.0']` with the path to your EEGLAB directory.
4. Change `savePath` to where you want to save your data.

### Analysis

1. Change `max_iter` to your desired maximum number of iterations.
2. Change `cond` to the condition you want to analyse.
3. Change `n` to the number of bumps you expect to be present.
4. If you want to change the visualisation behaviour (which can be costly), then open `hierarchical_gibbs/visualize.m` and change the following variables:
    * `burnIn`: The number of samples from the first iterations of the sampler that will not be displayed.
    * `renderLag`: After how many iterations all figures are updated.
    * `thinLag`: If the value of `thinLag` is x, then only display every xth sample in the figures (important for histograms).

## Running the analysis

Simply run `main.m`. Once the run has completed, all sampled parameters will be saved to the `savePath` directory (with the run number appended to it).
