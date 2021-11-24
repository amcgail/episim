# episim

I wrote this collection of notebooks, and the `episim` module, to document my research as a part of the publication `Only as Safe as your Riskiest Contact` (forthcoming).
The module includes my quite inefficient code for simulating stochastic networked disease spread, but the final version of the paper made use of Ryan McGee's [SEIRS+](https://github.com/ryansmcgee/seirsplus) framework.

## Guide to the files

+ `epimodel/sampling.py` Contains all the algorithms I use in this study for sampling from a network for vaccination.
+ `simulation_helper.py` Contains some helpers: the default argsets (strat + R0_mean + VACCINATE_P + INITIAL_INFECT) and a clever iterator over limited combinations of arguments to the simulation,
+ `simulation_manager.py` All code for interfacing with SEIRS+: the simulation_manager class, which includes network generation, simulation, results loading and saving, and summarization.
+ `1a simulations*.ipynb` Runs the central argument-sets on three contexts: hs1, hs2, and the synthetic network. 
  These are separated into separate notebooks, but are essentially identical code, except for the bit that initializes the network.
+ `1b extending*.ipynb` Extends analyses to sampling.friendHighDegNormalErr and sampling.friendHighDegRandTopN at various levels of error.
  This helps to demonstrate how much the central method of the paper, sampling.friendHighDeg, is robust to errors in estimation of the extent of individuals' contacts' contact with others.
+ `2b critical stats*.ipynb` These notebooks develop and implement a method for visualizing the results of these many, many simulations.
+ `3 network diagram*.ipynb` These create diagrams of the networks used in these three contexts.
+ `4 example simulations.ipynb` Produces figures which summarize and give examples of the simulations.
+ `5 final results.ipynb` Produces a few paragraphs which show up in the text.

## Reproducing the results in this paper.

Ideally you should be able to just run these notebooks, although there are a few things to keep in mind:

+ You should have at least 16GB of RAM, and a lot of time. The `synthset 3000` simulation maxed out my desktop for 24h.
+ The simulations are stochastic, so the results you get will be slightly different from those reported in the paper.
+ My specific results are not included in this repository, because they are quite large. If you would like access, send an email to am2873@cornell.edu
+ Make sure the relevant manager.load_models('<NAME>') and manager.dump_models('<NAME>') match each other