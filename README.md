# episim

I wrote this collection of notebooks, and the `episim` module, to document my research as a part of the publication `Only as Safe as your Riskiest Contact` (forthcoming).
The module includes my quite inefficient code for simulating stochastic networked disease spread, but the final version of the paper made use of Ryan McGee's [SEIRS+](https://github.com/ryansmcgee/seirsplus) framework.

## Guide to the files

+ `epimodel/sampling.py`
+ `1a*.ipynb` Runs the central argument-sets on three contexts: hs1, hs2, and the synthetic network. 
  These are separated into separate notebooks, but are essentially identical code, except for the bit that initializes the network.
+ `1b*.ipynb` Extends analyses to `sampling.friendHighDegNormalErr` and `sampling.friendHighDegRandTopN` at various levels of error.
  This helps to demonstrate how much the central method of the paper, `sampling.friendHighDeg`, is robust to errors in estimation of the extent of individuals' contacts' contact with others.
+ `2b*.ipynb` These notebooks develop and implement a method for visualizing the results of these many, many simulations.
+ `3*.ipynb` These create diagrams of the networks used in these three contexts.

## Reproducing the results in this paper.

Ideally you should be able to just run these notebooks, although there are a few things to keep in mind:

+ You should have at least 16GB of RAM, and a lot of time. The `synthset 3000` simulation maxed out my desktop for 24h.