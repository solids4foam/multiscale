# OpenFOAM toolbox for RVE-based solid mechanics multiscale modelling
This toolbox contains code and test cases from the paper Wu K, Tukovic Z, Cardiff P, Ivankovic (2002) Hierarchical RVE-based multiscale modeling of nonlinear heterogeneous materials using the finite volume method, International Journal for Multiscale Computational Engineering, 10.1615/IntJMultCompEng.2022040214.

# How to install
This toolbox requires foam-extend-4.1 and solids4foam-v2.0: see [https://solids4foam.github.io](https://solids4foam.github.io) for solids4foam installation instructions. Once solids4foam-v2.0 is compiled with foam-extend-4.1, the multiscale toolbox can be compiled with the included build script:
```
> ./Allwmake
```

# Test cases
The test cases are located in the `tutorials` directory. To run the single-scale hole-in-a-plate case:
```
> cd tutorials/singleScale/incr_TL && solids4Foam
```
To run the single-scale hole-in-a-plate case:
```
> cd tutorials/multiscale/incr_TL/macro && solid4Foam
```

# How to cite
```
@article{Wu_2022
	author  = {Ke Wu and Željko Tukovic and Philip  Cardiff and Alojz Ivankovic},
	title   = {HIERARCHICAL RVE-BASED MULTISCALE MODELING OF NONLINEAR HETEROGENEOUS MATERIALS USING THE FINITE VOLUME METHOD},
	journal = {International Journal for Multiscale Computational Engineering},
	issn    = {1543-1649},
        doi     = {10.1615/IntJMultCompEng.2022040214},
	year    = {2022},
	volume  = {20},
	number  = {2},
	pages   = {83--103}
}
```

# Who can I contact about this toolbox?
* philip.cardiff@ucd.ie
* zeljko.tukovic@fsb.hr


# Acknowledgements
This work has emanated from research conducted with the financial support of the Science Foundation Ireland under Grant No. 16/RC/3872.
Additionally, P Cardiff gratefully acknowledges financial support from the Irish Research Council through the Laureate programme, Grant No. IRCLA/2017/45, and from Bekaert through the University Technology Centre (UTC Phases I and II) at University College Dublin.
