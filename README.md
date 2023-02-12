# Shared source code

The source code in this repository is an extension of [DAFoam discrete-adjoint optimization platform]([DAFoam: Discrete Adjoint with OpenFOAM for High-fidelity Multidisciplinary Design Optimization | DAFoam](https://dafoam.github.io/index.html)). Some of the codes are directly modified from the source file of DAFoam. It can only be used when it is compiled together with DAFoam. Contact the author for further instructions on compilation.

* `wall.C` contains a simple implementation of the wall-distance field calculation method proposed in the article.
* `DALaunderSharmaKE.*` is an implementation of the modified Launder-Sharma k-epsilon model that is suitable for topology optimization. 
* `DAkOmegaSST.*` is a modified version of the original k-Omega SST model in DAFoam. It uses the wall distance field calculated by the proposed method and it also has a penalty term related to Darcy's source term.