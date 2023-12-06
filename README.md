# ResORR - Reservoir Operations driven River Regulation

ResORR is a scalable, minimally parameterized model to estimate regulation of rivers due to upstream reservoir operations primarily using satellite observations of reservoir dynamcis. It is designed within the Reservoir Assessment Tool (RAT) framework, and can be easily initialized and run using an existing setup of [RAT 3.0](https://rat-satellitedams.readthedocs.io/en/latest/). 

Please refer to the **[documentation](https://resorr.readthedocs.io/en/latest/)** for a detailed description of the model and jupyter notebooks on how to run the model.
[![Documentation Status](https://readthedocs.org/projects/resorr/badge/?version=latest)](https://resorr.readthedocs.io/en/latest/?badge=latest)


## Installation
The ResORR model will be included in the future versions of RAT, until then, the model can be installed as a python module. If you are using `conda` or `mamba` for environment management, you can install the package using the following command:
```bash
git clone https://github.com/UW-SASWE/ResORR.git
conda develop ResORR/src/
```
Please note that you would need to have `conda-build` in your environment to be able to use the `conda develop` command. If you do not have `conda-build` installed, you can install it using the following command:
```bash
conda install conda-build
```

