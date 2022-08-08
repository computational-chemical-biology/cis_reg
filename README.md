# Rhodosporidium cis-regulatory elements discovery 
Companion material to the [](publication).

## Installation

```
conda env create -f environment.yml
conda activate cis_reg
Rscript install_packages.R
```

## Usage
 
```
conda activate cis_reg
jupyter notebook
```

## Update

If a new library needs to be installed, don't forget to update the environment.yml file

```
conda env export | grep -v "^prefix: " > environment.yml
```

### License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details

