# MultiGWAS Table of Contents
   * [MultiGWAS Table of Contents](#multigwas-table-of-contents)
   * [Installation](#installation)
      * [General steps to install multiGWAS on a Linux system](#general-steps-to-install-multigwas-on-a-linux-system)
      * [Specific instructions to install multiGWAS on a Linux Ubuntu](#specific-instructions-to-install-multigwas-on-a-linux-ubuntu)
         * [Install external software:](#install-external-software)
         * [Install MultiGWAS tool:](#install-multigwas-tool)
   * [Executing the MultiGWAS tool](#executing-the-multigwas-tool)
      * [Using the command line interface:](#using-the-command-line-interface)
      * [Using the graphical user interface:](#using-the-graphical-user-interface)
   * [Running the examples](#running-the-examples)
   * [General usage](#general-usage)
   * [Genomic data formats](#genomic-data-formats)
      * [Genotype:](#genotype)
      * [Phenotype](#phenotype)
   * [Configuration file](#configuration-file)
   * [Considerations](#considerations)
      * [Implementation](#implementation)
      * [Number of SNPs in Manhattan and QQ plots](#number-of-snps-in-manhattan-and-qq-plots)
      * [Correction for multiple testing](#correction-for-multiple-testing)

## Full installation steps:
We provide two ready-to-use full installations: one for Ubuntu 18.xx and another for Ubuntu 20.xx.  For both installations, open a linux terminal, clone the MultiGWAS repository, change to the multiGWAS folder, and run the install script:

### Ubuntu 18.xx
  ```
	git clone https://github.com/agrosavia-bioinformatics/multiGWAS-ubuntu18.git
	cd multiGWAS-ubuntu18
	sh INSTALL.sh 
  ```

### Ubuntu 20.xx
  ```
	git clone https://github.com/agrosavia-bioinformatics/multiGWAS-ubuntu20.git
	cd multiGWAS-ubuntu20
	sh INSTALL.sh 
  ```


