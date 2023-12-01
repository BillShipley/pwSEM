# pwSEM
An R package to impliment my new piecewise SEM function called pwSEM. This supports and combination linear, generalized linear, generalized linear mixed models as well as additive, generalized additive, generalized additive mixed models.  The ability to use (generalized) additive (mixed) models to test the d-separation claims is also supported.

Installation

First, make sure that you have the latest version of R on your computer, or at least a very recent version. To install the package from GitHub, type the following commands in your R or RStudio console:

# install package dependencies first
install.packages("ggm")
install.packages("gamm4")
install.packages("mgcv")

# install and load devtools to be able to install packages from GitHub with install_github
install.packages("devtools")
library(devtools)

# install pwSEM from Bill's GitHub
install_github("BillShipley/pwSEM")
library(pwSEM)
?pwSEM

If everything worked, the last command should have opened the help file for the pwSEM function.
If things do not work
Step 1

Read the error messages and make sure all packages dependencies are installed and loaded, especially package lavaan and ggm. If a message says that a package could not be loaded, try installing it manually by typing:

# manually installing a dependencies
install.packages("packagename")

until all packages are installed. A current bug in install_github on Windows prevents the installation of package dependencies of dependencies (lavaan and ggm).
Step 2

Although this should not be required, for certain packages with compiled code, Rtools and MiKTeX need to be installed on Windows to be able to build source packages. For Mac users, Xcode is required and can be installed through the apple store. Here is a more detailed list of prerequisites for building source packages for Windows, Mac and Linux.
Step 3

In case something goes wrong with the package installation and the previous instructions do not work, here is the code for every function in the package. The code can be copied and pasted in the R console to get the definition of all functions. However, in this case, help files won't be available.
