# *Sherlock-Genome*
## *Sherlock-Genome*: A R Shiny App for Genomic Analysis and Visualization

### Welcome!

Welcome to *Sherlock-Genome*: A R Shiny App for Genomic Analysis and Visualization! This tool was developed to enable researchers to review and explore project results all in one place, conduct a wide range of genomic analyses, and generate several different types of visualizations. 

### Documentation, Data Requirements, and Demos
There is documentation available [here](https://github.com/xtmgah/Sherlock-Genome/wiki). This documentation highlights each module in the app and their functionality, in addition to information on uploading a user's own data.

If a user wishes to upload their own data, there is a [Data Requirements spreadsheet](https://docs.google.com/spreadsheets/d/1XjoYzG1mQQiw0Dqm4shtXPPnJDAmLY4GPHsi1xuazKg/edit#gid=2016382736) that explains each file required for each module, how to generate the files, and what column names are necessary. These data requirements can also be found in the Documentation module of *Sherlock-Genome* under Data Requirement Info.

There are also a few video tutorials available to get you started:

1. Downloading the Genomic Data and moving it into the application correctly.
2. Access the Data Requirements Info tab in the Documentation module.
3. Uploading user project data.
3. Demonstration of each module.

### Get Started Using *Sherlock-Genome*

To access and use the *Sherlock-Genome* R Shiny App, follow the steps below:

#### Install the SKIT package

*Sherlock-Genome* requires several packages to run. These are installed if needed and loaded automatically upon opening the app. However, one package, SKIT, must be installed on its own outside of the app. SKIT is a method that can be used in association testing. To download and install the SKIT package, follow the steps below.

1. Go to the [page](https://dceg.cancer.gov/tools/analysis/skit) that includes the SKIT package file for download.
2. Download the tar.gz file.
3. Open R or RStudio.
3. Install the SKIT package-
- Option 1: Use the command: `install.packages('/path/to/skit-0.0.2.tar.gz', repos=NULL, type='source')` to install the package. **Be sure to set up the path to the tar.gz file.**

- Option 2: Click 'Install' under the Packages tab in RStudio, select 'Install from Package Archive File', and then browse for the tar.gz file you just downloaded. Click 'Install'.

#### Run *Sherlock-Genome* on your local machine

1. Download the application repository for the app from GitHub. The repository link can be found [here](https://github.com/xtmgah/Sherlock-Genome.git).
2. Decompress the file.
3. Download the 'Genomic Data' folder from [Google Drive](https://drive.google.com/file/d/11q90cjMoiBzLMTyw7tVuXtYb8Y3T3n8t/view?usp=drive_link) to add into the application.
4. Locate the 'Genomic Data' folder you just downloaded, unzip it, and move it into the www/ folder in the application.
5. Open the ui.R file, which can be found in the main Sherlock_Genome directory. 
6. Check the working directory by using the `getwd()` command. The working directory should be something like: `/path/to/Sherlock_Genome`. If it is not, close R or RStudio and open ui.R again. Or set, the working directory yourself, using the `setwd('/path/to/Sherlock_Genome')`.
7. Run the app using the `runApp()` command. Or, if using RStudio, you can click the 'Run App' button.

