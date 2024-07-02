# *Sherlock-Genome*
## *Sherlock-Genome*: A R Shiny App for Genomic Analysis and Visualization

### Welcome!

Welcome to *Sherlock-Genome*: A R Shiny App for Genomic Analysis and Visualization! This tool was developed to enable researchers to review and explore project results all in one place, conduct a wide range of genomic analyses, and generate several different types of visualizations. Below is a visual of the many modules and functions available in the app:

![Visual of the capabilites of Sherlock-Genome](https://github.com/xtmgah/Sherlock-Genome/blob/master/Documentation_screenshots/Sherlock-Genome_main_figure.png)

### Documentation, Data Requirements, and Demos
There is documentation available [here](https://github.com/xtmgah/Sherlock-Genome/wiki/Documentation-for-Sherlock%E2%80%90Genome). This documentation highlights each module in the app and their functionality, in addition to information on uploading a user's own data.

If a user wishes to upload their own data, there is a [Data Requirements spreadsheet](https://docs.google.com/spreadsheets/d/15CiRPx3A5unRMVemROYR632x7A9dpcT-/edit?usp=sharing&ouid=113024141265719270420&rtpof=true&sd=true) that explains each file required for each module, how to generate the files, and what column names are necessary. These data requirements can also be found in the Documentation module of *Sherlock-Genome* under Data Requirement Info.

There are also a few video tutorials available to get you started. Click any of the links below to take you to the video.

1. [Downloading the Genomic Data and moving it into the application correctly.](https://drive.google.com/file/d/1C3xqRqdCx2D1iQQt6s33N9gFsMTlxaPS/view?usp=sharing)
2. [Access the Data Requirements Info tab in the Documentation module.](https://drive.google.com/file/d/1-D5mue2lD180qFViIDjS6mIt6YJP9xmr/view?usp=sharing)
3. [Uploading user project data.](https://drive.google.com/file/d/1VgEskRW7Uwa-ZWPlTlWsKQy_JdsCN26i/view?usp=sharing)
3. [Demonstration of each module.](https://drive.google.com/file/d/1F42IULr_DydC4DshPM9JFPPXAdEHezpO/view?usp=sharing)

### Get Started Using *Sherlock-Genome*

To access and use the *Sherlock-Genome* R Shiny App, follow the steps below. **Users must have R version of at least 4.0.0 installed.**

#### Install the SKIT package

*Sherlock-Genome* requires several packages to run. These are installed if needed and loaded automatically upon opening the app. However, one package, SKIT, must be installed on its own outside of the app. SKIT is a method that can be used in association testing. To download and install the SKIT package, follow the steps below.

1. Go to the [page](https://dceg.cancer.gov/tools/analysis/skit) that includes the SKIT package file for download.
2. Download the tar.gz file.
3. Open R or RStudio.
3. Install the SKIT package-
- Option 1: Use the command: `install.packages('/path/to/skit-0.0.2.tar.gz', repos=NULL, type='source')` to install the package. **Be sure to set up the path to the tar.gz file.**

- Option 2: Click 'Install' under the Packages tab in RStudio, select 'Install from Package Archive File', and then browse for the tar.gz file you just downloaded. Click 'Install'.

#### Run *Sherlock-Genome* on your local machine

1. Download the application repository for the app from GitHub. The repository link can be found [here](https://github.com/xtmgah/Sherlock-Genome.git). To download the repository, click on the link provided, then click on the green `Code` dropdown near the top right of the screen. Select 'Download ZIP' from the dropdown. Th app zip file will download (most likely to your Downloads folder, unless you have downloads set to be saved elsewhere.)
2. Go to the folder you downloaded the app zip file to. Decompress the file.
3. Download the 'Genomic Data' folder from [Google Drive](https://drive.google.com/file/d/1tZ6aPA5LvVDnIbrxUrnUkrukySNQ5GGr/view?usp=sharing) to add into the application.
4. Locate the 'Genomic Data' folder you just downloaded, unzip it, and move it into the folder named 'www/' in the application.
5. Open the ui.R file using R or RStudio, which can be found in the main Sherlock-Genome directory. 
6. Check the working directory by using the `getwd()` command. The working directory should be something like: `/path/to/Sherlock-Genome-master`. If it is not, close R or RStudio and open ui.R again. Or set, the working directory yourself, using the `setwd('/path/to/Sherlock-Genome-master')`.
7. Run the app using the `runApp()` command. Or, if using RStudio, you can click the 'Run App' button.

A FEW ADDITIONAL NOTES: 
1. If you are having trouble getting the app to run initially due to package installation errors, try running lines 6-54 in the ui.R file before running the app.
2. We suggest opening the app in a browser to open in a browser window after initializing the app to run- sometimes the figures are too big to show up in just the app window.
