Readme
================
S Chatterjee
2023-04-15

# BFDAMs: Bayesian Functional Data Analysis over Dependent Regions and Its Application for Identification of Differentially Methylated Regions.

# Methodological Pipeline

<img src="C:/Users/suvchat/OneDrive - Indiana University/Current_Projects/BFDA/BFDAMs/Picture1.jpg" width="700" align="center">

## Usage

BFDAM1(): This function fits a Bayesian functional data analysis model
by assuming the contiguous sequence of windows/regions are independent
to each other while accounting for dependency of observations inside a
particular window/region.

BFDAM2(): This function fits a Bayesian functional data analysis model
by assuming the contiguous sequence of windows/regions are dependent to
each other as well as accounting for dependency of observations inside a
particular window/region.

## Arguments

Both functions BFDAM1() and BFDAM2() have the same set of arguments.
These functions have 8 arguments. For example, to run the functions
under their default settings the user can use the following;

``` r
Ex: fit = fit.BFDAM1(beta.df, chr.pos, nCpGs = 100,
                     iter=20000, burn = 1000, seed = 1234,
                     control.group = 'Normal',
                     case.group = 'Tumor')

Ex: fit = fit.BFDAM2(beta.df, chr.pos, nCpGs = 100, 
                     iter = 20000, burn = 1000, seed = 1234,
                     control.group = 'Normal',
                     case.group = 'Tumor')
```

### The table below details the required arguments:

| Parameter     | Default  | Description                                                                                                          |
|:--------------|:--------:|:---------------------------------------------------------------------------------------------------------------------|
| beta.df       |          | A data frame object that will contain the methylation data where rows are CpG sites and the columns are the samples. |
| chr.pos       |          | A data frame object that will contain the chromosome number and position of the CpG sites.                           |
| nCpGs         |   100    | Number of CpG sites to be included per genomic region/window.                                                        |
| iter          |  20000   | Number of MCMC iterations to be run                                                                                  |
| burn          |   1000   | Number of burn-ins                                                                                                   |
| seed          |   1234   | A random number for reproducibilty of results                                                                        |
| control.group | ‘Normal’ | Label for controls in beta.df data frame                                                                             |
| case.group    | ‘Tumor’  | Label for cases in beta.df data frame                                                                                |

## Values

### Both the functions BFDAM1 and BFDAM2 will return a data frame with the following values summarized in the table below.

| Object    | Description                                                                                                                                                     |
|:----------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| chr       | Chromosome number of the genomic region/window modeled.                                                                                                         |
| window    | Genomic region/window number for the respective contiguous set of CpGs.                                                                                         |
| pos.start | The starting genomic position of the modeled genomic window/region.                                                                                             |
| pos.end   | The ending genomic position of the modeled genomic window/region.                                                                                               |
| status    | The predicted differential methylation status for the respective genomic window/region modeled. Each window/region can be predicted as either a DMR or Non-DMR. |

## Working Example

Since in your paper we recommend using the Bayesian Functional data
dependent model (BFDAM2) for detecting differentially methylated
regions, here we provide a minimum working example (MWE) for users to
implement BFDAM2. In this demo representation we will restrict ourselves
to a lower number of iterations and burn-ins to save some time in
reproducing this MWE.

**A brief description of the example data**

The example data is a subset of the real data used in the paper and can
be found here. This example data contains the methylation values on
chromosome 12 surrounding (+- 1Mb) the popular lung cancer gene KRAS.

Below we show the required input data for using the functions.

**Guide for loading in required data**

For all the analysis shown below the user needs to input the following
example data: From the “Data” folder please load the example data
(“LUAD2015_ExampleData.RData”) into your R environment. All the results
and analysis shown below is based on this example data set. Note: Once
the example data is loaded into the R environment the object will be
denoted as a list of 2 elements with the displayed name “data”.

**Guide for loading in helper functions**

Both of our developed models BFDAM1 and BFDAM2 uses a set of helper
functions to run. The user needs to source in the helper function file
(“BFDAMs_UtilityFns.R”) into their R environment before running the
analysis shown below. For plotting purposes (as shown below) the user
needs to source in the function “BFDAMs_PlotFn.R”. The above two R files
can be found under the “R” folder.

**Guide for loading in the main functions**

Both of our developed models BFDAM1 (“BFDA_independent_main.R”) and
BFDAM2 (“BFDA_dependent_main.R”) can be found under the “R” folder.
Since in this MWE we are only implementing the BFDAM2 model (dependent
method) the user needs to source in the file (“BFDA_dependent_main.R”)
into their R environment before running this demo analysis.

**Input Data (Snapshot of methylation beta-values data)**

``` r
load('C:\\Users\\suvchat\\Documents\\BFDA_GitHub\\LUAD2015_ExampleData.RData')
beta.df = data$beta.df
beta.df[1:5,1:8]
```

    ##            GSM1632880_Tumor GSM1632881_Normal GSM1632882_Tumor GSM1632883_Tumor
    ## cg01389585       0.39451706         0.6009137       0.32448934       0.37254102
    ## cg12642725       0.87778511         0.8162177       0.67998402       0.85591416
    ## cg19205533       0.28720751         0.4551641       0.09996816       0.36298645
    ## cg08830758       0.15903439         0.1938489       0.08301977       0.21019356
    ## cg24979348       0.01993344         0.0152749       0.01365316       0.01624346
    ##            GSM1632884_Normal GSM1632885_Tumor GSM1632886_Tumor GSM1632887_Tumor
    ## cg01389585        0.46949142       0.43230172       0.51259288       0.48507750
    ## cg12642725        0.79901909       0.80612768       0.82314826       0.77993175
    ## cg19205533        0.33971963       0.29181200       0.38693604       0.22205678
    ## cg08830758        0.13571027       0.18929850       0.18164961       0.12048252
    ## cg24979348        0.01438445       0.01261574       0.01088544       0.01509158

**Input Data (Snapshot of CpG site location data)**

``` r
load('C:\\Users\\suvchat\\Documents\\BFDA_GitHub\\LUAD2015_ExampleData.RData')
chr.pos = data$chr.pos
chr.pos[1:5,]
```

    ##              chr      pos
    ## cg01389585 chr12 15359440
    ## cg12642725 chr12 15371727
    ## cg19205533 chr12 15373987
    ## cg08830758 chr12 15374175
    ## cg24979348 chr12 15374303

**Running BFDAM2 (Dependent Bayesian FDA model) function on example
data**

``` r
load('C:\\Users\\suvchat\\Documents\\BFDA_GitHub\\LUAD2015_ExampleData.RData')
source('C:\\Users\\suvchat\\Documents\\BFDA_GitHub\\BFDAMs_UtilityFns.R')
source('C:\\Users\\suvchat\\Documents\\BFDA_GitHub\\BFDAM_dependent_main.R')
beta.df = data$beta.df
chr.pos = data$chr.pos
BFDAM2 = fit.BFDAM2(beta.df = beta.df, chr.pos = chr.pos, 
                    nCpGs = 100, iter = 500, burn = 50, 
                    seed = 1234, control.group = 'Normal',
                    case.group = 'Tumor')
```

    ## [1] "Running Genomic Window 1 of Chromosome 12"
    ## [1] "Running Genomic Window 2 of Chromosome 12"
    ## [1] "Running Genomic Window 3 of Chromosome 12"
    ## [1] "Running Genomic Window 4 of Chromosome 12"
    ## [1] "Running Genomic Window 5 of Chromosome 12"
    ## [1] "Running Genomic Window 6 of Chromosome 12"
    ## [1] "Running Genomic Window 7 of Chromosome 12"
    ## [1] "Running Genomic Window 8 of Chromosome 12"
    ## [1] "Running Genomic Window 9 of Chromosome 12"
    ## [1] "Running Genomic Window 10 of Chromosome 12"
    ## [1] "Running Genomic Window 11 of Chromosome 12"
    ## [1] "Running Genomic Window 12 of Chromosome 12"
    ## [1] "Running Genomic Window 13 of Chromosome 12"
    ## [1] "Running Genomic Window 14 of Chromosome 12"
    ## [1] "Running Genomic Window 15 of Chromosome 12"
    ## [1] "Running Genomic Window 16 of Chromosome 12"
    ## [1] "Running Genomic Window 17 of Chromosome 12"
    ## [1] "Running Genomic Window 18 of Chromosome 12"
    ## [1] "Running Genomic Window 19 of Chromosome 12"

**Visualizing the fitted curves for genomic window/region 13 in
Chromosome 12 of example data**

To aid in visualizing the fitted functional mean curves for any given
chromosome and genomic region/window we also provide the users with the
function “plot.fittedCurves()”. Below we demonstrate the use of this
function and plot the functional curves for the 2 groups. The blue and
red lines represent the fitted functional mean methylation curves for
healthy and cancer subjects respectively.

``` r
source('C:\\Users\\suvchat\\Documents\\BFDA_GitHub\\BFDAMs_PlotFn.R')
print(plot.fittedCurves(beta.df = beta.df, chr.pos = chr.pos, 
                        nCpGs = 100, iter = 500, burn = 50, 
                        seed = 1234, control.group = 'Normal',
                        case.group = 'Tumor',chr=12, window.plot=13))
```

![](BFDAM_readme_GitHub_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

**Understanding the output format**

The output will be a list item where each element of this list will be
the outputs pertaining to respective chromosome. So, let’s say your
methylation data spans over 22 chromosomes then the output will be a
list of 22 elements where each element will be the data frame of results
pertaining to the respective chromosome. In our example data since we
have only considered chromosome 12 our output is a list of length 1.

**Making sense of the output**

The output table as shown below depicts the genomic regions/windows that
were analyzed and their differential methylation status is depicted in
the status column. Windows/regions which are denoted as ‘DMR’ are
regions where the mean methylation rates were found to be significantly
different between the cancer and healthy subjects.

``` r
BFDAM2[[1]]
```

    ##    chr window pos.start  pos.end  status
    ## 1   12      1  15359440 16064533 Non-DMR
    ## 2   12      2  16064569 18797084     DMR
    ## 3   12      3  18824440 19836352 Non-DMR
    ## 4   12      4  19864189 21590544     DMR
    ## 5   12      5  21590632 22196840 Non-DMR
    ## 6   12      6  22198364 23858674     DMR
    ## 7   12      7  23909887 25055518     DMR
    ## 8   12      8  25055676 25538177     DMR
    ## 9   12      9  25538221 26288686 Non-DMR
    ## 10  12     10  26306247 27167667 Non-DMR
    ## 11  12     11  27170664 27863646     DMR
    ## 12  12     12  27863705 29302016 Non-DMR
    ## 13  12     13  29302035 30167296 Non-DMR
    ## 14  12     14  30251537 31254007 Non-DMR
    ## 15  12     15  31254036 31812558     DMR
    ## 16  12     16  31812580 32555628     DMR
    ## 17  12     17  32571851 34261006 Non-DMR
    ## 18  12     18  34261116 34496342 Non-DMR
    ## 19  12     19  34496768 34846306     DMR

## Contributions

If you find small bugs, larger issues, or have suggestions, please email
the maintainers at <suvchat@iu.edu> or <shrabanti.chowdhury@mssm.edu>.
Contributions (via pull requests or otherwise) are welcome.
