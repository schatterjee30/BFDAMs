
# Bayesian Functional Data Analysis over Dependent Regions and Its Application for Identification of Differentially Methylated Regions

This repository describes the implementation of our developed Bayesian
functional data analysis model BFDAM2 for identifying differentially
methylated regions in 450K array DNA methylation data. This repository
also contains all the required functions and files for reproducing the
results found in our manuscript.

# Methodological Pipeline

![alt text](https://github.com/schatterjee30/BFDAMs/blob/main/Images/Pipeline.jpg)

## Usage

BFDAM1(): This function fits a Bayesian functional data analysis model
by assuming the contiguous sequence of windows/regions are independent
to each other while accounting for the dependency of observations inside a
particular window/region.

BFDAM2(): This function fits a Bayesian functional data analysis model
by assuming the contiguous sequence of windows/regions are dependent to
each other as well as accounting for the dependency of observations inside a
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
| beta.df       |          | A data frame object that will contain the methylation data where rows are CpG sites and columns are the samples. |
| chr.pos       |          | A data frame object that will contain the chromosome number and position of the CpG sites.                           |
| nCpGs         |   100    | Number of CpG sites to be included per genomic region/window.                                                        |
| iter          |  20000   | Number of MCMC iterations to be run                                                                                  |
| burn          |   1000   | Number of burn-ins                                                                                                   |
| seed          |   1234   | A random number for reproducibility of results                                                                        |
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

Since in our paper, we recommend using the Bayesian Functional data
dependent model (BFDAM2) for detecting differentially methylated
regions, here we provide a minimum working example (MWE) for users to
implement BFDAM2. In this demo representation, we will restrict ourselves
to a lower number of iterations and burn-ins to save some time in
reproducing this MWE.

**A brief description of the example data**

The example data is a subset of the real data used in the paper and can
be found here. This example data contains the methylation values on
chromosome 12 spanning 10000 bases up and downstream from the popular
lung cancer gene KRAS.

**Guide for loading in required data**

For all the analysis shown below the user needs to input the following
example data: From the “Data” folder please load the example data
(“LUAD2015_ExampleData.RData”) into your R environment. All the results
and analysis shown below is based on this example data set. Note: Once
the example data is loaded into the R environment the object will be
denoted as a list of 2 elements with the displayed name “data”.

**Guide for loading in plotting function**

For plotting purposes (as shown below) the user needs to source in the
function “BFDAMs_PlotFn.R”. This R file can be found under the “R”
folder.

**Guide for loading in the main functions**

Both of our developed models BFDAM1 (“BFDA_independent_main.R”) and
BFDAM2 (“BFDA_dependent_main.R”) can be found under the “R” folder.
Since in this MWE we are only implementing the BFDAM2 model (dependent
method) the user needs to source in the file (“BFDA_dependent_main.R”)
into their R environment before running this demo analysis.

**Input Data (Snapshot of methylation beta-values data)**

``` r
beta.df = data$beta.df
beta.df[1:5,1:8]
```

    ##            GSM1632880_Tumor GSM1632881_Normal GSM1632882_Tumor GSM1632883_Tumor GSM1632884_Normal
    ## cg01389585       0.39451706         0.6009137       0.32448934       0.37254102        0.46949142
    ## cg12642725       0.87778511         0.8162177       0.67998402       0.85591416        0.79901909
    ## cg19205533       0.28720751         0.4551641       0.09996816       0.36298645        0.33971963
    ## cg08830758       0.15903439         0.1938489       0.08301977       0.21019356        0.13571027
    ## cg24979348       0.01993344         0.0152749       0.01365316       0.01624346        0.01438445
    ##            GSM1632885_Tumor GSM1632886_Tumor GSM1632887_Tumor
    ## cg01389585       0.43230172       0.51259288       0.48507750
    ## cg12642725       0.80612768       0.82314826       0.77993175
    ## cg19205533       0.29181200       0.38693604       0.22205678
    ## cg08830758       0.18929850       0.18164961       0.12048252
    ## cg24979348       0.01261574       0.01088544       0.01509158

**Input Data (Snapshot of CpG site location data)**

``` r
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

**Visualizing the fitted curves for genomic window/region 4 in
Chromosome 12 of example data**

To aid in visualizing the fitted functional mean methylation curves for
any given chromosome and genomic region/window we also provide the users
with the function “plot.fittedCurves()”. Below we demonstrate the use of
this function and plot the functional curves for the 2 groups. The blue
and red lines represent the fitted functional mean methylation curves
for healthy and cancer subjects respectively. Note that the y-axis
mentions “M-values” because our developed models logit transforms the
methylation beta-values to M-values before performing any downstream
statistical analysis.

``` r
print(plot.fittedCurves(beta.df = beta.df, chr.pos = chr.pos, 
                        nCpGs = 100, iter = 500, burn = 50, 
                        seed = 1234, control.group = 'Normal',
                        case.group = 'Tumor',chr=12, window.plot=4))
```

![alt text](https://github.com/schatterjee30/BFDAMs/blob/main/Images/fitting_plot.png)

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
output = BFDAM2[[1]]
output[, 1:5]
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

**Extracting required CpG sites for downstream biological analysis**

One of the key features of DMR analysis is to get a list of CpGs which
are differentially methylated in given region and then to further
investigate the charecteristics of these CpGs from a biological
perspective. If the user is interested in performing downstream
functional genomics analysis then the user can get the list of CpGs
which belong to a particular genomic window/region that was deemed as
DMR/Non-DMR by our function. For Ex. in our example data analysis
genomic window/region 4 was found to be a DMR by our dependent model
BFDAM2. If the user wants to extract the CpGs belonging to region 4 to
perform some downstream analysis then the user can do the following.

``` r
region.CpGs = output[output$chr==12 & output$window==4, ]
region.CpGs$CpGs
```

    ## [1] "'cg06539938', 'cg10104336', 'cg18456512', 'cg00697301', 'cg07238058', 'cg16436164', 'cg24035484', 'cg14184954', 'cg25401612', 'cg27546797', 'cg16316755', 'cg03763001', 'cg20255933', 'cg25756003', 'cg12486762', 'cg24827381', 'cg03815913', 'cg05405094', 'cg21427902', 'cg15735030', 'cg07274333', 'cg11931798', 'cg07424132', 'cg27099262', 'cg23032200', 'cg27302255', 'cg23655651', 'cg17460855', 'cg00063703', 'cg17067190', 'cg09060449', 'cg01047904', 'cg23367478', 'cg17535647', 'cg05660670', 'cg26387689', 'cg22927788', 'cg23752086', 'cg24975564', 'cg04012266', 'cg04101351', 'cg22355517', 'cg25912911', 'cg11416338', 'cg23015991', 'cg26571814', 'cg18639524', 'cg12136731', 'cg13446110', 'cg02631767', 'cg21913301', 'cg22999327', 'cg06026769', 'cg19264056', 'cg22056480', 'cg24336686', 'cg01153442', 'cg04724646', 'cg23371800', 'cg20994699', 'cg13063900', 'cg24509104', 'cg16169129', 'cg04091078', 'cg18109798', 'cg14525610', 'cg02769624', 'cg25885914', 'cg14628803', 'cg05717050', 'cg10143449', 'cg24755515', 'cg03355286', 'cg06947614', 'cg04250181', 'cg17572313', 'cg07396272', 'cg03895880', 'cg00995065', 'cg14997702', 'cg04579507', 'cg10832076', 'cg16923485', 'cg19659215', 'cg11704114', 'cg15583072', 'cg19610177', 'cg12127290', 'cg14191244', 'cg09860935', 'cg05987389', 'cg14622193', 'cg08551027', 'cg01055462', 'cg21962791', 'cg09738429', 'cg23515255', 'cg06157277', 'cg03618715', 'cg00912639'"

The above output gives the list of 100 CpGs that were found in
window/region 4 of chromosome 12 which were determined by our model to
be differentially methylated.

## Contributions

If you find small bugs, larger issues, or have suggestions, please email
the maintainers at <suvchat@iu.edu> or <shrabanti.chowdhury@mssm.edu>.
Contributions (via pull requests or otherwise) are welcome.
