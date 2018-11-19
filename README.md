#####################################################################
#  Data Exploration in  Ecotoxicogenomics
#  - how to work with toxiconomical data? 
#####################################################################

#Purpose of my report is to give an introduction to the field of toxico- and ecotoxicogenomics, 
#to explore what kind of data used in this field and how the data can be accessed.

#The initial goal of this script was to explore 
##what kind of data is available in different databases
##how this data can be imported and read
##how the data can be used by using entry level bioinformatic and statistical tools in R.

#First I worked with the database CEBS and look into how you can access datafiles directly from the websites using URL 
#instead of downloading files to your working directory.
#However during the process I found studying the microarrays interesting and that is why I decided to focus on learning more
#about the microarrays and CEL-files and the codes. 


#####################################################################
##### CEBS database - how to import data from a new database?  ######
#####################################################################

#Getting data straight from the database from the web using data.table package and 
#function fread() 
#CEBS file orientation:
#KEYS FILE = has definitions for all the abbreviations used for assay name 
#column headers in the zip and sample files.
#SAMPLE FILE = has 20 records in a tab delineated text file and gives a preview 
#of the raw data contained in the zip file. 
#ZIP FILE = has all the raw data as a tab delineated text file

??data.table
require(data.table)
library(data.table)

CEBS.sample <- fread("ftp://anonftp.niehs.nih.gov/ntp-cebs/datatype/REPRODUCTIVE/NTP/REPRODUCTIVE_NTP_SAMPLE.txt")
head(CEBS.sample)
CEBS.keys <- fread("ftp://anonftp.niehs.nih.gov/ntp-cebs/datatype/REPRODUCTIVE/NTP/keys.txt")
View(CEBS.keys)

immunology.data <- fread("ftp://anonftp.niehs.nih.gov/ntp-cebs/datatype/IMMUNOLOGY/NTP/IMMUNOLOGY_NTP_SAMPLE.txt")
View(immunology.data)
immunology.keys <- fread("ftp://anonftp.niehs.nih.gov/ntp-cebs/datatype/IMMUNOLOGY/NTP/keys.txt")
View(immunology.keys)



#####################################################################
######   TG-GATE Open Database - microarray data exploration  #######
#####################################################################
setwd("C:/R files/DATA/celfiles")

#Open TG-GATEs is a public toxicogenomics database. It consists of 170 compounds.
#TG-GATE commonly stores the data using Affymetrix GeneChips and stores data in CEL-files. 
#Affymetrix is a product of Thermo Fisher Scientific company and makes quartz chips for analysis of DNA Microarrays called GeneChip® arrays. 
#Affymetrix's GeneChip arrays scan for the presence of particular genes in a biological sample. Within this area, 
#Affymetrix is focused on oligonucleotide microarrays.  

#Data from TG-GATE: (https://toxico.nibiohn.go.jp/open-tggates/english/compound_search/screen3/compound?compound_id=00137&design_name=Rat230_2&organ_id=ORGA0010&test_type=in+vitro)
#substance=ethanol 
#eperiment= rat in vitro 2h, 8h and 24h; control, low, middle and high dose groups, two samples of each group

#NOTE! For the full run of the script download the cel.files from the link above or from github: https://github.com/hakkinenmiia/project1


### IMPORT .CEL FILE DATA ###
 
library(affy)

#specify the path on your computer where the folder that contains the CEL-files is located
celpath = "C:/R files/DATA/celfiles"

#import CEL files containing raw probe-level data into an R AffyBatch object
data <- ReadAffy(celfile.path=celpath)

dim(data)
class(data)

#How to retrieve intensities from affy?
#These methods extract the intensities of all probes (both PM=perfect match and MM=mismatch probes) from the AffyBatch
#Extract the matrix of expression values, which is in some sense 'the contents of a .cel file'.
expr = exprs(data)
#Function exprs= access the expression and error measurements of assay data stored in an object
#exprs returns a (usually large) matrix of expression values
#you can also use  
int = intensity(data)
#Both exprs and intensity do the same thing (note you got identical results!). 
#They extract the matrix of expression values, which is in some sense 'the contents of a .cel file'.


#To retrieve a matrix with PM (perfect match) probe intensities use the pm() method, 
#e.g. the command below will retrieve the PM intensities of the first 5 rows of the data set.
pm(data)[1:5,] 
#you are reading multiple cel-files simultaniously so you find file names (which are basically different samples)
#as column names and probe numbers as row names. Note that probe name is not same as probe ID.
#pm() reorders the rows while exprs() and intensity() retrieve the intensities of the probes in the original order (as they occur in the CEL files)

#Apart from the expression data itself, microarray data sets need to include information about the samples that were hybridized to the arrays, e.g. sample labels. 
#AffyBatch objects have several characteristics. One of them is called phenoData
ph= data@phenoData
ph
ph$sample
ph@data
ph@varMetadata

#Microarray data sets should include information on the probes. AffyBatches have a slot called featureData, 
#a data frame that contains labels for the probes.
feat=data@featureData
feat #An object of class 'AnnotatedDataFrame': none <- as we can see, there is no information about the probe annotation
feat@data

#Microarray data sets should also include information on the experiment.
exp= data@experimentData
exp #as we can see, there is no experiment annotations, this is not uncommon on datasets from public data repositories. 

#retrieve name of the CDT file associated with the array data
cdfName(data)

#retrieve names of the probe set IDs on the arrays:
featureNames(data)
length(featureNames(data)) #number of probe sets on the arrays 

pm(data)[1:2,]

ph@data$sample[1:24]

### QUALITY CONTROL ### 
#Microarray data requires quality control. Micro arrays contain a lot of noise due different amounts of RNA synthesized 
#imperfections on the array surface and synthesis of the probe and from differences in hybridization conditions.
#This noise should be removed in order to make biologically meaningfull conclusions about the data.
#In this script I only try out three of the steps of the quality control:
#add annotation to phenodata, create a picture of micro array data and normalization of the data.

#Adding annotation to pheno data and give samples informative names:
ph@data[ ,1] = c("control1_2h","control2_2h","low1_2h","low2_2h","mid1_2h","mid2_2h","high1_2h","high2_2h",
                 "control1_8h","control2_8h","low1_8h","low2_8h","mid1_8h","mid2_8h","high1_8h","high2_8h",
                 "control1_24h","control2_24h","low1_24h","low2_24h","mid1_24h","mid2_24h","high1_24h","high2_24h")
ph@data

#Create a micro array picture:
#this code creates individual jpg. pictures of each microarray to your working directory folder 
for (i in 1:24)
{
  name = paste("image",i,".jpg",sep="")
  jpeg(name)
  image(data[,i],main=ph@data$sample[i])
  dev.off()
}

#you can also create collective image of all the arrays. 
#Because we have 24 microarrays we use 4x6 dimentions 
op = par(mfrow = c(4,6))
for (i in 1:24){image(data[,i],main=ph@data$sample[i])}

#Normalization of the data:
#The standard method for normalization is RMA. 
#RMA is one of the few normalization methods that only uses the PM probes

data.rma <- rma(data)
data.matrix <- exprs(data.rma)

head(data.matrix)
#(row names are probe IDs and columns different samples)

#What does Background correction do for the dataset? 
#Background correnction to correct for spatial variation within individual arrays, 
#log transformation to improve the distribution of the data, quantile normalization 
#to correct for variation between the arrays and probe normalization to correct for 
#variation within probe sets.

### Further steps ###
#I did not perform the entire process of the quality control of the data, merely tried out
#some of the codes and aspects of it. 
#To study the data further would include reading and retrieving intensities from different probes,
#create visual plots and histograms to assess the quality of the data and compare raw and normalized data.
#However any meaningful comparison of the intensities would require metadata of the probe IDs to see which genes
#we are comparing. Doing this blindsided (without IDs) would not provide any meaningful results.
#Unfortunately this IDs were not available nor to be found online for this dataset. 

#Sources:
#Datatable: https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html 
#Description of TG-GATE and CEL file data: https://dbarchive.biosciencedbc.jp/en/open-tggates/data-2.html
#What is affymetrix genechip microarray:  data https://www.vsni.co.uk/software/genstat/htmlhelp/marray/AffymetrixChips.htm 
#Open microarray data files: https://www.vsni.co.uk/software/genstat/htmlhelp/marray/OpenMicroarrayDataFiles.htm
#How to read CEL.files: https://wiki.bits.vib.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor#Looking_at_the_data
#                       https://wiki.bits.vib.be/index.php/How_to_retrieve_intensities_using_affy
#Normalization of the data: https://wiki.bits.vib.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor#Looking_at_the_data
#Microarray picture: https://wiki.bits.vib.be/index.php/How_to_create_microarray_pictures


