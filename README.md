# Group-9
Part 1
Q1.
Downloading gene_expression.tsv
```{r}
download.file("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv",destfile="gene_expression.tsv")
```

Downloading growth_data.csv
```{r}
download.file("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv",destfile="growth_data.csv")
```
#Load the data into R
```{r}
gene_expression <- read.table("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv")
```
#Set the gene identifiers as row names
```{r}
gene_expression <- read.table("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv",stringsAsFactors = FALSE, header=TRUE, row.names=1)
gene_expression
```
#Display the first six rows
```{r}
head(gene_expression)
```

Q2.
#Calculate the mean of the columns for each gene
#Show the first six rows
```{r}
gene_expression$Mean <- rowMeans(gene_expression)
head(gene_expression)
```

Q3
#Order the data frame by the mean column in descending order
#Extract the top 10 genes with the highest mean expression and display the top 10 genes
```{r}
gene_expression_sorted <- gene_expression[order(-gene_expression$Mean),]
top_10_genes <- head(gene_expression_sorted,10)
top_10_genes
```

Q4.
#Count the number of genes where the mean expression is less than 10 and display the results
```{r}
num_genes_below_10 <- sum(gene_expression$Mean < 10)
num_genes_below_10
```

Q5.
#Create a histogram of the mean values
```{r}
hist(gene_expression$Mean)
```

Q6.
#Download and load the dataset in R
Downloading growth_data.csv
```{r}
download.file("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv",destfile="growth_data.csv")
```
Following chunk is used to import the csv file to R object and display the column names.
#Import the CSV file and display the column names
```{r}
growth_data <- read.csv("growth_data.csv")
colnames(growth_data)
```

Q7.
#Check the circumference columns are numeric and if not, convert them

```{r}
growth_data$Circumf_2005_cm <- as.numeric(gsub("[^0-9.]","",growth_data$Circumf_2005_cm))
growth_data$Circumf_2020_cm <- as.numeric(gsub("[^0-9,]","",growth_data$Circumf_2020_cm))
```
#Subset data by site
```{r}
control_site <- subset(growth_data, Site == "northeast")
treatment_site <- subset(growth_data, Site == "southwest")
```

#Calculating mean and standard deviation for each site and year
#Display the summary table

```{r}
stats <- data.frame(
  Site = c("Control", "Control", "Treatment", "Treatment"),
  Year = c("2005", "2020", "2005", "2020"),
  Mean = c(mean(control_site$Circumf_2005_cm, na.rm = TRUE),
           mean(control_site$Circumf_2020_cm, na.rm = TRUE),
           mean(treatment_site$Circumf_2005_cm, na.rm = TRUE),
           mean(treatment_site$Circumf_2020_cm, na.rm = TRUE)),
  SD = c(sd(control_site$Circumf_2005_cm, na.rm = TRUE),
         sd(control_site$Circumf_2020_cm, na.rm = TRUE),
         sd(treatment_site$Circumf_2005_cm, na.rm = TRUE),
         sd(treatment_site$Circumf_2020_cm, na.rm = TRUE))
)
print(stats)
```

#Includes in variable (growth_data)
```{r}
growth_data
```

Q8.
Made box plot of tree circumference at the start and end of the study at control and treatment site.
#Create the boxplot

```{r}
boxplot(Circumference ~ Site * Time, data = growth_data_long, main = "Tree Circumference at Start and End of Study by Site", xlab = "Site and Time", ylab = "Circumference (cm)", col = c("lightblue", "lightgreen"), las = 2, # Makes x-axis labels vertical 
names = c("Control\nStart (2005)", "Control\nEnd (2020)", "Treatment\nStart (2005)", "Treatment\nEnd (2020)"))
```

Q9.
#Ensure the circumference columns are numeric

```{r}
growth_data$Circumf_2010_cm <- as.numeric(growth_data$Circumf_2010_cm)
growth_data$Circumf_2020_cm <- as.numeric(growth_data$Circumf_2020_cm)
```
#Calculated the growth over the last 10 years (2020-2010)
```{r}
growth_data$Growth_10yrs <- growth_data$Circumf_2020_cm - growth_data$Circumf_2010_cm
```
#Calculated the mean growth at each site and display the result
```{r}
mean_growth <- aggregate(Growth_10yrs ~ Site, data = growth_data, FUN = mean, na.rm = TRUE)
print(mean_growth)
```

Q10.
#Perform a two sample t-test to compare growth between sites and display the t-test result
```{r}
t_test_result <- t.test(Growth_10yrs ~ Site, data = growth_data)
print(t_test_result)
```

Assignment 4_Part2

Organism - E.coli

## downloading the E coli gene
```{r, Download_Seq}
library("R.utils")

URL="http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
download.file(URL,destfile="ecoli_cds.fa.gz")
gunzip("ecoli_cds.fa.gz")
list.files()
```

organism - Cryobacterium sp. SO1 (GCA_004210215)

## downloading the Cryobacterium gene
```{r}
library("R.utils")

URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_98_collection/cryobacterium_sp_so1_gca_004210215/cds/Cryobacterium_sp_so1_gca_004210215.ASM421021v1.cds.all.fa.gz"
download.file(URL,destfile="Cryobacterium_cds.fa.gz")
gunzip("Cryobacterium_cds.fa.gz")
list.files()
```

Organism - E coli


```{r}
library("seqinr")
```

#loading the sequence for E coli
```{r}
cds <- seqinr::read.fasta("ecoli_cds.fa")
str(head(cds))
```

#Number of coding sequences
```{r}
length(cds)
```

#Total length of sequence
```{r}
head(summary(cds))
len <- as.numeric(summary(cds)[,1])
```

```{r}
sum(len)
```

#mean and median coding sequence length of E coli
```{r}
mean(len)
```

```{r}
median(len)
```

#boxplot of coding sequence length in E coli
```{r}
boxplot(len,ylab="sequence length (bp)")
```
