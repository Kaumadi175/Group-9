# Group-9
Part 1
Q1.
#Define the URL for the raw data
#Download the file
```{r}
url <- "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv"
download.file(url, destfile = "gene_expression.tsv")
```

#Read data and make identifiers as row names
#Display the first six rows
```{r}
gene_data <- read.delim("gene_expression.tsv", header = TRUE, row.names = 1)
head(gene_data)
```

Q2.
#Calculate the row wise mean of the sample
#Display the first six rows with the new column
```{r}
gene_data$Mean_Expression <- rowMeans(gene_data)
head(gene_data)
```

Q3
```{r}
sorted_data <- gene_data[order(gene_data$Mean_Expression, decreasing = TRUE), ]
top_10_genes <- head(sorted_data, 10)
top_10_genes
```

Q4.
#Count the number of genes where the mean expression is less than 10 and display the results
#Filter genes with mean expression < 10 and display the result
```{r}
num_genes_below_10 <- sum(gene_data$Mean_Expression < 10)
num_genes_below_10
```

Q5.
#Histogram plot of the mean value
```{r}
hist(
  gene_data$Mean_Expression,
  breaks = 20,                      # Adjust the number of bins
  col = "skyblue",                  # Set bar color
  border = "black",                 # Set border color
  main = "Histogram of Mean Gene Expression",  # Title
  xlab = "Mean Expression Value",   # X-axis label
  ylab = "Frequency"                # Y-axis label
)
```

Q6.
#Define the URL for the raw CSV data
#Read the CSV file into a R object and display the column names
```{r}
url <- "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv"
growth_data <- read.csv(url, header = TRUE)
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
#Sequence composition
```{r}
GC(cds[[1]])
```
[1] 0.5151515

#number of each nucleotide present in E coli
```{r}
count(cds[[1]],1)
count(cds[[1]],2)
count(cds[[1]],3)
```
#Unlist first three sequences from the cds list
```{r}
length(unlist(cds[1:3]))
```
```{r}
GC(unlist(cds[1:3]))
```

#bar plots for nucleotide in E coli
```{r}
dna <- unlist(cds)
dna_composition <- count(dna,1)
```

```{r}
barplot(dna_composition,xlab="nucleotides",ylab="frequency", main="E coli CDS composition")
```

```{r}
translate(cds[[1]])
prot <- lapply(cds, translate)
```

```{r}
aa <- unique(prot[[2]])
aa <- aa[aa != "*"]
length(aa)
```

```{r}
count(prot[[1]],wordsize=1,alphabet=aa)
count(prot[[2]],wordsize=1,alphabet=aa)
```

#bar plots for amino acids in E coli
```{r}
aminoacids <- unlist(cds)
aminoacids_composition <- count(aminoacids,1)
barplot(aminoacids_composition,xlab="aminoacids",ylab="frequency", main="E coli CDS composition")
```

#codon usage for a single protein sequence
```{r}
uco(cds[[2]])

uco(cds[[2]],index="rscu")

uco(cds[[2]],index="rscu",as.data.frame=TRUE)
```

#Sequences over- or under-represented in an Ecoli
```{r}
prots <- unlist(prot)

mycounts <- count(prots,wordsize=3,alphabet=aa)

str(mycounts)

head(mycounts)

myfreq <- count(prots,wordsize=3,alphabet=aa,freq=TRUE)
```



#Cryobacterium sp organism 
```{r}
cds <- seqinr::read.fasta("Cryobacterium_cds.fa")
str(head(cds))

length(cds)

head(summary(cds))
len <- as.numeric(summary(cds)[,1])
sum(len)
mean(len)
median(len)
boxplot(len,ylab="sequence length (bp)")

GC(cds[[1]])

count(cds[[1]],1)

count(cds[[1]],2)

count(cds[[1]],3)

length(unlist(cds[1:3]))

dna <- unlist(cds)

dna_composition <- count(dna,1)
barplot(dna_composition,xlab="nucleotides",ylab="frequency", main="Cryobacterium CDS composition")

translate(cds[[1]])
prot <- lapply(cds, translate)

aa <- unique(prot[[2]])
aa <- aa[aa != "*"]
length(aa)

count(prot[[1]],wordsize=1,alphabet=aa)

count(prot[[2]],wordsize=1,alphabet=aa)

aminoacids <- unlist(cds)
aminoacids_composition <- count(aminoacids,1)
barplot(aminoacids_composition,xlab="aminoacids",ylab="frequency", main="Cryobacterium CDS composition")

uco(cds[[2]])

uco(cds[[2]],index="rscu")

uco(cds[[2]],index="rscu",as.data.frame=TRUE)

prots <- unlist(prot)
mycounts <- count(prots,wordsize=3,alphabet=aa)

str(mycounts)

head(mycounts)

myfreq <- count(prots,wordsize=3,alphabet=aa,freq=TRUE)
```
