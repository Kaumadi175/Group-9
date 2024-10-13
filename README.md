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
