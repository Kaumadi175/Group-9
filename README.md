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
