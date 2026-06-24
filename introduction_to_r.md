---
title: Concepts and Techniques
subtitle: "Introduction to R"
layout: home
nav_order: 2
description: "Intensive bioinformatics course covering RNA-seq, ATAC-seq, CUT&RUN, and single-cell genomics."
---

# INTRODUCTION TO R

Throughout the Genomics of Gene Regulation course, we will be working with R using the RStudio environment. This R Markdown file serves as a brief introduction for students who are not familiar with R and as a refresher for students who have used this programming language before. R is a programming language as well as a statistical and graphics package. The ability to use R for statistical computing combined with its modularity for generating informative, compelling, and publication-ready plots/graphics make this programming language widely used in the fields of genomics and computational biology. The data bootcamps in this course will take you through the use of many R functions, explaining each step as we analyze the genomics datasets before us.

R has a standard interface which you can feel free to use. However, for this course we will be scripting/coding in R using the RStudio interface. RStudio has several advantages over the base R interface. These include syntax highlighting (e.g., becomes useful when there is a mistake in your line of code or when you are trying to track how many brackets you end to use to ), multiple windows for plots and R objects, scripts and command history.

# BECOMING FAMILIAR WITH RSTUDIO

R starts automatically when you open RStudio. The bottom left window is the console, where R runs. When using the console you have to run each command one-by-one. This can become confusing if you have to create a long workflow here you need to go back to previous commands that you ran before or save the commands so that you can share them with other researchers. For this reason we use R scrips.

To begin we are going to open a script file in RStudio by clicking File -> New File-> R Script. The script is just a set of commands that you can save in a file and open later to know exactly how, for example, you analyzed a dataset. Once you open the script, you can edit it in the script window in RStudio. There are different versions of script files in R. For example this text is being written in R Markdown which is essentially a script file that allows for better annotation of code and separation of code chunks to make following a workflow easier. To create an R Markdown file click File -> New File-> R Markdown. We encourage you to use R Markdown files as much as possible in your research as they are easy to follow and can be used to generate high quality reports that can be shared with other researchers. This makes it easier for your peers to reproduce your analysis if the need to. To learn more about R Markdown follow this link [Rmarkdown](https://rmarkdown.rstudio.com/lesson-1.html).

In the command line, you can type your command after the ">" symbol and then hit enter. When working with R Markdown you can type the command in the code chunk and run the chunk by clicking the run button on the top right of th chunk. If you are working on an R script you can run a command by selecting it(or just clicking on the line) and clicking run on the top right of the script window.

```r
#in the code chunk you can comment by using the "#" symbol
#anything that starts with a "#" symbol will not be executed

x=5 #here i am creating a new variable x that is equal to 5
x+3 # doing a simple calculation
```

# R DATA TYPES

R has a lot of different data types, many of which we will encounter in this course. Two of the most basic data types are vectors and data frames. These are fundamental to the logic with which R operates.

## Vectors

A vector is a collection of numbers, characters, logical statemtns, etc. You can create a vector in R by using the combine function **c()**. This function will put all objects inside the parantheses in a vector.

```r
#below we are creating three different classes of vectors
a = c(1, 4, 0, -5) #numeric vector
a
b = c("five", "nine", "six") #character vector
b
c = c(TRUE, FALSE, T, F) # logical vector
c
```

You can do operations with vectors in R.
```r
d = a*2
d
e = a+d
e
length(e) #the length function gives you the number of elements, or objects, in a vector
d[2] #here i am extracting the second element of the vector d
f = c(a,d) # here I am using the combine function to combine two vectors into a new vector
f
```

## Data Frames

A data frame is another important data type in R. It is similar to a matrix but it has name that identify the columns. In a data frame a column can be numeric, alphanumeric, or a logical indicator (TRUE or FALSE). Let's create a data frame and explore some operations.

```r
gene = c("TFAP2A", "GAPDH", "HK1", "MSX1")
normalized_expression = c(0.5, 6, 4, 0.8)
type = c("Transcription Factor", "Housekeeping", "Housekeeping", "Transcription Factor")

mydata = data.frame(gene, normalized_expression, type)
mydata

names(mydata) #this command will give you the column names of your data frame
dim(mydata) #this command will give you the dimensions of your data frame in rows x columns
head(mydata) #if you had a very long data frame, the head function will show you the first lines
tail(mydata) #if you had a very long data frame, the head function will show you the last lines
head(mydata, n=2) #if you add the n argument, you can specify how many lines you want to show
```

Let's say you want to know the classes of each of the columns in your data frame. For this you would use the "str()" function.

```r
str(mydata)
```

Notice how the "gene" and "type" columns are character vectors, whereas the "normalized_expression" column is a numeric vector. This is actually very important to check every time you are working with a data frame. This is because sometimes R will not execute a function if the class of the object is not what the function expects.

If you are running this using the R Markdown file, take a moment to look at the top right tab. Notice how all the vectors, variables, and data frames that you specified are saved there. Now to move on to some functions. Let's say that I want to know the mean expression level of the 4 genes in my imaginary cell. Then I want to know the mean expression of the housekeeping genes. This brings us to an important operand, the dollar sign "$". You can use the dollar sign to subset specific columns from your data frame. It is very useful.

```r
mydata$gene
mydata$normalized_expression

f = mydata$type #notice how when i save a column to a variable, it is saved as a vector
class(f)
```

OK, I want to find the mean expression of the 4 genes.

```r
mean(mydata$normalized_expression)
sd(mydata$normalized_expression) #this function gives you the standard deviation
```

To calculate the mean expression of the housekeeping genes we will use the subset function. The subset function allows you to subset a data frame based on a specification.

```r
housekeeping_data_frame = subset(mydata, type == "Housekeeping")
housekeeping_data_frame
```

Notice the structure of the function. I essentially told R to subset the data frame and only return the lines where the "type" column is equal to "Housekeeping". The double equal sign "==" is a logical operand and it is saying that type has to be exactly equal to "Housekeeping". In this case, if you wrote "housekeeping" with a lowercase "h", then the command would give you an error.

```r
mean(housekeeping_data_frame$normalized_expression)
```

Say that I want to subset based on a numeric specification. In the case below the function will be returning lines where the "normalized_expression" value is less than 1.

```r
subset(mydata, normalized_expression < 1)
subset(mydata, normalized_expression <= 0.5) #"<=" is less than or equal to
```

As functions in R become more complicated and involve more arguments, it can be difficult to remember what each argument is specifying. One useful trick is to use the "?" symbol followed by the name of the function. When you do this, a description of the function as well as each argumnent will appear in the bottom right panel of RStudio.

```r
?subset
```

If you are not sure about the exact name of the function but know it roughly, you can use the double question mark symbol "??". This will output a list of possible functions, vignettes, and help pages where you can look for your specific function.

Finally, one last import aspect of using R is knowing the working directory where you are working from. This becomes especially important when you have to save your script or when you have to import a data frame. To get the working directory where you are currently in use the "getwd()" function. To set the worlking directory use the "setwd()" function.

```r
getwd() #this gives me the current working directory, i want to change it
setwd("/Users/fjodormerkuri/Desktop/brazil_course/")
```

I changed the working directory to a folder that I created for the course on my desktop. Now if I save a script it will be automatically saved to this folder. Also I try to upload something, R will look specifically in this folder. If my file is not found there, then it will give me an error. whenever you change the working directory on R Markdown, you are not changing it for the R environment. This is not something we'll get into but just know that whenever you change teh working directory on R Markdown, also copy and paste the command on the console and run it there too.

If you are ever stuck with something, here are some links to a few nice R resources. Often just googling the problem or error will lead to a solution as it is likely that a lot of other people had the same problem.

[Quick R](http://www.statmethods.net/) - nice wiki that overs a bunch of base R functions

[Data Camp](https://www.datacamp.com/courses/introduction-to-r) - online interactive R tutorial

[R reference Card from CRAN](http://cran.r-project.org/doc/contrib/Short-refcard.pdf) 





