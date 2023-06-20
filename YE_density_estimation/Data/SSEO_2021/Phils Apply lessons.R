# Construct a 5x6 matrix
X <- rnorm(30, nrow=5, ncol=6)

# Sum the values of each column with `apply()`
.....(X, 2, sum)

X<-matrix(rnorm(30),nrow=5,ncol=6)
X.apply<-apply(X,MARGIN=2,FUN=sum)	#Margin=2 columns, 1=rows
sum(X[,1])		#same as doing this..., summing each column

#The lapply() Function
#You want to apply a given function to every element of a list and obtain a list as a result. When you execute ?lapply, you see that the syntax looks like the apply() function.
#The difference is that:
#    It can be used for other objects like dataframes, lists or vectors; and
#    The output returned is a list (which explains the “l” in the function name), which has the same number of elements as the object passed to it.

A<-matrix(rnorm(30),nrow=5,ncol=6)
B<-matrix(runif(30),nrow=5,ncol=6)
C<-matrix(rbeta(30, shape1=0.5, shape2=0.5),nrow=5,ncol=6)

MyList<-list(A,B,C)

lapply(MyList,"[",,2)	#Grabs the second column from each matrix and makes
				#a new list, each with a vector of that column
lapply(MyList,"[",1,2)	#select a single element from each list, in this
				#case 1st row, second column...
lapply(MyList,"[",3,)	#grab the third row from each item in list

#The sapply() function works like lapply(), but it tries to simplify the output to the most elementary 
#data structure that is possible. And indeed, sapply() is a ‘wrapper’ function for lapply().
#An example may help to understand this: let’s say that you want to repeat the extraction operation 
#of a single element as in the last example, but now take the first element of the second row 
#(indexes 2 and 1) for each matrix.

#Applying the lapply() function would give us a list unless you pass simplify=FALSE as a 
#parameter to sapply(). Then, a list will be returned. See how it works in the code chunk below:

# Return a list with `lapply()`
lapply(MyList,"[", 2, 1 )

# Return a vector with `sapply()`
sapply(MyList,"[", 2, 1 )

# Return a list with `sapply()`
sapply(MyList,"[", 2, 1, simplify=F)

# Return a vector with `unlist()`
unlist(lapply(MyList,"[", 2, 1 ))

Z<-sapply(MyList,"[",1,1)
Z<-rep(Z,c(3,1,2))		#REPEATS each element the amount of times in the c vector...so 1st element 3x, 2nd element once, 3rd element twice...

#The mapply() Function
#The mapply() function stands for ‘multivariate’ apply. Its purpose is to be able to vectorize arguments 
#to a function that is not usually accepting vectors as arguments.
#In short, mapply() applies a Function to Multiple List or multiple Vector Arguments.
#Let’s look at a mapply() example where you create a 4 x 4 matrix with a call to the rep() function repeatedly:

Q1 <- matrix(nrow=4,ncol=4,c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4)))
Q2<-mapply(rep,1:4,4)		#vectorize the action of the function rep()
Q2

B_means<-apply(B,2,mean)	#mean of columns; knows column bc of MARGIN=2
mean(B[,2])

B_sdev<-apply(B,2,sd)

B_trans1<-sweep(B,2,B_means,"-")	#calcualtes the difference between the means (B_means) and original data points, B
#This means: “take the elements of the columns of the dataset MyPoints, and subtract the mean, 
#dataPoints_means, from each of them”.

B_trans2<-sweep(B_trans1,2,B_sdev,"/") #divide all those differences by the sd vector...

B_trans3<-sweep(sweep(B,2,apply(B,2,mean),"-"),2,B_sdev,"/") #same same all nested together

##the AGGREGATE () FUNCTION
library(stats)

DProgr<-seq(1,100, by=1)
Qty<-round(rnorm(100,mean=50,sd=10),0.1)				#?rnorm
Del<-as.logical(rbinom(100,size=1, prob=0.5))
DepPC<-as.character(rep(c("A","B","C","D"),25))

EX<-as.data.frame(cbind(DepPC,DProgr,Qty,Del))
str(EX)
EX$DProgr<-as.numeric(EX$DProgr)
EX$Qty<-as.integer(EX$Qty)
EX$DepPC<-as.factor(EX$DepPC)
sapply(EX, class)

#Here, you are interested in knowing where the product sells best in which department, for example. 
#that’s why you should regroup the data by department, summing up the sales, Qty, for each 
#department DepPC with the help of the aggregate() function:

aggregate(EX$Qty, by=EX["DepPC"], FUN=sum)

# Plot the sales per department
ggplot(data=aggregate(EX$Qty,by=EX["DepPC"],FUN=sum), aes(x=DepPC, y=x)) +
  geom_point()+
  ggtitle("Sales per department - All")







