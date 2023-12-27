#!/usr/bin/bash  
### The line above is where we indicate which program will be used for the script
### Notice that comments are denoted with the '#' symbol

### This is a function declared using 'function':
### Function to calculate nth digit of triangular number series, as well as printing 
function triangular(){
    tri=0
    i=1
    while [ $i -le $1 ]     ### Here's an example of a while loop
    do
	tri=$((tri+i))     ### Any arithmetic expressions should be enclosed in double parentheses, or evaulated using let or expr
	echo $tri >> $2    ### We can use I/O redirection to store our values into a file
	i=$((i+1))
    done
    echo $tri
}

### And this is a function declared without 'function':
### Function to calculate nth digit of Fibonacci sequence
fibonacci(){
    a=1
    b=1
    fn=$b
    echo $a >> $2
    for ((i=1;i<$1;i++))    ### Here's a C-style for loop
    do
	echo $b >> $2
	fn=$((a+b))
	a=$b
	b=$fn
    done
    echo $a
}


### Actual 
if [ $# -lt 3 ] ### Always good housekeeping to make sure if we require 3 arguments, we get 3 arguments (this is a minimum check)
then
    echo "You need to provide 3 arguments: 'name', 'N', 'output'"
    exit 1
fi
for i in $@   ### And here's a sample-style for loop
do
    echo $i
done
rm -f $3

if [ $1 == "Fibonacci" ] || [ $1 == "fibonacci" ]  ### Here's a conditional statement
then
    echo Fibonacci
    fibonacci $2 $3 ### Notice that the global positional variables in our script are different in number compared to those in functions
else
    echo Triangular
    triangular $2 $3
fi