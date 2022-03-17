import sys
import numpy

"""The following uses Python to challenge you to create an algorithm for finding
matches between a set of aligned strings. Minimal familiarity with Python is 
necessary, notably list and Numpy array slicing. 
"""

"""Problem 1.

Let X be a list of M binary strings (over the alphabet { 0, 1 }) each of length 
N. 

For integer 0<=i<=N we define an ith prefix sort as a lexicographic sort 
(here 0 precedes 1) of the set of ith prefixes: { x[:i] | x in X }.
Similarly an ith reverse prefix sort is a lexicographic sort of the set of
ith prefixes after each prefix is reversed.

Let A be an Mx(N+1) matrix such that for all 0<=i<M, 0<=j<=N, A[i,j] is the 
index in X of the ith string ordered by jth reverse prefix. To break ties 
(equal prefixes) the ordering of the strings in X is used. 

Complete code for the following function that computes A for a given X.

Here X is a Python list of Python strings. 
To represent A we use a 2D Numpy integer array.

Example:

>>> X = getRandomX() #This is in the challenge1UnitTest.py file
>>> X
['110', '000', '001', '010', '100', '001', '100'] #Binary strings, M=7 and N=3
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>> 

Hint:
Column j (0 < j <= N) of the matrix can be constructed from column j-1 and the 
symbol in each sequence at index j-1.  

Question 1: In terms of M and N what is the asymptotic cost of your algorithm?
    O(MN)
"""

def constructReversePrefixSortMatrix(X):
    """ Referred to the code presented by Prof. Paten, I've written detailed notes and made revisions to demonstrate my understanding of the code """
    # Creates the Mx(N+1) matrix
    A = numpy.empty(shape=[len(X), 1 if len(X) == 0 else len(X[0])+1 ], dtype=int)
    # Code to write - you're free to define extra functions
    # (inline or outside of this function) if you like.
    A[:, 0] = list(range(0, len(X)))  # sets the first column of the A array as the range of number from 0 to
    if len(X) == 0:  # sets the length with this if/else loop
        length = 0
    else:
        length = len(X[0])
    for i in range(length):  # goes in a range of the length we determined earlier
        first = []; second = []
        for j in A[:, i]:  # goes till the end of the column
            if X[j][i] == "0":  # checks to see if the value at the given index is a 0
                first.append(j)  # appends the index if a 0 is identified
            else:
                second.append(j)  # appends the index if a 1 is identified
        currRange = first + second  # sets the range/column
        A[:, i+1] = currRange  # sets the column to the determined range
    return A

"""Problem 2:

Following on from the previous problem, let Y be the MxN matrix such that for
all 0 <= i < M, 0 <= j < N, Y[i,j] = X[A[i,j]][j].

Complete the following to construct Y for X.

Hint: You can either use your solution to constructReversePrefixSortMatrix()
or adapt the code from that algorithm to create Y without using
constructReversePrefixSortMatrix().

Question 2: In terms of M and N what is the asymptotic cost of your algorithm?
    O(MN)
"""
def constructYFromX(X):
    #Creates the MxN matrix
    Y = numpy.empty(shape=[len(X), 0 if len(X) == 0 else len(X[0]) ], dtype=int)
    #Code to write - you're free to define extra functions
    #(inline or outside of this function) if you like.
    A = constructReversePrefixSortMatrix(X)  # saves the
    for i in range(0, len(X)):  # iterates through the number of values in X
        for j in range(0, len(X[0])):  # iterates through the length of a value within X
            Y[i, j] = X[A[i, j]][j]  # uses the X and A to set the values within Y
    return Y

"""Problem 3.

Y is a transformation of X. Complete the following to construct X from Y,
returning X as a list of strings as defined in problem 1.
Hint: This is the inverse of X to Y, but the code may look very similar.

Question 3a: In terms of M and N what is the asymptotic cost of your algorithm?
    O(max(M * N * long matches)
    
Question 3b: What could you use the transformation of Y for?
Hint: consider the BWT.
    The transformation of Y can be used to determine the the A which can then be used to recreate the X

Question 3c: Can you come up with a more efficient data structure for storing Y?
    Use of compression algorithms can help with more efficient data construction for storing Y. 
    Compression algorithms, such as the huffman algorithm, can be a more efficient data structure for storing Y. 
"""
def constructXFromY(Y):
    """ Referred to the code presented by Prof. Paten, I've written detailed notes and made revisions to demonstrate my understanding of the code """
    # This problem attempts to reconstruct X using the Y we calculated in the previous problem.
    # The best approach, and the approach used in this problem is working backwArds.
    # The problem recreates the values within "A" and then is used to create X which is an array
    # The X array is then formatted to match the X we are trying to produce
    # Creates the MxN matrix
    x, y = Y.shape
    X = numpy.empty(shape=[x, 0 if x == 0 else y], dtype=int)
    # Code to write - you're free to define extra functions
    # (inline or outside of this function) if you like.
    currRange = range(0, x)  # sets the range to from 0 to the number of columns
    for i in range(y):  # iterates through all the columns
        for j in range(x):  # iterates through all the rows
            X[currRange[j], i] = Y[j, i]  # sets the value in X based off the array we created in the second problem
        first = []; second = []  # initializes lists to hold the indexes containing 0's and 1's
        for k in range(x):  # iterates through all the rows
            if Y[k, i] == 0:  # checks to see the value in that index is 0
                first.append(currRange[k])  # appends the index to zeros index
            else:
                second.append(currRange[k])  # appends the index to the ones index
        currRange = first + second  # adds the indexes together to determine the final range
    formattedX = list(map(lambda i: "".join(map(str, i)), X))  # formats the X array into the appropriate manner
    return formattedX

"""Problem 4.

Define the common suffix of two strings to be the maximum length suffix shared
by both strings, e.g. for "10110" and "10010" the common suffix is "10" because
both end with "10" but not both "110" or both "010".

Let D be a Mx(N+1) Numpy integer array such that for all 1<=i<M, 1<=j<=N,
D[i,j] is the length of the common suffix between the substrings X[A[i,j]][:j]
and X[A[i-1,j]][:j].

Complete code for the following function that computes D for a given A.

Example:

>>> X = getRandomX()
>>> X
['110', '000', '001', '010', '100', '001', '100']
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>> D = constructCommonSuffixMatrix(A, X)
>>> D
array([[0, 0, 0, 0],
       [0, 1, 2, 2],
       [0, 1, 2, 3],
       [0, 1, 1, 1],
       [0, 0, 2, 2],
       [0, 1, 0, 0],
       [0, 1, 1, 3]])

Hints:

As before, column j (0 < j <= N) of the matrix can be constructed from column j-1
and the symbol in each sequence at index j-1.

For an efficient algorithm consider that the length of the common suffix
between X[A[i,j]][:j] and X[A[i-k,j]][:j], for all 0<k<=i is
min(D[i-k+1,j], D[i-k+2,j], ..., D[i,j]).

Question 4: In terms of M and N what is the asymptotic cost of your algorithm?
    O(MN)
"""
def constructCommonSuffixMatrix(A, X):
    """ Referred to the code presented by Prof. Paten, I've written detailed notes and made revisions to show my understanding of the code """
    # Creates the Mx(N+1) D matrix
    D = numpy.zeros(shape=A.shape, dtype=int)
    x, y = A.shape
    # Code to write - you're free to define extra functions
    # (inline or outside of this function) if you like.
    D[:, 0] = 0 ** y  # sets the first column to all 0's
    if x == 0:  # this if/else loop determines the length for the first for loop
        length = 0
    else:
        length = len(X[0])
    for i in range(length):  # goes in range of the length we determined
        first = []; second = []  # initializes the lists to get the indexes of 0's and 1's
        p = 0; q = 0
        for j in range(0, x):  # goes in range of the number the number of rows
            p = min(p, D[j, i]+1)  # saves the min value between the max and index val as p
            q = min(q, D[j, i]+1)  # saves the min value between the max and index val as q
            index = X[A[j,i]][i]  # saves the index within X
            if index == "0":  # checks to see if the value at the given index contains a 0
                first.append(p)  # appends the suffix length to the zeros list
                p = sys.maxsize  # resets the value of p as the maxsize, which means that it fetches the largest value a variable of a data type can store
            elif index == "1":  # checks to see if the value at the given index contains a 0
                second.append(q)  # appends the suffix length to the ones list
                q = sys.maxsize  # resets the value of q as the maxsize, which means that it fetches the largest value a variable of a data type can store
        currRange = first + second  # updates the range
        D[:, i+1] = currRange  # sets the values in the column to the range
    return D

"""Problem 5.

For a pair of strings X[x], X[y], a long match ending at j is a common substring
of X[x] and X[y] that ends at j (so that X[x][j] != X[y][j] or j == N) that is longer
than a threshold 'minLength'. E.g. for strings "0010100" and "1110111" and length
threshold 2 (or 3) there is a long match "101" ending at 5.

The following algorithm enumerates for all long matches between all substrings of
X, except for simplicity those long matches that are not terminated at
the end of the strings.

Question 5a: What is the asymptotic cost of the algorithm in terms of M, N and the
number of long matches?
    O(MN*longMatches)
    
Question 5b: Can you see any major time efficiencies that could be gained by
refactoring?
    Refactoring the code to have less references and iterations through b and c can help with time efficiencies. 
    
Question 5c: Can you see any major space efficiencies that could be gained by
refactoring?
    If you avoid saving the arrays into separate arrays, you can help with saving space within the code. Rather, it 
    would be space efficient to refer to it, or make references to it. 
    
Question 5d: Can you imagine alternative algorithms to compute such matches?,
if so, what would be the asymptotic cost and space usage?
    I cannot think of another algorithm so compute such matches. However, in the case that such an algorithm exists, 
    I would expect the asymptotic cost to be O(MN) and the space usage to be less than that of the existing algorithm. 
"""
def getLongMatches(X, minLength):
    assert minLength > 0

    A = constructReversePrefixSortMatrix(X)
    D = constructCommonSuffixMatrix(A, X)

    #For each column, in ascending order of column index
    for j in range(1, 0 if len(X) == 0 else len(X[0])):
        #Working arrays used to store indices of strings containing long matches
        #b is an array of strings that have a '0' at position j
        #c is an array of strings that have a '1' at position j
        #When reporting long matches we'll report all pairs of indices in b X c,
        #as these are the long matches that end at j.
        b, c = [], []

        #Iterate over the aligned symbols in column j in reverse prefix order
        for i in range(len(X)):
            #For each string in the order check if there is a long match between
            #it and the previous string.
            #If there isn't a long match then this implies that there can
            #be no long matches ending at j between sequences indices in A[:i,j]
            #and sequence indices in A[i:,j], thus we report all long matches
            #found so far and empty the arrays storing long matches.
            if D[i,j] < minLength:
                for x in b:
                    for y in c:
                        #The yield keyword converts the function into a
                        #generator - alternatively we could just to append to
                        #a list and return the list

                        #We return the match as tuple of two sequence
                        #indices (ordered by order in X) and coordinate at which
                        #the match ends
                        yield (x, y, j) if x < y else (y, x, j)
                b, c = [], []

            #Partition the sequences by if they have '0' or '1' at position j.
            if X[A[i,j]][j] == '0':
                b.append(A[i,j])
            else:
                c.append(A[i,j])

        #Report any leftover long matches for the column
        for x in b:
            for y in c:
                yield (x, y, j) if x < y else (y, x, j)