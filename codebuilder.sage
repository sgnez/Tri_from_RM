import copy
################################################################################
#                  INSTRUCTIONS FOR MODIFICATION OF THIS FILE                  #
################################################################################

#Modify this file as you wish, go to the folder containing this file, and run
# sage codebuilder.sage && mv codebuilder.sage.py codebuilder.py
# or 
# sagealone codebuilder.sage && mv codebuilder.sage.py codebuilder.py
#Reset the kernel for the modifications to reflect

################################################################################
#                             VARIABLE DEFINITIONS                             #
################################################################################

code = ''
for i in range(1,21):
    code+='x'+str(i)+'= var(\'x'+str(i)+'\')'
    if i<20:
        code+='\n'
exec(code)

code = 'variables = ['
for i in range(1,21):
    code+='x'+str(i)
    if i<20:
        code+=', '
code +=']'

exec(code)

################################################################################
#          FUNCTIONS THAT MAKE THE GENERATOR MATRIX FROM POLYNOMIALS           #
################################################################################

#Recursive function that efficiently generates the subspaces
def Rec_Gen_Maker(poly,num_var,depth,arr,G):
    if depth==num_var:
        if int(poly)%2==1:
            G.append(([1]+arr).copy())
        return
    for i in range(2):
        arr[depth]=i
        Rec_Gen_Maker(poly.subs(variables[depth]==i),num_var,depth+1,arr,G)

#Returns the generator matrix
def Generator_Matrix(poly,num_var):
    arr=[0]*num_var
    G=[]
    Rec_Gen_Maker(poly, num_var, 0,arr,G)
    G=matrix(G)
    G=G.transpose()
    return matrix(Integers(2),G)

################################################################################
#                    OUTPUTS THE WEIGHT ENUMERATOR FUNCTION                    #
################################################################################

def Weight_Enumerator(A):
    cols=A.ncols()
    rows=A.nrows()
    Prof=[0]*(cols+1)
    num = 0
    for i in range(2^rows):

        dig = num.digits(2)
        num = num +1
        dig = vector(dig + [0]*(rows-len(dig)))
        Prof[(dig*A).hamming_weight()]+=1
    return Prof

################################################################################
#                    CHECKS IF A SUBSPACE IS TRIORTHOGONAL                     #
################################################################################

def Is_Tri(X):
    S = 0
    cs = X.ncols()
    for a in X:
        for b in X:
            for c in X:
                for i in range(cs):
                    S+=a[i]*b[i]*c[i]
                if (S%2)==1:
                    return False
    return True

################################################################################
#         THIS PART OF THE CODE GENERATES THE ODD AND EVEN DESCENDANTS         #
#          SEE THE MAIN JUPYTER NOTE BOOK FOR INSTRUCTIONS ON HOW TO           #
#                             USE THESE FUNCTIONS.                             #
################################################################################

#Recursive function that generates the triorthogonal matrices 
def Rec_Tri_Mat_Maker(G,c,k,depth,chosen,A,B,func,arg):
    if chosen>1:
        if A.rank()<chosen:
            return False
    if depth==c:
        if chosen==k:            
            X=A.augment(B).echelon_form()
            X=X[:,k:]
#------------------------------------------------------------------------------#
# TODO: At this point the matrix X is a triorthogonal even descendants of a    #
#       the subspace G. We can do anything with X through the function f       #
#  We generically choose to check for distance                                 #
#  This part should be called through the Generate_Even_Desc(G,k,arg) function #
#  or the Generate_Odd_Desc(G,k,arg) function the input is passed over from    #
#  arg there                                                                   #
#  If returned "True", the search immediately ends
            return func(X,arg)
#------------------------------------------------------------------------------#
        return False
    
    if chosen<k:
        if type(A)==type([]):
            A=matrix(G.column(depth)).transpose()
            if Rec_Tri_Mat_Maker(G,c,k,depth+1,chosen+1,A,B,func,arg):
                return True
            A=[]
        else:
            A=A.augment(G.column(depth))
            if Rec_Tri_Mat_Maker(G,c,k,depth+1,chosen+1,A,B,func,arg):
                return True
            A=A.delete_columns([A.ncols()-1])
    
    if type(B)==type([]):
        B=matrix(G.column(depth)).transpose()
        if Rec_Tri_Mat_Maker(G,c,k,depth+1,chosen,A,B,func,arg):
            return True
        B=[]
    else:
        B=B.augment(G.column(depth))
        if Rec_Tri_Mat_Maker(G,c,k,depth+1,chosen,A,B,func,arg):
            return True
        B=B.delete_columns([B.ncols()-1])
    return False


#This function generates all even descendants of G with k logical variables, and runs the function 
# func on them with inputs func(X,arg), where X is the triorthogonal matrix
def Generate_Even_Desc(G,k,func,arg):
    c=G.ncols()
    A=[]
    B=[]
    return Rec_Tri_Mat_Maker(G,c,k,0,0,A,B,func,arg)

#This function generates all odd descendants of G with k logical variables, and runs the function 
# func on them with inputs func(X,arg), where X is the triorthogonal matrix
def Generate_Odd_Desc(G,k,func,arg):
    c=G.ncols()
    r=G.nrows()
    for marker in range(c):
        Desc = copy.copy(G)
        for i in range(1,r):
            if Desc[i,marker]==1:
                Desc[i,:]+=Desc[0,:]
        rws = range(1,r)
        cls = [i for i in range(c) if i not in [marker]]
        Desc=Desc[rws,cls]
        #At this point, one row and one column is removed from Desc, it is sufficient to consider
        #  odd descendants of Desc
        if Generate_Even_Desc(Desc,k,func,arg):
            return True
    return False

################################################################################
#        FUNCTIONS FOR CHECKING THE DISTANCE OF A TRIORTHOGONAL MATRIX         #
################################################################################

#Recursive function for checking the distance
def Rec_Weight_Maker(G,k,s,d,chosen,marker,Vec):
    if marker+d-chosen>s:
        return True
    if chosen==d:
        col=G*Vec
        if not(col[0:k].is_zero()) and col[k:].is_zero() == True:
            return False
        return True
    for i in range(marker,s):
        Vec[i]=1
        if not(Rec_Weight_Maker(G,k,s,d,chosen+1,i+1,Vec)):
            return False
        Vec[i]=0
    return True

#The following function returns true if the distance of Desc with k logical 
#qubits is larger than d
def Is_Dist_Larger(Desc, d):
    s=Desc.ncols()
    k=sum(Desc.columns()).hamming_weight()
    Vec=vector(Integers(2),[0]*s)
    return Rec_Weight_Maker(Desc,k,s,d,0,0,Vec)














