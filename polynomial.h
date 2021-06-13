////////////////////////////////////////////////////////////////////////////////
//                             POLYNOMIAL LIBRARY                             //
//             DEFINES THE 'POLYNOMIAL' CLASS USED FOR STORING AND            //
//                      MANIPULATING BINARY POLYNOMIALS.                      //
////////////////////////////////////////////////////////////////////////////////

class Polynomial
{
	public:

////////////////////////////////////////////////////////////////////////////////
//  The constructor only sets num_terms = 0
		Polynomial();

////////////////////////////////////////////////////////////////////////////////
//  Unsigned char variable, the polynomial terms are stored as 
//    binary numbers. For example the term x1*x3 is stored as 0b00000101.
		uchar poly[265];

////////////////////////////////////////////////////////////////////////////////
//	Stores the number of terms in the polynomial.
		int num_terms;

////////////////////////////////////////////////////////////////////////////////
//	Adds a term in the format of uchar (see above). Do not feed repeated terms.
		void add_term(uchar t);

////////////////////////////////////////////////////////////////////////////////
//  Prints the polynomial
		void print();

////////////////////////////////////////////////////////////////////////////////
//	Prints the polynomial in the binary format.
		void print_bin();

////////////////////////////////////////////////////////////////////////////////
//  Sets the num_term=0, and removes all auxiliary buffers.
//  See below for instruction
		void clear();

////////////////////////////////////////////////////////////////////////////////
//  Returns the truth table of the polynomial as bits of a unsigned long number.
//  This function ONLY works for polynomials of first 6 variables, i.e., x1,x2,
//  ...,x6. The truth table is the list all evaluations of the polynomial on 
//  the binary digits of numbers 0 to 63 as values of variables.
		ulong truth_table();

////////////////////////////////////////////////////////////////////////////////
//	Equal operator for the polynomials. Does not copy the buffers though.
		void operator=(Polynomial p);

////////////////////////////////////////////////////////////////////////////////
//  Check if the polynomials are equal to each other UP TO A PERMUTATION OF 
//  VARIABLES.
    bool operator==(Polynomial p);

////////////////////////////////////////////////////////////////////////////////
//  Sets buffers. The polynomial class needs three more instances of itself to
//  to operate the simplification process. It is not memory and time efficient 
//  to generate these buffers each time, so I feed the same buffer to different
//  polynomials. This function is used to assign these buffers.  
		void set_buffers(Polynomial buffers[3]);

////////////////////////////////////////////////////////////////////////////////
//  This function permutes the variables such that x1 appears the most in the 
//  terms, then x2, then x3, and so on. It also sorts the the terms in the poly 
//  array.
		void sort();

////////////////////////////////////////////////////////////////////////////////
//	It performs the change of variables and affine transformations to reduce the
//  number of terms in the polynomial. The basic steps of an affine 
//  transformation are the following operations:
//  1. plus_one: changes a variable x_i to x_i+1
//  2. transposition: changes a variable x_i to x_i + x_j
//  The function takes two variables wait and random_jumps.                     
//  At each step, the function performs random_jumps many random transpositions
//  and next perform a greedy minimization of the number of terms in the 
//  polynomial using transpositions and plus_ones until the number of terms is
//  constant. At the end of each step, the number of terms might decrease. If 
//  weight does not decrease for "wait" many consecutive steps, the function 
//  returns the polynomial with the minimum number of terms derived in any of 
//  the steps. 
		void quick_simplify(int wait,int random_jumps);
	private:
//  The following contains the number of repetition of each variable.		
		int profile[8];
//  The following is quick_simplify minus sort
		void simplify(int wait,int random_steps);
//  Performs all possible transpositions and accepts them if the number of 
//  terms decreases.
		bool transpositions();
//  Performs all possible plus_ones and accepts them if the number of 
//  terms decreases.
		bool plus_ones();
//  Makes the transposition and reports the extra terms into buffer 		
		void transposition(int var_source,int var_target);
//  Makes the plus_one and reports the extra terms into buffer 		
		void plus_one(int var);
//  Mixes the self and buffer into replica 
		void mix_buffer_to_replica();
//  Swap bits number var1 and var2 of num and returns the results
		uchar swap_bit(uchar num,int var1,int var2);
//  Permutes the bits of num according to perm
		uchar permute_bits(uchar num, int perm[8]);
//  Permutes the polynomial p (keeping the profile fixed) and checks if it 
//  matches self. Returns true if it matches. It is used in == operator.
		bool rec_permutation_maker(int perm[8],Polynomial p,int lev);
//  Makes the profile array
		void profile_maker();
//  Permutes the variables of the polynomial such that the profile array is 
//  decreasing.
		void sort_variables();
//  Sort the terms in poly
		void sort_terms();
//  Pointers to the three buffers
		Polynomial *replica,*buffer,*second_replica;
};


using namespace std;




Polynomial::Polynomial()
{
	num_terms=0;
}




void Polynomial::operator=(Polynomial p)
{
	num_terms=p.num_terms;
	for(int i=0;i<num_terms;i++)
		poly[i]=p.poly[i];
}




bool Polynomial::operator==(Polynomial p)
{
	if(num_terms!=p.num_terms)
		return false;
	p.profile_maker();
	profile_maker();
	for(int i=0;i<8;i++)
		if(profile[i]!=p.profile[i])
			return false;
	int perm[8];
	for(int i=0;i<8;i++)
		perm[i]=-1;
	return rec_permutation_maker(perm,p,0);
}




void Polynomial::set_buffers(Polynomial buffers[3])
{
  replica=&buffers[0];//replica1;
  buffer=&buffers[1];//buffer1;
  second_replica=&buffers[2];//second_replica1;
}




void Polynomial::clear()
{
	num_terms=0;
	replica=NULL;
	buffer=NULL;
	second_replica=NULL;
}




uchar Polynomial::permute_bits(uchar num, int perm[8])
{
	uchar out=0;
	for(int i=0;i<8;i++)
		if(num&((uchar)1<<i))
			out+=(uchar)1<<(perm[i]);
	return out;
}




bool Polynomial::rec_permutation_maker(int perm[8],Polynomial p,int lev)
{
	if(lev==8){
		bool flag;
		uchar c;
		for(int j,i=0;i<num_terms;i++){
			flag=false;
			c=permute_bits(poly[i],perm);
			for(j=0;j<num_terms;j++)
				if(p.poly[j]==c)
					flag=true;
			if(!flag)
				return false;
		}
		return true;
	}
	if(profile[lev]==0){
		perm[lev]=lev;
		return rec_permutation_maker(perm,p,lev+1);
	}
	for(int i=0;i<8;i++)
		if(profile[lev]==profile[i])
			if(perm[i]==-1){
				perm[i]=lev;
				if(rec_permutation_maker(perm,p,lev+1))
					return true;
				perm[i]=-1;
			}
	return false;
}




uchar Polynomial::swap_bit(uchar num, int var1,int var2)
{
	uchar output;
	bool f1=false,f2=false;
	if(num&((uchar)1<<var1))
		f1=true;
	if(num&((uchar)1<<var2))
		f2=true;

	if(f1)
		output=num|((uchar)1<<var2);
	if(!f1)
		output=num&(~((uchar)1<<var2));

	if(f2)
		output=output|((uchar)1<<var1);
	if(!f2)
		output=output&(~((uchar)1<<var1));

	return output;
}




void Polynomial::profile_maker()
{
	int i,j,c;
  for(i=0;i<8;i++)
    profile[i]=0;
  for(i=0;i<num_terms;i++){
    c=poly[i];
    for(j=0;j<8;j++){
      profile[j]+=c&((uchar)1);
      c=c>>1;
    }
  }
	return ;
}




void Polynomial::sort_variables()
{
	int i,j;
	int c;
	bool flag=true;
	while(flag){
		flag=false;
		for(i=0;i<7;i++)
			if(profile[i]<profile[i+1]){
				flag=true;
				for(j=0;j<num_terms;j++)
					poly[j]=swap_bit(poly[j],i,i+1);
				c=profile[i];
				profile[i]=profile[i+1];
				profile[i+1]=c;
			}
	}
}




void Polynomial::sort_terms()
{
	int i,j;
	bool flag=true;
	uchar swap_c;
	while(flag){
		flag=false;
		for(i=0;i<num_terms-1;i++)
			if(poly[i]>poly[i+1]){
				flag=true;
				swap_c=poly[i];
				poly[i]=poly[i+1];
				poly[i+1]=swap_c;
			}
	}
}




void Polynomial::sort()
{
	profile_maker();
	sort_variables();
	sort_terms();
}

void Polynomial::quick_simplify(int wait,int random_jumps)
{
	if(replica==NULL){
		cout<<"ERROR: buffer in quick simplify is not set";
		abort();
	}
	simplify(wait,random_jumps);
	sort();
}




void Polynomial::simplify(int wait,int random_steps)
{
	int min_terms=255;
	int steps=0;
	int a,b;
	*second_replica=*this;
	while(steps<wait){
		steps++;
		for(int i=0;i<random_steps;i++){
			a=rand()%8;
			b=rand()%8;
			if(a!=b){
				transposition(a,b);
				*this=(*replica);
			}
		}
		while(plus_ones()||transpositions());
		if(num_terms<min_terms){
			min_terms=num_terms;
			steps=0;
			*second_replica=*this;
		}
	}
	*this=*second_replica;
	return ;
}





ulong Polynomial::truth_table()
{
	ulong a=0;
	for(int i=0;i<num_terms;i++)
		a=a^monomial_truth(poly[i]);
	return a;
}





void Polynomial::add_term(uchar t)
{
	poly[num_terms]=t;
	num_terms++;
	return ;
}




void Polynomial::print_bin()
{
  for(int i=0;i<num_terms;i++)
    bin_representation(poly[i]);
  cout<<"\n";
  return ;
}




void Polynomial::print()
{
  bool flag2=false;
  for(int i=0;i<num_terms;i++){
    if(flag2)
      cout<<" + ";
    bool flag=false;
    uchar s=poly[i];
    for(int j =0;j<8;j++){
      if(s%2){
        if(flag)
          cout<<"*";
        cout<<"x"<<j+1;
        flag=true;
      }
      s=s>>1;
    }
    flag2=true;
  }
  return ;
}




void Polynomial::mix_buffer_to_replica()
{
	bool flag;
	(*replica).num_terms=0;
	for(int i=0;i<num_terms;i++){
		flag=true;
		for(int j=0;j<(*buffer).num_terms;j++)
			if((*buffer).poly[j]==poly[i]){
				(*buffer).poly[j]=0;
				flag=!flag;
			}
		if(flag){
			(*replica).poly[(*replica).num_terms]=poly[i];
			(*replica).num_terms++;
		}
	}
	for(int i=0;i<(*buffer).num_terms;i++){
		
		if((*buffer).poly[i]){
			flag=true;
			for(int j=i+1;j<(*buffer).num_terms;j++)
				if((*buffer).poly[j]==(*buffer).poly[i]){
					(*buffer).poly[j]=0;
					flag=!flag;
				}
			if(flag){
				(*replica).poly[(*replica).num_terms]=(*buffer).poly[i];
        (*replica).num_terms++;
			}
			(*buffer).poly[i]=0;
		}
	}
	return ;
}




void Polynomial::plus_one(int var)
{
	(*buffer).num_terms=0;
	for(int i=0;i<num_terms;i++)
		if(poly[i]&((uchar)1<<var)){
			(*buffer).poly[(*buffer).num_terms]=poly[i]^((uchar)1<<var);
			(*buffer).num_terms++;
		}
	mix_buffer_to_replica();
	return ;	
}





void Polynomial::transposition(int var_source, int var_target)
{
	(*buffer).num_terms=0;
	for(int i=0;i<num_terms;i++)
		if(poly[i]& ((uchar)1<<var_source)){
			(*buffer).poly[(*buffer).num_terms]=poly[i]^((uchar)1<<var_source)|((uchar)1<<var_target);
			(*buffer).num_terms++;
		}
	mix_buffer_to_replica();
	return ;	
}





bool Polynomial::plus_ones()
{
  bool flag=false;
  for(int var=0;var<8;var++){
    plus_one(var);
    if((*replica).num_terms<num_terms){
			*this=(*replica);
      flag=true;
    }
  }
  return flag ;
}





bool Polynomial::transpositions()
{
  bool flag=false;;
  for(int var_source=0;var_source<8;var_source++)
    for(int var_target=0;var_target<8;var_target++)
      if(var_target!=var_source){
        transposition(var_source,var_target);
        if( (*replica).num_terms < num_terms){
					*this=(*replica);
          flag=true;
        }
      }
  return flag;
}
