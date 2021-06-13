////////////////////////////////////////////////////////////////////////////////
//                              THE MAIN C++ CODE                             //
//        FAST AND PARALLELIZED CODE FOR FINDING REPRESENTATIVES OF ALL       //
//      LOW WEIGHT REED-MULLER POLYNOMIALS OF 8 VARIABLES AND OF DEGREE 4     //
////////////////////////////////////////////////////////////////////////////////
#include "stdio.h"
#include "stdlib.h"
#include <time.h>
#include <pthread.h>
#include <chrono>
#include <iostream>
#include <fstream>
using namespace std;
using namespace chrono;
#define ulong unsigned long 
#define uchar unsigned char


int MAX_NUM_POLYS = 10000;   //Maximum number of polynomials in each step. Increase if necessary.
int NUM_THREADS = 8;         //Number of threads. 
int trigger_wait=1;          //Wait value of the first level simplification process. (see
														 //  quick_simplify in the Polynomial class.)
int trigger_random_jumps=0;  //Random_jump value of the first level simplification process.
int weight=36;               //Hamming weight of the code
int num_bases=2;            //Number of base polynomials (see the paper and jupyter notebook) 
                             //  for the definition of the base polynomials. 

class Polynomial;

// Input data struction of the threads/
struct thread_data {
	 int weight;
	 int  num_bases;
	 Polynomial *Base1,*Base2;
   Polynomial *poly_list;
	 int num_polys;
};

////////////////////////////////////////////////////////////////////////////////
//                     DEFINITION OF SOME BASIC FUNCTIONS                     //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Monomials correspond to uchar numbers. This code supports up to 8 different
//   variables. For example, the monomial x4*x2*x1 is coded as 0b00001011
// The functions bin_represntation, bin_represntation, bin_represntation print 
//   the binary code of numbers i.e. bin_represntation(0b00110101)
//   prints 00110101
// The function monomial_truth takes a monomial and returns its truth table. 
//   (See the Polynomial class for the definition of the truth table)
// The truth table is only supported for monomials of 6 variables, i.e., from 
//   0b00000000 to 0b00111111.
// The truth table is simply all evaluations of the monomial. It is a list of 
//   64 numbers, and is encoded in a ulong number.
// You can print the truth table with bin_represntation. For example the truth 
//   table of x2*x4 is simply bin_represntation(monomial_truth(0b00001010)). 
//   (the variables start from x1)
// The function hamming_weight returns the hamming weight of a ulong table.
// In this code we generally ignore the constant term. It will change the 
//   weight form w-> 64-w.
     

//List of all second order monomials
uchar second_order_monomials[21];
ulong second_order_table[21];
//List of all one bit ulong numbers. Having them in memory makes the hamming weight 
//  calculation slightly faster.
ulong one_bit[64];

template <typename T>
void bin_representation(T t)
{
	for(int i=8*sizeof(t)-1;i>=0;i--)
    cout<<(int)( (t>>i)%2);
  cout<<" ";
  return ;
}

ulong monomial_truth(uchar monomial)
{
	ulong table=0;
	for(uchar i=0;i<64;i++)
		if ((~((~monomial)|i))==0)
			table=table | (((ulong) 1)<<i);
	return table;
}


int hamming_weight(ulong n)
{
	int i=0;
	int weight=0;
	//This way of writing this code might seem weird, but this runs much faster than  
	// a single for loop!
	while(i^8){
		i++;
		weight+=
		  (n&one_bit[0])+
		 ((n&one_bit[1])>>1)+
		 ((n&one_bit[2])>>2)+
		 ((n&one_bit[3])>>3)+
		 ((n&one_bit[4])>>4)+
		 ((n&one_bit[5])>>5)+
		 ((n&one_bit[6])>>6)+
		 ((n&one_bit[7])>>7);
		n>>=8;
	}
	return weight;
}
	
//initializes second_order_monomials
void init()
{
	int indx=0;
	for(int i =0;i<6;i++){
		second_order_monomials[indx]=((uchar)1)<<i;
		indx++;
	}
	for(int i=0;i<6;i++)
		for(int j=i+1;j<6;j++){
			second_order_monomials[indx]=(((uchar)1)<<i) ^(((uchar)1)<<j);
			indx++;
		}
	for(int i=0;i<21;i++)
		second_order_table[i]=monomial_truth(second_order_monomials[i]);

	for(int i=0;i<64;i++)
		one_bit[i]=1<<i;
  
	return ;	
}


////////////////////////////////////////////////////////////////////////////////
//                        LOADING THE POLYNOMIAL CLASS                        //
//           PLEASE SEE THE polynomial.h LIBRARY FOR A COMPREHENSIVE          //
//                  EXPLANATION OF FUNCTIONS AND PROPERTIES.                  //
////////////////////////////////////////////////////////////////////////////////

#include "polynomial.h"

//SOME of the public objects of the Polynomial class:

//  int num_terms:
//  Stores the number of terms in the polynomial.

//  void add_term(uchar t):
//  Adds a term in the format of uchar (see above). Do not feed repeated terms.

//  void print():
//  Prints the polynomial

//  void clear():
//  Sets the num_term=0, and removes all auxiliary buffers.
//  See below for instruction

//  ulong truth_table():
//  Returns the truth table of the polynomial as bits of a unsigned long number.
//  This function ONLY works for polynomials of first 6 variables, i.e., x1,x2,
//  ...,x6. The truth table is the list all evaluations of the polynomial on 
//  the binary digits of numbers 0 to 63 as values of variables.

//  void operator=(Polynomial p):
//  Equal operator for the polynomials. Does not copy the buffers though.

//  bool operator==(Polynomial p):
//  Check if the polynomials are equal to each other UP TO A PERMUTATION OF 
//  VARIABLES.

//  void set_buffers(Polynomial buffers[3]):
//  Sets buffers. The polynomial class needs three more instances of itself to
//  to operate the simplification process. It is not memory and time efficient 
//  to generate these buffers each time, so I feed the same buffer to different
//  polynomials. This function is used to assign these buffers.  

//  void sort():
//  This function permutes the variables such that x1 appears the most in the 
//  terms, then x2, then x3, and so on. It also sorts the the terms in the poly 
//  array.

//  void quick_simplify(int wait,int random_jumps);
//  It performs the change of variables and affine transformations to reduce the
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

////////////////////////////////////////////////////////////////////////////////
//      THE FOLLOWING TWO FUNCTIONS ARE THE MAIN PARTS OF THE CODE THAT       //
//     TEST ALL POSSIBLE COMBINATIONS OF THE TWO BODY TERMS TO ADD TO THE     //
//      POLYNOMIAL AND RETURN THE SIMPLIFIED POLYNOMIALS WITH THE CORRECT     //
//                    HAMMING WEIGHT AS A POLYNOMIAL ARRAY                    //
////////////////////////////////////////////////////////////////////////////////

//The main recursive funciton. It will be called from the generate_poly_list 
//  function below. This function encodes the polynomial truth table in ulong
//  variable for faster processing.
void rec(ulong table,int lev,int code,int &target,Polynomial &poly, Polynomial &Base1, Polynomial &Base2, Polynomial *poly_list, int *num_polys)
{
	if(lev==21){
		int ham_w = hamming_weight(table);
		if(ham_w==target || ham_w==64-target){
			poly.num_terms=0;
			for(int i=0;i<Base1.num_terms;i++)
				poly.add_term(((uchar)1)^(Base1.poly[i]<<2));	
			for(int i=0;i<Base2.num_terms;i++)
				poly.add_term(((uchar)2)^(Base2.poly[i]<<2));	

			int i=21;
			while(code){	
				i--;
				if(code&(uchar)1)
					poly.add_term(((uchar)1)^((uchar)2)^(second_order_monomials[i]<<2));
				code=code>>1;
			}
			if(ham_w==64-target)
	      poly.add_term(((uchar)1)^((uchar)2));
			poly.quick_simplify(trigger_wait,trigger_random_jumps);
			for(int i=0;i<*num_polys;i++)
				if(poly_list[i]==poly)
					return ;
			poly_list[*num_polys]=poly;
			(*num_polys)++;
			if((*num_polys)>MAX_NUM_POLYS){
				printf("\nERROR: Too many polynomials found\n");
				abort();
			}
		}
		return ;
	}
	rec(table                        ,lev+1,code<<1    ,target,poly,Base1,Base2,poly_list,num_polys);
	rec(table^second_order_table[lev],lev+1,(code<<1)^1,target,poly,Base1,Base2,poly_list,num_polys);
	return ;
}

//This function takes two base polynomials, and places all possible polynomials with that base and
//  given weight into the poly_list array. The num_polys will reflect the total number of 
//  polynomials in the list and will be modified as we call this function. This function does one 
//  quick level of polynomial simplification and do not add repreated polynomials. Lastly, it can 
//  be called on the same poly_list over and over again and it will just add extra polynomials that 
//  it finds to the list. 

void generate_poly_list(Polynomial &Base1,Polynomial &Base2, int weight,Polynomial *poly_list, int &num_polys)
{
	Polynomial poly, buffers[3];
	poly.set_buffers(buffers);
	int target=weight-hamming_weight(Base1.truth_table())-hamming_weight(Base2.truth_table());
	rec(Base1.truth_table()^Base2.truth_table(),0,0,target,poly,Base1,Base2,poly_list,&num_polys);
}

////////////////////////////////////////////////////////////////////////////////
//                  EXTRA LEVELS OF POLYNOMIAL SIMPLIFICATION                 //
//      THE FOLLOWING FUNCTIONS TAKE A POLYNOMIAL LIST AND SIMPLIFY THEM      //
//      AS MUCH AS POSSIBLE AND REMOVE THE AFFINE EQUIVALENT POLYNOMIALS.     //
////////////////////////////////////////////////////////////////////////////////

//The following function is usually called from the simplify_poly_list function.
void shorten_poly_list(int wait, int random_jumps, Polynomial *poly_list, int &num_polys)
{
	Polynomial buffers[3];
	bool *active_poly=new bool[MAX_NUM_POLYS];
	for(int i=0;i<num_polys;i++)
		active_poly[i]=true;
	for(int i=0;i<num_polys;i++){
		poly_list[i].set_buffers(buffers);
		poly_list[i].quick_simplify(wait,random_jumps);
	}
	for(int i=0;i<num_polys;i++)
		if(active_poly[i])
			for(int j=i+1;j<num_polys;j++)
				if(poly_list[i]==poly_list[j])
					active_poly[j]=false;
	int new_num_polys=0;
	for(int i=0;i<num_polys;i++)
		if(active_poly[i]){
			if(i>new_num_polys)
				poly_list[new_num_polys]=poly_list[i];
			new_num_polys++;
		}
	num_polys=new_num_polys;
}

void simplify_poly_list(Polynomial *poly_list, int &num_polys)
{
	shorten_poly_list(10,3,poly_list,num_polys);
	shorten_poly_list(50,5,poly_list,num_polys);
	for(int i=0;i<100;i++)
		shorten_poly_list(50,5,poly_list,num_polys);
	for(int i=0;i<10;i++)
		shorten_poly_list(100,10,poly_list,num_polys);
	for(int i=0;i<100;i++)
		shorten_poly_list(50,5,poly_list,num_polys);
	for(int i=0;i<10;i++)
		shorten_poly_list(100,10,poly_list,num_polys);
}

////////////////////////////////////////////////////////////////////////////////
//               THREAD ASSIGNMENT FUNCTIONS FOR PARALLELIZING.               //
////////////////////////////////////////////////////////////////////////////////

//Main thread function.
void *thread_function(void *var) 
{
	struct thread_data *data;
	data = (struct thread_data *) var;
	Polynomial Base1, Base2;
	for(int i=0;i<data->num_bases;i++){
		Base1.clear();
		Base2.clear();
		Base1=data->Base1[i];
		Base2=data->Base2[i];
		generate_poly_list(Base1,Base2,data->weight,data->poly_list,data->num_polys);
	}
	shorten_poly_list(10,3,data->poly_list,data->num_polys);
	shorten_poly_list(20,4,data->poly_list,data->num_polys);
  pthread_exit(NULL);
}

//This function mixes the polynomial lists that outputs of each thread
void mix_poly_lists(Polynomial *final_poly_list,int &final_num_polys, thread_data *data)
{
	int c=0;
	for(int i=0;i<NUM_THREADS;i++)
		for(int j=0;j<data[i].num_polys;j++){
			final_poly_list[c]=data[i].poly_list[j];
			c++;
		}
	simplify_poly_list(final_poly_list,final_num_polys);
	return ;
}


thread_data *read_file_init_thread_data()
{
  thread_data *data;
	ifstream file;
	file.open(".//poly_finder_instructions.txt");
	file>>weight;
	file>>num_bases;
	file>>NUM_THREADS;
	file>>trigger_wait;
	file>>trigger_random_jumps;
	file>>MAX_NUM_POLYS;
	data=new thread_data[NUM_THREADS];
  Polynomial *Base1 = new Polynomial[num_bases];
  Polynomial *Base2 = new Polynomial[num_bases];
	int bn;
	unsigned int n;
	for(int i=0;i<num_bases;i++){
		Base1[i].clear();
		Base2[i].clear();
		file>>bn;
		for(int j=0;j<bn;j++){
			file>>n;
			Base1[i].add_term((uchar)n);
		}
		file>>bn;
		for(int j=0;j<bn;j++){
			file>>n;
			Base2[i].add_term((uchar)n);
		}
	}
	file.close();

  int c;
  for(int i=0;i<NUM_THREADS;i++){
    data[i].num_bases=0;
    data[i].weight=weight;
    data[i].poly_list = new Polynomial[MAX_NUM_POLYS];
  }
  c=0;
  for(int i=0;i<num_bases;i++){
    data[c].num_bases++;
    c++;
    c%=NUM_THREADS;
  }
  for(int i=0;i<NUM_THREADS;i++){
    data[i].Base1=new Polynomial[data[i].num_bases];
    data[i].Base2=new Polynomial[data[i].num_bases];
  }

  c=0;
  for(int i=0;i<NUM_THREADS;i++)
    for(int j=0;j<data[i].num_bases;j++){
      data[i].Base1[j]=Base1[c];
      data[i].Base2[j]=Base2[c];
      c++;
    }
	return data;
}


int main()
{
	init();
	cout<<"Reading the instruction file and polynomials...";
  thread_data *data = read_file_init_thread_data();
  cout<<"Done!"<<endl;
  Polynomial **poly_list = new Polynomial*[NUM_THREADS];
  for(int i = 0; i < NUM_THREADS; ++i)
    poly_list[i] = new Polynomial[MAX_NUM_POLYS];
  int *num_polys=new int[NUM_THREADS];
	
	pthread_t threads[NUM_THREADS];
  pthread_attr_t attr;
  void *status;

  // Initialize and set thread joinable
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	cout<<"Finding all polynomials with weight "<<weight<<" ... "<<endl;
	cout<<"Starting timer!"<<endl;
	high_resolution_clock::time_point start = high_resolution_clock::now();
	int rc;
  for(int i = 0; i < NUM_THREADS; i++ ) {
     printf("Creating thread %d \n",i);
     rc = pthread_create(&threads[i], &attr, thread_function, (void *)&data[i]);
     if (rc) {
        printf("Error:unable to create thread.");
        exit(-1);
     }
  }
  pthread_attr_destroy(&attr);
  for(int i = 0; i < NUM_THREADS; i++ ) {
    rc = pthread_join(threads[i], &status);
    if (rc) {
      printf("Error:unable to join thread.");
      exit(-1);
    }
  }
	cout<<"All threads done. All polynomials found."<<endl;
	int tot_num=0;
	for(int i=0;i<NUM_THREADS;i++)
		tot_num+=data[i].num_polys;
	cout<<"Total number of polynomials after one level of pruning: "<<tot_num<<endl;

  high_resolution_clock::time_point stop = high_resolution_clock::now();
  duration<double> duration = duration_cast<microseconds>(stop - start);
	cout<< "Time: "<< duration.count() << " seconds" << endl;
	cout<< "Simplifying polynomials, removing equivalent polynomials ... ";
	int final_num_polys=0;
	for(int i=0;i<NUM_THREADS;i++)
		final_num_polys+=data[i].num_polys;
	Polynomial *final_poly_list=new Polynomial[final_num_polys];
		
	mix_poly_lists(final_poly_list,final_num_polys,data);

	cout<<"Done!"<<endl;

	stop = high_resolution_clock::now();
	duration = duration_cast<microseconds>(stop - start);
  
  cout<< "Total time elapsed: "<< duration.count() << " seconds" << endl;	
	cout<< "Number of polynomial representatives: " << final_num_polys<<endl;
	cout<< "List of representatives: "<<endl;
	cout<<"[";
	for(int i=0;i<final_num_polys;i++){
		if(i)
			cout<<",\n";
    final_poly_list[i].print();
  }
	cout<<"]";
	return 0;
}



