# include <iostream>
# include <fstream>
# include <cstdlib>
# include <cmath>
# include <ctime>
using namespace std;
extern int seed;
extern int Max;

inline int m(int i,int j,int N) {
	return (j-1)*N+i-1;
}

inline double abs2 (double i) {
	return i>=0 ? i : -i;
}

inline double normif(double *a,int N) {
	double t=0;
	for(int i=0;i<N;i++) {
		if(t<abs2(a[i])) t=abs2(a[i]);
	}
	return t;
}

inline double dot(double *a,double *b,int N){
	double r=0;
	for(int i=0;i<N;i++) r+=a[i]*b[i];
	return r;
}

inline double norm2(double *a,int N) {               //a的模的平方！！！！
	double r=0;
	for(int i=0;i<N;i++) r+=a[i]*a[i];
	return r;
}

inline double sign(double i) {
	srand((unsigned)seed);
	seed++;
	if(i==0) return -1.0+2.0*(rand()%2);
	return i>0 ? 1.0 :-1.0;
}

class element_b_m{
	public:
		int n,is_zero;
		double* a;
		element_b_m():n(0),is_zero(0),a(NULL){}
		void initial(int i){
			n=i;
			is_zero=0;
		}
		element_b_m& operator=(element_b_m x){
			n=x.n;
			is_zero=x.is_zero;
			for(int i=0;i<n*n;i++) a[i]=x.a[i];
			return *this;
		}
		void zeros(){
			for(int j=0;j<n*n;j++) a[j]=0;
		}
};

class block_matrix{
	public:
		int M,N,n;
		element_b_m *A;
		block_matrix():M(0),N(0),n(0),A(NULL){}
		void initial(int a,int b,int c){
			M=a;
			N=b;
			n=c;
			A=new element_b_m[M*N];
			for(int i=0;i<M*N;i++){
				(A+i)->initial(n);
			}
		}
		block_matrix& operator=(block_matrix x){
			M=x.M;
			N=x.N;
			n=x.n;
			for(int i=1;i<=M;i++){
				for(int j=1;j<=N;j++) A[m(i,j,M)]=x.A[m(i,j,M)];
			}
			return *this;
		}
};

class element_b_v{
	public:
		int n,is_zero;
		double *b;
		element_b_v():n(0),is_zero(0),b(NULL){};
		void initial(int i){
			n=i;
			is_zero=0;
			b=new double[n];
			for(int j=0;j<i;j++) b[j]=0;
		}
		void zeros(){
			for(int i=0;i<n;i++) b[i]=0;
			return;
		}
		double normif() {
			double t=0;
			for(int i=0;i<n;i++) {
				if(t<abs2(b[i])) t=abs2(b[i]);
			}
			return t;
		}
		double norm2() {               //a的模的平方！！！！
			double r=0;
			for(int i=0;i<n;i++) r+=b[i]*b[i];
			return r;
		}
		void rand2(){
			for(int k=0;k<n;k++) b[k]=rand()%10;
			is_zero=1;
			seed++;
			return;
		}
		element_b_v& operator=(element_b_v x){
			n=x.n;
			is_zero=x.is_zero;
			for(int i=0;i<n;i++) b[i]=x.b[i];
			return *this;
		}
};

class block_vector{
	public:
		int m,n;
		element_b_v *b;
		block_vector():m(0),n(0),b(NULL){}
		void initial(int i,int j){
			m=i;
			n=j;
			b=new element_b_v[m];
			for(int k=0;k<m;k++){
				(b+k)->initial(n);
			}
		}
		double normif() {
			double t=0;
			for(int i=0;i<m;i++) {
				if(t<b[i].normif()) t=b[i].normif();
			}
			return t;
		}
		double norm2() {               //a的模的平方！！！！
			double r=0;
			for(int i=0;i<m;i++) r+=b[i].norm2();
			return r;
		}
		block_vector& operator=(block_vector x){
			m=x.m;
			n=x.n;
			for(int i=0;i<m;i++) b[i]=x.b[i];
			return *this;
		}
};

element_b_v operator+(element_b_v& a,element_b_v& b);

element_b_v operator-(element_b_v& a,element_b_v& b);

element_b_v operator*(element_b_m& A,element_b_v& b);

element_b_v operator*(double k,element_b_v& b);

block_vector operator+(block_vector& a,block_vector& b);

block_vector operator-(block_vector& a,block_vector& b);

block_vector operator*(double k,block_vector& b);

double dot(block_vector& a,block_vector& b) ;

void rand1(double* a,int N) ;

void rand2(double* a,int N) ;

double condition_num_if(double *A,int N) ;

void gauss (double *A,double *b,int N) ;

void gauss_pca (double *A,double *b,int N) ;

void sqrt_method(double *A,double *b,int N) ;

void improved_sqrt_method (double *A,double *b,int N) ;

void reflect(double *w,double *a,int N) ;

void householder(double *a,int N) ;

void qr(double *A,int M,int N) ;

void qrsolve(double *A,double *b,int N) ;

double ls(double *A,double *b,double *x,int M,int N) ;

void b_jacobi(block_matrix& A,block_vector& b,block_vector& r) ;

void b_g_s(block_matrix& A,block_vector& b,block_vector& r) ;

void b_sor(block_matrix& A,block_vector& b,block_vector& r,double w) ;

void cg(block_matrix& A,block_vector& b,block_vector& r) ;
