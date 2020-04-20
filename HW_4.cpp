# include "n_algebra.h"

double h=.05;
double pi=3.1415926;
element_b_m T;

double f(double x,double y){
	//return 0;
	return sin(x*y);
}

double u_boundary(double x,double y){
	//return 0;
	return x*x+y*y;
}

void T0(int n){
	T.initial(n);
	T.a=new double[n*n];
	T.n=n;
	T.is_zero=1;
	for(int i=1;i<n;i++){
		T.a[m(i,i,n)]=1+h*h/4.0;
		T.a[m(i,i+1,n)]=-1.0/4.0;
		T.a[m(i+1,i,n)]=-1.0/4.0;		
	}
	T.a[m(n,n,n)]=1+h*h/4.0;
	return;
}

void initial(block_matrix& A,block_vector& b){
	int n=A.M;

	for(int i=1;i<n;i++){
		A.A[m(i,i,n)].a=new double[n*n];
		A.A[m(i,i+1,n)].a=new double[n*n];
		A.A[m(i+1,i,n)].a=new double[n*n];
		A.A[m(i,i,n)].is_zero=1;
		A.A[m(i+1,i,n)].is_zero=1;
		A.A[m(i,i+1,n)].is_zero=1;
		for(int j=1;j<n+1;j++){
			A.A[m(i,i+1,n)].a[m(j,j,n)]=-1.0/4.0;
			A.A[m(i+1,i,n)].a[m(j,j,n)]=-1.0/4.0;
		}
		A.A[m(i,i,n)]=T;
	}
	A.A[m(n,n,n)].a=new double[n*n];
	A.A[m(n,n,n)]=T;   

	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++) b.b[i].b[j]=f((j+1)*h,(i+1)*h)*h*h/4.0;
		b.b[i].b[0]+=u_boundary(0,(i+1)*h)/4.0;
		b.b[i].b[n-1]+=u_boundary(1,(i+1)*h)/4.0;
		b.b[i].is_zero=1;
	}
	for(int i=0;i<n;i++){
		b.b[0].b[i]+=u_boundary((i+1)*h,0)/4.0;
		b.b[n-1].b[i]+=u_boundary((i+1)*h,1)/4.0;
	}
	return;
}

int main(){
	double N_t=1/h;
	int N=N_t;
	block_matrix A;
	block_vector b,r1,r2;
	A.initial(N-1,N-1,N-1);	b.initial(N-1,N-1);	r1.initial(N-1,N-1);r2.initial(N-1,N-1);
	T0(N-1);
	initial(A,b);
	cg(A,b,r1);
	b_sor(A,b,r2,.9);
	r1=r1-r2;
	cout<<"the two method have a L2 difference of "<<r1.norm2()<<endl;
	return 0;
}
