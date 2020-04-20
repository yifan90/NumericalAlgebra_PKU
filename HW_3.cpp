# include "n_algebra.h"

double h=.1;
int Max=4000;
element_b_m T;

double f(double x,double y){
	return 0;
	return 2*(x*(x-1)+y*(y-1));
}

double u(double x,double y){
	return 0;
	return x*(x-1)*y*(y-1);
}

void T0(int n){
	T.initial(n);
	T.n=n;
	T.is_zero=1;
	for(int i=1;i<n;i++){
			T.a[m(i,i,n)]=4;
			T.a[m(i,i+1,n)]=-1;
			T.a[m(i+1,i,n)]=-1;		
	}
	T.a[m(n,n,n)]=4;
    return;
}

void initial0(block_matrix& A,block_vector& b){
	int n=A.M;
	
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			b.b[i].b[j]=f((j+1)*h,(i+1)*h)*h*h;
		}
		b.b[i].is_zero=1;
	}

	for(int i=1;i<n;i++){
		A.A[m(i,i,n)]=T;
		A.A[m(i+1,i,n)].is_zero=1;
		A.A[m(i,i+1,n)].is_zero=1;
		for(int j=1;j<n+1;j++){
			A.A[m(i,i+1,n)].a[m(j,j,n)]=-1;
			A.A[m(i+1,i,n)].a[m(j,j,n)]=-1;
		}
	}
	A.A[m(n,n,n)]=T;   
	return;
}

void initial(block_matrix& A){
	int n=A.M;
	for(int i=1;i<=n;i++) A.A[m(i,i,n)]=T;
	return;
}

void b_jacobi(block_matrix& A,block_vector& b,block_vector& r,block_vector& e){
	int N=A.M;
	int n=A.n;
	int step=0;
	block_vector x1,x2,t;
	element_b_v s;
	x1.initial(N,n);x2.initial(N,n);t.initial(N,n);s.initial(n);
	//for(int i=0;i<N;i++) x2.b[i].rand2();
	x2=r;
	do{
		x1=x2;
		x2=b;
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				if(j==i) continue;
				if(A.A[m(i+1,j+1,N)].is_zero==0) continue;
				s=A.A[m(i+1,j+1,N)]*x1.b[j];
				x2.b[i]=x2.b[i]-s;
			}
			gauss_pca(A.A[m(i+1,i+1,N)].a,x2.b[i].b,n);
		}
		t=x2-x1;
		step++;
		initial(A);
	}while(step<Max && t.norm2()>pow(10,-12));

	initial(A);
	t=b;
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			if(A.A[m(i+1,j+1,N)].is_zero==0) continue;
			s=A.A[m(i+1,j+1,N)]*x2.b[j];
			t.b[i]=t.b[i]-s;
		}
	}
	cout<<"using block jacobi L2_error is "<<sqrt(t.norm2())<<endl;

	if(step==Max) cout<<"fails to converge after "<<Max<<" times"<<endl;
	else cout<<"totally iterates "<<step<<" times"<<endl;
	
	x1=r-e;
	t=x2-e;
	cout<<"the convergence speed R is about "<<-log(t.norm2()/x1.norm2())/step<<endl;
	r=x2;
	return;
}

void b_g_s(block_matrix& A,block_vector& b,block_vector& r,block_vector& e){
	int N=A.M;
	int n=A.n;
	int step=0,end=0;
	block_vector x1,x2,t;
	element_b_v s;
	x1.initial(N,n);x2.initial(N,n);t.initial(N,n);	s.initial(n);
	//for(int i=0;i<N;i++) x2.b[i].rand2();
	x2=r;
	do{
		x1=x2;
		x2=b;
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				if(j==i) continue;
				if(j>i){
                    if(A.A[m(i+1,j+1,N)].is_zero==0) continue;
					s=A.A[m(i+1,j+1,N)]*x1.b[j];
					x2.b[i]=x2.b[i]-s;
				}
				else{
                    if(A.A[m(i+1,j+1,N)].is_zero==0) continue;
                    s=A.A[m(i+1,j+1,N)]*x2.b[j];
					x2.b[i]=x2.b[i]-s;
				}
			}
			gauss_pca(A.A[m(i+1,i+1,N)].a,x2.b[i].b,n);
		}
		t=x2-x1;
		step++;
		initial(A);
	}while(step<Max && t.norm2()>pow(10,-12));
	
	initial(A);
	t=b;
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			if(A.A[m(i+1,j+1,N)].is_zero==0) continue;
			s=A.A[m(i+1,j+1,N)]*x2.b[j];
			t.b[i]=t.b[i]-s;
		}
	}
	cout<<"using block Gauss-Seidel L2_error is "<<sqrt(t.norm2())<<endl;

	if(step==Max) cout<<"fails to converge after "<<Max<<" times"<<endl;
	else cout<<"toltally iterates "<<step<<" times"<<endl;
    x1=r-e;
	t=x2-e;
	cout<<"the convergence speed R is about "<<-log(t.norm2()/x1.norm2())/step<<endl;
	r=x2;
	
	return;
}

void b_sor(block_matrix& A,block_vector& b,block_vector& r,block_vector& e,double w){
	int N=A.M;
	int n=A.n;
	int step=0;
	block_vector x1,x2,t;
	element_b_v s;
	x1.initial(N,n);
	x2.initial(N,n);
	t.initial(N,n);
	s.initial(n);
	//for(int i=0;i<N;i++) x2.b[i].rand2();
	x2=r;
	do{
		x1=x2;
		x2=b;
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				if(j==i) continue;
				if(j>i){
                    if(A.A[m(i+1,j+1,N)].is_zero==0) continue;
					s=A.A[m(i+1,j+1,N)]*x1.b[j];
					x2.b[i]=x2.b[i]-s;
				}
				else{
                    if(A.A[m(i+1,j+1,N)].is_zero==0) continue;
                    s=A.A[m(i+1,j+1,N)]*x2.b[j];
					x2.b[i]=x2.b[i]-s;
				}
			}
			x2.b[i]=w*x2.b[i];
			s=A.A[m(i+1,i+1,N)]*x1.b[i];
			s=(1-w)*s;
			x2.b[i]=x2.b[i]+s;
			gauss_pca(A.A[m(i+1,i+1,N)].a,x2.b[i].b,n);
		}
		t=x2-x1;
		step++;
		initial(A);
	}while(step<Max && t.norm2()>pow(10,-12));
    
	initial(A);
	t=b;
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			if(A.A[m(i+1,j+1,N)].is_zero==0) continue;
			s=A.A[m(i+1,j+1,N)]*x2.b[j];
			t.b[i]=t.b[i]-s;
		}
	}
	cout<<"using block SOR L2_error is "<<sqrt(t.norm2())<<" the w is "<<w<<endl;

	if(step==Max) cout<<"fails to converge after "<<Max<<" times"<<endl;
	else cout<<"totally iterates "<<step<<" times"<<endl;

    x1=r-e;
	t=x2-e;
	cout<<"the convergence speed R is about "<<-log(t.norm2()/x1.norm2())/step<<endl;
	r=x2;
	
	return;
}

/* 
double w_opt(block_matrix& A,block_vector& b,block_vector& r,block_vector& e){
	double t=1.2,s=.1,wopt;
	int i=1,n,m=Max;
	block_vector r0;
	r0.initial(r.m,r.n);
	r0=r;
	for(;t<2;t+=s){
		initial(A);
		r=r0;
		n=b_sor(A,b,r,e,t);
		if(n<m){
			m=n;
			wopt=t;
		}
		cout<<"w="<<t<<" take "<<n<<"steps"<<endl;
	}
	return wopt;
}*/

int main(){
	double N_t=1/h;
	int N=N_t;
	block_matrix A;
	block_vector b,r1,r2,r3,exact;
	element_b_v s;
	A.initial(N-1,N-1,N-1);b.initial(N-1,N-1);r1.initial(N-1,N-1);r2.initial(N-1,N-1);r3.initial(N-1,N-1);exact.initial(N-1,N-1);s.initial(N-1);
    T0(N-1); 
	initial0(A,b);

	for(int i=0;i<N-1;i++){
		for(int j=0;j<N-1;j++){
			exact.b[i].b[j]=u((j+1)*h,(i+1)*h);
		}
		exact.b[i].is_zero=1;
		r1.b[i].ones();
	}
	r2=r1;
	r3=r1;
	initial(A);

	for(int i=0;i<N-1;i++){
		for(int j=0;j<N-1;j++){
			if(A.A[m(i+1,j+1,N-1)].is_zero==0) continue;
			s=A.A[m(i+1,j+1,N-1)]*exact.b[j];
			b.b[i]=b.b[i]-s;
		}
	}
	cout<<"using the exact solution to the equation, we have a L2 error of "<<sqrt(b.norm2())<<endl<<endl;

	initial(A);
	b_jacobi(A,b,r1,exact);
	r1=r1-exact;
	cout<<"L2 error with the exact solution "<<sqrt(r1.norm2())<<endl;
	cout<<"infinite norm error with the exact solution "<<r1.normif()<<endl<<endl;
	
	initial(A);
	b_g_s(A,b,r2,exact);
    r2=r2-exact;
    cout<<"L2 error with the exact solution "<<sqrt(r2.norm2())<<endl;
	cout<<"infinite norm error with the exact solution "<<r1.normif()<<endl<<endl;

	initial(A);
	b_sor(A,b,r3,exact,.9);
    r3=r3-exact;
    cout<<"L2 error with the exact solution "<<sqrt(r3.norm2())<<endl;
	cout<<"infinite norm error with the exact solution "<<r1.normif()<<endl<<endl;

	initial(A);
	//cout<<"the best choice of w is about "<<w_opt(A,b,r3,exact)<<endl;
	return 0;
}
