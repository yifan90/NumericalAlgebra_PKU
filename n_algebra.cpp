# include "n_algebra.h"

int seed=(int) (time(0)%100);
int Max=2000;

double dot(block_vector& a,block_vector& b){
	int N=a.m;
	int n=a.n;
	double r=0;
	for(int i=0;i<N;i++) r+=dot(a.b[i].b,b.b[i].b,n);
	return r;
}

void rand1(double* a,int N) {
	double temp=0;
	srand((unsigned)seed);
	for (int i=0;i<N;i++){
		a[i]=(rand()%200)*(-1.0+2.0*(rand()%2));
		temp+=abs(a[i]);
	}
	for(int i=0;i<N;i++) a[i]/=temp;
	seed++;
	return;
}

void rand2(double* a,int N) 
{
	srand((unsigned)seed);
	for (int i=0;i<N;i++){
		a[i]=(rand()%100)*(-1.0+2.0*(rand()%2));
	}
	seed++;
	return;
}

double condition_num_if(double *A,int N)
{
	int l=0;
	double t1=0,t2=0,t3=0,M=0;
	double *A1=new double[N*N],*A2=new double[N*N];

	for(int i=0;i<N*N;i++) A1[i]=A[i];
	for(int i=1;i<N+1;i++) {
		for(int j=1;j<N+1;j++) A2[m(i,j,N)]=A[m(j,i,N)];
	}

	for(int i=1;i<N+1;i++) {
		double temp=0;
		for (int j=1;j<N+1;j++) temp+=abs(A[m(i,j,N)]);
		if(M<temp) M=temp;
	}
	t1=M;

	for(;l<10;l++) {
		int k=1,j,num=0;
		double n_z,*x=new double[N],*v=new double[N],*z=new double[N],*x2=new double[N];
		rand1(x,N);

		while(k==1 && num<4*N) {
			double t=0;
			for(int i=0;i<N*N;i++) A1[i]=A[i];
			for(int i=1;i<N+1;i++) {
				for(int j=1;j<N+1;j++) A2[m(i,j,N)]=A[m(j,i,N)];
			}
		    for(int i=0;i<N;i++) x2[i]=x[i];

			gauss_pca(A2,x,N);
			t3=0;
			for(int i=0;i<N;i++) {
				t3+=abs(x[i]);
				v[i]=sign(x[i]);
			}
			gauss_pca(A1,v,N);
			j=0;
			n_z=0;
			for(int i=0;i<N;i++) {
				if(n_z<abs(v[i])) {
					n_z=abs(v[i]);
					j=i;
				}
			}
			for(int i=0;i<N;i++) t+=x2[i]*v[i];
			if(n_z<=t) {
				k=0;
			}
			else {
				for(int i=0;i<N;i++) x[i]=0;
				x[j]=1;
				num++;
			}
		}
		if(t3>t2) t2=t3;
		delete []x; delete []v; delete []z; delete []x2;
	}
	delete []A1;
	delete []A2;
	return t1*t2;
}

void gauss (double *A,double *b,int N) {
	for(int k=1;k<N;k++) {
		for(int i=k+1;i<N+1;i++) {
			A[m(i,k,N)]/=A[m(k,k,N)];
			b[i-1]-=b[k-1]*A[m(i,k,N)];
			for(int j=k+1;j<N+1;j++) A[m(i,j,N)]-=A[m(k,j,N)]*A[m(i,k,N)];
		}
	}
	for(int k=N;k>0;k--) {
		b[k-1]/=A[m(k,k,N)];
		for(int i=k-1;i>0;i--) b[i-1]-=b[k-1]*A[m(i,k,N)];
	}
}

void gauss_pca (double *A,double *b,int N) {
	int ma;
	double t;
	for(int k=1;k<N;k++) { 
		ma=k;
		t=A[m(k,k,N)];
		for(int i=k+1;i<N+1;i++) {
			if(abs2(A[m(i,k,N)])>abs2(t)) {
				ma=i;
				t=A[m(i,k,N)];
			}
		}
		if(ma!=k) {
			for(int j=1;j<N+1;j++) {
				t=A[m(ma,j,N)];
				A[m(ma,j,N)]=A[m(k,j,N)];
				A[m(k,j,N)]=t;
			}
			t=b[ma-1];
			b[ma-1]=b[k-1];
			b[k-1]=t;			
		}
		for(int i=k+1;i<N+1;i++) {
			A[m(i,k,N)]/=A[m(k,k,N)];
			for(int j=k+1;j<N+1;j++) A[m(i,j,N)]-=A[m(k,j,N)]*A[m(i,k,N)];
		}
	}
	for(int i=1;i<N+1;i++) {
		for(int j=i+1;j<N+1;j++) b[j-1]-=b[i-1]*A[m(j,i,N)];
	}
	for(int k=N;k>0;k--) {
		b[k-1]/=A[m(k,k,N)];
		for(int i=k-1;i>0;i--) b[i-1]-=b[k-1]*A[m(i,k,N)];
	}
}

void sqrt_method(double *A,double *b,int N){
	for(int k=1;k<N+1;k++) {
		A[m(k,k,N)]=sqrt(A[m(k,k,N)]);
		for(int i=k+1;i<N+1;i++) A[m(i,k,N)]/=A[m(k,k,N)];
		for(int j=k+1;j<N+1;j++) {
			for(int i=j;i<N+1;i++) A[m(i,j,N)]-=A[m(i,k,N)]*A[m(j,k,N)];
		}
	}
	for(int i=1;i<N+1;i++) {
		b[i-1]/=A[m(i,i,N)];
		for(int j=i+1;j<N+1;j++) b[j-1]-=b[i-1]*A[m(j,i,N)];
	}
	for(int k=N;k>0;k--) {
		b[k-1]/=A[m(k,k,N)];
		for(int i=k-1;i>0;i--) b[i-1]-=b[k-1]*A[m(i,k,N)];
	}
}

void improved_sqrt_method (double *A,double *b,int N){
	double *v=new double[N];
	double t=0;
	v[0]=1;
	for(int j=1;j<N+1;j++) {
		for(int i=1;i<j;i++) v[i]=A[m(j,i,N)]*A[m(i,i,N)];
                t=0;
		for(int k=1;k<j;k++) t+=v[k]*A[m(j,k,N)]; 
		A[m(j,j,N)]-=t;
		for(int i=j+1;i<N+1;i++) {
			t=0;
			for(int k=1;k<j;k++) t+=v[k]*A[m(i,k,N)];
			A[m(i,j,N)]=(A[m(i,j,N)]-t)/A[m(j,j,N)];
		}
	}
	for(int i=1;i<N+1;i++) {
		for(int j=i+1;j<N+1;j++) b[j-1]-=b[i-1]*A[m(j,i,N)];
	}
	for(int i=0;i<N;i++) b[i]/=A[m(i+1,i+1,N)];
	for(int k=N;k>0;k--) {
		for(int i=k-1;i>0;i--) b[i-1]-=b[k-1]*A[m(k,i,N)];
	}
}

void reflect(double *w,double *a,int N) {
	double k=0;
	k=dot(w,a,N)/norm2(w,N);
	for(int i=0;i<N;i++) a[i]-=2*k*w[i];
	return;
}

void householder(double *a,int N){
	double n=0;
	n=abs2(a[0])-sqrt(norm2(a,N));
	if(n==0) return;
	n=a[0]-sqrt(norm2(a,N));
	a[0]=sqrt(norm2(a,N));
	for(int i=1;i<N;i++) a[i]/=n;
	return;
}

void qr(double *A,int M,int N){
	int f=0,cl=1,s;   
	while(cl<=N && (cl<M || f+1<cl)){
		double temp=0;
		while(norm2(A+m(f+2,cl,M),M-f-1)==0) cl++;
		householder(A+m(f+1,cl,M),M-f);
		temp=A[m(f+1,cl,M)];
		A[m(f+1,cl,M)]=1;
		for(int i=cl+1;i<=N;i++) reflect(A+m(f+1,cl,M),A+m(f+1,i,M),M-f);
		A[m(f+1,cl,M)]=temp;
		f++;
		cl++;
	}
}

void qrsolve(double *A,double *b,int N){            //only for A is non-singular;
	qr(A,N,N);
	for(int i=0;i<N-1;i++){
		double temp=0;
		if(norm2(A+m(i+2,i+1,N),N-i-1)==0) continue;
		temp=A[m(i+1,i+1,N)];
		A[m(i+1,i+1,N)]=1;
		reflect(A+m(i+1,i+1,N),b+i,N-i);
		A[m(i+1,i+1,N)]=temp;
	}
	for(int i=N;i>0;i--){
		b[i-1]/=A[m(i,i,N)];
		for(int j=i-1;j>0;j--) b[j-1]-=b[i-1]*A[m(j,i,N)];
	}
}

double ls(double *A,double *b,double *x,int M,int N){
	int c=0,k=0;
	int *r=new int[N];
	for(int i=0;i<N;i++) x[i]=0;
	qr(A,M,N);
	for(int i=0;i<N;i++){
		double temp=0;
		if(i==M-1){
			if(c<M-1) ;
			else{
				if(A[m(M,M,M)]!=0) {
					r[c]=M;
					c++;
					break;
				}
				else break;
			}
		}
		temp=A[m(c+1,i+1,M)];
		if(temp==0){
			x[i]=0;
			continue;
		}
		else{
			if(norm2(A+m(c+2,i+1,M),M-c-1)==0) ;
			else{
				A[m(c+1,i+1,M)]=1;
				reflect(A+m(c+1,i+1,M),b+c,M-c);
				A[m(c+1,i+1,M)]=temp;
			}
            r[c]=i+1;
			c++;
		}
	}
	for(int i=c-1;i>=0;i--){
		x[r[i]-1]=b[r[i]-1]/A[m(i+1,r[i],M)];
		for(int j=i-1;j>=0;j--) b[r[j]-1]-=x[r[i]-1]*A[m(j+1,r[i],M)];
	}
	delete []r;
	return sqrt(norm2(b+c+1,M-c-1));
}

element_b_v operator+(element_b_v& a,element_b_v& b){
	element_b_v r;
	int n=a.n;
	r.initial(n);
	for(int i=0;i<n;i++){
		r.b[i]=a.b[i]+b.b[i];
	}
	return r;
}

element_b_v operator-(element_b_v& a,element_b_v& b){
	element_b_v r;
	int n=a.n;
	r.initial(n);
	for(int i=0;i<n;i++) {
		r.b[i]=a.b[i]-b.b[i];
	}
	return r;
}

element_b_v operator*(element_b_m& A,element_b_v& b){
	element_b_v r;
	int n=b.n;
	if(A.n!=n) cout<<"error"<<endl;
	r.initial(n);
	for(int i=0;i<n;i++){
		r.b[i]=0;
		for(int j=0;j<n;j++){
			r.b[i]+=A.a[m(i+1,j+1,n)]*b.b[j];
		}
	}
	return r;
}

block_vector operator+(block_vector& a,block_vector& b){
	int N=a.m;
	int n=a.n;
	block_vector r;
	r.initial(N,n);
	for(int i=0;i<N;i++) r.b[i]=a.b[i]+b.b[i];
	return r;
}

block_vector operator-(block_vector& a,block_vector& b){
	int N=a.m;
	int n=a.n;
	block_vector r;
	r.initial(N,n);
	for(int i=0;i<N;i++) r.b[i]=a.b[i]-b.b[i];
	return r;
}

element_b_v operator*(double k,element_b_v& b){
	int n=b.n;
	element_b_v r;
	r.initial(n);
	for(int i=0;i<n;i++) r.b[i]=k*b.b[i];
	return r;
}

block_vector operator*(double k,block_vector& b){
	int m=b.m;
	int n=b.n;
	block_vector r;
	r.initial(m,n);
	for(int i=0;i<m;i++) r.b[i]=k*b.b[i];
	return r;
}

void b_jacobi(block_matrix& A,block_vector& b,block_vector& r){
	int N=A.M;
	int n=A.n;
	int step=0;
	block_vector x1,x2,t;
	element_b_v s;
	element_b_m temp;
	x1.initial(N,n);x2.initial(N,n);t.initial(N,n);s.initial(n);temp.initial(n);
	temp.a=new double[n*n];
	for(int i=0;i<N;i++) x2.b[i].rand2();
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
			temp=A.A[m(i+1,i+1,N)];
			gauss_pca(temp.a,x2.b[i].b,n);
		}
		t=x2-x1;
		step++;
	}while(step<Max && t.norm2()>pow(10,-12));

	r=x2;
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
	return;
}

void b_g_s(block_matrix& A,block_vector& b,block_vector& r){
	int N=A.M;
	int n=A.n;
	int step=0,end=0;
	block_vector x1,x2,t;
	element_b_v s;
	element_b_m temp;
	x1.initial(N,n);x2.initial(N,n);t.initial(N,n);	s.initial(n); temp.initial(n);
	temp.a=new double[n*n];
	for(int i=0;i<N;i++) x2.b[i].rand2();
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
			temp=A.A[m(i+1,i+1,N)];
			gauss_pca(temp.a,x2.b[i].b,n);
		}
		t=x2-x1;
		step++;
	}while(step<Max && t.norm2()>pow(10,-12));
	
	r=x2;
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
   	return;
}

void b_sor(block_matrix& A,block_vector& b,block_vector& r,double w){
	int N=A.M;
	int n=A.n;
	int step=0;
	block_vector x1,x2,t;
	element_b_v s;
	element_b_m temp;
	x1.initial(N,n); x2.initial(N,n); t.initial(N,n); s.initial(n); temp.initial(n);
	temp.a=new double[n*n];
	for(int i=0;i<N;i++) x2.b[i].rand2();
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
			temp=A.A[m(i+1,i+1,N)];
			gauss_pca(temp.a,x2.b[i].b,n);
		}
		t=x2-x1;
		step++;
	}while(step<Max && t.norm2()>pow(10,-12));
    
	r=x2;
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
	return;
}

void cg(block_matrix& A,block_vector& b,block_vector& r){
	int N=A.M;
	int n=A.n;
	int k=0;
	double w1,w2,alpha,beta;
	block_vector x1,t1,t2,s;
	element_b_v temp;
	x1.initial(N,n); t1.initial(N,n); t2.initial(N,n);s.initial(N,n);temp.initial(n);
	
	for(int i=0;i<N;i++) x1.b[i].rand2();
	t1=b;
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			if(A.A[m(i+1,j+1,N)].is_zero==0) continue;
			temp=A.A[m(i+1,j+1,N)]*x1.b[j];
			t1.b[i]=t1.b[i]-temp;
		}
	}
	w1=dot(t1,t1);

	while(k<n*n && w1>pow(10,-12)){
		k+=1;
		if(k==1)s=t1;
		else{
			beta=w1/w2;
			t2=beta*s;
			s=t1+t2;
		}
		for(int i=0;i<N;i++){
			t2.b[i].zeros();
			for(int j=0;j<N;j++){
				if(A.A[m(i+1,j+1,N)].is_zero==0) continue;
				temp=A.A[m(i+1,j+1,N)]*s.b[j];
				t2.b[i]=t2.b[i]+temp;
			}
		}
		alpha=w1/dot(s,t2);
		s=alpha*s;
		x1=x1+s;
		t2=alpha*t2;
		t1=t1-t2;
		w2=w1;
		w1=dot(t1,t1);
	}
	r=x1;
	cout<<"using conjecture-gradient-method, the L2_error is "<<sqrt(w1)<<endl;
	if(k==n*n) cout<<"fails to converge after "<<n*n<<" times"<<endl;
	else cout<<"totally iterates "<<k<<" times"<<endl;
	return;
}
