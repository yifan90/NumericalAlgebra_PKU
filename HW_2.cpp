# include "n_algebra.h"

void initial1(double *A,double *b,int N) {
	for(int i=0;i<N*N;i++) A[i]=0;
	for(int i=1;i<N;i++) {
		A[m(i,i,N)]=6;
		A[m(i+1,i,N)]=8;
		A[m(i,i+1,N)]=1;
		b[i]=15;
	}
	A[m(N,N,N)]=6;
	b[0]=7;
	b[N-1]=14;
}

void initial2(double *A,int N){
	for(int i=0;i<N*N;i++) A[i]=0;
	for(int i=1;i<N;i++) {
		A[m(i,i,N)]=10;
		A[m(i+1,i,N)]=1;
		A[m(i,i+1,N)]=1;
		//b[i]=12;
	}
	A[m(N,N,N)]=10;
	//b[0]=11;
	//b[N-1]=11;	
	//rand2(b,N);
}

void initial3(double *A,double *b,int N){
	for(int i=1;i<N+1;i++) {
		b[i-1]=0;
		for(int j=1;j<N+1;j++) {
			b[i-1]+=1.0/(i+j-1);
			A[m(i,j,N)]=1.0/(i+j-1);
		}
	}
}

void q1(){
	int N=100;
	double *A=new double[N*N];
	double *b=new double[N];
	ofstream s1("q1.txt");
	
	initial1(A,b,N);
	gauss(A,b,N);
	s1<<"gauss "<<N<<endl;
	for(int i=0;i<N;i++) s1<<b[i]<<endl;
	
	initial1(A,b,N);
	gauss_pca(A,b,N);
    s1<<endl<<"gauss_pca "<<N<<endl;
	for(int i=0;i<N;i++) s1<<b[i]<<endl;
	
    initial1(A,b,N);
	qrsolve(A,b,N);
	s1<<endl<<"qr "<<N<<endl;
	for(int i=0;i<N;i++) s1<<b[i]<<endl;
	
	delete []A;
	delete []b;
	return;
}

void q2(){
	int N=100;
	double *A=new double[N*N];
	double *x=new double[N];
	double *b=new double[N];
	ofstream s1("q2.txt");
	
	rand2(x,N);
	initial2(A,N);
	for(int i=0;i<N;i++){
		b[i]=0;
		for(int j=1;j<=N;j++) b[i]+=x[j-1]*A[m(i+1,j,N)];
	}
	gauss(A,b,N);
	s1<<"gauss "<<N<<endl;
	for(int i=0;i<N;i++) {
		s1<<b[i]<<endl;
		b[i]-=x[i];
	}
	s1<<"norm(x-x')="<<sqrt(norm2(b,N))<<endl;

	initial2(A,N);
	for(int i=0;i<N;i++){
		b[i]=0;
		for(int j=1;j<=N;j++) b[i]+=x[j-1]*A[m(i+1,j,N)];
	}
	gauss_pca(A,b,N);
    s1<<endl<<"gauss_pca "<<N<<endl;
    for(int i=0;i<N;i++) {
		s1<<b[i]<<endl;
		b[i]-=x[i];
	}
	s1<<"norm(x-x')="<<sqrt(norm2(b,N))<<endl;
	
	initial2(A,N);
	for(int i=0;i<N;i++){
		b[i]=0;
		for(int j=1;j<=N;j++) b[i]+=x[j-1]*A[m(i+1,j,N)];
	}
	sqrt_method(A,b,N);
	s1<<endl<<"sqrt_method "<<N<<endl;
    for(int i=0;i<N;i++) {
		s1<<b[i]<<endl;
		b[i]-=x[i];
	}
	s1<<"norm(x-x')="<<sqrt(norm2(b,N))<<endl;
	
	initial2(A,N);
	for(int i=0;i<N;i++){
		b[i]=0;
		for(int j=1;j<=N;j++) b[i]+=x[j-1]*A[m(i+1,j,N)];
	}
	improved_sqrt_method(A,b,N);
	s1<<endl<<"improved_sqrt_method "<<N<<endl;
	for(int i=0;i<N;i++) {
		s1<<b[i]<<endl;
		b[i]-=x[i];
	}
	s1<<"norm(x-x')="<<sqrt(norm2(b,N))<<endl;

    initial2(A,N);
	for(int i=0;i<N;i++){
		b[i]=0;
		for(int j=1;j<=N;j++) b[i]+=x[j-1]*A[m(i+1,j,N)];
	}
	qrsolve(A,b,N);
	s1<<endl<<"qr "<<N<<endl;
	for(int i=0;i<N;i++) {
		s1<<b[i]<<endl;
		b[i]-=x[i];
	}
	s1<<"norm(x-x')="<<sqrt(norm2(b,N))<<endl;
	
	delete []A;	delete []b; delete []x;
	return;
}

void q3(){
	int N=20;
	double *A=new double[N*N];
	double *b=new double[N];
	double *x=new double[N];
	ofstream s1("q3.txt");
	for(int i=0;i<N;i++) x[i]=1;

	initial3(A,b,N);
	gauss(A,b,N);
	s1<<"gauss "<<N<<endl;
	for(int i=0;i<N;i++) {
		s1<<b[i]<<endl;
		b[i]-=x[i];
	}
	s1<<"norm(x-x')="<<sqrt(norm2(b,N))<<endl;

	initial3(A,b,N);
	gauss_pca(A,b,N);
    s1<<endl<<"gauss_pca "<<N<<endl;
    for(int i=0;i<N;i++) {
		s1<<b[i]<<endl;
		b[i]-=x[i];
	}
	s1<<"norm(x-x')="<<sqrt(norm2(b,N))<<endl;
	
	initial3(A,b,N);
	sqrt_method(A,b,N);
	s1<<endl<<"sqrt_method "<<N<<endl;
    for(int i=0;i<N;i++) {
		s1<<b[i]<<endl;
		b[i]-=x[i];
	}
	s1<<"norm(x-x')="<<sqrt(norm2(b,N))<<endl;
	
	initial3(A,b,N);
	improved_sqrt_method(A,b,N);
	s1<<endl<<"improved_sqrt_method "<<N<<endl;
	for(int i=0;i<N;i++) {
		s1<<b[i]<<endl;
		b[i]-=x[i];
	}
	s1<<"norm(x-x')="<<sqrt(norm2(b,N))<<endl;

    initial3(A,b,N);
	qrsolve(A,b,N);
	s1<<endl<<"qr "<<N<<endl;
	for(int i=0;i<N;i++) {
		s1<<b[i]<<endl;
		b[i]-=x[i];
	}
	s1<<"norm(x-x')="<<sqrt(norm2(b,N))<<endl;
	
	delete []A;	delete []b; delete []x;
	return;
}

void q4(){
	double A[21],y2[7],x[3]={0};
	double y[7]={1,.8125,.75,1,1.3125,1.75,2.3125};
	double t[7]={-1,-.75,-.5,0,.25,.5,.75};
	double b=0;
	ofstream s("q4.txt");
	for(int i=0;i<7;i++) {
		A[m(i+1,1,7)]=t[i]*t[i];
        A[m(i+1,2,7)]=t[i];
        A[m(i+1,3,7)]=1;
		y2[i]=y[i];
	}
	b=ls(A,y,x,7,3);
	s<<"min{norm(Ax-b)}="<<b<<endl;
	s<<"with a="<<x[0]<<" b="<<x[1]<<" c="<<x[2]<<endl;
	s<<"and now Ax-b="<<endl;
	for(int i=0;i<7;i++) s<<x[0]*t[i]*t[i]+x[1]*t[i]+x[2]-y2[i]<<endl;
	s.close();
	return;
}

void q5(){
	double A[308],A2[308],y[28],y2[28],x[11]={0};
	double b=0;
	ifstream in("data.txt");
	ofstream s("q5.txt");
	for(int i=0;i<28;i++) {
		in>>y[i];
		y2[i]=y[i];
		for(int k=1;k<9;k++) in>>A[m(i+1,k,28)];
        in>>A[m(i+1,11,28)]>>A[m(i+1,9,28)]>>A[m(i+1,10,28)];
	}
	for(int i=0;i<308;i++) A2[i]=A[i];
	b=ls(A,y,x,28,11);
	s<<"min{norm(Ax-b)}="<<b<<endl;
	s<<"x="<<endl;
	for(int k=0;k<11;k++) s<<x[k]<<endl;
	s<<"and now Ax-b="<<endl;
	for(int j=0;j<28;j++) {
		b=0;
		for(int k=0;k<11;k++) b+=(x[k]*A2[m(j+1,k+1,28)]);
		s<<b-y2[j]<<endl;
	}
	in.close(); s.close();
	return;
}

int main(){
	q1();
	q2();
	q3();
	q4();
	q5();
	return 0;
}
