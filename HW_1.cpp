# include "n_algebra.h"

void hilbert(double *A,int N)
{
	for(int i=1;i<N+1;i++) {
		for(int j=1;j<N+1;j++) {
			A[m(i,j,N)]=1.0/(i+j-1);
		}
	}
	return;
}

void initial(double *A,int N)
{
	for(int i=1;i<=N;i++){
		A[m(i,i,N)]=1;
		A[m(i,N,N)]=1;
		for(int j=1;j<i;j++) A[m(i,j,N)]=-1;
		for(int j=i+1;j<N;j++) A[m(i,j,N)]=0;
	}	
	return;
}

void q1(){
	double t;
	ofstream s1("q1.txt");
	for(int i=5;i<31;i++) {
		double *A=new double[i*i];
		hilbert(A,i);
		t=condition_num_if(A,i);
		s1<<i<<' '<<t<<' '<<t/(pow(1+sqrt(2),4*i)/sqrt(i))<<endl;
		delete []A;
	}
    return;
}

void q2(){
    double t;
	ofstream s2("q2.txt");
    
	for(int i=5;i<31;i++) {
		double t,n1,n2,n3,n4,*A=new double[i*i],*A2=new double[i*i],*x0=new double[i],*x1=new double[i],*x2=new double[i];		
		s2<<i<<endl;
		initial(A,i);
		initial(A2,i);
		//hilbert(A,i);
		//hilbert(A2,i);

		t=condition_num_if(A,i);
		rand1(x0,i);
		s2<<"x0=";
		for(int k=0;k<i;k++) s2<<x0[k]<<' ';
		s2<<endl;
		for(int k=0;k<i;k++) {
			x1[k]=0;
			for(int j=0;j<i;j++) x1[k]+=(x0[j]*A[m(k+1,j+1,i)]);
		}
		s2<<"x1=Ax0=";
		for(int k=0;k<i;k++) {
			s2<<x1[k]<<' ';
			x2[k]=x1[k];
		}
		s2<<endl;
		gauss_pca(A,x2,i);
		s2<<"solve Ax'=x1"<<endl;
		s2<<"x'=";
                for(int k=0;k<i;k++) s2<<x2[k]<<' ';
		s2<<endl;
		n1=normif(x0,i);
		n2=normif(x1,i);
		s2<<"x0-x'=";
		for(int k=0;k<i;k++) {x0[k]-=x2[k]; s2<<x0[k]<<' ';}
		s2<<endl;
		n3=normif(x0,i);
		s2<<"x1-Ax'=";
		for(int k=0;k<i;k++) {
			for(int j=0;j<i;j++) x1[k]-=x2[j]*A2[m(k+1,j+1,i)];
			s2<<x1[k]<<' ';
		}
		s2<<endl;
		n4=normif(x1,i);
		s2<<"k="<<t<<" n(x0-x')/n(x0)="<<n3/n1<<' '<<"k*n(x1-Ax')/n(x1)="<<(t*n4)/n2<<endl;
	    delete []A; delete []A2; delete []x0; delete []x1; delete []x2;
	}
	return;
}

int main () {
	q1();
	q2();
	return 0;
}
