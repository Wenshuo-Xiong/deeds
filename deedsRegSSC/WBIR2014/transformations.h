
/* several functions to interpolate deformations and
 calculate Jacobians */
template <typename TypeI>

void interp3(TypeI* interp,TypeI* input,float* x1,float* y1,float* z1,int m,int n,int o,int m2,int n2,int o2,bool flag)	   //PASS 跟2013的一样
{
	for(int k=0;k<o;k++){
		for(int j=0;j<n;j++){
			for(int i=0;i<m;i++){
				int x=floor(x1[i+j*m+k*m*n]); int y=floor(y1[i+j*m+k*m*n]);  int z=floor(z1[i+j*m+k*m*n]); 
				float dx=x1[i+j*m+k*m*n]-x; float dy=y1[i+j*m+k*m*n]-y; float dz=z1[i+j*m+k*m*n]-z;
				
				if(flag){
					x+=j; y+=i; z+=k;
				}
				interp[i+j*m+k*m*n]=(1.0-dx)*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
				(1.0-dx)*dy*(1.0-dz)*input[min(max(y+1,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
				dx*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
				(1.0-dx)*(1.0-dy)*dz*input[min(max(y,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
				dx*dy*(1.0-dz)*input[min(max(y+1,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
				(1.0-dx)*dy*dz*input[min(max(y+1,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
				dx*(1.0-dy)*dz*input[min(max(y,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
				dx*dy*dz*input[min(max(y+1,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2];
			}
		}
	}
}

void filter1(float* imagein,float* imageout,int m,int n,int o,float* filter,int length,int dim)			 //PASS 一样
{
	int i,j,k,f;
	int i1,j1,k1;
	int hw=(length-1)/2;
	
	for(i=0;i<(m*n*o);i++){
		imageout[i]=0.0;
	}
	
	for(k=0;k<o;k++){
		for(j=0;j<n;j++){
			for(i=0;i<m;i++){
				for(f=0;f<length;f++){
					//replicate-padding
					if(dim==1)
						imageout[i+j*m+k*m*n]+=filter[f]*imagein[max(min(i+f-hw,m-1),0)+j*m+k*m*n]; 
					if(dim==2)
						imageout[i+j*m+k*m*n]+=filter[f]*imagein[i+max(min(j+f-hw,n-1),0)*m+k*m*n]; 
					if(dim==3)
						imageout[i+j*m+k*m*n]+=filter[f]*imagein[i+j*m+max(min(k+f-hw,o-1),0)*m*n]; 
				}
			}
		}
	}
}

void volfilter(float* imagein,int m,int n,int o,int length,float sigma)		//PASS	在Matlab中计算MIND时用过，用于1.简便计算邻域的灰度差（利用卷积）2.高斯滤波
{
	
	int hw=(length-1)/2;
	int i,j,f;
	float hsum=0;
	float* filter=new float[length];
	for(i=0;i<length;i++){
		filter[i]=exp(-pow((i-hw),2)/(2*pow(sigma,2)));
		hsum=hsum+filter[i];
	}
	for(i=0;i<length;i++){
		filter[i]=filter[i]/hsum;				//这两步在Matlab中就是用fspecial创建一个高斯滤波（一维）窗 。
	}
	float* image1=new float[m*n*o];
    for(i=0;i<m*n*o;i++){
        image1[i]=imagein[i];
    }
    filter1(image1,imagein,m,n,o,filter,length,1);
	filter1(imagein,image1,m,n,o,filter,length,2);
	filter1(image1,imagein,m,n,o,filter,length,3);			//迭代三次（每一维单独滤一遍）

	delete []image1;
	delete []filter;	
	
}

void upsampleDeformations2scale(float* u1,float* v1,float* w1,float* u0,float* v0,float* w0,int m,int n,int o,int m2,int n2,int o2)

{
			  //u1,v1,w1为输出形变场    u0,v0,w0为输入形变场    scale为上采样的比例（3，2，1）
	
	float scale_m=(float)m/(float)m2;
	float scale_n=(float)n/(float)n2;
	float scale_o=(float)o/(float)o2;
	
	float* x1=new float[m*n*o];
	float* y1=new float[m*n*o];
	float* z1=new float[m*n*o];
	for(int k=0;k<o;k++){
		for(int j=0;j<n;j++){
			for(int i=0;i<m;i++){
				x1[i+j*m+k*m*n]=j/scale_n;
				y1[i+j*m+k*m*n]=i/scale_m;
				z1[i+j*m+k*m*n]=k/scale_o;
			}
		}
	}
	
	interp3(u1,u0,x1,y1,z1,m,n,o,m2,n2,o2,false);		  //在flag为false的情况下，就仅仅是个插值函数，按比例上采样插值（为了减少计算量）
	interp3(v1,v0,x1,y1,z1,m,n,o,m2,n2,o2,false);
	interp3(w1,w0,x1,y1,z1,m,n,o,m2,n2,o2,false);
	
	delete []x1;
	delete []y1;
	delete []z1;
	
	for(int i=0;i<m*n*o;i++){					 //最后乘了一个比例系数？？
		u1[i]*=scale_n;
		v1[i]*=scale_m;
		w1[i]*=scale_o;
	}
    
}

template <typename TypeW>

void warpImage(TypeW* warped,TypeW* im1,float* u1,float* v1,float* w1)	  //PASS  最简单的x+u(x)
{
    int m=image_m; int n=image_n; int o=image_o; int sz=m*n*o;
    interp3(warped,im1,u1,v1,w1,m,n,o,m,n,o,true);
}

float jacobian(float* u1,float* v1,float* w1,int m,int n,int o,int factor)	 //跟MST一样，还没看
{
	
	float factor1=1.0/(float)factor;
	float jmean=0.0;
	float jstd=0.0;
	int i;
	float grad[3]={-0.5,0.0,0.5};
	float* Jac=new float[m*n*o];
	
	float* J11=new float[m*n*o];
	float* J12=new float[m*n*o];
	float* J13=new float[m*n*o];
	float* J21=new float[m*n*o];
	float* J22=new float[m*n*o];
	float* J23=new float[m*n*o];
	float* J31=new float[m*n*o];
	float* J32=new float[m*n*o];
	float* J33=new float[m*n*o];
	
	for(i=0;i<(m*n*o);i++){
		J11[i]=0.0;
		J12[i]=0.0;
		J13[i]=0.0;
		J21[i]=0.0;
		J22[i]=0.0;
		J23[i]=0.0;
		J31[i]=0.0;
		J32[i]=0.0;
		J33[i]=0.0;
	}
	
	float neg=0; float Jmin=1; float Jmax=1; float J;
	float count=0; float frac;
	
	filter1(u1,J11,m,n,o,grad,3,2);
	filter1(u1,J12,m,n,o,grad,3,1);
	filter1(u1,J13,m,n,o,grad,3,3);
	
	filter1(v1,J21,m,n,o,grad,3,2);
	filter1(v1,J22,m,n,o,grad,3,1);
	filter1(v1,J23,m,n,o,grad,3,3);
	
	filter1(w1,J31,m,n,o,grad,3,2);
	filter1(w1,J32,m,n,o,grad,3,1);
	filter1(w1,J33,m,n,o,grad,3,3);
	
	for(i=0;i<(m*n*o);i++){
		J11[i]*=factor1;
		J12[i]*=factor1;
		J13[i]*=factor1;
		J21[i]*=factor1;
		J22[i]*=factor1;
		J23[i]*=factor1;
		J31[i]*=factor1;
		J32[i]*=factor1;
		J33[i]*=factor1;
	}	
	
	for(i=0;i<(m*n*o);i++){
		J11[i]+=1.0;
		J22[i]+=1.0;
		J33[i]+=1.0;
	}
	for(i=0;i<(m*n*o);i++){
		J=J11[i]*(J22[i]*J33[i]-J23[i]*J32[i])-J21[i]*(J12[i]*J33[i]-J13[i]*J32[i])+J31[i]*(J12[i]*J23[i]-J13[i]*J22[i]);
		jmean+=J;
		if(J>Jmax)
			Jmax=J;
		if(J<Jmin)
			Jmin=J;
		if(J<0)
			neg++;
		count++;
		Jac[i]=J;
	}
	jmean/=(m*n*o);
	for(int i=0;i<m*n*o;i++){
		jstd+=pow(Jac[i]-jmean,2.0);
	}
	jstd/=(m*n*o-1);
	jstd=sqrt(jstd);
	frac=neg/count;
	cout<<"std(J)="<<round(jstd*100)/100.0;
	//cout<<"Range: ["<<Jmin<<", "<<Jmax<<"]round(jmean*100)/100.0<<
    cout<<" (J<0)="<<round(frac*1e7)/100.0<<"e-7  ";
	delete []Jac;
	
	
	delete []J11;
	delete []J12;
	delete []J13;
	delete []J21;
	delete []J22;
	delete []J23;
	delete []J31;
	delete []J32;
	delete []J33;
	
	return jstd;
	
	
}




