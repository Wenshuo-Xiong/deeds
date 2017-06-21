/* several functions to interpolate and symmetrise deformations (as well as make them diffeomorphic)
 calculates Jacobian and harmonic Energy       几个用于插值和对称（同时也微分同胚）的函数                              */

template <typename TypeI>

void interp3(TypeI* interp , TypeI* input , float* x1 , float* y1 , float* z1 , int m , int n , int o , int m2 , int n2 , int o2 , bool flag)	   //3维图像的双线性插值法
				
					//注：只实现了其中一维的插值，对于三个坐标要用三次此程序

			  //参数：interp代表插值后灰度序列，input代表原灰度序列，x1,y1,z1代表变换后的图像在原图像的xyz坐标网格序列 （所以在上采样时要先把内部按比例填充取值）


{

	for(int k=0;k<o;k++){
		for(int j=0;j<n;j++){
			for(int i=0;i<m;i++){

				int x=floor(x1[i+j*m+k*m*n]); int y=floor(y1[i+j*m+k*m*n]);  int z=floor(z1[i+j*m+k*m*n]); 	  //当用一维数组存储3维数据时，计算都是[i+j*m+k*m*n];
				float dx=x1[i+j*m+k*m*n]-x; float dy=y1[i+j*m+k*m*n]-y; float dz=z1[i+j*m+k*m*n]-z;			  //x,y,z为原图的整数项，dx,dy,dz为原图的小数项
				
				if(flag){
					x+=j; y+=i; z+=k;
				}
				//双线性插值法在3维图像中的应用：(1-△x)(1-△y)*f(i,j)+△x(1-△y)*f(i,j+1)+(1-△x)△y*f(i+1,j)+△x△y*f(i+1,j+1)

//在flag为1时，下式可改写为：interp[i+j*m+k*m*n]=input[(x+y*m+z*m*n)+(i+j*m+k*m*n)]--->V1[i+j*m+k*m*n]=-u[(x+y*m+z*m*n)+(i+j*m+k*m*n)]--->V1(x)=-u(0+x){假设还不需要用插值}
					//当再一次循环时，V2(x)=-u(x+V1(x))=-u(x+ -u(x))

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

void filter1(float* imagein,float* imageout,int m,int n,int o,float* filter,int length,int dim)	   //又是把3维滤波分为三个维单独滤波（4个for：复杂度突破天际）
{
	int i,j,k,f;
	int i1,j1,k1;
	int hw=(length-1)/2;

	//imagein: 输入图像 ；  imageout：输出图像； m，n，o：图像维度  ；filter：滤波窗   ；length：    ；dim：维度
	
	for(i=0;i<(m*n*o);i++){
		imageout[i]=0.0;
	}
	
	for(k=0;k<o;k++){
		for(j=0;j<n;j++){
			for(i=0;i<m;i++){
				for(f=0;f<length;f++){
					//replicate-padding	 重复填充				   //max（min）的嵌套使用：可以保证第1个数卡在第2，3个数的定义域内
					if(dim==1)
						imageout[i+j*m+k*m*n]+=filter[f]*imagein[max(min(i+f-hw,m-1),0)+j*m+k*m*n]; 	   //对于每一个像素点i+j*m+k*m*n:用一个条形滤波器filter来滤每一维度
					if(dim==2)																			   //边缘保证不越界（边缘处理很疏忽）总共下来三维滤波要用三次。
						imageout[i+j*m+k*m*n]+=filter[f]*imagein[i+max(min(j+f-hw,n-1),0)*m+k*m*n]; 
					if(dim==3)
						imageout[i+j*m+k*m*n]+=filter[f]*imagein[i+j*m+max(min(k+f-hw,o-1),0)*m*n]; 
				}
			}
		}
	}
}


void upsampleDeformations2(float* u1,float* v1,float* w1,float* u0,float* v0,float* w0,int m,int n,int o,int m2,int n2,int o2)//上采样：即插值（避免重复计算，但损失了精度）
{
	 //u1,v1,w1:上采样后偏差； u0,v0,w0:上采样前偏差； m,n,o:后图的尺寸； m2,n2,o2:前图的尺寸
	
	float scale_m=(float)m/(float)m2;			 //计算采样比例
	float scale_n=(float)n/(float)n2;
	float scale_o=(float)o/(float)o2;
	
	float* x1=new float[m*n*o];				   //x1,y1,z1 大小为后图尺寸的坐标序列 
	float* y1=new float[m*n*o];
	float* z1=new float[m*n*o];

	for(int k=0;k<o;k++){					  //坐标变换
		for(int j=0;j<n;j++){
			for(int i=0;i<m;i++){

				x1[i+j*m+k*m*n]=j/scale_n;
				y1[i+j*m+k*m*n]=i/scale_m;
				z1[i+j*m+k*m*n]=k/scale_o;
			}
		}
	}
	
	interp3(u1,u0,x1,y1,z1,m,n,o,m2,n2,o2,false);
	interp3(v1,v0,x1,y1,z1,m,n,o,m2,n2,o2,false);
	interp3(w1,w0,x1,y1,z1,m,n,o,m2,n2,o2,false);
	
	delete []x1;
	delete []y1;
	delete []z1;

		
}


void fastInverse(float* ui,float* vi,float* wi,float* u,float* v,float* w,int m,int n,int o)  //利用论文思路快速计算逆矩阵用于之后的对称变换的计算
{
	//输入参数：ui,vi,wi:逆变换； u,v,w:原图灰度； m,n,o：图像维度

	float* un=new float[m*n*o];
	float* vn=new float[m*n*o];
	float* wn=new float[m*n*o];
	float* uin=new float[m*n*o];
	float* vin=new float[m*n*o];
	float* win=new float[m*n*o];

	for(int i=0;i<m*n*o;i++){
		uin[i]=-u[i];
		vin[i]=-v[i];
		win[i]=-w[i];
		ui[i]=0.0;
		vi[i]=0.0;
		wi[i]=0.0;
	}
	for(int it=0;it<10;it++)			//Vn(x)=-u[x+Vn-1(x)]迭代更新：u(x)为前向形变场，v(x)为反向形变场。
	{
		interp3(un,uin,ui,vi,wi,m,n,o,m,n,o,true);					//ui,vi,wi为论文中的Vn-1(x)用于输入计算;un,vn,wn为输出的坐标序列Vn(x);uin,vin,win为u(x)映射
		interp3(vn,vin,ui,vi,wi,m,n,o,m,n,o,true);
		interp3(wn,win,ui,vi,wi,m,n,o,m,n,o,true);
		for(int i=0;i<m*n*o;i++){
			ui[i]=un[i];
			vi[i]=vn[i];
			wi[i]=wn[i];

		}
	}
	delete []uin;
	delete []vin;	
	delete []win;
	delete []un;
	delete []vn;	
	delete []wn;
}


void combineDeformation(float* u3,float* v3,float* w3,  float* u1,float* v1,float* w1,  float* u2,float* v2,float* w2,int m,int n,int o)
{
	//u3,v3,w3用于接收最后的形变场参数	u1,v1,w1是传入参量   u2,v2,w2是

	float* uc=new float[m*n*o];
	float* vc=new float[m*n*o];
	float* wc=new float[m*n*o];

	interp3(uc,u1,  u2,v2,w2,m,n,o,m,n,o,true);
	interp3(vc,v1,  u2,v2,w2,m,n,o,m,n,o,true);
	interp3(wc,w1,  u2,v2,w2,m,n,o,m,n,o,true);

	for(int i=0;i<m*n*o;i++){
		u3[i]=uc[i]+u2[i];
		v3[i]=vc[i]+v2[i];
		w3[i]=wc[i]+w2[i];

	}
	delete []uc;
	delete []vc;
	delete []wc;

}

void diffeomorphic(float* u1,float* v1,float *w1,int m,int n,int o,int expsteps,int factor)
{
	float* u1b=new float[m*n*o];
	float* v1b=new float[m*n*o];
	float* w1b=new float[m*n*o];

	float* u2=new float[m*n*o];
	float* v2=new float[m*n*o];
	float* w2=new float[m*n*o];
	
	float factor1=1.0/(float)factor;
	float coeff=1.0/(float)pow(2.0,expsteps);
	for(int i=0;i<m*n*o;i++){
		u1b[i]=coeff*u1[i]*factor1;
		v1b[i]=coeff*v1[i]*factor1;
		w1b[i]=coeff*w1[i]*factor1;

	}
	for(int it=0;it<expsteps;it++){
		combineDeformation(u2,v2,w2,u1b,v1b,w1b,u1b,v1b,w1b,m,n,o);
		for(int i=0;i<m*n*o;i++){
			u1b[i]=u2[i];
			v1b[i]=v2[i];
			w1b[i]=w2[i];
		}
	}
	for(int i=0;i<m*n*o;i++){
		u1[i]=u2[i]*(float)factor;
		v1[i]=v2[i]*(float)factor;
		w1[i]=w2[i]*(float)factor;

	}
	delete []u2;
	delete []v2;
	delete []w2;
	delete []u1b;
	delete []v1b;
	delete []w1b;

}


void symmetricMapping(float* u,float* v,float* w,float* u2,float* v2,float* w2,int m,int n,int o,int factor)
{
	float factor1=1.0/(float)factor;
	float* usym=new float[m*n*o];		//前向：对于每一维的形变场都要生成一个对称序列
	float* vsym=new float[m*n*o];
	float* wsym=new float[m*n*o];

	float* usym2=new float[m*n*o];		 //后向：对于每一维的形变场都要生成一个对称序列
	float* vsym2=new float[m*n*o];
	float* wsym2=new float[m*n*o];

	for(int i=0;i<m*n*o;i++){
		usym[i]=u[i]*(factor1*0.5);		 //u,v,w是前向形变场，乘一个因数存入sym中
		vsym[i]=v[i]*(factor1*0.5);
		wsym[i]=w[i]*(factor1*0.5);
		usym2[i]=u2[i]*(factor1*0.5);
		vsym2[i]=v2[i]*(factor1*0.5);
		wsym2[i]=w2[i]*(factor1*0.5);

	}
	float* ui=new float[m*n*o];
	float* vi=new float[m*n*o];
	float* wi=new float[m*n*o];
	float* u2i=new float[m*n*o];
	float* v2i=new float[m*n*o];
	float* w2i=new float[m*n*o];
	fastInverse(ui,vi,wi,usym,vsym,wsym,m,n,o);
	fastInverse(u2i,v2i,w2i,usym2,vsym2,wsym2,m,n,o);	   //到此，前向变换的半逆和后向变换的半逆都已算出（即[A(0.5)]-1 和[B(0.5)]-1）
	
	combineDeformation(u,v,w,u2i,v2i,w2i,usym,vsym,wsym,m,n,o);		//A(0.5)o B(0.5)-1	 -->As
	combineDeformation(u2,v2,w2,ui,vi,wi,usym2,vsym2,wsym2,m,n,o);	//B(0.5)o A(0.5)-1	 -->Bs

//	combineDeformation(u,v,w,usym,vsym,wsym,u2i,v2i,w2i,m,n,o); //the other way around
//	combineDeformation(u2,v2,w2,usym2,vsym2,wsym2,ui,vi,wi,m,n,o);
	for(int i=0;i<m*n*o;i++)									  //这个factor到底是干什么的，算之前除了进来，算之后又乘了回去
	{															//	（factor是论文中的搜索半径{8，6，4，2...}
		u[i]*=(float)factor;
		v[i]*=(float)factor;
		w[i]*=(float)factor;
		u2[i]*=(float)factor;
		v2[i]*=(float)factor;	
		w2[i]*=(float)factor;
	}
	
	delete []ui;
	delete []vi;
	delete []wi;
	delete []u2i;
	delete []v2i;
	delete []w2i;
	delete []usym;
	delete []vsym;
	delete []wsym;
	delete []usym2;
	delete []vsym2;
	delete []wsym2;
	
}

void consistentMapping(float* u,float* v,float* w,float* u2,float* v2,float* w2,int m,int n,int o,int factor)  //另一种保持逆一致性的方法
{
	float factor1=1.0/(float)factor;
	float* us=new float[m*n*o];
	float* vs=new float[m*n*o];
	float* ws=new float[m*n*o];
	float* us2=new float[m*n*o];
	float* vs2=new float[m*n*o];
	float* ws2=new float[m*n*o];
    
	for(int i=0;i<m*n*o;i++){
		us[i]=u[i]*factor1; vs[i]=v[i]*factor1;	ws[i]=w[i]*factor1;
		us2[i]=u2[i]*factor1; vs2[i]=v2[i]*factor1;	ws2[i]=w2[i]*factor1;
	}
    
    for(int it=0;it<10;it++){
        interp3(u,us2,us,vs,ws,m,n,o,m,n,o,true);
        interp3(v,vs2,us,vs,ws,m,n,o,m,n,o,true);
        interp3(w,ws2,us,vs,ws,m,n,o,m,n,o,true);
        for(int i=0;i<m*n*o;i++){
            u[i]=0.5*us[i]-0.5*u[i];
            v[i]=0.5*vs[i]-0.5*v[i];
            w[i]=0.5*ws[i]-0.5*w[i];
            
        }
        interp3(u2,us,us2,vs2,ws2,m,n,o,m,n,o,true);
        interp3(v2,vs,us2,vs2,ws2,m,n,o,m,n,o,true);
        interp3(w2,ws,us2,vs2,ws2,m,n,o,m,n,o,true);
        for(int i=0;i<m*n*o;i++){
            u2[i]=0.5*us2[i]-0.5*u2[i];
            v2[i]=0.5*vs2[i]-0.5*v2[i];
            w2[i]=0.5*ws2[i]-0.5*w2[i];
        }
        
        for(int i=0;i<m*n*o;i++){
            us[i]=u[i]; vs[i]=v[i]; ws[i]=w[i];
            us2[i]=u2[i]; vs2[i]=v2[i]; ws2[i]=w2[i];
        }
        
    }
    
    
	for(int i=0;i<m*n*o;i++){
		u[i]*=(float)factor;
		v[i]*=(float)factor;
		w[i]*=(float)factor;
		u2[i]*=(float)factor;
		v2[i]*=(float)factor;
		w2[i]*=(float)factor;
	}
	
    
    delete us; delete vs; delete ws;
    delete us2; delete vs2; delete ws2;
}

float harmonicEnergy(float* u,float* v,float* w,int m,int n,int o){
	int sz=m*n*o;
	float grad[3]={-0.5,0.0,0.5};
	float* output=new float[sz];
	
	float energy=0.0;
	filter1(u,output,m,n,o,grad,3,1);
	for(int i=0;i<sz;i++){
		energy+=pow(output[i],2.0);
	}
	filter1(u,output,m,n,o,grad,3,2);
	for(int i=0;i<sz;i++){
		energy+=pow(output[i],2.0);
	}
	filter1(u,output,m,n,o,grad,3,3);
	for(int i=0;i<sz;i++){
		energy+=pow(output[i],2.0);
	}
	filter1(v,output,m,n,o,grad,3,1);
	for(int i=0;i<sz;i++){
		energy+=pow(output[i],2.0);
	}
	filter1(v,output,m,n,o,grad,3,2);
	for(int i=0;i<sz;i++){
		energy+=pow(output[i],2.0);
	}
	filter1(v,output,m,n,o,grad,3,3);
	for(int i=0;i<sz;i++){
		energy+=pow(output[i],2.0);
	}
	filter1(w,output,m,n,o,grad,3,1);
	for(int i=0;i<sz;i++){
		energy+=pow(output[i],2.0);
	}
	filter1(w,output,m,n,o,grad,3,2);
	for(int i=0;i<sz;i++){
		energy+=pow(output[i],2.0);
	}
	filter1(w,output,m,n,o,grad,3,3);
	for(int i=0;i<sz;i++){
		energy+=pow(output[i],2.0);
	}
	energy/=(float)(sz);
	return energy;
}




float jacobian(float* u1,float* v1,float* w1,int m,int n,int o,int factor){
	
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
	
	filter1(u1,J11,m,n,o,grad,3,2);			//输入图像:u1; 输出:J11; m,n,o; 滤波器:grad; 长度:3
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
	cout<<"Jacobian of deformations| Mean (std): "<<round(jmean*1000)/1000.0<<" ("<<round(jstd*1000)/1000.0<<")\n";
	cout<<"Range: ["<<Jmin<<", "<<Jmax<<"] Negative fraction: "<<frac<<"\n";
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




