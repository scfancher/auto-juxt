#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "mersenne_twister.hpp"

#define sn 2
#define N 1000
#define nn 100
#define sl 1.0
#define fl 1000.0
#define nl 100
#define D .1
#define bl 50

double SyMat(int vecM[N*N], int vecsx[N], int vecsy[N], int vecsz[N], const int ns, const double lam, const int NS)
{
	double *M = new double[ns*ns];
	double *Minv = new double[ns*ns];
	int i,j,k,tempM;
	double R,Sy;

	for(i = 0; i < ns; i++)
	{
		for(j = 0; j < ns; j++)
		{
			if(i==j)
			{
				tempM = vecM[i*NS+j];
				for(k = ns; k < NS; k++)
				{
					tempM += vecM[i*NS+k];
				}
				M[i*ns+j] = 1.0+tempM*lam*lam/4.0;
				Minv[i*ns+j] = 1.0;
			}
			else
			{
				M[i*ns+j] = vecM[i*NS+j]*lam*lam/4.0;
				Minv[i*ns+j] = 0.0;
			}
		}
	}

	for(i = ns-1; i > 0; i--)
	{
		for(j = i-1; j >= 0; j--)
		{
			if(M[j*ns+i]!=0.0)
			{
				R = M[j*ns+i]/M[i*ns+i];
				M[j*ns+i] = 0.0;
				for(k = 0; k < ns; k++)
				{
					Minv[j*ns+k] -= Minv[i*ns+k]*R;
					if(k<i)
					{
						M[j*ns+k] -= M[i*ns+k]*R;
					}
				}
			}
		}
	}

	for(i = 0; i < ns; i++)
	{
		Minv[i] /= M[0];
	}

	Sy = 0;
	for(i = 0; i < ns; i++)
	{
		for(j = 0; j < ns; j++)
		{
			if(i==j)
			{
				Sy += Minv[i]*Minv[j]/2.0;
			}
			else
			{
				R = sqrt((vecsx[i]-vecsx[j])*(vecsx[i]-vecsx[j])+(vecsy[i]-vecsy[j])*(vecsy[i]-vecsy[j])+(vecsz[i]-vecsz[j])*(vecsz[i]-vecsz[j]));
				Sy += Minv[i]*Minv[j]/(2.0*R);
			}
		}
	}

	delete[] M;
	delete[] Minv;

	return Sy;
}

double Lfun(double r, double lam)
{
	double L;

	if(r==0.0)
	{
		L = (2.0*lam*lam)*(1.0-(1.0+1.0/lam)*exp(-1/lam));
	}
	else if(r<1.0)
	{
		L = (2.0*lam*lam*lam/r)*(r/lam-(1.0+1.0/lam)*exp(-1/lam)*sinh(r/lam));
	}
	else
	{
		L = (2.0*lam*lam*lam/r)*exp(-r/lam)*((1.0/lam)*cosh(1.0/lam)-sinh(1.0/lam));
	}

	return L;
}

double Sz(double vecsx[N], double vecsy[N], double vecsz[N], const int ns, double lam)
{
	double Ss,S1,S2,S3,r1,r2,r12,L1,lfact;
	int i,j;
	S1 = 0.0;
	S2 = 0.0;
	S3 = 0.0;
	
	for(i = 0; i < ns-1; i++)
	{
		r1 = sqrt((vecsx[i]-vecsx[0])*(vecsx[i]-vecsx[0])+(vecsy[i]-vecsy[0])*(vecsy[i]-vecsy[0])+(vecsz[i]-vecsz[0])*(vecsz[i]-vecsz[0]));
		L1 = Lfun(r1,lam);

		for(j = i+1; j < ns; j++)
		{
			r2 = sqrt((vecsx[j]-vecsx[0])*(vecsx[j]-vecsx[0])+(vecsy[j]-vecsy[0])*(vecsy[j]-vecsy[0])+(vecsz[j]-vecsz[0])*(vecsz[j]-vecsz[0]));
			r12 = sqrt((vecsx[i]-vecsx[j])*(vecsx[i]-vecsx[j])+(vecsy[i]-vecsy[j])*(vecsy[i]-vecsy[j])+(vecsz[i]-vecsz[j])*(vecsz[i]-vecsz[j]));
			if(r12<1.0)
			{
				lfact = 1.0;
			}
			else
			{
				lfact = 1.0/r12;
			}
			S1 += L1*Lfun(r2,lam)*lfact;
		}
		
		S2 += L1*L1/2.0;
		S3 += L1;
	}
	r1 = sqrt((vecsx[ns-1]-vecsx[0])*(vecsx[ns-1]-vecsx[0])+(vecsy[ns-1]-vecsy[0])*(vecsy[ns-1]-vecsy[0])+(vecsz[ns-1]-vecsz[0])*(vecsz[ns-1]-vecsz[0]));
	L1 = Lfun(r1,lam);
	S2 += L1*L1/2.0;
	S3 += L1;	

	Ss = (S1+S2)/(S3*S3);

	return Ss;
}

double AveNN(double vecsx[N], double vecsy[N], double vecsz[N], const int ns)
{
	double *nn = new double[ns];
	double ann,d;
	int i,j;
	ann = 0;

	d = sqrt((vecsx[0]-vecsx[1])*(vecsx[0]-vecsx[1])+(vecsy[0]-vecsy[1])*(vecsy[0]-vecsy[1])+(vecsz[0]-vecsz[1])*(vecsz[0]-vecsz[1]));
	nn[0] = d;
	nn[1] = d;

	if(ns>2)
	{
		for(i = 2; i < ns; i++)
		{
			d = sqrt((vecsx[0]-vecsx[i])*(vecsx[0]-vecsx[i])+(vecsy[0]-vecsy[i])*(vecsy[0]-vecsy[i])+(vecsz[0]-vecsz[i])*(vecsz[0]-vecsz[i]));
			if(d<nn[0])
			{
				nn[0] = d;
			}
			nn[i] = d;
		}
		for(i = 1; i < ns-1; i++)
		{
			for(j = i+1; j < ns; j++)
			{
				d = sqrt((vecsx[i]-vecsx[j])*(vecsx[i]-vecsx[j])+(vecsy[i]-vecsy[j])*(vecsy[i]-vecsy[j])+(vecsz[i]-vecsz[j])*(vecsz[i]-vecsz[j]));
				if(d<nn[i])
				{
					nn[i] = d;
				}
				if(d<nn[j])
				{
					nn[j] = d;
				}
			}
		}
	}
	
	for(i = 0; i < ns; i++)
	{
		ann += nn[i]/((double)ns);
	}

	delete[] nn;	

	return ann;
}

int main(void)
{
	MTRand rg(1);
	FILE *fp1;
	fp1 = fopen("NLdata.txt","w");
	double Juxt3,ParaO3,ParaNN3,dx,dy,dz,tempS,lamoa;
	double *dvecsx3 = new double[N];
	double *dvecsy3 = new double[N];
	double *dvecsz3 = new double[N];
	double *tempx = new double[N];
	double *tempy = new double[N];
	double *tempz = new double[N];
	int *ivecsx3 = new int[N];
	int *ivecsy3 = new int[N];
	int *ivecsz3 = new int[N];
	int *posx = new int[6*N];
	int *posy = new int[6*N];
	int *posz = new int[6*N];
	int *M3 = new int[N*N];
	int i,j,k,l,count,r1,nnd,minnd,veccnt,dsn;
	bool check,check2,check3,check4;

	ivecsx3[0] = 0;
	ivecsy3[0] = 0;
	ivecsz3[0] = 0;
	M3[0] = 0;

	if(N>1)
	{
		for(i = 1; i < N; i++)
		{
			count = 0;
			for(j = 0; j < i; j++)
			{
				veccnt = 0;
				if(ivecsx3[j] >= 0)
				{
					posx[count+veccnt] = ivecsx3[j]+2;
					posy[count+veccnt] = ivecsy3[j];
					posz[count+veccnt] = ivecsz3[j];
					veccnt += 1;
				}
				if(ivecsy3[j] >= 0)
				{
					posx[count+veccnt] = ivecsx3[j];
					posy[count+veccnt] = ivecsy3[j]+2;
					posz[count+veccnt] = ivecsz3[j];
					veccnt += 1;
				}
				if(ivecsz3[j] >= 0)
				{
					posx[count+veccnt] = ivecsx3[j];
					posy[count+veccnt] = ivecsy3[j];
					posz[count+veccnt] = ivecsz3[j]+2;
					veccnt += 1;
				}
				if(ivecsx3[j] <= 0)
				{
					posx[count+veccnt] = ivecsx3[j]-2;
					posy[count+veccnt] = ivecsy3[j];
					posz[count+veccnt] = ivecsz3[j];
					veccnt += 1;
				}
				if(ivecsy3[j] <= 0)
				{
					posx[count+veccnt] = ivecsx3[j];
					posy[count+veccnt] = ivecsy3[j]-2;
					posz[count+veccnt] = ivecsz3[j];
					veccnt += 1;
				}
				if(ivecsz3[j] <= 0)
				{
					posx[count+veccnt] = ivecsx3[j];
					posy[count+veccnt] = ivecsy3[j];
					posz[count+veccnt] = ivecsz3[j]-2;
					veccnt += 1;
				}
				k = 0;
				check = (k < veccnt);
				while(check)
				{
					l = j+1;
					check2 = (l < i);
					check3 = true;
					while(check2)
					{
						if(ivecsx3[l]==posx[count+k] && ivecsy3[l]==posy[count+k] && ivecsz3[l]==posz[count+k])
						{
							check2 = false;
							check3 = false;			
						}
						else
						{
							l += 1;
							check2 = (l < i);
						}
					}
					if(check3)
					{
						if(k==0)
						{
							nnd = posx[count+k]*posx[count+k]+posy[count+k]*posy[count+k]+posz[count+k]*posz[count+k];
							k += 1;
							check = (k < veccnt);
						}
						else
						{
							if(posx[count+k]*posx[count+k]+posy[count+k]*posy[count+k]+posz[count+k]*posz[count+k] < nnd)
							{
								for(l = k; l < veccnt; l++)
								{
									posx[count+l-k] = posx[count+l];
									posy[count+l-k] = posy[count+l];
									posz[count+l-k] = posz[count+l];
								}
								veccnt -= k;
								k = 1;
								check = (k < veccnt);
							}
							else if(posx[count+k]*posx[count+k]+posy[count+k]*posy[count+k]+posz[count+k]*posz[count+k]==nnd)
							{
								k += 1;
								check = (k < veccnt);
							}
							else
							{
								for(l = k+1; l < veccnt; l++)
								{
									posx[count+l-1] = posx[count+l];
									posy[count+l-1] = posy[count+l];
									posz[count+l-1] = posz[count+l];
								}
								veccnt -= 1;
								check = (k < veccnt);
							}
						}
					}
					else
					{
						for(l = k+1; l < veccnt; l++)
						{
							posx[count+l-1] = posx[count+l];
							posy[count+l-1] = posy[count+l];
							posz[count+l-1] = posz[count+l];
						}
						veccnt -= 1;
						check = (k < veccnt);
					}
				}
				if(veccnt > 0)
				{
					if(count==0)
					{
						count = veccnt;
						minnd = nnd;
					}
					else
					{
						if(nnd < minnd)
						{
							for(k = 0; k < veccnt; k++)
							{
								posx[k] = posx[count+k];
								posy[k] = posy[count+k];
								posz[k] = posz[count+k];
							}
							count = veccnt;
							minnd = nnd;
						}
						else if(nnd==minnd)
						{
							count += veccnt;
						}
					}
				}
			}
			r1 = (int)rg.randExc(count);
			ivecsx3[i] = posx[r1];
			ivecsy3[i] = posy[r1];
			ivecsz3[i] = posz[r1];

			M3[i*N+i] = 0;
			for(j = 0; j < i; j++)
			{
				if((ivecsx3[i]-ivecsx3[j])*(ivecsx3[i]-ivecsx3[j])+(ivecsy3[i]-ivecsy3[j])*(ivecsy3[i]-ivecsy3[j])+(ivecsz3[i]-ivecsz3[j])*(ivecsz3[i]-ivecsz3[j])==4)
				{
					M3[i*N+i] += 1;
					M3[j*N+j] += 1;
					M3[i*N+j] = -1;
					M3[j*N+i] = -1;
				}
				else
				{
					M3[i*N+j] = 0;
					M3[j*N+i] = 0;
				}
			}
		}
	}

	for(k = 0; k < nl; k++)
	{
		lamoa = pow(10.0,log10(sl)+log10(fl/sl)*((double) k)/((double)(nl-1)));
		if (lamoa<((double) k)+sl)
		{
			lamoa = ((double) k)+1.0;
		}

		printf("%2.1f\n",lamoa);

		for(j = 0; j < N; j++)
		{
			dvecsx3[j] = (double)ivecsx3[j];
			dvecsy3[j] = (double)ivecsy3[j];
			dvecsz3[j] = (double)ivecsz3[j];
		}

		for (l = 0; l < nn; ++l)
		{
			i = (int) pow(10.0,log10((double)(sn-1))+log10(((double)(N-1))/((double)(sn-1)))*((double) l)/((double)(nn-1)));
			if (i<l+sn-1)
			{
				i = l+sn-1;
			}

			printf("%d\n",i);
			Juxt3 = SyMat(M3, ivecsx3, ivecsy3, ivecsz3, i+1, lamoa, N);
			ParaO3 = Sz(dvecsx3, dvecsy3, dvecsz3, i+1, lamoa);

			check = true;
			count = 0;
			while(check)
			{
				r1 = (int)rg.randExc(i);
				dx = rg.randNorm(0.0,D);
				dy = rg.randNorm(0.0,D);
				dz = rg.randNorm(0.0,D);
				for(j = 0; j < N; j++)
				{
					if(r1==0 && j>0)
					{
						tempx[j] = dvecsx3[j]-dx;
						tempy[j] = dvecsy3[j]-dy;
						tempz[j] = dvecsz3[j]-dz;
					}
					else if(r1>0 && j==r1)
					{
						tempx[j] = dvecsx3[j]+dx;
						tempy[j] = dvecsy3[j]+dy;
						tempz[j] = dvecsz3[j]+dz;
					}
					else
					{
						tempx[j] = dvecsx3[j];
						tempy[j] = dvecsy3[j];
						tempz[j] = dvecsz3[j];
					}
				}
				tempS = Sz(tempx, tempy, tempz, i+1, lamoa);
				if(tempS<ParaO3)
				{
					for(j = 0; j < N; j++)
					{
						dvecsx3[j] = tempx[j];
						dvecsy3[j] = tempy[j];
						dvecsz3[j] = tempz[j];
					}
					ParaO3 = tempS;
					count = 0;
				}
				else if(tempS==ParaO3)
				{
					for(j = 0; j < N; j++)
					{
						dvecsx3[j] = tempx[j];
						dvecsy3[j] = tempy[j];
						dvecsz3[j] = tempz[j];
					}
					count += 1;
					if(count>=bl)
					{
						check = !check;
						ParaNN3 = AveNN(dvecsx3, dvecsy3, dvecsz3, i+1);
					}
				}
				else
				{
					count += 1;
					if(count>=bl)
					{
						check = !check;
						ParaNN3 = AveNN(dvecsx3, dvecsy3, dvecsz3, i+1);
					}
				}
			}

			fprintf(fp1,"%3.1f\t%d\t%8.8e\t%8.8e\t%8.8e\n",lamoa,i+1,ParaNN3,Juxt3,ParaO3);

		}
		printf("\n");
	}

	fclose(fp1);

	system("pause");
	
	return 0;
}