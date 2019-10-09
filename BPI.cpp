
/*****************************Back-projection Stacking and imaging***********************************/
//Time: Feb. 09, 2019
//Place: Stanford, Mitchell, Rm. 359

//Change for just 1 slice

#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "math.h"
#include "string.h"

#define len1 23
#define len2 13
#define MAXREC 2500


void **alloc2(size_t n1,size_t n2, size_t size);
void free2(void **p);


int main()
{	

/////////PART 1
/***************************Define the arguments*******************************/

	FILE *fp_control;           //Control file pointer
	FILE *fp_info_geophone;	    //the information of geophone pointer
	FILE *fp_t_travel;			//travel time file pointer

	FILE *fp_data;		    	//seismic data pointer

	FILE *fp_image;		    	//Final imaging result file pointer

	FILE *fp_vz_test;			//Examine whether the recorder data are well read in

	int num_geo;              //the total number of geophones
	int num_curr;             //current number of geophones
	int x0,x1;
	int y0,y1;
	int z0,z1;                //The start and the end of the grid searching space
	int dx_trav;
	int dx,dy,dz;             //3D interval
	int fre;                  //sampling frequency
	int daytime;              //how many hours for a day

	char info_geophone[400];           //Name of vz component recorded data file
	char t_travelfile[400];     //Name of vz component recorded data file

	float  **data;              /*store all the recorded data read in: 2D*/

	int twin;				/*time window for event location: s*/
	int nt;

	int nx,ny,nz;                /*horizontal and vertical searching grid points*/
	float cx,cy,cz;
	float dist;	       	/*Distance between the searching grid point and the location of the geophone*/

	int travmod_x;
	int travmod_z;	/*Dimensions of traveltime model*/

	int table_x,table_z;		/*Dimensions of traveltime table*/
        float  *t_travel;           /*store the travel time read in*/

	int n_point;
	float part;
	float time,t1,t2;

	int  t_counter;           /*store the travel time counter*/
	
	float  *shift;        /*store the time-shifted vz data, the volume should be the same as data_in*/
	
	float  *stacking;	      /*1D array, store the stacking trace for each grid point*/

	int T_count;             /*counter for time window*/

	float **image;           /*store the first part of the final imaging result*/
	float image01;           /*store the intermediate imaging result*/

	int i,j,k,h;
	int ix,iy,iz;

	float x;
	float y;                //information read in from info_geophone

	int counter=0;
	char name_geo[len1];

	float gx[MAXREC];
	float gy[MAXREC];        //current geophone (in the (20)?-th day)

	char name[len2]=".037.new.dat";     //??????????????????????????????

/////////PART 2
/***************************Read the control file******************************/

		if((fp_control=fopen("work.prn","r"))==NULL)
		{
			return 0;
		}
		
		fscanf(fp_control,"%d",&num_geo);
		fscanf(fp_control,"%d",&num_curr);
		fscanf(fp_control,"%d%d",&x0,&x1);
		fscanf(fp_control,"%d%d",&y0,&y1);
		fscanf(fp_control,"%d%d",&z0,&z1);
		fscanf(fp_control,"%d%d%d",&dx,&dy,&dz);
        fscanf(fp_control,"%d",&daytime);
		fscanf(fp_control,"%d",&fre);
		fscanf(fp_control,"%d%d",&travmod_x,&travmod_z);
		fscanf(fp_control,"%d",&dx_trav);
		fscanf(fp_control,"%d",&twin);

		fscanf(fp_control,"%s",info_geophone);
		fscanf(fp_control,"%s",t_travelfile);

		fclose (fp_control);

//		printf("%s",name);

		nt=daytime*60*60*fre;   //How many sampling points in one day

		nx=(x1-x0)/dx+1;
		ny=(y1-y0)/dy+1;
		nz=(z1-z0)/dz+1;

 	        T_count=twin*fre;       //How many sampling points in a timewindow(s)

//		printf("%d\t%d\t%d\n",nx,ny,nz);
//              printf("T:%d\tT_count:%d\tnt:%d\ttravmod_z:%d\n",twin,T_count,nt,travmod_z);
/////////PART 3
/***************************Read in the geophone information file******************************/
// We store the location coordinate in gx[] & gy[];
// We store the recorded data in data[num_curr][nt];

		if((fp_info_geophone=fopen(info_geophone,"r"))==NULL)
		{
			return 0;
		}

		data=(float **)alloc2(nt,num_curr,sizeof(float));

		
/*		if((fp_vz_test=fopen("vz_test.dat","ab+"))==NULL)
		{
			return 0;
		}
*/


/////////Loop for different geophones////////////////
		for(i=0;i<num_geo;i++)
		{	
//			printf("%d\t%d\n",i,nt);
		
			for(j=0;j<len1;j++)
			{
				name_geo[j]='\0';
			}


			fscanf(fp_info_geophone,"%s%f%f",name_geo,&x,&y);
//			printf("%s\t%f\t%f\n",name_geo,x,y);
			strcat(name_geo,name);

//			printf("%s\n",name_geo);

			if((fp_data=fopen(name_geo,"rb"))!=NULL)
			{	
			//	printf("%d\n",counter);
				gx[counter]=x;
				gy[counter]=y;
				
				printf("%d\t%f\t%f\t%s\n",counter,gx[counter],gy[counter],name_geo);
				fread(data[counter],sizeof(float),nt,fp_data);

//				for(j=0;j<50;j++)
//					printf("%f\n",data[counter][j]);

				counter++;

				fclose(fp_data);
			}

//		    fclose(fp_vz_test);

		}
			
		fclose(fp_info_geophone);
//    		printf("T:%d\tT_count:%d\tnt:%d\tdx_trav:%d\n",twin,T_count,nt,dx_trav);

///////////////END OF LOOP//////////////////////

/////////PART 4: 
//////////////////////////READ IN THE TRAVELTIME TABLE//////////////////////////////	

	if((fp_t_travel=fopen(t_travelfile,"rb"))==NULL)
	{
		return 0;
	}

	t_travel=(float *)malloc(travmod_x*travmod_z*sizeof(float));    
	fread(t_travel,sizeof(float),travmod_x*travmod_z,fp_t_travel);//I need to read in all the traveltime data, 
								 //although the first line of the data shows the travel time from the first point (0, 0) to the surface, useless.
	fclose(fp_t_travel);

/////////PART 5: 
//////////////////////////LOOP FOR EVERY SEARCHING GRID POINT///////////////////////////////
//our searching grid point coordinate is (cx,cy,cz) 

		shift=(float *)malloc(nt*sizeof(float));
		stacking=(float *)malloc(nt*sizeof(float));


//printf("%d\t%d\n",twin,T_count);
//printf("%d\n",nx*ny*(nz-1)*60*60/twin);
		image=(float **)alloc2(nx*ny*(nz-1),daytime*60*60/twin,sizeof(float));

//printf("%d\t%d\n",T,T_count);
		if((fp_image=fopen("image.dat","wb"))==NULL)
		{
			return 0;
   	        }

		for(i=0;i<(daytime*60*60/twin);i++)
			for(j=0;j<nx*ny*(nz-1);j++)
			{
				image[i][j]=0.0;
			}


	//   fwrite(image[0],sizeof(float),nx*ny*(nz-1)*60*60/T,fp_image);///////////////////////
//		printf("%d\t%d\t%d\n",nx,ny,nz);


//		for(ix=0;ix<nx;ix++)
		for(iy=0;iy<ny;iy++)	
		{
			cy=y0+iy*dy;

			for(ix=0;ix<nx;ix++)
//			for(iy=0;iy<ny;iy++)
			{
				cx=x0+ix*dx;
//				cy=y0+iy*dy;

				for(iz=1;iz<nz;iz++)
			//	for(iz=0;iz<nz-1;iz++)
				{
					cz=z0+iz*dz;
					
					for(k=0;k<nt;k++)
					{
						stacking[k]=0.0;
					}
					////////////LOOP FOR EACH TRACE//////////////
					for(counter=0;counter<num_curr;counter++)
					{
						dist=sqrt((gx[counter]-cx)*(gx[counter]-cx)+(gy[counter]-cy)*(gy[counter]-cy));

						n_point=(int)(dist/dx_trav);
						printf("dist: %f\t dx_trav: %d\t travmod_z: %d\n",dist,dx_trav,travmod_z);

						part=(dist-n_point*dx_trav)/(dx_trav*1.0);
						
						t1=0.0;
						t2=0.0;
						t1=t_travel[iz*(dz/dx_trav)*travmod_z+n_point];
						t2=t_travel[iz*(dz/dx_trav)*travmod_z+n_point+1];

						time=part*(t2-t1)+t1;

						t_counter=(int) round(time*fre);

						printf("%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%d\n",iy,ix,iz,cx,cy,cz,gx[counter],gy[counter],dist,n_point,t1,t2,t_counter);

						for(h=0;h<nt;h++)
						{
							shift[h]=0.0;
						}

						for(k=0;k<(nt-t_counter);k++)
						{
				    		shift[k]=data[counter][k+t_counter];	
			//	    		shift[k]=data[counter][k+t_counter]*pow(10,2);				

			//				printf("%d\t%f\n",k,shift[k]);
						}
		
						for(k=0;k<nt;k++)
						{
							stacking[k]+=shift[k];
			//				printf("%d\t%f\n",k,stacking[k]);
						}
					}
				//////////////////////////////////////////////////
					
					for(i=0;i<daytime*60*60/twin;i++)//Change to the night time 
					{	
						image01=0.0;

						for(j=0;j<fre*twin;j++)
						{
							if(fabs(stacking[j+i*fre*twin])>fabs(image01))
							{
								image01=fabs(stacking[j+i*fre*twin]);
					//			printf("%d\t%d\t%d\t%d\t%d\t%f\t%d\n",ix,iy,iz,i,j,image01,T);
							}
						}

						printf("%d\t%d\t%d\t%d\t%f\t\n",iy,ix,iz,i,image01);

							image[i][iy*nx*(nz-1)+ix*(nz-1)+(iz-1)]=image01;

					}

				}

			}

		}

	   printf("%d\t%d\n",twin,T_count);

	   fwrite(image[0],sizeof(float),nx*ny*(nz-1)*daytime*60*60/twin,fp_image);///////////////////////

	   fclose(fp_image);


}

/**************************************分配二维动态内存*************************************/
void **alloc2(size_t n1,size_t n2, size_t size)    /*n1表示纵向网格点数，n2表示横向网格点数*/
{
	size_t i2;
	void **p;
	 if((p=(void **)malloc(n2*sizeof(void *)))==NULL)
	    return NULL;
	 if((p[0]=(void *)malloc(n1*n2*size))==NULL)
	 {
		  free(p);
		 return NULL;
	 }

	   for(i2=0;i2<n2;i2++)
		   p[i2]=(char *)p[0]+i2*n1*size;
	   return p;
}

/*************************************释放二维动态内存***************************************/
void free2(void **p)
{
  free(p[0]);
  free(p);
}

