// #include <stdio.h>
// #include <math.h>
// #include <string.h>
#include "gl_functions.h"

#include <cstdio> 
#include <cmath>
#include <cstring>
#include <Kokkos_Core.hpp>
#include "env_kokkos.h"

extern int nface, nspt, natm, nchr;
extern int **extr_v;						//[3][nspt]
extern int **extr_f;						//[2][nface]
extern int **face, **face_copy;				//[3][nface]
extern double **vert, **snrm;				//[3][nspt];
extern double *tr_xyz, *tr_q;				//[3][nface]
extern double *tr_area, *bvct, *xvct;		//[nface];
extern double **atmpos;						//[3][natm/nchr];
extern double *atmrad, *atmchr, *chrpos;	//[natm/nchr]; 
// extern double *dev_tr_xyz, *dev_tr_q, *dev_tr_area, *dev_bvct;
extern double **tr_xyz2D, **tr_q2D;

/* function computing the area of a triangle given vertices coodinates */
double triangle_area(double v[3][3]) {
    int i;
    double a[3],b[3],c[3],aa,bb,cc,ss,area;
    for (i=0;i<=2;i++){
        a[i]=v[i][0]-v[i][1];
        b[i]=v[i][0]-v[i][2];
        c[i]=v[i][1]-v[i][2];
    }
    aa=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    bb=sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
    cc=sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
    ss=0.5*(aa+bb+cc);
    area=sqrt(ss*(ss-aa)*(ss-bb)*(ss-cc));
    return(area);
}

/* function read in molecule information */
void readin(char fname[16], char density[16]) {
    FILE *fp;
    char c;
	char fpath[256];
	char fname_tp[256];

	FILE *fpw; //yang
	char c1,c2,c3,buff[256]; //yang
	int iii,jjj,count; //yang

    int i,j,k,i1,i2,i3,j1,j2,j3,ii,jj,kk,namelength=4,nfacenew,ichanged;
    double den,prob_rds,a1,a2,a3,b1,b2,b3,a_norm,r0_norm,v0_norm;
	double r0[3],v0[3],v[3][3],r[3][3];
	int    idx[3],jface[3],iface[3];
    double rs,rd[3],pot=0.0,sum=0.0,pot_temp=0.0;
    double temp_x,temp_q,tchg,tpos[3],dist_local,area_local;
	double cos_theta,G0,tp1,G1,r_s[3];
	double xx[3],yy[3];



	/*count atom lines*/
	sprintf(fpath,"../test_proteins/");
	sprintf(fname_tp, "%s%s.pqr",fpath,fname);
   	fp=fopen(fname_tp,"r");
   	count = 0;
	while(fgets(buff,256,fp)) {
		if (strstr(buff,"ATOM")!=NULL) {
			count++;	
		}
	}
	printf("count is:%d\n",count); //yang
	fclose(fp);


	fp=fopen(fname_tp,"r");
	sprintf(fname_tp, "%s%s.xyzr",fpath,fname);
   	fpw=fopen(fname_tp,"w");

   	nchr = count;
   	natm = count;
	if ((atmrad=(double *) (Kokkos::kokkos_malloc(natm*sizeof(double))))==NULL) {
		printf("error in allcating atmrad");
	}
	atmpos=Make2DDoubleArray(3,natm,(char*) "atmpos");
	if ((atmchr=(double *) (Kokkos::kokkos_malloc(nchr*sizeof(double))))==NULL) {
		printf("error in allcating atmchr");
	}
	if ((chrpos=(double *) (Kokkos::kokkos_malloc(3*nchr*sizeof(double))))==NULL) { 
		printf("error in allcating chrpos");
	}


	for (i=0;i<count;i++) {
		fscanf(fp,"%s %d %s %s %d %lf %lf %lf %lf %lf",&c1,&iii,&c2,&c3,&jjj,&a1,&a2,&a3,&b1,&b2);
		sprintf(buff,"%.3f\t\t %.3f\t\t %.3f\t\t %.4f\n",a1,a2,a3,b2);
		fputs(buff,fpw);
		chrpos[3*i]=a1;
		chrpos[3*i+1]=a2;
		chrpos[3*i+2]=a3;
		atmchr[i]=b1;

		atmpos[0][i]=a1;
		atmpos[1][i]=a2;
		atmpos[2][i]=a3;		
		atmrad[i]=b2;
	}

	fclose(fp);
	fclose(fpw);
	printf("finish reading pqr and writing xyzr file...\n");



	/*read in vertices*/
	sprintf(fpath,"../test_proteins/");
	sprintf(fname_tp,"./msms -if %s%s.xyzr -prob 1.4 -dens %s -of %s%s ",fpath,fname,density,fpath,fname);
	printf("%s\n",fname_tp);

	system(fname_tp);

	/* read in vert */
	sprintf(fname_tp, "%s%s.vert",fpath,fname);	//or use "strcat(fname_tp,".vert")"
	printf("%s\n",fname_tp);

	/* open the file and read through the first two rows*/
	fp=fopen(fname_tp,"r");
	for (i=1;i<=2;i++){
        while (c=getc(fp)!='\n'){
        }
    }
    fscanf(fp,"%d %d %lf %lf ",&nspt,&natm,&den,&prob_rds);
    printf("nspt=%d, natm=%d, den=%lf, prob=%lf\n", nspt,natm,den,prob_rds);
	/*allocate variables for vertices file*/
	extr_v=Make2DIntArray(3,nspt,(char*) "extr_v");
	vert=Make2DDoubleArray(3,nspt,(char*) "vert");
	snrm=Make2DDoubleArray(3,nspt,(char*) "snrm");

    for (i=0;i<=nspt-1;i++){
        fscanf(fp,"%lf %lf %lf %lf %lf %lf %d %d %d",&a1,&a2,&a3,&b1,&b2,&b3,&i1,&i2,&i3);

		/*radial projection to improve msms accuracy, ONLY FOR SPHERE!!!!!!!!
        a_norm=sqrt(a1*a1+a2*a2+a3*a3);
		b1=a1/a_norm;
		a1=b1*rds;
		b2=a2/a_norm;
		a2=b2*rds;
		b3=a3/a_norm;
		a3=b3*rds;*/

		vert[0][i]=a1;
        vert[1][i]=a2;
        vert[2][i]=a3;
        snrm[0][i]=b1;
        snrm[1][i]=b2;
        snrm[2][i]=b3;
        extr_v[0][i]=i1;
        extr_v[1][i]=i2;
        extr_v[2][i]=i3;
    }
    fclose(fp);
	printf("finish reading vertices file...\n");

	/* read in faces */

	sprintf(fname_tp, "%s%s.face",fpath,fname);
   	fp=fopen(fname_tp,"r");
	for (i=1;i<=2;i++){
        while (c=getc(fp)!='\n'){
        }
    }
    fscanf(fp,"%d %d %lf %lf ",&nface,&natm,&den,&prob_rds);
    printf("nface=%d, natm=%d, den=%lf, prob=%lf\n", nface,natm,den,prob_rds);

	extr_f=Make2DIntArray(2,nface,(char*) "extr_f");
	face=Make2DIntArray(3,nface,(char*) "face");


    for (i=0;i<=nface-1;i++){
        fscanf(fp,"%d %d %d %d %d",&j1,&j2,&j3,&i1,&i2);
        face[0][i]=j1;
        face[1][i]=j2;
        face[2][i]=j3;
        extr_f[0][i]=i1;
        extr_f[1][i]=i2;
    }
    fclose(fp);
	printf("finish reading face file...\n");

	// /*read atom coodinates and radius */
	// sprintf(fname_tp, "%s%s.xyzr",fpath,fname);
	// fp=fopen(fname_tp,"r");

	// if ((atmrad=(double *) malloc(natm*sizeof(double)))==NULL) {
	// // if ((atmrad=(double *) (Kokkos::kokkos_malloc(natm*sizeof(double))))==NULL) {
	// 	printf("error in allcating atmrad");
	// }
	// atmpos=Make2DDoubleArray(3,natm,(char*) "atmpos");

	// printf("natm is:%d\n",natm); //yang
	// for (i=0;i<=natm-1;i++){
	// 	fscanf(fp,"%lf %lf %lf %lf ",&a1,&a2,&a3,&b1);
	// 	atmpos[0][i]=a1;
	// 	atmpos[1][i]=a2;
	// 	atmpos[2][i]=a3;
	// 	atmrad[i]=b1;
    // }
	// fclose(fp);
	// printf("finish reading position file...\n");

	// /*read charge coodinates and charge */
	// sprintf(fname_tp, "%s%s.pqr",fpath,fname);
	// fp=fopen(fname_tp,"r");

	// nchr=natm;
	// printf("nchr is:%d\n",nchr); //yang
	// // if ((atmchr=(double *) malloc(nchr*sizeof(double)))==NULL){
	// if ((atmchr=(double *) (Kokkos::kokkos_malloc(nchr*sizeof(double))))==NULL) {
	// 	printf("error in allcating atmchr");
	// }
	// // if ((chrpos=(double *) malloc(3*nchr*sizeof(double)))==NULL){
	// if ((chrpos=(double *) (Kokkos::kokkos_malloc(3*nchr*sizeof(double))))==NULL) { 
	// 	printf("error in allcating chrpos");
	// }

	// for (i=0;i<=nchr-1;i++){
	// 	// fscanf(fp,"%lf %lf %lf %lf ",&a1,&a2,&a3,&b1);
	// 	fscanf(fp,"%s %d %s %s %d %lf %lf %lf %lf %lf",&c1,&iii,&c2,&c3,&jjj,&a1,&a2,&a3,&b1,&b2);
	// 	chrpos[3*i]=a1;
	// 	chrpos[3*i+1]=a2;
	// 	chrpos[3*i+2]=a3;
	// 	atmchr[i]=b1;
    // }
	// fclose(fp);
	// printf("finish reading charge file...\n");

	/* delele triangles with extreme small areas and closed to each other */
	nfacenew=nface;
	face_copy=Make2DIntArray(3,nface,(char*) "face_copy");
	for (i=0;i<3; i++) memcpy(face_copy[i],face[i],nface*sizeof(int));

	for (i=0;i<nface;i++){
		for (ii=0;ii<3;ii++){
			iface[ii]=face[ii][i];
			xx[ii]=0.0;
		}
        for (ii=0;ii<3;ii++){
			for (kk=0;kk<3;kk++){
				r[kk][ii]=vert[kk][iface[ii]-1];
				xx[kk]=xx[kk]+1.0/3.0*r[kk][ii];
			}
		}

        area_local=triangle_area(r);

		if (area_local<1e-5) {
			printf("The %d th triangle has small area:%e\n",i,area_local);
			goto exit;
		}

		for (j=i-10;(j>=0 && j<i);j++){
			for(jj=0;jj<3;jj++) {
				jface[jj]=face[jj][j];
				yy[jj]=0.0;
			}
			for(jj=0;jj<3;jj++){
				for (kk=0;kk<3;kk++){
					r[kk][jj]=vert[kk][jface[jj]-1];
					yy[kk]=yy[kk]+1.0/3.0*r[kk][jj];
				}
			}
			dist_local=0.0;
			for (jj=0;jj<3;jj++) dist_local+=(xx[jj]-yy[jj])*(xx[jj]-yy[jj]);
			dist_local=sqrt(dist_local);
			if (dist_local<1e-6) {
				printf("particles %d and %d are too close: %e\n", i,j,dist_local);
				goto exit;
			}
		}

		continue;
exit:	ichanged=nface-nfacenew;
		for (ii=i-ichanged;ii<nface;ii++) {
			for (jj=0;jj<3;jj++){
				face_copy[jj][ii]=face_copy[jj][ii+1];
			}
		}
        nfacenew=nfacenew-1;
		// printf("nface in readin exit is %d\n",nface);
	}

	printf("%d faces are deleted\n",nface-nfacenew);
	nface=nfacenew;
	// printf("nface after=nfacenew is %d\n",nface);



	for(i = 0; i < 3; i++) free(face[i]);
	free(face);

	face=Make2DIntArray(3,nface,(char*) "face msms");
	for (i=0; i<nface; i++){
		for (j=0; j<3; j++) face[j][i]=face_copy[j][i];
	}
	for(i = 0; i<3; i++) free(face_copy[i]);
	free(face_copy);

//tr_xyz: The position of the particles on surface
//tr_q:	  The normail direction at the particle location
//tr_area: the triangular area of each element
	// tr_xyz=(double *) calloc(3*nface, sizeof(double));
	// tr_q=(double *) calloc(3*nface, sizeof(double));
	// tr_area=(double *) calloc(nface, sizeof(double));
	// bvct=(double *) calloc(2*nface, sizeof(double));
	tr_xyz2D=Make2DDoubleArray(3,nface,"tr_xyz2D");
    tr_q2D=Make2DDoubleArray(3,nface,"tr_q2D");

	tr_xyz=(double *) (Kokkos::kokkos_malloc(3*nface * sizeof(double)));
	tr_q=(double *) (Kokkos::kokkos_malloc(3*nface * sizeof(double)));
	tr_area=(double *) (Kokkos::kokkos_malloc(nface * sizeof(double)));
	bvct=(double *) (Kokkos::kokkos_malloc(2*nface * sizeof(double)));
  	
	// ViewVectorDouble tr_xyz("tr_xyz", 3*nface);
  	// ViewVectorType::HostMirror h_y = Kokkos::create_mirror_view( tr_xyz );
  	// ViewVectorType::HostMirror h_x = Kokkos::create_mirror_view( tr_q);
  	// ViewMatrixType::HostMirror h_A = Kokkos::create_mirror_view( tr_area);
	// ViewMatrixType::HostMirror h_A = Kokkos::create_mirror_view( bvct);

	// Kokkos::View<double**,Kokkos::CudaUVMSpace> tr_xyz2D ("tr_xyz2D", 3,nface);
	// Kokkos::View<double**,Kokkos::CudaUVMSpace> tr_q2D ("tr_q2D", 3,nface);


    for (i=0;i<nface;i++){
        for (j=0;j<=2;j++){
            idx[j]=face[j][i];
        }
        for (j=0;j<=2;j++){
            r0[j]=0;
            v0[j]=0;
            for (k=0;k<=2;k++){
				r0[j]=r0[j]+vert[j][idx[k]-1]/3.0;
                v0[j]=v0[j]+snrm[j][idx[k]-1]/3.0;
				r[j][k]=vert[j][idx[k]-1];
				v[j][k]=snrm[j][idx[k]-1];
            }

		}
		v0_norm=sqrt(v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]);
		for (k=0;k<=2;k++){
			v0[k]=v0[k]/v0_norm;
		}

		/* radial projection for sphere only!!!
		r0_norm=sqrt(r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2]);
		for (k=0;k<=2;k++){
			r0[k]=r0[k]/r0_norm*rds;
		}*/

        for (j=0;j<=2;j++){
            tr_xyz[3*i+j]=r0[j];
            tr_q[3*i+j]=v0[j];
            tr_xyz2D[j][i] = r0[j];
            tr_q2D[j][i] = v0[j];

            // tr_xyz(3*i+j)=r0[j];
            // tr_q(3*i+j)=v0[j];
            // tr_xyz2D[j][i] = r0[j];
            // tr_q2D[j][i] = v0[j]; 
            // tr_xyz2D(j,i) = r0[j];
            // tr_q2D(j,i) = v0[j];            
        }
        tr_area[i]=triangle_area(r);
        sum=sum+tr_area[i];

		// // for kokkos host and device:
        // for (j=0;j<3*nface;j++){
        // 	dev_tr_xyz(j)=tr_xyz[j];
		// 	dev_tr_q(j)=tr_q[j];
        // }

        // for (j=0;j<nface;j++){
		// 	dev_tr_area(j)=tr_area[j];
        // }

        // for (j=0;j<2*nface;j++){
		// 	dev_bvct(j)=bvct[j];
        // }

	}
    printf("total area = %f\n",sum);
    
}
