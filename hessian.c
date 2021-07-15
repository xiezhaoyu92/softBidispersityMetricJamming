//
//  hessian.c
//  softBidispersityMetricJamming
//
//  Created by Zhaoyu Xie on 6/6/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#include "hessian.h"
#include "surface.h"
#include "particle.h"
#include "solveHessian.h"
#include <stdlib.h>
#include <math.h>

#ifndef PI
#define PI 3.141592653589793
#endif
#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

int overlapQ(particle *p1, particle *p2) {
    if(abs(p1->coord[0]-p2->coord[0])>1||abs(p1->coord[1]-p2->coord[1])>1||abs(p1->coord[2]-p2->coord[2])>1)
        return FALSE;
    
    double r1 = p1->rr;
    double r2 = p2->rr;
    double s1 = sin(p1->th);
    double s2 = sin(p2->th);
    double c1 = cos(p1->th);
    double c2 = cos(p2->th);
    double c12 = cos(p1->ph-p2->ph);
    double s12 = sin(p1->ph-p2->ph);
    double R = sqrt(r1*r1+r2*r2-2*r1*r2*c1*c2-2*r1*r2*s1*s2*c12);
    
    if(R<p1->rad+p2->rad)
        return TRUE;
    else
        return FALSE;
}

double tesseralHarmonicsZDthth(int n, int m, double th, double ph) {
    if(m!=0||n>10) {
        printf("Error: harmonics doesn't exist for (%d, %d)!\n", n, m);
        exit(1);
    }
    switch(n) {
        case 0:
            return 0;
            break;
        case 1:
            return -1.0/2*sqrt(3/PI)*cos(th);
            break;
        case 2:
            return -3.0/2*sqrt(5/PI)*cos(2*th);
            break;
        case 3:
            return -3.0/16*sqrt(7/PI)*(cos(th)+15*cos(3*th));
            break;
        case 4:
            return -5.0/8*sqrt(9/PI)*(cos(2*th)+7*cos(4*th));
            break;
        case 5:
            return -15.0/256*sqrt(11/PI)*(2*cos(th)+21*cos(3*th)+105*cos(5*th));
            break;
        case 6:
            return -21.0/256*sqrt(13/PI)*(5*cos(2*th)+24*cos(4*th)+99*cos(6*th));
            break;
        case 7:
            return -7.0/2048*sqrt(15/PI)*(25*cos(th)+243*cos(3*th)+825*cos(5*th)+3003*cos(7*th));
            break;
        case 8:
            return -9.0/1024*sqrt(17/PI)*(35*cos(2*th)+154*cos(4*th)+429*cos(6*th)+1430*cos(8*th));
            break;
        case 9:
            return -45.0/65536*sqrt(19/PI)*(98*cos(th)+924*cos(3*th)+2860*cos(5*th)+7007*cos(7*th)+21879*cos(9*th));
            break;
        case 10:
            return -55.0/65536*sqrt(21/PI)*(294*cos(2*th)+1248*cos(4*th)+3159*cos(6*th)+7072*cos(8*th)+20995*cos(10*th));
            break;
        default:
            exit(1);
    }
}
double tesseralHarmonicsZDthph(int n, int m, double th, double ph) {
    if(m!=0||n>10) {
        printf("Error: harmonics doesn't exist for (%d, %d)!\n", n, m);
        exit(1);
    }
    return 0;
}
double tesseralHarmonicsZDphph(int n, int m, double th, double ph) {
    if(m!=0||n>10) {
        printf("Error: harmonics doesn't exist for (%d, %d)!\n", n, m);
        exit(1);
    }
    return 0;
}

double VththS(int count, double nm[][3], particle *p1, particle *p2, double epsilon) {
    double r1 = p1->rr;
    double r2 = p2->rr;
    double s1 = sin(p1->th);
    double s2 = sin(p2->th);
    double c1 = cos(p1->th);
    double c2 = cos(p2->th);
    double c12 = cos(p1->ph-p2->ph);
    double s12 = sin(p1->ph-p2->ph);
    double R = sqrt(r1*r1+r2*r2-2*r1*r2*c1*c2-2*r1*r2*s1*s2*c12);
    
    double drth1 = 0;
    double drthth1 = 0;
    for(int i=0; i<count; i++) {
        drth1 += nm[i][2]*tesseralHarmonicsZDth((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
        drthth1 += nm[i][2]*tesseralHarmonicsZDthth((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
    }
    double Rth1 = (r1*drth1-(drth1*s1+r1*c1)*r2*s2*c12-(drth1*c1-r1*s1)*r2*c2)/R;
    double Rthth1 = (drth1*drth1+r1*drthth1-(drthth1*s1+drth1*c1+drth1*c1-r1*s1)*r2*s2*c12-(drthth1*c1-drth1*s1-drth1*s1-r1*c1)*r2*c2-Rth1*Rth1)/R;
    
    double sigma = p1->rad+p2->rad;
    return epsilon/sigma/sigma*Rth1*Rth1-epsilon/sigma*(1-R/sigma)*Rthth1;
}
double VphphS(int count, double nm[][3], particle *p1, particle *p2, double epsilon) {
    double r1 = p1->rr;
    double r2 = p2->rr;
    double s1 = sin(p1->th);
    double s2 = sin(p2->th);
    double c1 = cos(p1->th);
    double c2 = cos(p2->th);
    double c12 = cos(p1->ph-p2->ph);
    double s12 = sin(p1->ph-p2->ph);
    double R = sqrt(r1*r1+r2*r2-2*r1*r2*c1*c2-2*r1*r2*s1*s2*c12);
    
    double drph1 = 0;
    double drphph1 = 0;
    for(int i=0; i<count; i++) {
        drph1 += nm[i][2]*tesseralHarmonicsZDph((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
        drphph1 += nm[i][2]*tesseralHarmonicsZDphph((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
    }
    double Rph1 = (r1*drph1-(drph1*c12-r1*s12)*r2*s1*s2-drph1*r2*c1*c2)/R;
    double Rphph1 = (drph1*drph1+r1*drphph1-(drphph1*c12-drph1*s12-drph1*s12-r1*c12)*r2*s1*s2-drphph1*r2*c1*c2-Rph1*Rph1)/R;
    
    double sigma = p1->rad+p2->rad;
    return epsilon/sigma/sigma*Rph1*Rph1-epsilon/sigma*(1-R/sigma)*Rphph1;
}
double VthphS(int count, double nm[][3], particle *p1, particle *p2, double epsilon) {
    double r1 = p1->rr;
    double r2 = p2->rr;
    double s1 = sin(p1->th);
    double s2 = sin(p2->th);
    double c1 = cos(p1->th);
    double c2 = cos(p2->th);
    double c12 = cos(p1->ph-p2->ph);
    double s12 = sin(p1->ph-p2->ph);
    double R = sqrt(r1*r1+r2*r2-2*r1*r2*c1*c2-2*r1*r2*s1*s2*c12);
    
    double drth1 = 0;
    double drph1 = 0;
    double drthph1 = 0;
    for(int i=0; i<count; i++) {
        drth1 += nm[i][2]*tesseralHarmonicsZDth((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
        drph1 += nm[i][2]*tesseralHarmonicsZDph((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
        drthph1 += nm[i][2]*tesseralHarmonicsZDthph((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
    }
    double Rth1 = (r1*drth1-(drth1*s1+r1*c1)*r2*s2*c12-(drth1*c1-r1*s1)*r2*c2)/R;
    double Rph1 = (r1*drph1-(drph1*c12-r1*s12)*r2*s1*s2-drph1*r2*c1*c2)/R;
    double Rthph1 = (drph1*drth1+r1*drthph1-(drthph1*s1+drph1*c1)*r2*s2*c12+(drth1*s1+r1*c1)*r2*s2*s12-(drthph1*c1-drph1*s1)*r2*c2-Rth1*Rph1)/R;
    
    double sigma = p1->rad+p2->rad;
    return epsilon/sigma/sigma*Rth1*Rph1-epsilon/sigma*(1-R/sigma)*Rthph1;
}
double VththP(int count, double nm[][3], particle *p1, particle *p2, double epsilon) {
    double r1 = p1->rr;
    double r2 = p2->rr;
    double s1 = sin(p1->th);
    double s2 = sin(p2->th);
    double c1 = cos(p1->th);
    double c2 = cos(p2->th);
    double c12 = cos(p1->ph-p2->ph);
    double s12 = sin(p1->ph-p2->ph);
    double R = sqrt(r1*r1+r2*r2-2*r1*r2*c1*c2-2*r1*r2*s1*s2*c12);

    double drth1 = 0;
    double drth2 = 0;
    for(int i=0; i<count; i++) {
        drth1 += nm[i][2]*tesseralHarmonicsZDth((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
        drth2 += nm[i][2]*tesseralHarmonicsZDth((int)(nm[i][0]), (int)(nm[i][1]), p2->th, p2->ph);
    }
    double Rth1 = (r1*drth1-(drth1*s1+r1*c1)*r2*s2*c12-(drth1*c1-r1*s1)*r2*c2)/R;
    double Rth2 = (r2*drth2-(drth2*s2+r2*c2)*r1*s1*c12-(drth2*c2-r2*s2)*r1*c1)/R;
    double Rth1th2 = (-(drth1*s1+r1*c1)*(drth2*s2+r2*c2)*c12-(drth1*c1-r1*s1)*(drth2*c2-r2*s2)-Rth1*Rth2)/R;
    
    double sigma = p1->rad+p2->rad;
    return epsilon/sigma/sigma*Rth1*Rth2-epsilon/sigma*(1-R/sigma)*Rth1th2;
}
double VthphP(int count, double nm[][3], particle *p1, particle *p2, double epsilon) {
    double r1 = p1->rr;
    double r2 = p2->rr;
    double s1 = sin(p1->th);
    double s2 = sin(p2->th);
    double c1 = cos(p1->th);
    double c2 = cos(p2->th);
    double c12 = cos(p1->ph-p2->ph);
    double s12 = sin(p1->ph-p2->ph);
    double R = sqrt(r1*r1+r2*r2-2*r1*r2*c1*c2-2*r1*r2*s1*s2*c12);
    
    double drth1 = 0;
    double drph2 = 0;
    for(int i=0; i<count; i++) {
        drth1 += nm[i][2]*tesseralHarmonicsZDth((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
        drph2 += nm[i][2]*tesseralHarmonicsZDph((int)(nm[i][0]), (int)(nm[i][1]), p2->th, p2->ph);
    }
    double Rth1 = (r1*drth1-(drth1*s1+r1*c1)*r2*s2*c12-(drth1*c1-r1*s1)*r2*c2)/R;
    double Rph2 = (r2*drph2-(drph2*c12+r2*s12)*r1*s2*s1-drph2*r1*c2*c1)/R;
    double Rth1ph2 = (-(drth1*s1+r1*c1)*s2*(drph2*c12+r2*s12)-(drth1*c1-r1*s1)*c2*drph2-Rth1*Rph2)/R;
    
    double sigma = p1->rad+p2->rad;
    return epsilon/sigma/sigma*Rth1*Rph2-epsilon/sigma*(1-R/sigma)*Rth1ph2;
}
double VphphP(int count, double nm[][3], particle *p1, particle *p2, double epsilon) {
    double r1 = p1->rr;
    double r2 = p2->rr;
    double s1 = sin(p1->th);
    double s2 = sin(p2->th);
    double c1 = cos(p1->th);
    double c2 = cos(p2->th);
    double c12 = cos(p1->ph-p2->ph);
    double s12 = sin(p1->ph-p2->ph);
    double R = sqrt(r1*r1+r2*r2-2*r1*r2*c1*c2-2*r1*r2*s1*s2*c12);

    double drph1 = 0;
    double drph2 = 0;
    for(int i=0; i<count; i++) {
        drph1 += nm[i][2]*tesseralHarmonicsZDph((int)(nm[i][0]), (int)(nm[i][1]), p1->th, p1->ph);
        drph2 += nm[i][2]*tesseralHarmonicsZDph((int)(nm[i][0]), (int)(nm[i][1]), p2->th, p2->ph);
    }
    double Rph1 = (r1*drph1-(drph1*c12-r1*s12)*r2*s1*s2-drph1*r2*c1*c2)/R;
    double Rph2 = (r2*drph2-(drph2*c12+r2*s12)*r1*s2*s1-drph2*r1*c2*c1)/R;
    double Rph1ph2 = (-(drph1*r2-r1*drph2)*s1*s2*s12-(r1*r2+drph1*drph2)*s1*s2*c12-drph1*drph2*c1*c2-Rph1*Rph2)/R;
    
    double sigma = p1->rad+p2->rad;
    return epsilon/sigma/sigma*Rph1*Rph2-epsilon/sigma*(1-R/sigma)*Rph1ph2;
}

double secondTotalDth(int m, int count, double nm[][3], int np, particle p[], double epsilon){
    double result = 0;
    for(int n=0; n<np; n++) {
        if(n==m||!overlapQ(&(p[m]), &(p[n])))
            continue;
        result += VththS(count, nm, &(p[m]), &(p[n]), epsilon);
    }
    return result;
}
double secondTotalDph(int m, int count, double nm[][3], int np, particle p[], double epsilon){
    double result = 0;
    for(int n=0; n<np; n++) {
        if(n==m||!overlapQ(&(p[m]), &(p[n])))
            continue;
        result += VphphS(count, nm, &(p[m]), &(p[n]), epsilon);
    }
    return result;
}
double secondPartialDSinglethph(int m, int count, double nm[][3], int np, particle p[], double epsilon){
    double result = 0;
    for(int n=0; n<np; n++) {
        if(n==m||!overlapQ(&(p[m]), &(p[n])))
            continue;
        result += VthphS(count, nm, &(p[m]), &(p[n]), epsilon);
    }
    return result;
}

int indexColumn(int np, int row, int column) {
    return column*2*np+row;
}
void hessianMatrix(int count, double nm[][3], int np, particle p[], double epsilon, double *hessian){
    for(int m=0; m<np; m++) {
        hessian[indexColumn(np, 2*m, 2*m)] = secondTotalDth(m, count, nm, np, p, epsilon);
        hessian[indexColumn(np, 2*m+1, 2*m+1)] = secondTotalDph(m, count, nm, np, p, epsilon);
        hessian[indexColumn(np, 2*m, 2*m+1)] = secondPartialDSinglethph(m, count, nm, np, p, epsilon);
        hessian[indexColumn(np, 2*m+1, 2*m)] = hessian[indexColumn(np, 2*m, 2*m+1)];
        for(int n=m+1; n<np; n++) {
            if(!overlapQ(&(p[m]),&(p[n])))
                continue;
            hessian[indexColumn(np, 2*m, 2*n)] = VththP(count, nm, &(p[m]), &(p[n]), epsilon);
            hessian[indexColumn(np, 2*m, 2*n+1)] = VthphP(count, nm, &(p[m]), &(p[n]), epsilon);
            hessian[indexColumn(np, 2*m+1, 2*n)] = VthphP(count, nm, &(p[n]), &(p[m]), epsilon);
            hessian[indexColumn(np, 2*m+1, 2*n+1)] = VphphP(count, nm, &(p[m]), &(p[n]), epsilon);
            hessian[indexColumn(np, 2*n, 2*m)] = hessian[indexColumn(np, 2*m, 2*n)];
            hessian[indexColumn(np, 2*n+1, 2*m)] = hessian[indexColumn(np, 2*m, 2*n+1)];
            hessian[indexColumn(np, 2*n, 2*m+1)] = hessian[indexColumn(np, 2*m+1, 2*n)];
            hessian[indexColumn(np, 2*n+1, 2*m+1)] = hessian[indexColumn(np, 2*m+1, 2*n+1)];
        }
    }
}

void solveHessian(int count, double nm[][3], int np, particle p[], double epsilon, int *info) {
    int nn=0;
    int pi[np];
    for(int i=0; i<np; i++)
        for(int j=i+1; j<np; j++)
            if(overlapQ(&(p[i]), &(p[j]))){
                if(nn==0){
                    pi[0]=i;
                    pi[1]=j;
                    nn=2;
                    continue;
                }
                int k;
                for(k=0; k<nn; k++)
                    if(pi[k]==i)
                        break;
                if(k==nn){
                    pi[k]=i;
                    nn++;
                }
                for(k=0; k<nn; k++)
                    if(pi[k]==j)
                        break;
                if(k==nn){
                    pi[k]=j;
                    nn++;
                }
            }
    particle pp[nn];
    for(int i=0; i<nn; i++)
        pp[i] = p[pi[i]];
    
    double *hessian = malloc(2*nn*2*nn*sizeof(double));
    for(int i=0; i<2*nn*2*nn; i++)
        hessian[i] = 0;
    hessianMatrix(count, nm, nn, pp, epsilon, hessian);
    
    for(int j=0; j<np; j++){
        p[j].forceTh = 0;
        p[j].forcePh = 0;
    }
    addForces(count, nm, p, np, epsilon);
    double forces[2*nn];
    for(int i=0; i<nn; i++) {
        forces[2*i] = p[pi[i]].forceTh;
        forces[2*i+1] = p[pi[i]].forcePh;
    }
    
    double oldForces[2*nn];
    for(int i=0; i<2*nn; i++)
        oldForces[i]=forces[i];
    double *oldHessian = malloc(2*nn*2*nn*sizeof(double));
    for(int i=0; i<2*nn*2*nn; i++)
        oldHessian[i]=hessian[i];
    
    int n = 2*nn;
    int nrhs = 1;
    int lda = n;
    int ldb = n;
    int ipiv[n];
    dgesv_(&n,&nrhs,hessian,&lda,ipiv,forces,&ldb,info);
/*
    for(int i=0; i<2*nn; i++){
        double sum=0;
        for(int j=0; j<2*nn; j++)
            sum += oldHessian[i*2*nn+j]*forces[j];
        printf("%.15lf %.15lf\n", forces[i], sum-oldForces[i]);
     }
*/
    for(int j=0; j<np; j++) {
        p[j].forceTh = 0;
        p[j].forcePh = 0;
    }
    for(int i=0; i<nn; i++) {
        p[pi[i]].forceTh=forces[2*i];
        p[pi[i]].forcePh=forces[2*i+1];
    }
    
    free(oldHessian);
    free(hessian);
}
/*
int main(){
    double epsilon;
    double rad;
    int np;
    int count=0;
    
    FILE *eps;
    eps = fopen("epsilon.dat", "r");
    fscanf(eps, "%lf", &epsilon);
    fclose(eps);
    FILE *radFile=NULL;
    radFile = fopen("rad.dat","r");
    fscanf(radFile, "%lf", &rad);
    fclose(radFile);
    FILE *npFile=NULL;
    npFile = fopen("npts.dat","r");
    fscanf(npFile, "%i", &np);
    fclose(npFile);
    FILE *coeffFile=NULL;
    coeffFile = fopen("coeffOutTemp.dat","r");
    for(char c=getc(coeffFile); c!=EOF; c=getc(coeffFile))
        if(c=='\n')
            count++;
    count++;
    fclose(coeffFile);
    
    double nm[count][3];
    coeffFile = fopen("coeffOutTemp.dat","r");
    for(int i=0; i<count; i++)
        fscanf(coeffFile, "%lf %lf %lf", nm[i], nm[i]+1, nm[i]+2);
    fclose(coeffFile);
    
    particle p[np];
    FILE *configuration=NULL;
    FILE *radiusFile=NULL;
    configuration = fopen("configOutTemp.asc","r");
    radiusFile = fopen("radii.dat", "r");
    for (int j = 0; j < np; j++) {
        fscanf(configuration, "%lf %lf", &(p[j].th), &(p[j].ph));
        p[j].rr = radiusFromCoefficients(count, nm, p[j].th, p[j].ph);
        p[j].x = p[j].rr*sin(p[j].th)*cos(p[j].ph);
        p[j].y = p[j].rr*sin(p[j].th)*sin(p[j].ph);
        p[j].z = p[j].rr*cos(p[j].th);
        fscanf(radiusFile, "%lf", &(p[j].rad));
    }
    fclose(configuration);
    fclose(radiusFile);
    
    double L=0;//L represent one side boundary
    for(int i=0; i<np; i++) {
        if(fabs(p[i].x)>L)
            L = fabs(p[i].x);
        if(fabs(p[i].y)>L)
            L = fabs(p[i].y);
        if(fabs(p[i].z)>L)
            L = fabs(p[i].z);
    }
    L = 2*L;//enlarge the box to make sure the deformation is inside.
    
    int nPart[3];
    nPart[0] = (int)(L/(rad))+1;
    nPart[1] = (int)(L/(rad))+1;
    nPart[2] = (int)(L/(rad))+1;
    double partBoundX[nPart[0]];
    double partBoundY[nPart[1]];
    double partBoundZ[nPart[2]];
    for(int i=0; i<nPart[0]; i++)
        partBoundX[i] = -L+i*2*rad;
    for(int i=0; i<nPart[1]; i++)
        partBoundY[i] = -L+i*2*rad;
    for(int i=0; i<nPart[2]; i++)
        partBoundZ[i] = -L+i*2*rad;
    
    for (int j = 0; j < np; j++) {
        partitionOneParticle(&(p[j]), nPart,partBoundX, partBoundY, partBoundZ);
    }
    
    double *hessian = malloc(2*np*2*np*sizeof(double));
    for(int i=0; i<2*np*2*np; i++)
        hessian[i] = 0;
    hessianMatrix(count, nm, np, p, epsilon, hessian);
    
    addForces(count, nm, p, np, epsilon);
    double *oldForces = malloc(2*np*sizeof(double));
    for(int i=0; i<np; i++) {
        oldForces[2*i] = p[i].forceTh;
        oldForces[2*i+1] = p[i].forcePh;
    }
    
    int info;
    solveHessian(count, nm, np, p, epsilon, &info);
    
    double *forces = malloc(2*np*sizeof(double));
    for(int i=0; i<np; i++) {
        forces[2*i] = p[i].forceTh;
        forces[2*i+1] = p[i].forcePh;
    }
    for(int i=0; i<2*np; i++){
        double sum=0;
        for(int j=0; j<2*np; j++)
            sum += hessian[i*2*np+j]*forces[j];
        printf("%.15lf %.15lf\n", forces[i], sum-oldForces[i]);
    }
    
    free(hessian);
    free(oldForces);
    free(forces);
    
}
*/







