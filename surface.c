//
//  surface.c
//  softBidispersityMetricJamming
//
//  Created by Zhaoyu Xie on 4/8/20.
//  Copyright © 2020 Zhaoyu Xie. All rights reserved.
//

#ifndef PI
#define PI 3.141592653589793
#endif

#include "energy.h"
#include "particle.h"
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

//only has zero term for each order, where theta is the only parameter.
double tesseralHarmonicsZ(int n, int m, double th, double ph) {
    if(m!=0||n>10) {
        printf("Error: harmonics doesn't exist for (%d, %d)!\n", n, m);
        exit(1);
    }
    switch(n) {
        case 0:
            return 1.0/2/sqrt(PI);
            break;
        case 1:
            return 1.0/2*sqrt(3/PI)*cos(th);
            break;
        case 2:
            return 1.0/8*sqrt(5/PI)*(1+3*cos(2*th));
            break;
        case 3:
            return 1.0/16*sqrt(7/PI)*(3*cos(th)+5*cos(3*th));
            break;
        case 4:
            return 1.0/128*sqrt(9/PI)*(9+20*cos(2*th)+35*cos(4*th));
            break;
        case 5:
            return 1.0/256*sqrt(11/PI)*(30*cos(th)+35*cos(3*th)+63*cos(5*th));
            break;
        case 6:
            return 1.0/32*sqrt(13/PI)*(-5+105*pow(cos(th),2)-315*pow(cos(th),4)+231*pow(cos(th),6));
            break;
        case 7:
            return 1.0/32*sqrt(15/PI)*(-35*cos(th)+315*pow(cos(th),3)-693*pow(cos(th),5)+429*pow(cos(th),7));
            break;
        case 8:
            return 1.0/256*sqrt(17/PI)*(35-1260*pow(cos(th),2)+6930*pow(cos(th),4)-12012*pow(cos(th),6)+6435*pow(cos(th),8));
            break;
        case 9:
            return 1.0/65536*sqrt(19/PI)*(4410*cos(th)+4620*cos(3*th)+5148*cos(5*th)+6435*cos(7*th)+12155*cos(9*th));
            break;
        case 10:
            return 1.0/512*sqrt(21/PI)*(-63+3465*pow(cos(th),2)-30030*pow(cos(th),4)+90090*pow(cos(th),6)-109395*pow(cos(th),8)+46189*pow(cos(th),10));
            break;
        default:
            exit(1);
    }
}
double tesseralHarmonicsZDth(int n, int m, double th, double ph) {
    if(m!=0||n>10) {
        printf("Error: harmonics doesn't exist for (%d, %d)!\n", n, m);
        exit(1);
    }
    switch(n) {
        case 0:
            return 0;
            break;
        case 1:
            return -1.0/2*sqrt(3/PI)*sin(th);
            break;
        case 2:
            return -3.0/2*sqrt(5/PI)*cos(th)*sin(th);
            break;
        case 3:
            return -3.0/16*sqrt(7/PI)*(sin(th)+5*sin(3*th));
            break;
        case 4:
            return 5.0/4*sqrt(9/PI)*cos(th)*(3-7*pow(cos(th),2))*sin(th);
            break;
        case 5:
            return -15.0/256*sqrt(11/PI)*(2*sin(th)+7*sin(3*th)+21*sin(5*th));
            break;
        case 6:
            return -21.0/512*sqrt(13/PI)*(5*sin(2*th)+12*sin(4*th)+33*sin(6*th));
            break;
        case 7:
            return -7.0/2048*sqrt(15/PI)*(25*sin(th)+81*sin(3*th)+165*sin(5*th)+429*sin(7*th));
            break;
        case 8:
            return -9.0/4096*sqrt(17/PI)*(70*sin(2*th)+154*sin(4*th)+286*sin(6*th)+715*sin(8*th));
            break;
        case 9:
            return -45.0/65536*sqrt(19/PI)*(98*sin(th)+308*sin(3*th)+572*sin(5*th)+1001*sin(7*th)+2431*sin(9*th));
            break;
        case 10:
            return -55.0/131072*sqrt(21/PI)*(294*sin(2*th)+624*sin(4*th)+1053*sin(6*th)+1768*sin(8*th)+4199*sin(10*th));
            break;
        default:
            exit(1);
    }
}
double tesseralHarmonicsZDph(int n, int m, double th, double ph) {
    if(m!=0||n>10) {
        printf("Error: harmonics doesn't exist for (%d, %d)!\n", n, m);
        exit(1);
    }
    return 0;
}

double radiusFromCoefficients(int count, double nm[][3], double th, double ph){
    double rr = 0;
    for(int i=0; i<count; i++) {
        rr += nm[i][2]*tesseralHarmonicsZ((int)(nm[i][0]), (int)(nm[i][1]), th, ph);
    }
    return rr;
}

double areaInt(double th, void *params) {
    double *coeff = (double *)params;
    int count = (int)(coeff[0]);
    double rr = 0;
    double ph = 0;//the value of ph doesn't matter
    double drth = 0;

    for(int i=0; i<count; i++) {
        rr += coeff[3*i+3]*tesseralHarmonicsZ((int)(coeff[3*i+1]), (int)(coeff[3*i+2]), th, ph);
        drth += coeff[3*i+3]*tesseralHarmonicsZDth((int)(coeff[3*i+1]), (int)(coeff[3*i+2]), th, ph);
    }
    return rr*sqrt(sin(th)*sin(th)*(rr*rr+drth*drth));
}
double areaFromCoefficients(int count, double nm[][3]) {
    double coeff[3*count+1];
    coeff[0] = count;
    for(int i=0; i<count; i++)
        for(int j=0; j<3; j++)
            coeff[3*i+j+1] = nm[i][j];
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    F.function = &areaInt;
    F.params = coeff;
    gsl_integration_qags(&F, 0, PI, 0, 1e-7, 1000, w, &result, &error);
    //gsl_integration_qag(&F, 0, PI, 0, 1e-7, 1000, 6, w, &result, &error);
    gsl_integration_workspace_free (w);
    return 2*PI*result;
}

double areaIntD(double th, void *params) {
    double *coeff = (double *)params;
    int count = (int)(coeff[0]);
    int c = (int)(coeff[1]);//c starts from 0
    double rr = 0;
    double ph = 0;//the value of ph doesn't matter
    double drth = 0;
    
    for(int i=0; i<count; i++) {
        rr += coeff[3*i+4]*tesseralHarmonicsZ((int)(coeff[3*i+2]), (int)(coeff[3*i+3]), th, ph);
        drth += coeff[3*i+4]*tesseralHarmonicsZDth((int)(coeff[3*i+2]), (int)(coeff[3*i+3]), th, ph);
    }
    double Yc = tesseralHarmonicsZ((int)(coeff[3*c+2]), (int)(coeff[3*c+3]), th, ph);
    double dYcth = tesseralHarmonicsZDth((int)(coeff[3*c+2]), (int)(coeff[3*c+3]), th, ph);
    return Yc*sqrt((drth*drth+rr*rr)*sin(th)*sin(th))+rr*(drth*dYcth+rr*Yc)*sin(th)*sin(th)/sqrt((drth*drth+rr*rr)*sin(th)*sin(th));
}
double areaFromCoefficientsD(int count, double nm[][3], int c) {
    double coeff[3*count+2];
    coeff[0] = count;
    coeff[1] = c;//c starts from 0
    for(int i=0; i<count; i++)
        for(int j=0; j<3; j++)
            coeff[3*i+j+2] = nm[i][j];
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    F.function = &areaIntD;
    F.params = coeff;
    gsl_integration_qags(&F, 0, PI, 0, 1e-7, 1000, w, &result, &error);
    //gsl_integration_qag(&F, 0, PI, 0, 1e-7, 1000, 6, w, &result, &error);
/*
    gsl_set_error_handler_off();
    int status = gsl_integration_qags(&F, 0, PI, 0, 1e-7, 1000, w, &result, &error);
    if(status) {
        if(status == GSL_EROUND) {
            for(int j=0; j<count; j++)
                printf("%d, %d, %.15lf\n", (int)(nm[j][0]), (int)(nm[j][1]), nm[j][2]);
            exit(-1);
        }
        else {
            printf ("gsl_integration failed, gsl_errno=%d\n", status);
            exit(-1);
        }
    }
    gsl_set_error_handler (NULL);
*/
    gsl_integration_workspace_free (w);
    return 2*PI*result;
}

double volumeInt(double th, void *params) {
    double *coeff = (double *)params;
    int count = (int)(coeff[0]);
    double rr = 0;
    double ph = 0;//the value of ph doesn't matter
    double drth = 0;
    
    for(int i=0; i<count; i++) {
        rr += coeff[3*i+3]*tesseralHarmonicsZ((int)(coeff[3*i+1]), (int)(coeff[3*i+2]), th, ph);
        drth += coeff[3*i+3]*tesseralHarmonicsZDth((int)(coeff[3*i+1]), (int)(coeff[3*i+2]), th, ph);
    }
    return rr*rr*drth*sin(th)*sin(th)*cos(th)-rr*rr*rr*pow(sin(th),3);
}
double volumeFromCoefficients(int count, double nm[][3]) {
    double coeff[3*count+1];
    coeff[0] = count;
    for(int i=0; i<count; i++)
        for(int j=0; j<3; j++)
            coeff[3*i+j+1] = nm[i][j];
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    F.function = &volumeInt;
    F.params = coeff;
    gsl_integration_qags(&F, PI, 0, 0, 1e-7, 1000, w, &result, &error);
    //gsl_integration_qag(&F, PI, 0, 0, 1e-7, 1000, 6, w, &result, &error);
    gsl_integration_workspace_free (w);
    return PI*result;
}

double volumeIntD(double th, void *params) {
    double *coeff = (double *)params;
    int count = (int)(coeff[0]);
    int c = (int)(coeff[1]);//c starts from 0
    double rr = 0;
    double ph = 0;//the value of ph doesn't matter
    double drth = 0;
    
    for(int i=0; i<count; i++) {
        rr += coeff[3*i+4]*tesseralHarmonicsZ((int)(coeff[3*i+2]), (int)(coeff[3*i+3]), th, ph);
        drth += coeff[3*i+4]*tesseralHarmonicsZDth((int)(coeff[3*i+2]), (int)(coeff[3*i+3]), th, ph);
    }
    double Yc = tesseralHarmonicsZ((int)(coeff[3*c+2]), (int)(coeff[3*c+3]), th, ph);
    double dYcth = tesseralHarmonicsZDth((int)(coeff[3*c+2]), (int)(coeff[3*c+3]), th, ph);
    return (2*rr*drth*Yc+rr*rr*dYcth)*sin(th)*sin(th)*cos(th)-3*rr*rr*Yc*pow(sin(th),3);
}
double volumeFromCoefficientsD(int count, double nm[][3], int c) {
    double coeff[3*count+2];
    coeff[0] = count;
    coeff[1] = c;//c starts from 0
    for(int i=0; i<count; i++)
        for(int j=0; j<3; j++)
            coeff[3*i+j+2] = nm[i][j];
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    F.function = &volumeIntD;
    F.params = coeff;
    gsl_integration_qags(&F, PI, 0, 0, 1e-6, 1000, w, &result, &error);
    //gsl_integration_qag(&F, PI, 0, 0, 1e-7, 1000, 6, w, &result, &error);
    gsl_integration_workspace_free (w);
    return PI*result;
}

void volumeNormal(int count, double nm[][3], double nv[]) {
    for(int c=0; c<count; c++)
        nv[c] = volumeFromCoefficientsD(count, nm, c);
    double norm = 0;
    for(int c=0; c<count; c++)
        norm += nv[c]*nv[c];
    norm = sqrt(norm);
    for(int c=0; c<count; c++)
        nv[c] = nv[c]/norm;
}
double volumeNormalNorm(int count, double nm[][3]) {
    double nv[count];
    for(int c=0; c<count; c++)
        nv[c] = volumeFromCoefficientsD(count, nm, c);
    double norm = 0;
    for(int c=0; c<count; c++)
        norm += nv[c]*nv[c];
    norm = sqrt(norm);
    return norm;
}
/*
double areaFromCoefficients(int count, double nm[][3]){
    double area;
    FILE *coeffTemp = fopen("coeffTemp","w");
    for(int i=0; i<count; i++)
        fprintf(coeffTemp, "%d %d %.15lf\n", (int)(nm[i][0]), (int)(nm[i][1]), nm[i][2]);
    fclose(coeffTemp);
    system("/Applications/Mathematica.app/Contents/MacOS/MathKernel -script calculateArea.wl");
    FILE *areaFile = fopen("area", "r");
    fscanf(areaFile, "%lf", &area);
    fclose(areaFile);
    system("rm coeffTemp area");
    return area;
}

void areaFromCoefficientsDAll(int count, double nm[][3], double areaD[]){
    FILE *coeffTemp = fopen("coeffTemp","w");
    for(int i=0; i<count; i++)
        fprintf(coeffTemp, "%d %d %.15lf\n", (int)(nm[i][0]), (int)(nm[i][1]), nm[i][2]);
    fclose(coeffTemp);
    system("/Applications/Mathematica.app/Contents/MacOS/MathKernel -script calculateAreaD.wl");
    FILE *areaDFile = fopen("areaD", "r");
    for(int i=0; i<count; i++)
        fscanf(areaDFile, "%lf", areaD+i);
    fclose(areaDFile);
    system("rm coeffTemp areaD");
}
*/
/*
int main(){
    int count = 0;
    FILE *coeffFile=NULL;
    coeffFile = fopen("coeffOutEmono100.dat","r");
    if (coeffFile) {
        for(char c=getc(coeffFile); c!=EOF; c=getc(coeffFile))
            if(c=='\n')
                count++;
        count++;
        fclose(coeffFile);
    } else {
        printf("coeffFile pointer is null\n");
        exit(1);
    }
    
    double nm[count][3];
    coeffFile = fopen("coeffOutEmono100.dat","r");
    if (coeffFile) {
        for(int i=0; i<count; i++)
            fscanf(coeffFile, "%lf %lf %lf", nm[i], nm[i]+1, nm[i]+2);
        fclose(coeffFile);
    } else {
        printf("coeffFile pointer is null\n");
        exit(1);
    }
    
    for(int i=0; i<count; i++)
        printf("%d, %d, %lf\n", (int)(nm[i][0]), (int)(nm[i][1]), nm[i][2]);
    
    double nv[count];
    volumeNormal(count, nm, nv);
    for(int c=0; c<count; c++)
        printf("%.15lf\n", nv[c]);
    return 0;
}
*/



