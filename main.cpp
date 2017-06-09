//
//  main.cpp
//  软件工程
//
//  Created by fanhao on 2017/6/7.
//  Copyright © 2017年 fanhao. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstring>
#define EPS (1E-7) //-------改为EPS---------
using namespace std;
int m;
struct circle {
    double x,y,r;
    int p; //0圆 -1直线 1点 2不存在
}C[1000]={0};
struct Root{
    double x1,x2;
    int p;//0两个解 1一个解 2无解
};



Root quadratic(double a,double b,double c){
    Root x;
    if (abs(a)<EPS && abs(b)<EPS) {
        x.p=2; return x;
    }
    if (abs(a)<EPS) {
        x.x1=-c/b,x.p=1; return x;
    }
    double derta=b*b-4*a*c;
    if (derta<0){
        x.p=2; return x;
    }
    if (derta==0){
        x.p=1; x.x1=-b/a/2.0; return x;
    }
    x.p=0;
    derta=sqrt(derta);
    x.x1=(-b-derta)/a/2.0,x.x2=(-b+derta)/a/2.0;
    return x;
}


circle cal(circle a,circle b,circle c){//计算与圆a、b、c相切的圆
    circle z; z.x=0,z.y=0,z.r=0;
    if (c.p==-1){//三条均为直线
        z.x=0,z.y=0,z.r=1,z.p=0;
        return z;
    }
    if (b.p==-1){
        if (a.x+b.x == 0){//两条直线相对
            z.x=0,z.y=0,z.r=1;
            if (abs(c.x*c.x+c.y*c.y - (1+c.r)*(1+c.r))<EPS){
                z.p=0; return z;
            }
            else {
                z.p=-2; return z;
            }
        }
        else{//直线相交 圆心在对角线上 化为二元一次方程
            if ((a.x+b.x)*(a.y+b.y)==1){//y=x
                if (a.x+b.x == 1){//r=1-x
                    Root r=quadratic(1,2*(1+c.r-c.x-c.y),c.x*c.x+c.y*c.y-c.r*c.r-2*c.r-1);
                    if (r.p==2) {
                        z.p=2; return z;
                    }
                    else if (r.p==1){
                        z.x=r.x1,z.y=z.x,z.r=1-z.x,z.p=0;
                        return z;
                    }
                    else {
                        if (r.x1>=0 && r.x1<=1) z.x=r.x1,z.y=z.x,z.r=1-z.x,z.p=0;
                        else z.x=r.x2,z.y=z.x,z.r=1-z.x,z.p=0;
                        return z;
                    }
                }
                else{//r=x+1
                    Root r=quadratic(1,2*(-1-c.r-c.x-c.y),c.x*c.x+c.y*c.y-c.r*c.r-2*c.r-1);
                    if (r.p==2) {
                        z.p=2; return z;
                    }
                    else if (r.p==1){
                        z.x=r.x1,z.y=z.x,z.r=1+z.x,z.p=0;
                        return z;
                    }
                    else {
                        if (r.x1<=0 && r.x1>=-1) z.x=r.x1,z.y=z.x,z.r=1+z.x,z.p=0;
                        else z.x=r.x2,z.y=z.x,z.r=1+z.x,z.p=0;
                        return z;
                    }
                }
            }
            else {//y=-x
                if (a.x+b.x == 1){//r=1-x
                    Root r=quadratic(1,2*(1+c.r-c.x+c.y),c.x*c.x+c.y*c.y-c.r*c.r-2*c.r-1);
                    if (r.p==2) {
                        z.p=2; return z;
                    }
                    else if (r.p==1){
                        z.x=r.x1,z.y=-z.x,z.r=1-z.x,z.p=0;
                        return z;
                    }
                    else {
                        if (r.x1>=0 && r.x1<=1) z.x=r.x1,z.y=-z.x,z.r=1-z.x,z.p=0;
                        else z.x=r.x2,z.y=-z.x,z.r=1-z.x,z.p=0;
                        return z;
                    }
                }
                else{//r=1+x
                    Root r=quadratic(1,2*(-1-c.r-c.x+c.y),c.x*c.x+c.y*c.y-c.r*c.r-2*c.r-1);
                    if (r.p==2) {
                        z.p=2; return z;
                    }
                    else if (r.p==1){
                        z.x=r.x1,z.y=-z.x,z.r=1+z.x,z.p=0;
                        return z;
                    }
                    else {
                        if (r.x1<=0 && r.x1>=-1) z.x=r.x1,z.y=-z.x,z.r=1+z.x,z.p=0;
                        else z.x=r.x2,z.y=-z.x,z.r=1+z.x,z.p=0;
                        return z;
                    }
                }
            }
        }
    }
    if (a.p==-1){
        if (a.x!=0){
            if (b.y==c.y){
                z.p=2; return z;
            }
            if (a.x==1){
                double A=b.r*b.r-c.r*c.r+c.x*c.x-b.x*b.x+c.y*c.y-b.y*b.y;
                double B=(c.x-b.x-c.r+b.r)/(b.y-c.y);
                double C=(A/2-(c.r-b.r))/(c.y-b.y);
                double a1,a2,a3;
                a1=B*B;
                a2=2*B*C+2-2*b.x-2*B*b.y+2*b.r;
                a3=C*C-1-2*C*b.y-2*b.r+b.x*b.x+b.y*b.y-b.r*b.r;
                Root r=quadratic(a1, a2, a3);
                if (r.p==2) {
                    z.p=2; return z;
                }
                else if (r.p==1){
                    z.x=r.x1,z.y=B*z.x+C,z.r=1-z.x,z.p=0;
                    return z;
                }
                else {
                    z.x=r.x1,z.y=B*z.x+C,z.r=1-z.x,z.p=0;
                    if (z.x+z.r >1 || z.x-z.r<-1 || z.y+z.r >1 || z.y-z.r<-1)
                        z.x=r.x2,z.y=B*z.x+C,z.r=1-z.x,z.p=0;
                    return z;
                }
            }
            else {
                double A=b.r*b.r-c.r*c.r+c.x*c.x-b.x*b.x+c.y*c.y-b.y*b.y;
                double B=(c.x-b.x+c.r-b.r)/(b.y-c.y);
                double C=(A/2-(c.r-b.r))/(c.y-b.y);
                double a1,a2,a3;
                a1=B*B;
                a2=2*B*C-2-2*b.x-2*B*b.y-2*b.r;
                a3=C*C-1-2*C*b.y-2*b.r+b.x*b.x+b.y*b.y-b.r*b.r;
                Root r=quadratic(a1, a2, a3);
                if (r.p==2) {
                    z.p=2; return z;
                }
                else if (r.p==1){
                    z.x=r.x1,z.y=B*z.x+C,z.r=1+z.x,z.p=0;
                    return z;
                }
                else {
                    z.x=r.x1,z.y=B*z.x+C,z.r=1+z.x,z.p=0;
                    if (z.x+z.r >1 || z.x-z.r<-1 || z.y+z.r >1 || z.y-z.r<-1)
                        z.x=r.x2,z.y=B*z.x+C,z.r=1+z.x,z.p=0;
                    return z;
                }
            }
        }
        else {
            if (b.x==c.x){
                z.p=2; return z;
            }
            if (a.y==1){
                double A=b.r*b.r-c.r*c.r+c.x*c.x-b.x*b.x+c.y*c.y-b.y*b.y;
                double B=(c.y-b.y-c.r+b.r)/(b.x-c.x);
                double C=(A/2-(c.r-b.r))/(c.x-b.x);
                double a1,a2,a3;
                a1=B*B;
                a2=2*B*C+2-2*b.y-2*B*b.x+2*b.r;
                a3=C*C-1-2*C*b.x-2*b.r+b.x*b.x+b.y*b.y-b.r*b.r;
                Root r=quadratic(a1, a2, a3);
                if (r.p==2) {
                    z.p=2; return z;
                }
                else if (r.p==1){
                    z.y=r.x1,z.x=B*z.x+C,z.r=1-z.y,z.p=0;
                    return z;
                }
                else {
                    z.y=r.x1,z.x=B*z.y+C,z.r=1-z.y,z.p=0;
                    if (z.x+z.r >1 || z.x-z.r<-1 || z.y+z.r >1 || z.y-z.r<-1)
                        z.y=r.x2,z.x=B*z.y+C,z.r=1-z.y,z.p=0;
                    return z;
                }
            }
            else {
                double A=b.r*b.r-c.r*c.r+c.x*c.x-b.x*b.x+c.y*c.y-b.y*b.y;
                double B=(c.y-b.y+c.r-b.r)/(b.x-c.x);
                double C=(A/2-(c.r-b.r))/(c.x-b.x);
                double a1,a2,a3;
                a1=B*B;
                a2=2*B*C-2-2*b.y-2*B*b.x-2*b.r;
                a3=C*C-1-2*C*b.x-2*b.r+b.x*b.x+b.y*b.y-b.r*b.r;
                Root r=quadratic(a1, a2, a3);
                if (r.p==2) {
                    z.p=2; return z;
                }
                else if (r.p==1){
                    z.y=r.x1,z.x=B*z.x+C,z.r=1+z.y,z.p=0;
                    return z;
                }
                else {
                    z.y=r.x1,z.x=B*z.y+C,z.r=1+z.y,z.p=0;
                    if (z.x+z.r >1 || z.x-z.r<-1 || z.y+z.r >1 || z.y-z.r<-1)
                        z.y=r.x2,z.x=B*z.y+C,z.r=1+z.y,z.p=0;
                    return z;
                }
            }
        }
    }
    double r[4][4]={0},u[4]={0};
    r[1][1]=2*(b.x-a.x), r[1][2]=2*(b.y-a.y), r[1][3]=2*(b.r-a.r);
    r[2][1]=2*(c.x-a.x), r[2][2]=2*(c.y-a.y), r[2][3]=2*(c.r-a.r);
    r[3][1]=2*(c.x-b.x), r[3][2]=2*(c.y-b.y), r[3][3]=2*(c.r-b.r);
    
    u[1]=a.r*a.r-b.r*b.r + b.x*b.x-a.x*a.x + b.y*b.y-a.y*a.y;
    u[2]=a.r*a.r-c.r*c.r + c.x*c.x-a.x*a.x + c.y*c.y-a.y*a.y;
    u[3]=b.r*b.r-c.r*c.r + c.x*c.x-b.x*b.x + c.y*c.y-b.y*b.y;
    
    //--------------高斯消元
    double t[4]={0};
    for(int k=1;k<3;k++){//选主元
        double ab_max=-1; int max_ik=0;
        for(int i=k;i<=3;i++){
            if(abs(r[i][k])>ab_max){
                ab_max=abs(r[i][k]); max_ik=i;
            }
        }
        if(ab_max<EPS){//0矩阵情况
            z.p=2; return z;
        }
        else if(max_ik!=k){//是否是当前行，不是交换
            double temp;
            for(int j=1;j<=3;j++){
                temp=r[max_ik][j]; r[max_ik][j]=r[k][j]; r[k][j]=temp;
            }
            temp=u[max_ik]; u[max_ik]=u[k]; u[k]=temp;
        }
        //消元计算
        for(int i=k+1;i<=3;i++){
            r[i][k]/=r[k][k];
            for(int j=k+1;j<=3;j++) r[i][j]-=r[i][k]*r[k][j];
            u[i]-=r[i][k]*u[k];
        }
    }
    if (abs(r[3][3])<EPS) {
        z.p=2; return z;
    }
    else {//回代求解
        u[3]=u[3]/r[3][3];
        for(int i=2;i>0;i--){
            t[i]=u[i];
            for(int j=i+1;j<=3;j++) t[i]-=r[i][j]*t[j];
            t[i]/=r[i][i];
        }
    }
    z.p=0;
    z.x=t[1],z.y=t[2],z.r=t[3];
    return z;
}

bool isIntersect(circle a,circle b){
    if (a.p==-1){
        if (a.x!=0){
            if (abs(b.x-a.x)-b.r > -EPS) return 0;
            else return 1;
        }
        else {
            if (abs(b.y-a.y)-b.r > -EPS) return 0;
            else return 1;
        }
    }
    /*
    else if (a.p==1){
        if ((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y) >= b.r*b.r) return 0;
        else return 1;
        
    }*/
    else {
        if ((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y) - (a.r+b.r)*(a.r+b.r) > -EPS) return 0;
        else return 1;
    }
}

void solve(int num){
    circle maxx;
    maxx.x=0,maxx.y=0,maxx.r=0;
    
    for (int i=0;i<num;++i)
        for (int j=i+1;j<num;++j)
            for (int k=j+1;k<num;++k){
                circle ans=cal(C[i],C[j],C[k]);
                if (ans.p==2) continue;
                if (ans.x>1 || ans.x<-1 || ans.y>1 || ans.y<-1  || ans.r<=0) continue;
                
                int flag=1;
                for (int p=0;p<num;++p){
                    if (abs(ans.x-C[p].x)<EPS && abs(ans.y-C[p].y)<EPS &&
                        abs(ans.r-C[p].r)<EPS && ans.p==C[p].p) {//判断重复
                        flag=0; break;
                    }
                    if (isIntersect(C[p],ans)) {//判断是否相交
                        flag=0; break;
                    }
                    
                }
                if (flag && maxx.r<ans.r){
                    maxx=ans;
                }
            }
    C[num]=maxx;
}

int main() {
    cin>>m;
    memset(C, 0, sizeof(C));
    C[0].x=1,  C[0].y=0,  C[0].p=-1;
    C[1].x=0,  C[1].y=1,  C[1].p=-1;
    C[2].x=-1, C[2].y=0,  C[2].p=-1;
    C[3].x=0,  C[3].y=-1, C[3].p=-1;
    
    for (int i=4;i<m+4;++i) solve(i);
    
    double sum=0;
    for (int i=4;i<m+4;++i){
        printf("%d: The coordinate is (%.4lf,%.4lf).  r=%.4lf\n",i-3,C[i].x,C[i].y,C[i].r);
        sum+=C[i].r*C[i].r;
    }
    printf("r^2=%.4lf\n",sum);
    return 0;
}
