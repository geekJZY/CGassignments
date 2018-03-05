#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;


///Define Objects
class Quadrics{
public:
    int a02,a12,a22,a21,a00;
    double s0, s1, s2;
    Vec3f n0, n1, n2;
    Vec4f ColorH0, ColorH1, ColorH2, ColorH3;
    Vec3f pc;
    bool checkIn(Vec3f pos);
    double rayTracer(Vec3f pe, Vec3f npe);
    Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3f pc_, Vec4f ColorH0_, Vec4f ColorH1_, Vec4f ColorH2_, Vec4f ColorH3_);
    Vec3f nH(Vec3f Ph);
};

//Get the normal vector
Vec3f Quadrics::nH(Vec3f Ph){
    return (Ph - pc)/norm(Ph - pc);
}

Quadrics::Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3f pc_, Vec4f ColorH0_, Vec4f ColorH1_, Vec4f ColorH2_, Vec4f ColorH3_){
    a02 = a02_;
    a12 = a12_;
    a22 = a22_;
    a21 = a21_;
    a00 = a00_;
    s0 = s0_;
    s1 = s1_;
    s2 = s2_;
    pc = pc_;
    ColorH0 = ColorH0_;
    ColorH1 = ColorH1_;
    ColorH2 = ColorH2_;
    ColorH3 = ColorH3_;
    n0 = Vec3f(1,0,0);
    n1 = Vec3f(0,1,0);
    n2 = Vec3f(0,0,1);

}

bool Quadrics::checkIn(Vec3f pos){
    double F = pow((pos-pc).dot(n0)/s0,2)*a02+pow((pos-pc).dot(n1)/s1,2)*a12+pow((pos-pc).dot(n2)/s2,2)*a22+a00;
    if(F<0) return true;
    else return false;
}

double Quadrics::rayTracer(Vec3f pe, Vec3f npe){
    double A = a02*pow(n0.dot(npe)/s0,2)+a12*pow(n1.dot(npe)/s1,2)+a22*pow(n2.dot(npe)/s2,2);
    double B = a02*(2*n0.dot(npe))*(n0.dot(pe-pc))/(s0*s0)+
               a12*(2*n1.dot(npe))*(n1.dot(pe-pc))/(s1*s1)+
               a22*(2*n2.dot(npe))*(n2.dot(pe-pc))/(s2*s2)+
               a21*n2.dot(npe)/s2;
    double C = a02*pow(n0.dot(pe-pc)/s0,2) +
               a12*pow(n1.dot(pe-pc)/s1,2) +
               a22*pow(n2.dot(pe-pc)/s2,2) +
               a21*n2.dot(pe-pc)/s2 +
               a00;
//    cout << "A is " << A <<" B is " << B <<" C is " << C << endl;
    double delta = pow(B,2)-4*A*C;
    if(delta < 0) return -1;
    else return ((-B-sqrt(pow(B,2)-4*A*C))/(2*A));
}

class Planes{
public:
    Vec3f pc;
    Vec3f n2;
    Vec4f ColorH0, ColorH1, ColorH2, ColorH3;
    Planes(Vec3f n2_, Vec3f pc_, Vec4f ColorH0_, Vec4f ColorH1_, Vec4f ColorH2_, Vec4f ColorH3_);
    bool checkIn(Vec3f pos);
    Vec3f nH(void);
    double rayTracer(Vec3f pe, Vec3f npe);
};

Vec3f Planes::nH(void){
    return n2;
}

Planes::Planes(Vec3f n2_, Vec3f pc_, Vec4f ColorH0_, Vec4f ColorH1_, Vec4f ColorH2_, Vec4f ColorH3_){
    n2 = n2_;
    pc = pc_;
    ColorH0 = ColorH0_;
    ColorH1 = ColorH1_;
    ColorH2 = ColorH2_;
    ColorH3 = ColorH3_;
}

bool Planes::checkIn(Vec3f pos){
    double F = n2.dot(pos - pc);
    if(F<0) return true;
    else return false;
}

double Planes::rayTracer(Vec3f pe, Vec3f npe){
    return n2.dot(pc-pe)/n2.dot(npe);
}

///Define Light Source
class lightSource{
public:
    enum Mode
    {
      SPOT,
      DIRECTIONAL
    };

    Vec3f Pl, nLH_D; //light point OR light Direction
    Vec4f ColorL;   //color of light
    //params for color calculation
    double ks, kd, k0, kb;
    Mode lightMode;
    lightSource(Vec3f Pl, Vec3f nLH_, Vec4f ColorL, double ks, double kd, double k0, double kb, Mode lightMode);
    //H0 ambient, H1 diffuse, H2 specular
    Vec3b colorDeter(Vec3f Pe, Vec3f Ph, Vec3f nH, Vec4f ColorH0, Vec4f ColorH1, Vec4f ColorH2, Vec4f ColorH3);
};

lightSource::lightSource(Vec3f Pl_, Vec3f nLH_, Vec4f ColorL_, double ks_, double kd_, double k0_, double kb_, Mode lightMode_){
    Pl = Pl_;
    ColorL = ColorL_;
    ks = ks_;
    kd = kd_;
    k0 = k0_;
    kb = kb_;
    nLH_D = nLH_;
    lightMode = lightMode_;
}

Vec3b lightSource::colorDeter(Vec3f Pe, Vec3f Ph, Vec3f nH, Vec4f ColorH0, Vec4f ColorH1, Vec4f ColorH2, Vec4f ColorH3){
    double s, t, b; //the coefficient for specular, diffuse light and hale light.
    double cosRE, cosEH; //cos(angle between reflect light and eye)
    Vec3f rHE, nLH, nHE;
    Vec4f ColorPE;
    if(lightMode == SPOT) nLH = (Pl - Ph)/norm(Pl - Ph);
    else if(lightMode == DIRECTIONAL) nLH = nLH_D;
    nHE = (Ph - Pe)/norm(Ph - Pe);
    rHE = nLH - 2*(nLH.dot(nH))*nH;
    cosRE = rHE.dot(nHE);
    cosEH = -nHE.dot(nH);
    s = std::max((pow(cosRE,3)-0.9)/0.1,0.0);
    t = std::max(cosRE,0.0);
    b = std::max(((1-cosEH)-0.92)/0.08,0.0);
    if(lightMode == SPOT) ColorPE = ks*s*ColorH2.mul(ColorL) + kd*t*ColorH1.mul(ColorL)/(1+norm(Pl-Ph)) + k0*ColorH0 + b*kb*ColorH3;
    else if(lightMode == DIRECTIONAL) ColorPE = ks*s*ColorH2.mul(ColorL) + kd*t*ColorH1.mul(ColorL) + k0*ColorH0 + b*kb*ColorH3;

    //output and test
    static int i = 0;
    if((t>i/20.0)&&(t<(i+1)/20.0)){
        i ++;
        cout << "ColorPE (" << int(ColorPE[0]*255/ColorPE[3]) << ", " << int(ColorPE[1]*255/ColorPE[3])<<", " << int(ColorPE[2]*255/ColorPE[3]) <<")" << endl;
        cout << "t is " << int(t*20) << "  Color0 is " << k0*ColorH0 << " ColorD is " << kd*t*ColorH1.mul(ColorL)/(1+norm(Pl-Ph)) << endl;
    }

    return Vec3b(char(ColorPE[0]*255/ColorPE[3]),char(ColorPE[1]*255/ColorPE[3]),char(ColorPE[2]*255/ColorPE[3]));
}



Vec3f cameraPhGrab(Vec3f pe, Vec3f npe, double t){
    return (pe+npe*t);
}

int traceTimeComp(vector<double> traceVecTime){
    double timeMin = -1;
    int minIdx = -1;
    int idx = 0;
    while(traceVecTime[idx] <= 0){
        idx ++;
        if(idx == traceVecTime.size()) return -1;
    }
    timeMin = traceVecTime[idx];
    minIdx = idx;
    for(;idx < traceVecTime.size(); idx ++){
        if((traceVecTime[idx] < timeMin)&&(traceVecTime[idx] > 0)){
            timeMin = traceVecTime[idx];
            minIdx = idx;
        }
    }
    return minIdx;
}

int main(int argc, char *argv[]){

    vector<Scalar> color{Scalar(255,224,147),Scalar(93,66,255),Scalar(255,255,255),Scalar(230,180,80)};
    Mat img(Size(640,480),CV_8UC3,color[0]);
    double Sx = 16, Sy = 12;
    Vec3f pos;
    Vec3f Ph;
    int objectIdx;
    int xSample = 2, ySample = 2;
    vector<int> sample(5);
    double randX,randY;
    int M = img.cols, N = img.rows;

    //Graphics
    Quadrics sphere1(1,1,1,0,-1,5,5,5,Vec3f(0,0,6),Vec4f(53,36,155,255),Vec4f(113,86,255,255),Vec4f(1930,1660,2550,2550),Vec4f(0,0,0,255));
    Quadrics sphere2(1,1,1,0,-1,3,3,3,Vec3f(0,12,4),Vec4f(150,50,50,255),Vec4f(205,100,100,255),Vec4f(2550,2000,2000,2550),Vec4f(0,0,0,255));
    Planes plane1(Vec3f(0,0,1),Vec3f(0,0,0),Vec4f(100,50,50,255),Vec4f(220,170,170,255),Vec4f(2300,1800,1800,2550),Vec4f(10,10,10,255));
    Planes plane2(Vec3f(-1,0,0),Vec3f(12,0,0),Vec4f(50,100,50,255),Vec4f(170,220,170,255),Vec4f(1800,2300,1800,2550),Vec4f(10,10,10,255));
    Planes plane3(Vec3f(0,1,0),Vec3f(0,-9,0),Vec4f(50,50,100,255),Vec4f(170,170,220,255),Vec4f(1800,1800,2300,2550),Vec4f(10,10,10,255));
    lightSource lighter(Vec3f(-10,40,30), Vec3f(-0.5774,0.5774,0.5774),Vec4f(135,135,135,135),0.5,0.5,0.5,0.5, lightSource::SPOT);

    //camera
    Vec3f pe(-40,0,6);
    Vec3f vUp(-1,0,1);
    Vec3f v2(-20,0,1);
    Vec3f n0Vec(0,1,0);
    Vec3f n2 = v2/norm(v2);

    Vec3f n0 = n0Vec/norm(n0Vec);
    Vec3f n1 = n0.cross(n2);

    Vec3f pcMat = pe + v2;
    Vec3f pc(pcMat[0], pcMat[1], pcMat[2]);
    Vec3f th_P;

    Vec3i colorSum;
    double cameraX,cameraY; //the coordinate in the camera
    Vec3f npe;

    cout << n0 << endl << n1 << endl << n2 << endl;

    //checkIn
    if(sphere1.checkIn(pe)||sphere2.checkIn(pe)){
        cout << "fail in checking in " << endl;
        return 0;
    }

    for(int m = 0; m < M; m ++)
        for(int n = 0; n < N; n ++)
        {
            //with anti-aliasing
            randX = rand()%1000/1000.0;
            randY = rand()%1000/1000.0;
            colorSum = Vec3i(0,0,0);
            for(int xCnt = 0; xCnt < xSample; xCnt ++)
                for(int yCnt = 0; yCnt < ySample; yCnt ++)
                {
                    cameraX = (Sx*(m+(xCnt+randX)/xSample)/M-Sx/2.0);
                    cameraY = (Sy*(n+(yCnt+randY)/ySample)/N-Sy/2.0);
                    pos[0] = pc[0] + n1[0]*cameraY + n0[0]*cameraX;
                    pos[1] = pc[1] + n1[1]*cameraY + n0[1]*cameraX;
                    pos[2] = pc[2] + n1[2]*cameraY + n0[2]*cameraX;
                    npe = (pe - pos)/norm(pe - pos);

//                    if(fabs(npe.at<double>(0)+1)<1e-6){
//                        cout << "pos is " << pos << endl;
//                        cout << "npe is " << npe << endl;
//                        cout << "rayTracer is " << ellipsoid.rayTracer(pe,npe,n0,n1,n2) << endl;
//                        while(!getchar());
//                    }
                    //rayTracer of three objects
                    vector<double> timeRayVec;
                    timeRayVec.push_back(sphere1.rayTracer(pe,npe));
                    timeRayVec.push_back(sphere2.rayTracer(pe,npe));
                    timeRayVec.push_back(plane1.rayTracer(pe,npe));
                    timeRayVec.push_back(plane2.rayTracer(pe,npe));
                    timeRayVec.push_back(plane3.rayTracer(pe,npe));
                    objectIdx = traceTimeComp(timeRayVec);
                    if(objectIdx >= 0)
                        Ph = cameraPhGrab(pe, npe, timeRayVec[objectIdx]);

                    switch(objectIdx){
                        case -1 :
                            colorSum += Vec3b(50,50,50);
                            break;
                        case 0 :
                            colorSum += lighter.colorDeter(pe, Ph, sphere1.nH(Ph), sphere1.ColorH0, sphere1.ColorH1, sphere1.ColorH2, sphere1.ColorH3);
                            break;
                        case 1 :
                            colorSum += lighter.colorDeter(pe, Ph, sphere2.nH(Ph), sphere2.ColorH0, sphere2.ColorH1, sphere2.ColorH2, sphere2.ColorH3);
                            break;
                        case 2 :
                            colorSum += lighter.colorDeter(pe, Ph, plane1.nH(), plane1.ColorH0, plane1.ColorH1, plane1.ColorH2, plane1.ColorH3);
                            break;
                        case 3 :
                            colorSum += lighter.colorDeter(pe, Ph, plane2.nH(), plane2.ColorH0, plane2.ColorH1, plane2.ColorH2, plane2.ColorH3);
                            break;
                        case 4 :
                            colorSum += lighter.colorDeter(pe, Ph, plane3.nH(), plane3.ColorH0, plane3.ColorH1, plane3.ColorH2, plane3.ColorH3);
                            break;
                    }
                }
            img.at<Vec3b>(n,m) = Vec3b((colorSum)/(xSample*ySample));

        }

    flip(img, img, 1);
    imshow("lena", img);
    imwrite("img.png",img);
    resize(img,img,Size(160,120));
    imwrite("00.jpg",img);
    waitKey(0);
    return 0;
}



