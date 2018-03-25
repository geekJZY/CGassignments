#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

#define PI 3.1415927

using namespace std;
using namespace cv;


///Define Objects
class Quadrics{
private:
    Vec4f ColorH0, ColorH1, ColorH2;
    Mat textureMap;
public:
    enum SphereType{
        TEXTURE,
        SOLID_TEXTURE,
        SOLID_TEXTURE_PER,
        HOMO
    };
    SphereType sphereType;
    int a02,a12,a22,a21,a00;
    double s0, s1, s2;
    float x_solid, y_solid;
    Vec3d n0, n1, n2;
    Vec4f  ColorH3;
    Vec3d pc;
    Vec3d planeA, planeY, planeX, planeC, projPoint;
    bool checkIn(Vec3d pos);
    double rayTracer(Vec3d pe, Vec3d npe);
    Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3d pc_, Vec4f ColorH0_, Vec4f ColorH1_, Vec4f ColorH2_, Vec4f ColorH3_);
    Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3d pc_, Mat textureMap_, Vec4f ColorH3_);
    Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3d pc_, Mat textureMap_, Vec3d planeA_, Vec3d planeX_, Vec3d planeC_, Vec3d projPoint_, float x_solid, float y_solid, Vec4f ColorH3_, SphereType _sphereType = SOLID_TEXTURE);

    vector<Vec4f> TexColorDeter(Vec3d Ph);
    Vec3d nH(Vec3d Ph);
};

//Get the normal vector
Vec3d Quadrics::nH(Vec3d Ph){
    Vec3d vecH = 2*a02*n0.dot(Ph - pc)/(s0*s0)*n0 + 2*a12*n1.dot(Ph - pc)/(s1*s1)*n1 + 2*a22*n2.dot(Ph - pc)/(s2*s2)*n2 + a21/s2*n2;
    return vecH/norm(vecH);
}

Quadrics::Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3d pc_, Vec4f ColorH0_, Vec4f ColorH1_, Vec4f ColorH2_, Vec4f ColorH3_){
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
    n0 = Vec3d(1,0,0);
    n1 = Vec3d(0,1,0);
    n2 = Vec3d(0,0,1);
    sphereType = HOMO;
}

Quadrics::Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3d pc_, Mat textureMap_, Vec4f ColorH3_){
    a02 = a02_;
    a12 = a12_;
    a22 = a22_;
    a21 = a21_;
    a00 = a00_;
    s0 = s0_;
    s1 = s1_;
    s2 = s2_;
    pc = pc_;
    textureMap = textureMap_;
    ColorH3 = ColorH3_;
    n0 = Vec3d(1,0,0);
    n1 = Vec3d(0,1,0);
    n2 = Vec3d(0,0,1);
    sphereType = TEXTURE;
}

Quadrics::Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3d pc_, Mat textureMap_, Vec3d planeA_, Vec3d planeX_, Vec3d planeC_, Vec3d projPoint_, float x_solid_, float y_solid_, Vec4f ColorH3_, SphereType _sphereType){
    a02 = a02_;
    a12 = a12_;
    a22 = a22_;
    a21 = a21_;
    a00 = a00_;
    s0 = s0_;
    s1 = s1_;
    s2 = s2_;
    pc = pc_;
    textureMap = textureMap_;
    ColorH3 = ColorH3_;
    n0 = Vec3d(1,0,0);
    n1 = Vec3d(0,1,0);
    n2 = Vec3d(0,0,1);
    sphereType = _sphereType;
    planeA = planeA_;
    planeX = planeX_;
    planeY = planeA.cross(planeX);
    planeC = planeC_;
    x_solid = x_solid_;
    y_solid = y_solid_;
    projPoint = projPoint_;
}


bool Quadrics::checkIn(Vec3d pos){
    double F = pow((pos-pc).dot(n0)/s0,2)*a02+pow((pos-pc).dot(n1)/s1,2)*a12+pow((pos-pc).dot(n2)/s2,2)*a22+a00;
    if(F<0) return true;
    else return false;
}

double Quadrics::rayTracer(Vec3d pe, Vec3d npe){
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

vector<Vec4f> Quadrics::TexColorDeter(Vec3d Ph){
    vector<Vec4f> texColors(3);
    if(sphereType == HOMO){
        texColors[0] = ColorH0;
        texColors[1] = ColorH1;
        texColors[2] = ColorH2;
        return texColors;
    }
    else if(sphereType == TEXTURE){
        Vec3d npe;
        npe = (Ph-pc)/norm(Ph-pc);
        double x,y,z,theta,phi;
        int X,Y;
        x = npe.dot(n0);
        y = npe.dot(n1);
        z = npe.dot(n2);
        phi = acos(z);
        theta = acos(y/sin(phi));
        if(x<0) theta = 2*PI - theta;
        X = max(min(int(round((theta/(2*PI))*textureMap.cols)),textureMap.cols),0);
        Y = max(min(int(round((phi/(PI))*textureMap.rows)),textureMap.rows),0);

        texColors[1][0] = textureMap.at<Vec3b>(Y, X)[0];
        texColors[1][1] = textureMap.at<Vec3b>(Y, X)[1];
        texColors[1][2] = textureMap.at<Vec3b>(Y, X)[2];
        texColors[1][3] = 255;

        for(int i = 0; i < 3; i ++){
            texColors[0][i] = texColors[1][i]/5;
            texColors[2][i] = min(texColors[1][i], float(255.0));
            texColors[2][i] *= 10;
        }
        texColors[0][3] = 255;
        texColors[2][3] = 2550;
        return texColors;
    }
    else if(sphereType == SOLID_TEXTURE){
        Vec3d vpc;
        vpc = Ph-planeC;
        double x,y;
        int X,Y;
        x = vpc.dot(planeX);
        y = vpc.dot(planeY);
        X = x / x_solid * textureMap.cols;
        Y = y / y_solid * textureMap.rows;
        if(X < 0 || Y < 0 || X >= textureMap.cols || Y >= textureMap.rows){
            for(int i = 0; i < 3; i ++) texColors[i] = Vec4f(0,0,0,0);
            return texColors;
        }

        texColors[1][0] = textureMap.at<Vec3b>(Y, X)[0];
        texColors[1][1] = textureMap.at<Vec3b>(Y, X)[1];
        texColors[1][2] = textureMap.at<Vec3b>(Y, X)[2];
        texColors[1][3] = 255;

        for(int i = 0; i < 3; i ++){
            texColors[0][i] = texColors[1][i]/5;
            texColors[2][i] = min(texColors[1][i], float(255.0));
            texColors[2][i] *= 10;
        }
        texColors[0][3] = 255;
        texColors[2][3] = 2550;
        return texColors;
    }
    else if(sphereType == SOLID_TEXTURE_PER){
        Vec3d vph, pI;
        vph = Ph-projPoint;
        double t,x,y;
        int X,Y;
        t = planeA.dot(planeC - projPoint)/(planeA.dot(vph));
        if(t < 0 || t > 1){
            for(int i = 0; i < 3; i ++) texColors[i] = Vec4f(0,0,0,0);
            return texColors;
        }
        pI = projPoint + t * vph - planeC;
        x = pI.dot(planeX);
        y = pI.dot(planeY);
        X = x / x_solid * textureMap.cols;
        Y = y / y_solid * textureMap.rows;
        if(X < 0 || Y < 0 || X >= textureMap.cols || Y >= textureMap.rows){
            for(int i = 0; i < 3; i ++) texColors[i] = Vec4f(0,0,0,0);
            return texColors;
        }

        texColors[1][0] = textureMap.at<Vec3b>(Y, X)[0];
        texColors[1][1] = textureMap.at<Vec3b>(Y, X)[1];
        texColors[1][2] = textureMap.at<Vec3b>(Y, X)[2];
        texColors[1][3] = 255;

        for(int i = 0; i < 3; i ++){
            texColors[0][i] = texColors[1][i]/5;
            texColors[2][i] = min(texColors[1][i], float(255.0));
            texColors[2][i] *= 10;
        }
        texColors[0][3] = 255;
        texColors[2][3] = 2550;
        return texColors;
    }
}

class Planes{
private:
    Vec4f ColorH0, ColorH1, ColorH2;
    enum PlaneType{
        TEXTURE,
        HOMO
    };
    PlaneType planeType;
    Mat textureMap;
    Mat normalMap;
    float s0, s1; //params for texture map
public:
    Vec3d pc;
    Vec3d n1, n2, n0;
    Vec4f ColorH3;
    Planes(Vec3d n2_, Vec3d pc_, Vec4f ColorH0_, Vec4f ColorH1_, Vec4f ColorH2_, Vec4f ColorH3_);
    Planes(Vec3d n2_, Vec3d pc_, Vec3d n1_, Mat textureMap_, Mat normalMap_, float s0_, float s1_, Vec4f ColorH3_);
    vector<Vec4f> TexColorDeter(Vec3d ph);
    bool checkIn(Vec3d pos);
    Vec3d nH(void);
    Vec3d nH_Tex(Vec3d Ph);
    double rayTracer(Vec3d pe, Vec3d npe);
};

Vec3d Planes::nH_Tex(Vec3d Ph){
    Vec3d nH;
    double xH,yH,zH;
    double x = (Ph - pc).dot(n0)/s0;
    double y = (Ph - pc).dot(n1)/s1;
    int xI = min(max(int(round((x-floor(x))*normalMap.cols)),0),normalMap.cols);
    int yI = min(max(int(round((y-floor(y))*normalMap.rows)),0),normalMap.cols);
//        cout << "xI is " << xI << " yI is " << yI << endl;
    xH = normalMap.at<Vec3b>(yI, xI)[0]/255.0 - 0.5;
    yH = normalMap.at<Vec3b>(xI, yI)[1]/255.0 - 0.5;
    zH = normalMap.at<Vec3b>(xI, yI)[2]/255.0 - 0.5;

    nH = xH*n0 + yH*n1 + zH*n2;
    return nH;
}


Vec3d Planes::nH(void){
    return n2;
}

Planes::Planes(Vec3d n2_, Vec3d pc_, Vec4f ColorH0_, Vec4f ColorH1_, Vec4f ColorH2_, Vec4f ColorH3_){
    n2 = n2_;
    pc = pc_;
    ColorH0 = ColorH0_;
    ColorH1 = ColorH1_;
    ColorH2 = ColorH2_;
    ColorH3 = ColorH3_;
    planeType = HOMO;
}

Planes::Planes(Vec3d n2_, Vec3d pc_, Vec3d n1_, Mat textureMap_, Mat normalMap_, float s0_, float s1_, Vec4f ColorH3_){
    n1 = n1_;
    n2 = n2_;
    n0 = n1.cross(n2);
    pc = pc_;
    planeType = TEXTURE;
    ColorH3 = ColorH3_;
    textureMap = textureMap_;
    s0 = s0_;
    s1 = s1_;
    normalMap = normalMap_;
}

bool Planes::checkIn(Vec3d pos){
    double F = n2.dot(pos - pc);
    if(F<0) return true;
    else return false;
}

double Planes::rayTracer(Vec3d pe, Vec3d npe){
    return n2.dot(pc-pe)/n2.dot(npe);
}

vector<Vec4f> Planes::TexColorDeter(Vec3d ph = Vec3d(0,0,0)){
    vector<Vec4f> texColors(3);
    if(planeType == HOMO){
        texColors[0] = (ColorH0);
        texColors[1] = (ColorH1);
        texColors[2] = (ColorH2);
        return texColors;
    }
    else if(planeType == TEXTURE){
        double x = (ph - pc).dot(n0)/s0;
        double y = (ph - pc).dot(n1)/s1;
        int xI = round((x-floor(x))*textureMap.cols);
        int yI = round((y-floor(y))*textureMap.rows);
//        cout << "xI is " << xI << " yI is " << yI << endl;
        texColors[1][0] = textureMap.at<Vec3b>(yI, xI)[0];
        texColors[1][1] = textureMap.at<Vec3b>(yI, xI)[1];
        texColors[1][2] = textureMap.at<Vec3b>(yI, xI)[2];
        texColors[1][3] = 255;

        for(int i = 0; i < 3; i ++){
            texColors[0][i] = texColors[1][i]/5;
            texColors[2][i] = min(texColors[1][i], float(255.0));
            texColors[2][i] *= 10;
        }
        texColors[0][3] = 255;
        texColors[2][3] = 2550;
        return texColors;
    }
}

///Define Light Source
class lightSource{
public:
    enum Mode
    {
      SPOT,
      DIRECTIONAL,
      POINTLIGHT
    };

    Vec3d Pl, nLH_D; //light point OR light Direction
    double spotRange;
    double ks, kd;
    Vec4f ColorL;   //color of light
    //params for color calculation
    Mode lightMode;
    lightSource(Vec3d Pl, Vec3d nLH, double spotRange, double ks, double kd, Vec4f ColorL, Mode lightMode);
    //H0 ambient, H1 diffuse, H2 specular
    Vec4f colorSDcal(Vec3d Pe, Vec3d Ph, Vec3d nH, Vec4f ColorH0, Vec4f ColorH1, Vec4f ColorH2, Vec4f ColorH3);
};

lightSource::lightSource(Vec3d Pl_, Vec3d nLH_, double spotRange_, double ks_, double kd_, Vec4f ColorL_, Mode lightMode_){
    Pl = Pl_;
    ColorL = ColorL_;
    nLH_D = nLH_;
    ks = ks_;
    kd = kd_;
    spotRange = spotRange_;
    lightMode = lightMode_;
}

Vec4f lightSource::colorSDcal(Vec3d Pe, Vec3d Ph, Vec3d nH, Vec4f ColorH0, Vec4f ColorH1, Vec4f ColorH2, Vec4f ColorH3){
    double s, t; //the coefficient for specular, diffuse light and hale light.
    double cosRE; //cos(angle between reflect light and eye)
    Vec3d rHE, nLH, nHE;
    Vec4f ColorPE;
    if((lightMode == SPOT)||(lightMode == POINTLIGHT)) nLH = (Ph - Pl)/norm(Pl - Ph);
    else if(lightMode == DIRECTIONAL) nLH = nLH_D;
    nHE = (Pe - Ph)/norm(Pe - Ph);
    rHE = nLH - 2*(nLH.dot(nH))*nH;
    cosRE = rHE.dot(nHE);
    s = std::max((pow(cosRE,3)-0.9)/0.1,0.0);
    t = std::max(cosRE,0.0);
    if((lightMode == POINTLIGHT)||((lightMode == SPOT)&&(nLH.dot(nLH_D) > spotRange))) return (ks*s*ColorH2.mul(ColorL) + kd*t*ColorH1.mul(ColorL)/(1+norm(Pl-Ph)));
    else if(lightMode == DIRECTIONAL) return (ks*s*ColorH2.mul(ColorL) + kd*t*ColorH1.mul(ColorL)/100);
    else return Vec4f(0,0,0,0);
}


Vec3d cameraPhGrab(Vec3d pe, Vec3d npe, double t){
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


class colorShadowDeter{
public:
    enum ShadowMode{
        HARD_SHADOW,
        SOFT_SHADOW,
        DOWN_SHADOW
    };

    ShadowMode shadowMode;
    int objectIdx;
    Vec3d Ph;
    colorShadowDeter(ShadowMode shadowMode);
    void PhCal(Vec3d pe, Vec3d npe, vector<Planes> vecPlane, vector<Quadrics> vecQuad);
    Vec3b ColorDeter(Vec3d pe, double kb, double k0, vector<Planes> vecPlane, vector<Quadrics> vecQuad, vector<lightSource> lighters, vector<float> shadowFlag);
    vector<float> shadowRayTracing(vector<Planes> planes, vector<Quadrics> spheres, vector<lightSource> lighters);
};

colorShadowDeter::colorShadowDeter(ShadowMode shadowMode_){
    shadowMode = shadowMode_;
}

float shadowRayTrace(Quadrics quad, Vec3d pl, Vec3d Ph, Vec3d nLH_D, lightSource::Mode lightType, colorShadowDeter::ShadowMode shadowMode = colorShadowDeter::SOFT_SHADOW){
    double lenLH;
    Vec3d nlh;
    if((lightType == lightSource::POINTLIGHT)||(lightType == lightSource::SPOT)){
        lenLH = norm(Ph - pl);
        nlh = (Ph - pl)/lenLH;
    }
    else{
         nlh = -nLH_D;
         pl = Ph;
    }
    if(shadowMode == colorShadowDeter::DOWN_SHADOW){
        Ph -= quad.nH(Ph)*0.5;  //down for 0.5
    }

    double A = quad.a02*pow(quad.n0.dot(nlh)/quad.s0,2)+quad.a12*pow(quad.n1.dot(nlh)/quad.s1,2)+quad.a22*pow(quad.n2.dot(nlh)/quad.s2,2);
    double B = quad.a02*(2*quad.n0.dot(nlh))*(quad.n0.dot(pl-quad.pc))/(quad.s0*quad.s0)+
               quad.a12*(2*quad.n1.dot(nlh))*(quad.n1.dot(pl-quad.pc))/(quad.s1*quad.s1)+
               quad.a22*(2*quad.n2.dot(nlh))*(quad.n2.dot(pl-quad.pc))/(quad.s2*quad.s2)+
               quad.a21*quad.n2.dot(nlh)/quad.s2;
    double C = quad.a02*pow(quad.n0.dot(pl-quad.pc)/quad.s0,2) +
               quad.a12*pow(quad.n1.dot(pl-quad.pc)/quad.s1,2) +
               quad.a22*pow(quad.n2.dot(pl-quad.pc)/quad.s2,2) +
               quad.a21*quad.n2.dot(pl-quad.pc)/quad.s2 +
               quad.a00;
//    cout << "A is " << A <<" B is " << B <<" C is " << C << endl;
    double delta = pow(B,2)-4*A*C;
    if(delta < 0) return 0;
    delta = sqrt(delta);
    double tSmall = ((-B-delta)/(2*A));
    double tBig = ((-B+delta)/(2*A));
    if(shadowMode != colorShadowDeter::DOWN_SHADOW){
        if((lightType == lightSource::POINTLIGHT)||(lightType == lightSource::SPOT)){
            if((tSmall > 0)&&(tBig - 1e-4 < lenLH)) return (tBig - tSmall);
            else return 0;
        }
        else{
            if((tSmall + 1e-4 > 0)&&(tBig + 1e-4 > 0)) return (tBig - tSmall);
            else return 0;
        }
    }
    else{
        if((lightType == lightSource::POINTLIGHT)||(lightType == lightSource::SPOT)){
            if((tSmall > 0)&&(tBig < lenLH)) return (tBig - tSmall);
            else if((tSmall > 0)&&(tSmall < lenLH)&&(tBig >= lenLH)) return (lenLH - tSmall);
            else return 0;
        }
        else{
            if((tSmall + 1e-4 > 0)&&(tBig + 1e-4 > 0)) return (tBig - tSmall);
            if((tSmall + 1e-4 < 0)&&(tBig + 1e-4 > 0)) return (tBig);
            else return 0;
        }
    }
}

float shadowRayTrace(Planes plane, Vec3d pl, Vec3d Ph, Vec3d nLH_D, lightSource::Mode lightType, colorShadowDeter::ShadowMode shadowMode = colorShadowDeter::SOFT_SHADOW){
    double lenLH;
    Vec3d nlh;
    if((lightType == lightSource::POINTLIGHT)||(lightType == lightSource::SPOT)){
        lenLH = norm(Ph - pl);
        nlh = (Ph - pl)/lenLH;
    }
    else{
         nlh = -nLH_D;
         pl = Ph;
    }
    Ph -= plane.nH()*0.5;  //down for 0.5

    //plane ray trace
    double t = plane.n2.dot(plane.pc-pl)/plane.n2.dot(nlh);
    if((lightType == lightSource::POINTLIGHT)||(lightType == lightSource::SPOT)){
        if(t < lenLH) return (lenLH - t);
    }
    else{
        if(t > 0) return t;
    }
    return 0;
}

void colorShadowDeter::PhCal(Vec3d pe, Vec3d npe, vector<Planes> vecPlane, vector<Quadrics> vecQuad){
    vector<double> timeRayVec;

    for(int i = 0; i < vecQuad.size(); i ++)   timeRayVec.push_back(vecQuad[i].rayTracer(pe,npe));
    for(int i = 0; i < vecPlane.size(); i ++)  timeRayVec.push_back(vecPlane[i].rayTracer(pe,npe));
    objectIdx = traceTimeComp(timeRayVec);
    if(objectIdx >= 0) Ph = cameraPhGrab(pe, npe, timeRayVec[objectIdx]);
}

Vec3b colorShadowDeter::ColorDeter(Vec3d pe, double kb, double k0, vector<Planes> vecPlane, vector<Quadrics> vecQuad, vector<lightSource> lighters, vector<float> shadowFlag){
    Vec4f colorSum4f(0,0,0,0);

    double b; //the coefficient for hale light.
    double cosEH; //cos(angle between reflect light and eye)
    Vec3d rHE, nLH, nHE, nH;

    if(objectIdx < 0) return Vec3b(50,50,50);
    nHE = (pe - Ph)/norm(pe - Ph);
    nH = objectIdx < vecQuad.size()?vecQuad[objectIdx].nH(Ph):vecPlane[objectIdx - vecQuad.size()].nH();
    cosEH = nHE.dot(nH);
    b = std::max(((1-cosEH)-0.92)/0.08,0.0);

    if(objectIdx < vecQuad.size()){
        vector<Vec4f> textureClr = vecQuad[objectIdx].TexColorDeter(Ph);
        for(int i = 0; i < lighters.size(); i ++){
            if((shadowMode == HARD_SHADOW)&&(fabs(shadowFlag[i]) < 1e-2)) colorSum4f += lighters[i].colorSDcal(pe, Ph, nH, textureClr[0], textureClr[1], textureClr[2], vecQuad[objectIdx].ColorH3);
            else if(shadowMode == SOFT_SHADOW) colorSum4f += lighters[i].colorSDcal(pe, Ph, nH, textureClr[0], textureClr[1], textureClr[2], vecQuad[objectIdx].ColorH3)/(shadowFlag[i]+1);
            else if(shadowMode == DOWN_SHADOW) colorSum4f += lighters[i].colorSDcal(pe, Ph, nH, textureClr[0], textureClr[1], textureClr[2], vecQuad[objectIdx].ColorH3)/(std::max(double(shadowFlag[i]/0.5), 1.0));
            else if(shadowMode == DOWN_SHADOW) colorSum4f += lighters[i].colorSDcal(pe, Ph, nH, textureClr[0], textureClr[1], textureClr[2], vecQuad[objectIdx].ColorH3)/(std::max(double(shadowFlag[i]/0.5), 1.0));
        }
        colorSum4f += k0*textureClr[0] + b*kb*vecQuad[objectIdx].ColorH3;
    }
    else{
        vector<Vec4f> textureClr = vecPlane[objectIdx - vecQuad.size()].TexColorDeter(Ph);
//        for(int i = 0; i < textureClr.size(); i ++) cout << textureClr[i];
//        cout << endl;
        for(int i = 0; i < lighters.size(); i ++){
           if((shadowMode == HARD_SHADOW)&&(fabs(shadowFlag[i]) < 1e-2)) colorSum4f += lighters[i].colorSDcal(pe, Ph, vecPlane[objectIdx - vecQuad.size()].nH(), textureClr[0],
                                                textureClr[1], textureClr[2], vecPlane[objectIdx - vecQuad.size()].ColorH3);
           else if(shadowMode == SOFT_SHADOW) colorSum4f += lighters[i].colorSDcal(pe, Ph, vecPlane[objectIdx - vecQuad.size()].nH(), textureClr[0],
                                                textureClr[1], textureClr[2], vecPlane[objectIdx - vecQuad.size()].ColorH3)/(shadowFlag[i]+1);
           else if(shadowMode == DOWN_SHADOW) colorSum4f += lighters[i].colorSDcal(pe, Ph, vecPlane[objectIdx - vecQuad.size()].nH(), textureClr[0],
                                                textureClr[1], textureClr[2], vecPlane[objectIdx - vecQuad.size()].ColorH3)/( std::max(double(shadowFlag[i]/0.5), 1.0) );
        }
        colorSum4f += k0*textureClr[0] + b*kb*vecPlane[objectIdx - vecQuad.size()].ColorH3;
    }
    //output and test
//    if(norm(Ph-Vec3d(11,-9,1))<0.1){
//        cout << "Ph is " << Ph << endl;
//        cout << "nLH is " << nLH << " rHE is " << rHE << endl;
//        cout << "nHE is " << nHE << " rHE.dot(nHE) is " << rHE.dot(nHE) << endl;
//        cout << "ColorPE (" << int(ColorPE[0]*255/ColorPE[3]) << ", " << int(ColorPE[1]*255/ColorPE[3])<<", " << int(ColorPE[2]*255/ColorPE[3]) <<")" << endl;
//        cout << "t is " << t << "  s is " << s << " b is " << b << endl;
//    }

    return Vec3b(char(colorSum4f[0]*255/colorSum4f[3]),char(colorSum4f[1]*255/colorSum4f[3]),char(colorSum4f[2]*255/colorSum4f[3]));
}

vector<float> colorShadowDeter::shadowRayTracing(vector<Planes> planes,vector<Quadrics> spheres, vector<lightSource> lighters){
    vector<float> tracingTime(spheres.size());

    for(int i = 0; i < lighters.size(); i ++){
        tracingTime[i] = 0;
        if(shadowMode != DOWN_SHADOW){
            for(int j = 0; j < spheres.size(); j ++){
                tracingTime[i] += shadowRayTrace(spheres[j], lighters[i].Pl, Ph, lighters[i].nLH_D, lighters[i].lightMode);
            }
        }
        else{
            for(int j = 0; j < spheres.size(); j ++){
                tracingTime[i] += shadowRayTrace(spheres[j], lighters[i].Pl, Ph, lighters[i].nLH_D, lighters[i].lightMode, DOWN_SHADOW);
            }
            for(int j = 0; j < planes.size(); j ++){
                tracingTime[i] += shadowRayTrace(planes[j], lighters[i].Pl, Ph, lighters[i].nLH_D, lighters[i].lightMode, DOWN_SHADOW);
            }
        }
    }

    return tracingTime;
}

int main(int argc, char *argv[]){

    vector<Scalar> color{Scalar(255,224,147),Scalar(93,66,255),Scalar(255,255,255),Scalar(230,180,80)};
    Mat img(Size(640,480),CV_8UC3,color[0]);
    Mat texturePlane = imread("Texture.jpg");
    Mat sphereTexture = imread("ball_image_small.jpg");
    Mat sphereTexture2 = imread("ball_image_10_small.jpg");
    Mat normalTexture = imread("normal.jpg");
    double Sx = 16, Sy = 12;
    Vec3d pos;
    Vec3d Ph;
    int objectIdx;
    int xSample = 2, ySample = 2;
    double randX,randY;
    int M = img.cols, N = img.rows;

    //the efficient of background and outline
    double kb = 0.5, k0 = 0.5;
    colorShadowDeter CSDeter(colorShadowDeter::SOFT_SHADOW);
    vector<float> shadowFlag;

    //Graphics
    vector<Quadrics> spheres;
    vector<Planes> planes;
    vector<lightSource> lighters;

    spheres.push_back(Quadrics(1,1,1,0,-1,5,5,5,Vec3d(0,0,5),sphereTexture,Vec3d(0,-1,0),Vec3d(-1,0,0),Vec3d(3,-6,8),Vec3d(0,-12,5),6,6,Vec4f(0,0,0,255),Quadrics::SOLID_TEXTURE_PER));
    spheres.push_back(Quadrics(1,1,1,0,-1,3,3,3,Vec3d(0,12,4),sphereTexture2,Vec3d(0,-1,0),Vec3d(-1,0,0),Vec3d(2,8,6),Vec3d(0,4,4),4,4,Vec4f(0,0,0,255),Quadrics::SOLID_TEXTURE_PER));
    planes.push_back(Planes(Vec3d(0,0,1),Vec3d(0,0,0),Vec3d(0,1,0),texturePlane,normalTexture,10,10,Vec4f(10,10,10,255)));
//    planes.push_back(Planes(Vec3d(-1,0,0),Vec3d(12,0,0),Vec4f(50,100,50,255),Vec4f(170,220,170,255),Vec4f(160,210,160,255),Vec4f(10,10,10,255)));
//    planes.push_back(Planes(Vec3d(0,1,0),Vec3d(0,-9,0),Vec4f(50,50,100,255),Vec4f(170,170,220,255),Vec4f(1800,1800,2300,2550),Vec4f(10,10,10,255)));

    //area light
    lighters.push_back(lightSource(Vec3d(0,40,50), Vec3d(0.5773,-0.5773,-0.5773), 0.9,0.5,0.5,Vec4f(60,60,60,60), lightSource::POINTLIGHT));

    //camera
    Vec3d pe(40,25,6);
    Vec3d vUp(-1,0,1);
    Vec3d v2(20,10,1);
    Vec3d n0Vec(0,-1,0);
    Vec3d n2 = v2/norm(v2);

    Vec3d n0 = n0Vec/norm(n0Vec);
    Vec3d n1 = n0.cross(n2);

    Vec3d pcMat = pe + v2;
    Vec3d pc(pcMat[0], pcMat[1], pcMat[2]);
    Vec3d th_P;

    Vec3i colorSum;
    double cameraX,cameraY; //the coordinate in the camera
    Vec3d npe;

    cout << n0 << endl << n1 << endl << n2 << endl;

    //checkIn
    for(int i = 0; i < spheres.size(); i ++){
        if(spheres[i].checkIn(pe)){
            cout << "fail in checking in " << endl;
            return 0;
        }
    }
    for(int j = 0; j < planes.size(); j ++){
        if(planes[j].checkIn(pe)){
            cout << "fail in checking in " << endl;
            return 0;
        }
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

                    CSDeter.PhCal(pe, npe, planes, spheres);
                    if((xCnt == 0)&&(yCnt == 0)){
                        shadowFlag = CSDeter.shadowRayTracing(planes, spheres, lighters);
                    }
                    colorSum += CSDeter.ColorDeter(pe, kb, k0, planes, spheres, lighters, shadowFlag);
                }
            img.at<Vec3b>(n,m) = Vec3b((colorSum)/(xSample*ySample));

        }

//    flip(img, img, 1);
    imshow("lena", img);
    imwrite("img.png",img);
    resize(img,img,Size(160,120));
    imwrite("00.jpg",img);
    waitKey(0);
    return 0;
}

