#include "Quartics.h"

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
