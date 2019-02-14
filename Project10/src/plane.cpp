#include "plane.h"


Vec3d Planes::nH_Tex(Vec3d Ph){
    Vec3d nH;
    double xH,yH,zH;
    double x = (Ph - pc).dot(n0)/s0;
    double y = (Ph - pc).dot(n1)/s1;
    int xI = min(max(int(round((x-floor(x))*normalMap.cols)),0),normalMap.cols-1);
    int yI = min(max(int(round((y-floor(y))*normalMap.rows)),0),normalMap.rows-1);
//        cout << "xI is " << xI << " yI is " << yI << endl;
    xH = 2*normalMap.at<Vec3b>(yI, xI)[2]/255.0 - 1;
    yH = 2*normalMap.at<Vec3b>(yI, xI)[1]/255.0 - 1;
    zH = 2*normalMap.at<Vec3b>(yI, xI)[0]/255.0 - 1;

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

Planes::Planes(Vec3d n2_, Vec3d pc_, Vec3d n1_, int texture_, float nRR_, Mat textureMap_, Mat normalMap_, float s0_, float s1_, Vec4f ColorH3_){
    n1 = n1_;
    n2 = n2_;
    n0 = n1.cross(n2);
    pc = pc_;
    planeType = TEXTURE;
    ColorH3 = ColorH3_;
    textureMap = textureMap_;
    s0 = s0_;
    s1 = s1_;
    nRR = nRR_;
    normalMap = normalMap_;
    texture = texture_;
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
            texColors[2][i] *= 1;
        }
        texColors[0][3] = 255;
        texColors[2][3] = 255;
        return texColors;
    }
}
