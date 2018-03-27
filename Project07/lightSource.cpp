#include "lightSource.h"

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
