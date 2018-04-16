#include "lightSource.h"

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
