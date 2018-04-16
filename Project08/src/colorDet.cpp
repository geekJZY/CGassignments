#include "colorDet.h"

//The index of object is starting from 0
///Planes -> Spheres

int colorDet::traceTimeComp(vector<double> traceVecTime){
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

//initialize
colorDet::colorDet(vector<Planes> _planes, vector<Quadrics> _spheres, vector<Cube> _cubes, vector<lightSource> _lighters, Vec3d _Peye, Vec4f ColorH3_, double kb_, double k0_):
    planes(_planes),
    spheres(_spheres),
    cubes(_cubes),
    lighters(_lighters),
    Peye(_Peye),
    ColorH3(ColorH3_),
    kb(kb_),
    k0(k0_)
{
}

//give the incoming light and point, return the object hit and hit point
pair<int, Vec3d> colorDet::spaceTracer(Vec3d nep, Vec3d pe){
    vector<double> traceVecTime;
    int objIdx;
    for(int i = 0; i < planes.size(); i ++){
        traceVecTime.push_back(planes[i].rayTracer(pe, nep));
    }
    for(int i = 0; i < spheres.size(); i ++){
        traceVecTime.push_back(spheres[i].rayTracer(pe, nep));
    }
    for(int i = 0; i < cubes.size(); i ++){
        traceVecTime.push_back(cubes[i].rayTrace(pe, nep));
    }
    objIdx = traceTimeComp(traceVecTime); //rayTrace
    if (objIdx >= 0) return make_pair(objIdx, pe + nep * traceVecTime[objIdx]);
    else return make_pair(-1, -1);
}

void colorDet::projectColorInit(Mat projectImg_, Vec3d planeA_, Vec3d planeX_, Vec3d planeC_, float x_solid_, float y_solid_, Vec3d projP_){
    this -> projectImg = projectImg_;
    this -> planeA = planeA_;
    this -> planeX = planeX_;
    this -> planeC = planeC_;
    this -> x_solid = x_solid_;
    this -> y_solid = y_solid_;
    this -> projP = projP_;
    return;
}

///project image
Vec4f colorDet::projectColorDet(Vec3d ph){
    pair<int, Vec3d> rayTracing;

    //perspective projection
    Vec3d vph, pI;
    vph = ph-projP;
    Vec3b nhp = -vph/norm(vph);
    double t,x,y, t_other;
    int X,Y;
    t = planeA.dot(planeC - projP)/(planeA.dot(vph));
    if(t < 0 || t > 1){
        return Vec4f(0,0,0,0);
    }
    pI = projP + t * vph - planeC;
    x = pI.dot(planeX);
    Vec3b planeY = planeA.cross(planeX);
    y = pI.dot(planeY);
    X = x / x_solid * projectImg.cols;
    Y = y / y_solid * projectImg.rows;
    if(X < 0 || Y < 0 || X >= projectImg.cols || Y >= projectImg.rows){
        return Vec4f(0,0,0,0);
    }

    rayTracing = this->spaceTracer(nhp, ph);
    if(norm(rayTracing.second - ph) - t < -1e-3) return Vec4f(0,0,0,0);

    Vec4f projColor;
    projColor[0] = projectImg.at<Vec3b>(Y, X)[0];
    projColor[1] = projectImg.at<Vec3b>(Y, X)[1];
    projColor[2] = projectImg.at<Vec3b>(Y, X)[2];
    projColor[3] = 255;

    return projColor;
}

///I need add an infinite sphere to make rayTracing can always hit object
//give the object and point, return the color
Vec4f colorDet::colorReturn(int objFrom, int objInx, Vec3d ph, Vec3d nep, int tierCnt){
    //dividing by type of object
    if(tierCnt > 5) return Vec4f(100,100,100,200);
    if(objInx >= 0 && objInx < planes.size()){
        //plane
        if(planes[objInx].texture == TX_TOTAL_RL){
            Vec3d nr, nHE;
            pair<int, Vec3d> rayTracing;
            Vec3d normVec = planes[objInx].nH_Tex(ph);
            nr = nep - 2 * normVec.dot(nep) * normVec;
            rayTracing = this -> spaceTracer(nr, ph);
            //If rayTracing fail to hit object

            nHE = (Peye - ph)/norm(Peye - ph);
            double cosEH = nHE.dot(normVec);
            double b = std::max(((1-cosEH)-0.92)/0.08,0.0);
            return b*this -> kb*this -> ColorH3 + 0.8* this -> colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1);
        }
        else if(planes[objInx].texture == TX_HALF_RL_HALF_RR){
            Vec3d nr, nHE, vRR; // nr is reflection, vRR is refraction
            double nRR, delta, C, aRR, bRR;
            pair<int, Vec3d> rayTracing, rayRR;
            Vec3d normVec = planes[objInx].nH_Tex(ph);
            nr = nep - 2 * normVec.dot(nep) * normVec;
            rayTracing = this -> spaceTracer(nr, ph);
            nHE = (Peye - ph)/norm(Peye - ph);
            double cosEH = nHE.dot(normVec);
            double b = std::max(((1-cosEH)-0.92)/0.08,0.0);
            if(objFrom == objInx) nRR = 1/planes[objInx].nRR;
            else nRR = planes[objInx].nRR;
            //total reflection
            C = -nep.dot(normVec);
            delta = (C*C-1)/(nRR*nRR) + 1;
            if(delta < 0) //delta < 0 is reflection
                return b*this -> kb*this -> ColorH3 + 0.8* this -> colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1);
            aRR = -1/nRR;
            bRR = C/nRR - sqrt(delta);
            vRR = aRR * (-nep) + bRR * normVec;
            rayRR = this -> spaceTracer(vRR, ph);
            return b*this -> kb*this -> ColorH3 + 0.4* this -> colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1) + 0.4* this -> colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1);
        }
        else if(planes[objInx].texture == TX_OPAQUE){
            Vec4f colorSum4f(0,0,0,0);

            double b; //the coefficient for hale light.
            double cosEH; //cos(angle between reflect light and eye)
            Vec3d rHE, nLH, nHE, normVec;

            nHE = (Peye - ph)/norm(Peye - ph);
            normVec = planes[objInx].nH();
            cosEH = nHE.dot(normVec);
            b = std::max(((1-cosEH)-0.92)/0.08,0.0);

            vector<Vec4f> textureClr = planes[objInx].TexColorDeter(ph);

            for(int i = 0; i < lighters.size(); i ++){
                //perform rayTracing to the light and test if the hit point is the same as this point, if not, calculate the color and add it
                double dis_Hit2Obj;
                pair<int, Vec3d> lightTracing = this->spaceTracer((ph - lighters[i].Pl)/norm(ph - lighters[i].Pl), lighters[i].Pl);
                if(lightTracing.first == objInx) colorSum4f += lighters[i].colorSDcal(Peye, ph, normVec, textureClr[0],
                                                    textureClr[1], textureClr[2], this -> ColorH3);
            }
            colorSum4f += this -> k0*textureClr[0] + b*kb*this -> ColorH3;
//            cout << colorSum4f << endl;

            return colorSum4f;
        }
        else if(planes[objInx].texture == TX_HALF_OPAQUE_HALF_RL){

            Vec4f colorSum4f(0,0,0,0);

            double b; //the coefficient for hale light.
            double cosEH; //cos(angle between reflect light and eye)
            Vec3d rHE, nLH, nHE, normVec, nr;

            nHE = (Peye - ph)/norm(Peye - ph);
            normVec = planes[objInx].nH_Tex(ph);
//            cout << normVec << endl;
            cosEH = nHE.dot(normVec);
            b = std::max(((1-cosEH)-0.92)/0.08,0.0);

            vector<Vec4f> textureClr = planes[objInx].TexColorDeter(ph);

            for(int i = 0; i < lighters.size(); i ++){
                //perform rayTracing to the light and test if the hit point is the same as this point, if not, calculate the color and add it
                double dis_Hit2Obj;
                pair<int, Vec3d> lightTracing = this->spaceTracer((ph - lighters[i].Pl)/norm(ph - lighters[i].Pl), lighters[i].Pl);
                if(lightTracing.first == objInx) colorSum4f += lighters[i].colorSDcal(Peye, ph, normVec, textureClr[0],
                                                    textureClr[1], textureClr[2], this -> ColorH3);
            }
            colorSum4f += this -> k0*textureClr[0] + b*kb*this -> ColorH3;

            pair<int, Vec3d> rayTracing;
            nr = nep - 2 * normVec.dot(nep) * normVec;
            rayTracing = this -> spaceTracer(nr, ph);

            return 0.6*colorSum4f + 0.8* this ->colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1);
        }
    }
    else if(objInx >= planes.size() && objInx < planes.size() + spheres.size()){
        //sphere
        int sphereInx = objInx - planes.size();
        if(spheres[sphereInx].texture == TX_TOTAL_RL){
            Vec3d nr, nHE;
            pair<int, Vec3d> rayTracing;
            Vec3d normVec = spheres[sphereInx].nH(ph);
            nr = nep - 2 * normVec.dot(nep) * normVec;
            rayTracing = this -> spaceTracer(nr, ph);
            //If rayTracing fail to hit object

            nHE = (Peye - ph)/norm(Peye - ph);
            double cosEH = nHE.dot(normVec);
            double b = std::max(((1-cosEH)-0.92)/0.08,0.0);
            return b*this -> kb*this -> ColorH3 + 0.8* this -> colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1);
        }
        else if(spheres[sphereInx].texture == TX_HALF_RL_HALF_RR){
            Vec3d nr, nHE, vRR; // nr is reflection, vRR is refraction
            double nRR, delta, C, aRR, bRR;
            pair<int, Vec3d> rayTracing, rayRR;
            Vec3d normVec = spheres[sphereInx].nH(ph);
            nr = nep - 2 * normVec.dot(nep) * normVec;
            rayTracing = this -> spaceTracer(nr, ph);
            nHE = (Peye - ph)/norm(Peye - ph);
            double cosEH = nHE.dot(normVec);
            double b = std::max(((1-cosEH)-0.92)/0.08,0.0);
            if(objFrom == objInx) nRR = 1/spheres[sphereInx].nRR;
            else nRR = spheres[sphereInx].nRR;
            //total reflection
            C = -nep.dot(normVec);
            delta = (C*C-1)/(nRR*nRR) + 1;
            if(delta < 0) //delta < 0 is reflection
                return b*this -> kb*this -> ColorH3 + 0.8* this -> colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1);
            aRR = -1/nRR;
            bRR = C/nRR - sqrt(delta);
            vRR = aRR * (-nep) + bRR * normVec;
            rayRR = this -> spaceTracer(vRR, ph);
            return b*this -> kb*this -> ColorH3 + 0.4* this -> colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1) + 0.4* this -> colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1);
        }
        else if(spheres[sphereInx].texture == TX_OPAQUE){
            Vec4f colorSum4f(0,0,0,0);

            double b; //the coefficient for hale light.
            double cosEH; //cos(angle between reflect light and eye)
            Vec3d rHE, nLH, nHE, normVec;

            nHE = (Peye - ph)/norm(Peye - ph);
            normVec = spheres[sphereInx].nH(ph);
            cosEH = nHE.dot(normVec);
            b = std::max(((1-cosEH)-0.92)/0.08,0.0);

            vector<Vec4f> textureClr = spheres[sphereInx].TexColorDeter(ph);

            for(int i = 0; i < lighters.size(); i ++){
                //perform rayTracing to the light and test if the hit point is the same as this point, if not, calculate the color and add it
                double dis_Hit2Obj;
                pair<int, Vec3d> lightTracing = this->spaceTracer((ph - lighters[i].Pl)/norm(ph - lighters[i].Pl), lighters[i].Pl);
                if(lightTracing.first == objInx) colorSum4f += lighters[i].colorSDcal(Peye, ph, normVec, textureClr[0],
                                                    textureClr[1], textureClr[2], this -> ColorH3);
            }
            colorSum4f += this -> k0*textureClr[0] + b*kb*this -> ColorH3;

            return colorSum4f;
        }
        else if(spheres[sphereInx].texture == TX_HALF_OPAQUE_HALF_RL){

            Vec4f colorSum4f(0,0,0,0);

            double b; //the coefficient for hale light.
            double cosEH; //cos(angle between reflect light and eye)
            Vec3d rHE, nLH, nHE, normVec, nr;

            nHE = (Peye - ph)/norm(Peye - ph);
            normVec = spheres[sphereInx].nH(ph);
            cosEH = nHE.dot(normVec);
            b = std::max(((1-cosEH)-0.92)/0.08,0.0);

            vector<Vec4f> textureClr = spheres[sphereInx].TexColorDeter(ph);

            for(int i = 0; i < lighters.size(); i ++){
                //perform rayTracing to the light and test if the hit point is the same as this point, if not, calculate the color and add it
                double dis_Hit2Obj;
                pair<int, Vec3d> lightTracing = this->spaceTracer((ph - lighters[i].Pl)/norm(ph - lighters[i].Pl), lighters[i].Pl);
                if(lightTracing.first == objInx) colorSum4f += lighters[i].colorSDcal(Peye, ph, normVec, textureClr[0],
                                                    textureClr[1], textureClr[2], this -> ColorH3);
            }
            colorSum4f += this -> k0*textureClr[0] + b*kb*this -> ColorH3;

            pair<int, Vec3d> rayTracing;
            nr = nep - 2 * normVec.dot(nep) * normVec;
            rayTracing = this -> spaceTracer(nr, ph);

            return colorSum4f + 0.4* this ->colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1);
        }
        else if(spheres[sphereInx].texture == TX_INFINITE_SPHERE){
            vector<Vec4f> textureClr = spheres[sphereInx].TexColorDeter(ph);
            return textureClr[1];
        }
    }
    else if(objInx >= planes.size() + spheres.size() && objInx < planes.size() + spheres.size() + cubes.size()){
        //cubes
        int cubesInx = objInx - planes.size() - spheres.size();
        if(norm(cubes[cubesInx].pointSave-ph) > 1e-2) return Vec4f(0,0,0,1);
        if(cubes[cubesInx].texture == TX_TOTAL_RL){
            Vec3d nr, nHE;
            pair<int, Vec3d> rayTracing;
            Vec3d normVec = cubes[cubesInx].normalSave;
            nr = nep - 2 * normVec.dot(nep) * normVec;
            rayTracing = this -> spaceTracer(nr, ph);
            //If rayTracing fail to hit object

            nHE = (Peye - ph)/norm(Peye - ph);
            double cosEH = nHE.dot(normVec);
            double b = std::max(((1-cosEH)-0.92)/0.08,0.0);
            return b*this -> kb*this -> ColorH3 + 0.8* this -> colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1);
        }
        else if(cubes[cubesInx].texture == TX_HALF_RL_HALF_RR){
            Vec3d nr, nHE, vRR; // nr is reflection, vRR is refraction
            double nRR, delta, C, aRR, bRR;
            pair<int, Vec3d> rayTracing, rayRR;
            Vec3d normVec = cubes[cubesInx].normalSave;
            nr = nep - 2 * normVec.dot(nep) * normVec;
            rayTracing = this -> spaceTracer(nr, ph);
            nHE = (Peye - ph)/norm(Peye - ph);
            double cosEH = nHE.dot(normVec);
            double b = std::max(((1-cosEH)-0.92)/0.08,0.0);
            if(objFrom == objInx) nRR = 1/cubes[cubesInx].nRR;
            else nRR = cubes[cubesInx].nRR;
            //total reflection
            C = -nep.dot(normVec);
            delta = (C*C-1)/(nRR*nRR) + 1;
            if(delta < 0) //delta < 0 is reflection
                return b*this -> kb*this -> ColorH3 + 0.8* this -> colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1);
            aRR = -1/nRR;
            bRR = C/nRR - sqrt(delta);
            vRR = aRR * (-nep) + bRR * normVec;
            rayRR = this -> spaceTracer(vRR, ph);
            return b*this -> kb*this -> ColorH3 + 0.15* this -> colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1) + 0.8* this -> colorReturn(objInx, rayRR.first, rayRR.second, nr, tierCnt + 1);
        }
        else if(cubes[cubesInx].texture == TX_OPAQUE){
            Vec4f colorSum4f(0,0,0,0);

            double b; //the coefficient for hale light.
            double cosEH; //cos(angle between reflect light and eye)
            Vec3d rHE, nLH, nHE, normVec;

            nHE = (Peye - ph)/norm(Peye - ph);
            normVec = cubes[cubesInx].normalSave;
            cosEH = nHE.dot(normVec);
            b = std::max(((1-cosEH)-0.92)/0.08,0.0);

            vector<Vec4f> textureClr = cubes[cubesInx].colorSave;

            for(int i = 0; i < lighters.size(); i ++){
                //perform rayTracing to the light and test if the hit point is the same as this point, if not, calculate the color and add it
                double dis_Hit2Obj;
                pair<int, Vec3d> lightTracing = this->spaceTracer((ph - lighters[i].Pl)/norm(ph - lighters[i].Pl), lighters[i].Pl);
                if(lightTracing.first == objInx) colorSum4f += lighters[i].colorSDcal(Peye, ph, normVec, textureClr[0],
                                                    textureClr[1], textureClr[2], this -> ColorH3);
            }
            colorSum4f += this -> k0*textureClr[0] + b*kb*this -> ColorH3;

            return colorSum4f;
        }
        else if(cubes[cubesInx].texture == TX_HALF_OPAQUE_HALF_RL){

            Vec4f colorSum4f(0,0,0,0);

            double b; //the coefficient for hale light.
            double cosEH; //cos(angle between reflect light and eye)
            Vec3d rHE, nLH, nHE, normVec, nr;

            nHE = (Peye - ph)/norm(Peye - ph);
            normVec = cubes[cubesInx].normalSave;
            cosEH = nHE.dot(normVec);
            b = std::max(((1-cosEH)-0.92)/0.08,0.0);

            vector<Vec4f> textureClr = cubes[cubesInx].colorSave;

            for(int i = 0; i < lighters.size(); i ++){
                //perform rayTracing to the light and test if the hit point is the same as this point, if not, calculate the color and add it
                double dis_Hit2Obj;
                pair<int, Vec3d> lightTracing = this->spaceTracer((ph - lighters[i].Pl)/norm(ph - lighters[i].Pl), lighters[i].Pl);
                if(lightTracing.first == objInx) colorSum4f += lighters[i].colorSDcal(Peye, ph, normVec, textureClr[0],
                                                    textureClr[1], textureClr[2], this -> ColorH3);
            }
            colorSum4f += this -> k0*textureClr[0] + b*kb*this -> ColorH3;

            pair<int, Vec3d> rayTracing;
            nr = nep - 2 * normVec.dot(nep) * normVec;
            rayTracing = this -> spaceTracer(nr, ph);

            return colorSum4f + 0.8* this ->colorReturn(objInx, rayTracing.first, rayTracing.second, nr, tierCnt + 1);
        }
    }
}
