#include "colorUn.h"

colorUnite::colorUnite(void){
}

colorUnite::colorUnite(Vec3d Ph_, Vec3d norm_, Vec4f colorH0_, Vec4f colorH1_, Vec4f colorH2_, Vec4f colorH3_){
    Ph = Ph_;
    norm = norm_;
    ColorH0 = colorH0_;
    ColorH1 = colorH1_;
    ColorH2 = colorH2_;
    ColorH3 = colorH3_;
}
