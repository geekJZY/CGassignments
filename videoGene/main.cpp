
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>

using namespace std;
using namespace cv;

int main()
{

    VideoWriter writer("Video1.avi", CV_FOURCC('M', 'J', 'P', 'G'), 24.0, Size(1920, 1080));
    Mat frame;
    string dir = "/home/jzy/code_OpenCV/CGassignments/ziyujiang_finalvideo/1/";
    for(int i = 0; i < 120; i ++)
    {
        frame = imread(dir + "img" + to_string(i) + ".png");
        if(! frame.data){
            cout << "read image fail!" << endl;
            return 0;
        }
        writer << frame;
        imshow("video", frame);
        waitKey(10);
    }
}
