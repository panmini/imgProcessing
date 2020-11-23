#include "BMP.h"

int main()
{
    BMP bmp;
    bmp.readImage("hop.bmp");
    bmp.saveGrayScale("grayScale");
   // bmp.cropImage(0,0,374,234);
    bmp.finalisation();

    return 0;
}
