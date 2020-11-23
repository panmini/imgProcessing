#include<iostream>
#include<string>
#include<fstream>
#include <cstdint>
#include <cstring>
#include <math.h>
#include <algorithm>
#include "FileInfo.h"
#include "ImageInfo.h"

typedef unsigned char BYTE;
typedef unsigned short WORD;
typedef uint32_t DWORD;

class BMP
{
    public:
        BMP();

        FileInfo * fileinfo;
        ImageInfo * imageinfo;

        void readImage(string filename); //Reads Img from BMP File
        void saveGrayScale(string filename); //Convert bmp to intensity
        void saveGrayScale(BYTE * ,string filename);//Saves Intensity to same path as bmp.h (filename is the name of the file) 
        void saveImage(string filename);//Save the BMP file
        
        BYTE* zoom(int x1, int y1, int x2, int y2, int w, int h, BYTE* buffer1);
        void drawRectangle(BYTE* grayScale,int x1,int x2,int y1,int y2,int w,int h);

        BYTE* kMeans(BYTE* image, int w, int h, int N/*tag number*/); //Does kMeans operation 
        BYTE* erosion(BYTE* image, int width, int height); //Does Erosion operation
        BYTE* dilatation(BYTE* image, int width, int height);//Does Dilatation operation

        BYTE* label(BYTE* image, int W, int H); //Labelling
        int* distinction(BYTE* label, int w, int h);//Finds label counts and values
        void calcInvariantMoments(BYTE* image, int w, int h, double invariantMoments[7]);
        double calcNormalizedMoment(int p, int q, BYTE* image,int w, int h, double normalizeParameter00, double normalizeParameter01, double normalizeParameter10);
        double calcCentralMoment(int p, int q, BYTE* object, int width, int height, double normalizeParameter00, double normalizeParameter01 , double normalizeParameter10);
        double calcMoment(int p, int q, BYTE* object, int width, int height);
        double* recognation(BYTE* etiket, int w, int h, int* dis,double* moments,BYTE* grayScale);//Calculates objects' moment values.
        void matchMoments(double* moments, int count); //see if any of the moments matches given moments
	    
        void finalisation(); //Final function includes other func.
      
    private:


        BYTE * data; //Points RGB Image
        BYTE * grayData;//Points Intensity Image


        //BYTE * tempMatrix;

        int padding; //Width or Height must be dividible by 4 while operating img file conversions.
};