#include "BMP.h"

BMP::BMP()
{
    imageinfo = new ImageInfo;
    fileinfo = new FileInfo;
}
void BMP::readImage(string filename)
{
    ifstream file(filename , ios::binary); //opening file

    char filetemp[14];                     //store file informations
    char imagetemp[40];

    file.read(filetemp , sizeof(filetemp)); //read first 14 byte
    fileinfo->setFileHeader(filetemp);      //send allData
    file.read(imagetemp , sizeof(imagetemp));
    imageinfo->setImageHeader(imagetemp);
    fileinfo  ->  setOffSet(54);

    file.seekg(54);    //move cursor offset

    padding = 0;

    while((imageinfo->getWidth()*3+padding)%4 != 0) padding++; //calculate padding

    data = new BYTE[imageinfo -> getBiSize()];  //determine data size
    BYTE * pointerOfData = data;                //point data
    BYTE buffer[imageinfo->getWidth()*3];       //temprature memory

    for(int i=0;i<(int)(imageinfo->getBiSize()/sizeof(buffer));i++)
	{
        file.read((char*)buffer,sizeof(buffer));        //read first row image
        memcpy(pointerOfData,buffer,sizeof(buffer));    //copy buffer adresses
        pointerOfData += sizeof(buffer);                //move pointer
        file.read((char*)buffer,padding);               //move cursor
    }

    file.close();
}

void BMP::saveGrayScale(string filename)
{
    grayData = new BYTE[imageinfo->getBiSize()/3];
    BYTE * iterator = grayData;
    BYTE* pointerOfData = grayData;


    for(int i=0;i<(int)imageinfo->getBiSize();i+=3)
	{
        *iterator = BYTE((data[i]*0.21+data[i+1]*0.72+data[i+2]*0.07));
        iterator++;
    }

    ofstream file(filename + ".bmp" , ios::binary);
    file.write((char *)fileinfo->getAllHeader() , 14);
    file.write((char *)imageinfo->getAllInfo() , 40);


    unsigned int pixelNumber = 0;

    for(int i=0;i<(int)imageinfo->getBiSize()/3;i++)
	{
        file.write((char*)pointerOfData,1);
        file.write((char*)pointerOfData,1);
        file.write((char*)pointerOfData++,1);

        pixelNumber++;

        if(pixelNumber == imageinfo->getWidth())
		{
            BYTE pad = 0;
            for(int i=0;i<padding;i++) file.write((char*)&pad,1);
            pixelNumber = 0;
        }
    }

    file.close();
}
void BMP::saveGrayScale(BYTE * temp , string filename)
{
    BYTE * pointerOfData = temp;

    ofstream file(filename + ".bmp" , ios::binary);
    file.write((char *)fileinfo->getAllHeader() , 14);
    file.write((char *)imageinfo->getAllInfo() , 40);

    unsigned int pixelNumber = 0;

    for(int i=0;i<(int)imageinfo->getBiSize()/3;i++)
	{
        file.write((char*)pointerOfData,1);
        file.write((char*)pointerOfData,1);
        file.write((char*)pointerOfData++,1);

        pixelNumber++;

        if(pixelNumber == imageinfo->getWidth())
		{
            BYTE pad = 0;
            for(int i=0;i<padding;i++) file.write((char*)&pad,1);
            pixelNumber = 0;
        }
    }
    file.close();
}

void BMP::saveImage(string filename)
{
    ofstream file(filename + ".bmp" , ios::binary);
    file.write((char *)fileinfo->getAllHeader() , 14);
    file.write((char *)imageinfo->getAllInfo() , 40);

    BYTE * pointerOfData = data;
    unsigned int pixelNumber = 0;

    for(int i=0;i<(int)imageinfo->getBiSize();i++)
	{
        file.write((char*)pointerOfData++,1);

        pixelNumber++;
        if(pixelNumber >= imageinfo->getWidth() * 3)
		{
            BYTE pad = 0;
            for(int i=0;i<padding;i++) file.write((char*)&pad,1);
            pixelNumber = 0;
            i += padding;
        }
    }

    file.close();
}

BYTE* BMP::zoom(int x1, int y1, int x2, int y2, int w, int h, BYTE* buffer1) 
{
	if (x2 < x1)swap(x2, x1);
	if (y2 < y1)swap(y2, y1);

	int dX = x2 - x1;
	int dY = y2 - y1;
	BYTE* buffer2=new BYTE[dX * dY];
	int k = 0;

	for (int i = y1; i < y2; i++)
		for (int j = x1; j < x2; j++) 
		{
			if (j <= x2 && j >= x1) 
			{
				buffer2[k] = buffer1[i * w + j]; k++;
			}
		}
	return buffer2;
}

BYTE * BMP :: kMeans(BYTE* image, int w, int h, int N/*tag number*/)
{
   	int hist[256];//Store histogram
	for (int i = 0; i < 256; i++)
		hist[i] = 0;//Set to zero
	for (int i = 0; i < w*h; i++)
	{
		hist[image[i]] = hist[image[i]] + 1;//increment
	}
	//////////////////////////Tag number assignment
	int* t = new int[N];//Present tag value
	int* t_ = new int[N];//Next tag value
	int* color = new int[N];//Store intended tag value

	srand(time(NULL));
	for (int i = 0; i < N; i++)
	{
		t_[i] =rand()%256;
		color[i] = i *255 / (N - 1);
	}
	/////////////////////////bubble sort
	int i, j;
	bool swapped;
	for (i = 0; i < N - 1; i++)
	{
		swapped = false;
		for (j = 0; j < N - i - 1; j++)
		{
			if (t_[j] > t_[j + 1])
			{
				swap(t_[j], t_[j + 1]);
				swapped = true;
			}
		}

		// IF no two elements were swapped by inner loop, then break
		if (swapped == false)
			break;
	}
	/////////////////////////
	int tags[256];//Store tags and make modification on it
	bool check;//Check for present and next tag values
	/////////////////////////Variable for calculating center
	int* mx = new int[N];//mi*xi
	int* m = new int[N];//mi
	int* center = new int[N];//(mi1*xi1+mi2*xi2+......)/(mi1+mi2....)

	do
	{
		check = true;
		//////////////////////////Give new value
		for (int i = 0; i < N; i++)
		{
			t[i] = t_[i];
		}
		//////////////////////////Calculate Distance and set tags
		for (int i = 0; i < 256; i++)
		{
			int min = 0;
			for (int j = 1; j < N; j++)
			{
				if (abs(  t[j] - i ) < abs(  t[min] - i) )
				{
					min = j;
				}
			}
			tags[i] = t[min];
		}
		//////////////////////////Set variable to zero
		for (int i = 0; i < N; i++)
		{
			mx[i] = 0;
			m[i] = 0;
		}
		//////////////////////////Find values
		for (int i = 0; i < 256; i++)
		{
			for (int j = 0; j < N; j++)
			{
				if (tags[i] == t[j])
				{
					mx[j] += hist[i] * i;
					m[j] += hist[i];
				}
			}
		}
		//////////////////////////Calculate Center of Tags
		for (int i = 0; i < N; i++)
		{
			center[i] = mx[i] / (m[i] == 0 ? 1 : m[i]);
			t_[i] = center[i];

			if (t_[i] != t[i]) {
				check = false;
			}
		}
	} while (!check);//If true terminate the loop
	//////////////////////////Give initial value to tags with color
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (tags[i] == t[j])
			{
				tags[i] = color[j];
				break;
			}
		}
	}
	//////////////////////////Set image's new values
	for (int i = 0; i < w * h; i++)
	{
		image[i] = tags[image[i]];
	}
	//////////////////////////Return image
	return image;

}

BYTE * BMP :: erosion(BYTE* image, int width, int height) 
{
	BYTE* data = new BYTE[width*height];

	for (int i = 0; i < width*height; i++)data[i] = image[i];
	for (int i = 0; i < width; i++) 
	{
		for (int j = 0; j < height; j++) 
		{
			if (image[i + j * width] == 0) 
			{
				if		(i > 0			&& image[i - 1 +	width * j]		 == 255)	data[i + width * j] = 255;
				else if (j > 0			&& image[i +		width * (j - 1)] == 255)	data[i + width * j] = 255;
				else if (i + 1 < width	&& image[i + 1 +	width * j]		 == 255)	data[i + width * j] = 255;
				else if (j + 1 < height && image[i +		width * (j + 1)] == 255)	data[i + width * j] = 255;
			}
		}
	}

	return data;
}

BYTE * BMP :: dilatation(BYTE* image, int width, int height)
{
	BYTE* data = new BYTE[width*height];
	for (int i = 0; i < width*height; i++)data[i] = image[i];
	for (int i = 0; i < width; i++) 
	{
		for (int j = 0; j < height; j++) 
		{
			if(image[i + j * width] == 0) 
			{
				if		(i > 0			&& image[i - 1 + width * j]		  == 255)				data[(i - 1)+ width * j] = 0;
				if (j > 0			&& image[i	   + width * (j - 1)] == 255)				data[i + width * (j-1)] = 0;
				if (i + 1 < width	&& image[i + 1 + width * j]		  == 255)				data[(i + 1) + width * j] = 0;
				if (j + 1 < height && image[i	   + width * (j + 1)] == 255)				data[i + width * (j+1)] = 0;

			}
		}
	}
	return data;
}

BYTE* BMP::label(BYTE* image, int W, int H)
{

	BYTE* label = new BYTE[W * H];
	int min = 0;// minimum edge value of the 3x3 matrix
	int labelnum = 0;// label's value
	for (int i = 0; i < W * H; i++)
	{
		label[i] = 0;//set all labels to zero
	}

	int first = 0, last = 0;
	bool flag = false;

	do 
	{
		flag = false;
		for (int y = 1; y < H - 1; y++)
		{
			for (int x = 1; x < W - 1; x++)
			{
				if (image[x + y * W] == 0) 
				{
					first = label[x + y * W];
					min = 0;
					for (int j = -1; j <= 1; j++) {//calculate  minimum edge value of the 3x3 matrix
						for (int i = -1; i <= 1; i++) 
						{
							if (label[x + i + (y + j) * W] != 0 && min == 0) 
							{
								min = label[x + i + (y + j) * W];
							}

							else if (label[x + i + (y + j) * W] < min && label[x + i + (y + j) * W] != 0 && min != 0) 
							{
								min = label[x + i + (y + j) * W];
							}
						}
					}
					if (min != 0)  label[x + y * W] = min;
					else if (min == 0) 
					{
						labelnum++;
						label[x + y * W] = labelnum;
					}
					last = label[x + y * W];

					if (first != last) flag = true;
				}
			}
		}
	} while (flag);

	return label;
}

double BMP::calcMoment(int p, int q, BYTE* object, int width, int height)
{
	// calc (p+q)th moment of object
	double sum = 0;
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			sum += object[x + y * width] * pow(x, p) * pow(y, q);
		}
	}
	return sum;
}

double BMP::calcCentralMoment(int p, int q, BYTE* object, int width, int height, double normalizeParameter00, double normalizeParameter01 , double normalizeParameter10)
{
	// calculate central moment
	double sum = 0;
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			sum += object[x + y * width] * pow((x - (normalizeParameter10 / normalizeParameter00 )), p) * pow((y - (normalizeParameter01 / normalizeParameter00) ), q);
		}
	}
	return sum;
}

double BMP::calcNormalizedMoment(int p, int q, BYTE* image,int w, int h, double normalizeParameter00, double normalizeParameter01, double normalizeParameter10)
{
	// calculate normalized moments
	double temp = (((double)p + (double)q) / 2) + 1;
	return calcCentralMoment(p, q, image, w, h, normalizeParameter00, normalizeParameter01, normalizeParameter10) / pow(normalizeParameter00, temp);
}

void BMP::calcInvariantMoments(BYTE* image, int w, int h, double invariantMoments[7])
{
	// calculate 7 invariant moments


	double normalizeParameter00 = calcMoment(0, 0, image, w, h);
	double normalizeParameter01 = calcMoment(0, 1, image, w, h);
	double normalizeParameter10 = calcMoment(1, 0, image, w, h);

	double mom20 = calcNormalizedMoment(2, 0, image, w, h, normalizeParameter00, normalizeParameter01, normalizeParameter10), mom02 = calcNormalizedMoment(0, 2, image, w, h, normalizeParameter00, normalizeParameter01, normalizeParameter10);
	double mom11 = calcNormalizedMoment(1, 1, image, w, h, normalizeParameter00, normalizeParameter01, normalizeParameter10), mom30 = calcNormalizedMoment(3, 0, image, w, h, normalizeParameter00, normalizeParameter01, normalizeParameter10);
	double mom12 = calcNormalizedMoment(1, 2, image, w, h, normalizeParameter00, normalizeParameter01, normalizeParameter10), mom21 = calcNormalizedMoment(2, 1, image, w, h, normalizeParameter00, normalizeParameter01, normalizeParameter10);
	double mom03 = calcNormalizedMoment(0, 3, image, w, h, normalizeParameter00, normalizeParameter01, normalizeParameter10);

	invariantMoments[0] = mom20 + mom02;

	invariantMoments[1] = pow((mom20 - mom02), 2) + 4 * (pow(mom11, 2));

	invariantMoments[2] = pow((mom30 - 3 * mom12), 2) +
		pow((3 * mom21 - mom03), 2);

	invariantMoments[3] = pow((mom30 + mom12), 2) +
		pow((mom21 + mom03), 2);

	invariantMoments[4] = (mom30 - 3 * mom12) *
		(mom30 + mom12) *
		(pow(mom30 + mom12, 2) - 3 * pow(mom21 + mom03, 2)) +
		(3 * mom21 - mom03) * (mom21 + mom03) *
		(pow(3 * (mom30 + mom12), 2) - pow(mom21 + mom03, 2));

	invariantMoments[5] = (mom20 - mom02) * (pow(mom30 + mom12, 2) -
		pow(mom21 + mom03, 2)) + (4 * mom11* (mom30 + mom12) *
			mom21 + mom03);

	invariantMoments[6] = (3 * mom21 - mom03) * (mom30 + mom12) * (pow(mom30 + mom12, 2) -
		3 * pow(mom21 + mom03, 2)) - (mom30 - 3 * mom12) * (mom21 + mom03) *
		(3 * pow(mom30 + mom12, 2) - pow(mom21 + mom03, 2));


	for (int i = 0; i < 7; i++)
	{
		invariantMoments[i] = -1 * copysign(1.0, invariantMoments[i]) * log10(abs(invariantMoments[i]));
	}//Note that hu[0] is not comparable in magnitude as hu[6]. We can use use a log transform given below to bring them in the same range

	for (int i = 0; i < 7; i++)
		cout << invariantMoments[i] << " ";
	cout << "\n";

}

int* BMP::distinction(BYTE* label, int w, int h)
{
	BYTE* buffer = new BYTE[w * h];
	int* dis = new int[w * h + 1]; // farklı tag'leri tutmak icin
	for (int i = 0; i < w * h; i++)
	{
		buffer[i] = label[i];
		dis[i] = 0;
	}
	dis[w * h] = 0;

	std::sort(buffer, buffer + (w * h - 1)); //tagleri sırala k->b
	for (int i = 0; i < w * h - 1; i++) //dis'in son elemanı tag sayisini tutuyor
	{										//dis[w*h]-1 = nesne sayısı ilk eleman zemin (0) )
		if (buffer[i] != buffer[i + 1])
		{
			dis[dis[w * h]] = buffer[i];
			dis[w * h] += 1;

		}
	}

	for (int i = 0; i < dis[w * h]; i++)  //etiketler ve nesne print icin 
		cout << dis[i] << " ";
	cout << "\n\n" << dis[w * h]<<"\n\n";

	return dis;
}

double* BMP::recognation(BYTE* etiket, int w, int h, int* dis,double* moments,BYTE* grayScale)
{
	int n = dis[w * h ];
	//cout << n << "\n\n";

	//double* moments = new double[n*7];
	for (int x = 1; x < n; x++)
	{/////////////////////nesnelerin koordinatlarını bulma
		int x1 = -1, y1 = -1, x2 = -1, y2 = -1;
		bool flag = false;
		for (int j = 0; j < h; j++)
		{
			for (int i = 0; i < w; i++)
			{
				if ((int(etiket[i + j * w]) == dis[x]) && y1 == -1)
				{
					y1 = j;
				}
				else if ((int(etiket[w-i + (h-1-j) * w]) == dis[x]) && y2 == -1)
				{
					y2 = h - j;
				}

			}
		}
		for (int i = 0; i < w; i++)
		{
			for (int j = 0; j < h; j++)
			{
				if ((int(etiket[i + j * w]) == dis[x]) && x1 == -1)
				{
					x1 = i;
				}
				else if ((int(etiket[(w-i) + (h-1-j) * w]) == dis[x]) && x2 == -1)
				{
					x2 = w - i;
				}

			}
		}
		///////////////////////////////////////////////////////////////////
		//koordinatlara göre crop işlemi
		int dx = x2 - x1, dy = y2 - y1;
		//cout << x1 <<"-"<< x2 <<"-----"<< y1<<"-"<<y2 << endl;
		////////////////////////////////////////////////////////////////////
		//Etikete ve koordinat değerlerine göre cisimleri dikdörtgen içine almak
		drawRectangle(grayScale,x1,x2,y1,y2,w,h);
		////////////////////////////////////////////////////////////////////
		BYTE* obje = zoom(x1, y1, x2, y2, w, h, etiket);
		for (int i = 0; i < dx * dy; i++)
		{
			if ((int)obje[i] != dis[x])
			{
				obje[i] = 0;
			}
			else
			{
				obje[i] = 1;
			}
				//cout <<(int) obje[i];
		}
		double temp[7];
		calcInvariantMoments(obje, dx, dy, temp);
		for (int i = 0; i < 7; i++)
		{
			moments[x * 7 + i] = temp[i];
			//printf("%.12f  ", temp[i]);// temp[i] << " - ";
		}
		cout << "\n\n";
	}
	return moments;
}

void BMP::matchMoments(double* moments, int count)
{

	
	int dataNumber = 7;
	double epsilon = 0.3;

	//double database[7] = { 0.691898 , 3.42744  ,3.31272  ,3.58146  ,7.11494 ,- 2.23171 ,- 7.05517 //yıldız

	//};
	double database[49] = { 0.67661, 3.35843, 5.9718, 6.08609, 12.6059, 3.4112, 12.2338, //yıldız
							0.783572, 3.6836, 5.01136, 6.25786, -11.8925, -8.09966, -27.0249, //dikdörtgen
							0.69457, 3.38651, 5.25801, 5.82178, 11.1649, 3.32939, 11.627, //beşgen
							-1.03432, -1.92485, -2.98406, -2.76913, -5.48359, -3.47605, 3.30704,//nohut 
							0.745008, 2.18415, 3.79954, 4.82116, 9.02189, 2.95941, -10.0044,//badem 
							0.740741, 2.12562, 4.18878, 5.16142, 9.54202, -3.31631, -11.0005,//mercimek 
							0.792649, 3.25535, 4.29606, 6.34578, -11.731, -3.30692, -11.7971  //cubuk


	};

	int* dataCount = new int[dataNumber];

	for (int i = 0; i < dataNumber; i++)
		dataCount[i] = 0;

	string* dataName = new string[dataNumber];
	dataName[0]="yildiz";
	dataName[1]="dikdortgen";
	dataName[2]="besgen";
	dataName[3]="nohut";
	dataName[4]="badem";
	dataName[5]="mercimek";
	dataName[6]="cubuk";

	double sum = 0;
	for (int i = 1; i < count; i++)
	{
		for (int k = 0; k < dataNumber; k++)
		{
			sum = 0;
			for (int j = 0; j < 7; j++)
			{
				sum += abs(1/moments[i  * 7 + j] - 1/database[k*7+j]);
			}
			if (sum <= epsilon) {
				dataCount[k]++;
			}
		}
	}

	for (int i=0; i < dataNumber; i++) 
	{
		cout << dataName[i] << "---"<< dataCount[i] << "\n";
	}
}
void BMP::finalisation()
{
	BYTE* tempData = new BYTE[imageinfo->getWidth() * imageinfo->getHeight()];
    tempData = kMeans(grayData, imageinfo->getWidth()  , imageinfo->getHeight() , 2);
	saveGrayScale(tempData , "kmeans");

	tempData = erosion(tempData , imageinfo->getWidth() , imageinfo->getHeight() );
  	tempData = dilatation(tempData , imageinfo->getWidth() , imageinfo->getHeight() );
    saveGrayScale(tempData , "opening");

    
	tempData = dilatation(tempData , imageinfo->getWidth() , imageinfo->getHeight() );
    tempData = erosion(tempData , imageinfo->getWidth() , imageinfo->getHeight() );
    saveGrayScale(tempData , "closing");
 	

    BYTE* etiket=	label(tempData, imageinfo->getWidth() , imageinfo->getHeight());
	saveGrayScale(etiket , "labeling");


	int* dis = distinction(etiket, imageinfo->getWidth() , imageinfo->getHeight());
	double* moments = new double[(dis[imageinfo->getWidth()  * imageinfo->getHeight()]-1)*7];

	recognation(etiket, imageinfo->getWidth() , imageinfo->getHeight(), dis,moments, tempData);

	matchMoments(moments, dis[imageinfo->getWidth()  * imageinfo->getHeight()]);


}

void BMP::drawRectangle(BYTE* grayScale,int x1,int x2,int y1,int y2,int w,int h)
{
	for(int j=(y1)*w;j<(y2)*w;j+=w)
	{
		grayScale[j+x1-1]=0;
		grayScale[j+x2+1]=0;
	}

	for(int i=x1-1;i<x2+2	;i++)
	{
		grayScale[(y1)*w+i]=0;
		grayScale[(y2)*w+i]=0;
	}
	
	saveGrayScale(grayScale,"nesneler");
}
