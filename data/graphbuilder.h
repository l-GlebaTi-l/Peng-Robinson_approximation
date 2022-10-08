#include <windows.h>        /// only for struct BITMAPFILEHEADER and BITMAPINFOHEADER, which can change between different Windows/Linux versions. So, it VERSION CONTROL only.

long double stp(long double x);
long double stp(long double x,  int n);
long double sign(long double x);
long double module(long double x);
long double max3(long double *x);
long double min3(long double *x);
int n_in_char(const char* str);

int resolution_vysota;
int resolution_shirina;

typedef struct iXY {
    int x;
    int y;
} iXY, iPoint2, ScreenPoint;

long double stp(long double x)
{
    return x*x;
}

long double stp(long double x,  int n)
{
    long double res = 1;
    if(n > 0)
        for(int i = 0; i < n; i++){
            res *= x;
        }
    else
        for(int i = 0; i < -n; i++){
            res /= x;
        }
    return res;
}

long double sign(long double x)
{
    if(x > 0)
        return 1;
    else if(x < 0)
        return -1;
    return 0;
}

long double module(long double x)
{
    if(x < 0)
        return -x;
    return x;
}

long double max3(long double *x)
{
    if(x[0] > x[1] && x[0] > x[2])
        return x[0];
    else if(x[1] > x[2])
        return x[1];
    return x[2];
}

long double min3(long double *x)
{
    if(x[0] < x[1] && x[0] < x[2])
        return x[0];
    else if(x[1] < x[2])
        return x[1];
    return x[2];
}

int n_in_char(const char* str)
{
    int i = 0;
    while(str[i] != '\0')
        i++;
    return i;
}

unsigned char *loadBMP_32bit( const char *fname, int &mx, int &my )
{
    mx = my = -1;
    FILE *f = fopen( fname, "rb" );
    if( !f ) return NULL;
    BITMAPFILEHEADER bmpfhead;
    BITMAPINFOHEADER bmphead;

    size_t res;

    res = fread( &bmpfhead, 1, sizeof(BITMAPFILEHEADER), f );
    res += fread( &bmphead, 1, sizeof(BITMAPINFOHEADER), f );
    if( res != sizeof(BITMAPINFOHEADER) + sizeof(BITMAPFILEHEADER) ){
        fclose(f);
        return NULL;
    }

    fseek( f, 0, SEEK_END);
    int filesize = ftell(f);

    fseek( f, bmpfhead.bfOffBits, SEEK_SET);

    if(
//        bmphead.biSizeImage != filesize - 54 ||
        bmphead.biClrUsed != 0    ||
        bmphead.biPlanes   != 1    ||
       (bmphead.biSize!=40 && bmphead.biSize!=108 && bmphead.biSize!=124)||
//        bmpfhead.bfOffBits != 14+bmphead.biSize ||

        bmphead.biWidth <1 || bmphead.biWidth >10000 ||
        bmphead.biHeight<1 || bmphead.biHeight>10000 ||
        bmphead.biBitCount    != 32/* ||
        bmphead.biCompression !=  0*/
        )
    {
            fclose(f);
            return NULL;
    }
    mx = bmphead.biWidth;
    my = bmphead.biHeight;
    int mx4 = 4*mx;
    unsigned char *tmp_buf = new unsigned  char[mx4*my];
    res = fread( tmp_buf, 1, mx4*my, f);
    if( (int)res != mx4*my ) { delete []tmp_buf; fclose(f); return NULL; }

    fclose(f);

//    unsigned char *v = new unsigned char[mx*my*4];

//    v = tmp_buf;
//    unsigned char *ptr = v;
//    for(int y = my-1; y >= 0; y--) {
//        unsigned char *pRow = tmp_buf + mx4*y;
//        for(int x=0; x< mx; x++) {
//            *ptr++ = *pRow;
//            *ptr++ = *(pRow + 1);
//            *ptr++ = *(pRow + 2);
//            *ptr++ = *(pRow + 3);
//            pRow += 4;
//        }
//    }
//    delete []tmp_buf;
    return (tmp_buf);
}

unsigned char *loadBMP( const char *fname, int &mx, int &my )
{
    mx = my = -1;
    FILE *f = fopen( fname, "rb" );
    if( !f ) return NULL;
    BITMAPFILEHEADER bmpfhead;
    BITMAPINFOHEADER bmphead;

    size_t res;

    res = fread( &bmpfhead, 1, sizeof(BITMAPFILEHEADER), f );
    res += fread( &bmphead, 1, sizeof(BITMAPINFOHEADER), f );
    if( res != sizeof(BITMAPINFOHEADER) + sizeof(BITMAPFILEHEADER) ){
        fclose(f);
        return NULL;
    }

    fseek( f, 0, SEEK_END);
    int filesize = ftell(f);

    fseek( f, bmpfhead.bfOffBits, SEEK_SET);

    if(
//        bmphead.biSizeImage != filesize - 54 ||
        bmphead.biClrUsed != 0    ||
        bmphead.biPlanes   != 1    ||
       (bmphead.biSize!=40 && bmphead.biSize!=108 && bmphead.biSize!=124)||
        bmpfhead.bfOffBits != 14+bmphead.biSize ||

        bmphead.biWidth <1 || bmphead.biWidth >10000 ||
        bmphead.biHeight<1 || bmphead.biHeight>10000 ||
        bmphead.biBitCount    != 24 ||
        bmphead.biCompression !=  0
        )
    {
        if(bmphead.biBitCount == 32)
            return loadBMP_32bit(fname, mx, my);
        else{
            fclose(f);
            return NULL;
        }
    }
    mx = bmphead.biWidth;
    my = bmphead.biHeight;
    int mx3 = (3*mx+3) & (-4);
    unsigned char *tmp_buf = new unsigned  char[mx3*my];
    res = fread( tmp_buf, 1, mx3*my, f);
    if( (int)res != mx3*my ) { delete []tmp_buf; fclose(f); return NULL; }

    fclose(f);

    unsigned char *v = new unsigned char[mx*my*4];

    unsigned char *ptr = v;
    unsigned char *pRow = tmp_buf;
    for(int y = 0; y < my; y++) {
        for(int x = 0; x < mx; x++) {
            *ptr++ = *pRow;
            *ptr++ = *(pRow + 1);
            *ptr++ = *(pRow + 2);
            *ptr++ = 255;
            pRow += 3;
        }
    }
    delete (tmp_buf);
    return v;
}

int saveBMP( const char *fname, unsigned char *v, int mx, int my )	// ¬ каждом элементе упаковано все три RGB-байта
{
    BITMAPFILEHEADER bmpfhead;
    BITMAPINFOHEADER bmphead;

//	memset( &bf, 0, sizeof(bmpfhead) );
//	memset( &bh, 0, sizeof(bmphead) );

    int mx3 = (3*mx+3) & (-4);
	int filesize = 54 + my*mx3;

    bmpfhead.bfType = 19778;
    bmpfhead.bfSize = filesize;
    bmpfhead.bfReserved1 = 0;
    bmpfhead.bfReserved2 = 0;
    bmpfhead.bfOffBits = 54;

    bmphead.biSize = 40;
    bmphead.biWidth = mx;
    bmphead.biHeight = my;
    bmphead.biPlanes = 1;
    bmphead.biBitCount = 24;
    bmphead.biCompression = 0;
    bmphead.biSizeImage = 0;
    bmphead.biXPelsPerMeter = 3780;
    bmphead.biYPelsPerMeter = 3780;
    bmphead.biClrUsed = 0;
    bmphead.biClrImportant = 0;

	FILE *f = fopen( fname, "wb" );
	if( !f ) return -1;
//	size_t res;

	// пишем заголовок
	fwrite( &bmpfhead, 1, sizeof(BITMAPFILEHEADER), f );
	fwrite( &bmphead, 1, sizeof(BITMAPINFOHEADER), f );
//	res = fwrite( &bh, 1, sizeof(BMPheader), f );
//	if( res != sizeof(BMPheader) ) { fclose(f); return -1; }

	// приготовим временный буфер
	unsigned char *tmp_buf = new unsigned char[mx3*my];
	// ѕеренос данных (не забудем про RGB->BGR)
	unsigned char *ptr = v;
	for(int y = 0; y < my; y++) {
		unsigned char *pRow = tmp_buf + mx3*y;
		for(int x = 0; x < mx; x++) {
			*pRow = *ptr++;
			*(pRow + 1) = *ptr++;
			*(pRow + 2) = *ptr++;
			pRow+=3;
			ptr++;
		}
	}
//    unsigned char temp_buf = new unsigned char[mx*my*3];

	fwrite( tmp_buf, 1, mx3*my, f );
	fclose(f);
	delete []tmp_buf;
	return 0;	// OK
}

class GRAPH_TEXT {
protected:
    int img_x = 0;
    int img_y = 0;
    int bx;     //  ширина картинки буквы
    int by;     //  высота картинки буквы
    unsigned char* shr;   //  shrifty - изображение шрифтов
    iXY cp;
public:
    virtual int Convert(char a) = 0;
    virtual void say(int x, int y, double razm, const char* text, unsigned char* ipixels, int Ix, int Iy){  //  важность = -100
        if(y + razm > Iy)
            return;
        double r_x, r_y;
        int tbx, tby;
        r_x = bx/razm;
        r_y = by/razm;
        tbx = bx * razm/by;
        tby = by * razm/by;
        int k = 0;
        cp.x = x;
        cp.y = y;
        unsigned char* image;
        image = ipixels;
        int n = 0;
        while(text[k] != '\0'){
            if(text[k] != '\n' && text[k] != '\r'){
                n = Convert(text[k]);
                unsigned char* temp_t;
                int temp;
                if(n != -1)
                    for(int j = 0; j < tby; j++){
                        for(int i = 0; i < tbx; i++){
                            temp_t = shr + (((int)(j*r_y))*img_x + bx*n + ((int)(i*r_x)))*4;
                            if(cp.y + j >= 0 && cp.y + j < Iy && cp.x + i >= 0 && cp.x + i < Ix){
                                temp = (cp.y + j)*Ix + cp.x + i;
                                if(temp_t[3] != 0){
                                    image[temp*4] = temp_t[0];
                                    image[temp*4 + 1] = temp_t[1];
                                    image[temp*4 + 2] = temp_t[2];
                                }
                            }
                            else if(cp.x + i >= Ix && cp.y + j >= 0){
                                cp.y -= tby;
                                cp.x = x;
                            }
                            else break;
                        }
                    }
                cp.x += tbx;
            }
            else{
                cp.y -= tby;
                cp.x = x;
            }
            k++;
        }
    }
};

class TEXT_Eng : public GRAPH_TEXT {
public:
    TEXT_Eng(){
        unsigned char* t_shr;
        t_shr = loadBMP("data/english.bmp", img_x, img_y);
        if(t_shr == NULL)
            cout << "ERROR IN LOADING LANGUARGE!!!" << endl;

        for(int j = 0; j < img_y; j++){
            for(int i = 0; i < img_x; i++){
                t_shr[j*img_x*4 + i*4] = 255 - t_shr[j*img_x*4 + i*4 + 3];
                t_shr[j*img_x*4 + i*4 + 1] = 255 - t_shr[j*img_x*4 + i*4 + 3];
                t_shr[j*img_x*4 + i*4 + 2] = 255 - t_shr[j*img_x*4 + i*4 + 3];
//                t_shr[j*img_x*4 + i*4 + 3] = 0;
            }
        }

        shr = t_shr;
        by = img_y;
        bx = by;    //  img_x / 32
    }
    int Convert(char a){
        switch(a){
            case '0':
                return 0;
            case '1':
                return 1;
            case '2':
                return 2;
            case '3':
                return 3;
            case '4':
                return 4;
            case '5':
                return 5;
            case '6':
                return 6;
            case '7':
                return 7;
            case '8':
                return 8;
            case '9':
                return 9;
            case 'A':
            case 'a':
                return 10;
            case 'B':
            case 'b':
                return 11;
            case 'C':
            case 'c':
                return 12;
            case 'D':
            case 'd':
                return 13;
            case 'E':
            case 'e':
                return 14;
            case 'F':
            case 'f':
                return 15;
            case 'G':
            case 'g':
                return 16;
            case 'H':
            case 'h':
                return 17;
            case 'I':
            case 'i':
                return 18;
            case 'J':
            case 'j':
                return 19;
            case 'K':
            case 'k':
                return 20;
            case 'L':
            case 'l':
                return 21;
            case 'M':
            case 'm':
                return 22;
            case 'N':
            case 'n':
                return 23;
            case 'O':
            case 'o':
                return 24;
            case 'P':
            case 'p':
                return 25;
            case 'Q':
            case 'q':
                return 26;
            case 'R':
            case 'r':
                return 27;
            case 'S':
            case 's':
                return 28;
            case 'T':
            case 't':
                return 29;
            case 'U':
            case 'u':
                return 30;
            case 'V':
            case 'v':
                return 31;
            case 'W':
            case 'w':
                return 32;
            case 'X':
            case 'x':
                return 33;
            case 'Y':
            case 'y':
                return 34;
            case 'Z':
            case 'z':
                return 35;
            case '-':
                return 36;
            case '!':
                return 37;
            case '?':
                return 38;
            case '.':
                return 39;
            case ',':
                return 40;
            case ':':
                return 41;
            case ';':
                return 42;
            case '(':
                return 43;
            case ')':
                return 44;
            case '_':
                return 45;
            case '=':
                return 46;
            case '+':
                return 47;
            case '*':
                return 48;
            case '<':
                return 49;
            case '>':
                return 50;
            case '@':
                return 51;
            case '%':
                return 52;
            case '\'':
                return 53;
            case '"':
                return 54;
            default:
                return -1;
        }
    }
};

TEXT_Eng eng;

void Graph_get_ready(unsigned char* graph, int Ix, int Iy)
{
    int dx, dy;
    for(int j = 0; j < Iy; j++){    /// clearning the image
        for(int i = 0; i < Ix; i++){
            graph[j*Ix*4 + i*4    ] = 255;    /// B
            graph[j*Ix*4 + i*4 + 1] = 255;    /// G
            graph[j*Ix*4 + i*4 + 2] = 255;    /// R
            graph[j*Ix*4 + i*4 + 3] = 0;
        }
    }

    for(int i = 0; i < Ix; i++){
        graph[i*4    ] = 0;    /// B
        graph[i*4 + 1] = 0;    /// G
        graph[i*4 + 2] = 0;    /// R
        graph[i*4 + 3] = 0;
        graph[(Iy - 1)*Ix*4 + i*4    ] = 0;    /// B
        graph[(Iy - 1)*Ix*4 + i*4 + 1] = 0;    /// G
        graph[(Iy - 1)*Ix*4 + i*4 + 2] = 0;    /// R
        graph[(Iy - 1)*Ix*4 + i*4 + 3] = 0;
    }
    for(int j = 0; j < Iy; j++){    ///
        graph[j*Ix*4    ] = 0;    /// B
        graph[j*Ix*4 + 1] = 0;    /// G
        graph[j*Ix*4 + 2] = 0;    /// R
        graph[j*Ix*4 + 3] = 0;
        graph[j*Ix*4 + (Ix - 1)*4    ] = 0;    /// B
        graph[j*Ix*4 + (Ix - 1)*4 + 1] = 0;    /// G
        graph[j*Ix*4 + (Ix - 1)*4 + 2] = 0;    /// R
        graph[j*Ix*4 + (Ix - 1)*4 + 3] = 0;
    }

    {                /// расчерчивание графика полосами
        if(Ix >= Iy){
            dx = (int)((double)Ix / ((double)PARAM_GRAPH_SCALE * (double)Ix / (double)Iy));
            dy = (int)((double)Iy / ((double)PARAM_GRAPH_SCALE*1.5));
        }
        else{
            dx = (int)(Ix / (int)((double)PARAM_GRAPH_SCALE));
            dy = (int)((double)Iy*1.5 / ((double)PARAM_GRAPH_SCALE * (double)Iy / (double)Ix));
        }

        for(int i = 0; i < Ix; i += dx){
            for(int j = 0; j < Iy; j++){
                graph[j*Ix*4 + i*4] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 1] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 2] = PARAM_GRAPH_LINES_INVISIBLE;
            }
        }
        for(int j = 0; j < Iy; j += dy){
            for(int i = 0; i < Ix; i++){
                graph[j*Ix*4 + i*4] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 1] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 2] = PARAM_GRAPH_LINES_INVISIBLE;
            }
        }
    }
}

void Graph_print_error(unsigned char* graph, long double approximation_error, int Ix, int Iy)
{
    char string_error[126];
    if(approximation_error < 0.001)
        sprintf(string_error, " max error=%.5f%%", (float)(approximation_error*100));
    else if(approximation_error < 0.01)
        sprintf(string_error, " max error=%.3f%%", (float)(approximation_error*100));
    else if(approximation_error < 0.1)
        sprintf(string_error, " max error=%.2f%%", (float)(approximation_error*100));
    else
        sprintf(string_error, " max error=%.1f%%", (float)(approximation_error*100));
    eng.say(0, 0, PARAM_GRAPH_VALUE_SCALE*1.5, string_error, graph, Ix, Iy);
}

/// Graph_numeric() changes Ix, Iy to new!!!!!!!!
void Graph_numeric(unsigned char *&graph, int &Ix, int &Iy, long double P0, long double P1, long double W0, long double W1)    /// ”величиваем изображение дл€ накидывани€ размерностей.   примеру, по 120 пикселей с каждой стороны (чтобы нормально влезли числа, пока без увеличени€ на параметры изображени€)
{
    int dx, dy;
    if(Ix >= Iy){
        dx = (int)((double)Ix / ((double)PARAM_GRAPH_SCALE * (double)Ix / (double)Iy));
        dy = (int)((double)Iy / ((double)PARAM_GRAPH_SCALE*1.5));
    }
    else{
        dx = (int)(Ix / (int)((double)PARAM_GRAPH_SCALE));
        dy = (int)((double)Iy*1.5 / ((double)PARAM_GRAPH_SCALE * (double)Iy / (double)Ix));
    }

    int addx = PARAM_ADD_PIXELS_TO_GRAPH_LEFT, addy = PARAM_ADD_PIXELS_TO_GRAPH_DOWN;     /// на сколько увеличить изображение слева и снизу
    int dIx, dIy;
    dIx = Ix + addx;
    dIy = Iy + addy;
    unsigned char *diagram = new unsigned char[dIx*4*dIy];

    for(int j = 0; j < dIy; j++){    /// clearning the image
        for(int i = 0; i < dIx; i++){
            diagram[j*dIx*4 + i*4    ] = 255;    /// B
            diagram[j*dIx*4 + i*4 + 1] = 255;    /// G
            diagram[j*dIx*4 + i*4 + 2] = 255;    /// R
            diagram[j*dIx*4 + i*4 + 3] = 0;
        }
    }

    for(int j = 0; j < Iy; j++){
        for(int i = 0; i < Ix*4; i++){
            diagram[(j+addy)*dIx*4 + addx*4 + i] = graph[j*Ix*4 + i];
        }
    }

    char str[20];
    int n;
    double PT, WT, dP, dW; //  P temp, W temp
    PT = P0;
    WT = W0;
    dP = (P1 - P0) / ((double)Ix / (double)dx);
    dW = (W1 - W0) / ((double)Iy / (double)dy);

    for(int cX = addx; cX < dIx; cX += dx){
        if(PT == (int)PT)
            sprintf(str, "%.f", PT);
        else if(PT < 1.)
            sprintf(str, "%.2f", PT);
        else if(PT < 10.)
            sprintf(str, "%.1f", PT);
        else if(PT < 100.)
            sprintf(str, "%.1f", PT);
        else
            sprintf(str, "%.f", PT);
        eng.say(cX, addy - (PARAM_GRAPH_VALUE_SCALE*1.25), PARAM_GRAPH_VALUE_SCALE, str, diagram, dIx, dIy);
        PT += dP;
    }
    for(int cY = addy; cY < dIy; cY += dy){
        if(WT == 0)
            sprintf(str, "%.f", WT);
        else if(WT == (int)WT)
            sprintf(str, "%.f", WT);
        else if(WT < 0.01)
            sprintf(str, "%.5f", WT);
        else if(WT < 0.1)
            sprintf(str, "%.4f", WT);
        else if(WT < 1.)
            sprintf(str, "%.3f", WT);
        else if(WT < 10.)
            sprintf(str, "%.2f", WT);
        else if(WT < 100.)
            sprintf(str, "%.1f", WT);
        else
            sprintf(str, "%.f", WT);
        n = n_in_char(str);
        eng.say(addx - n*PARAM_GRAPH_VALUE_SCALE, cY, PARAM_GRAPH_VALUE_SCALE, str, diagram, dIx, dIy);
        WT += dW;
    }
    eng.say(dIx - PARAM_GRAPH_VALUE_SCALE*1.5*5, 0, PARAM_GRAPH_VALUE_SCALE*1.5, "P,bar", diagram, dIx, dIy);
    eng.say(PARAM_GRAPH_VALUE_SCALE*1.5*2 + PARAM_GRAPH_VALUE_SCALE*1.5*0.3, dIy - PARAM_GRAPH_VALUE_SCALE*1.5, PARAM_GRAPH_VALUE_SCALE*1.5, "W", diagram, dIx, dIy);

    Ix = dIx; Iy = dIy;
    unsigned char* temp_to_delete = graph;
    graph = diagram;
    delete(temp_to_delete);
}

void Build_Line(unsigned char *image, int Ix, int Iy, int koordinateX, int kFrom, int kTo, unsigned char colorR, unsigned char colorG, unsigned char colorB)
{
    if((kFrom > Iy-1 && kTo > Iy-1) || (kFrom < 0 && kTo < 0))
        return;
    if(kFrom > Iy-1) kFrom = Iy-1;
    if(kFrom < 0) kFrom = 0;
    if(kTo > Iy-1) kTo = Iy-1;
    if(kTo < 0) kTo = 0;
    if(kFrom > kTo){
        swap(kFrom, kTo);
    }
//    if((kFrom == kTo && kTo == 0) || (kFrom == kTo && kTo == Iy-1))
//        return;

    for(int j = kFrom; j <= kTo; j++){
        image[j*Ix*4 + (koordinateX-1)*4    ] = colorB;    /// B
        image[j*Ix*4 + (koordinateX-1)*4 + 1] = colorG;    /// G
        image[j*Ix*4 + (koordinateX-1)*4 + 2] = colorR;    /// R
        image[j*Ix*4 + (koordinateX-1)*4 + 3] = 0;
    }
}

void Build_Line(unsigned char *image, int Ix, int Iy, int koordinateX, int kFrom, int kTo, int LINE_type_param, int &LINE_type_counter, unsigned char colorR, unsigned char colorG, unsigned char colorB)
{
    if((kFrom > Iy-1 && kTo > Iy-1) || (kFrom < 0 && kTo < 0))
        return;
    if(kFrom > Iy-1) kFrom = Iy-1;
    if(kFrom < 0) kFrom = 0;
    if(kTo > Iy-1) kTo = Iy-1;
    if(kTo < 0) kTo = 0;
    if(kFrom > kTo){
        swap(kFrom, kTo);
    }
//    if((kFrom == kTo && kTo == 0) || (kFrom == kTo && kTo == Iy-1))
//        return;

    for(int j = kFrom; j <= kTo; j++){
        LINE_type_counter++;
        if((LINE_type_counter / LINE_type_param)%2 != 1){
            image[j*Ix*4 + (koordinateX-1)*4    ] = colorB;    /// B
            image[j*Ix*4 + (koordinateX-1)*4 + 1] = colorG;    /// G
            image[j*Ix*4 + (koordinateX-1)*4 + 2] = colorR;    /// R
            image[j*Ix*4 + (koordinateX-1)*4 + 3] = 0;
        }
    }
}

void Build_Line(unsigned char *image, int Ix, int Iy, int X0, int Y0, int X1, int Y1, unsigned char colorR, unsigned char colorG, unsigned char colorB)
{
    if((X0 < 0 && X1 < 0) || (X0 > Ix-1 && X1 > Ix-1) || (Y0 < 0 && Y1 < 0) || (Y0 > Iy-1 && Y1 > Iy-1)) return;


    long double k, b;

    if((module(X1 - X0)) >= (module(Y1 - Y0))){
        int koorY;
        k = ((long double)(Y1 - Y0))/((long double)(X1 - X0));
        b = (long double)(Y0) - (((long double)(Y1 - Y0))/((long double)(X1 - X0))) * (long double)(X0);
        for(int koorX = X0; koorX <= X1; koorX++){
            koorY = k*(long double)(koorX) + b;
            if(koorX >= 0 && koorX < Ix-1 && koorY >= 0 && koorY < Iy-1){
                image[koorY*Ix*4 + koorX*4    ] = colorB;
                image[koorY*Ix*4 + koorX*4 + 1] = colorG;
                image[koorY*Ix*4 + koorX*4 + 2] = colorR;
                image[koorY*Ix*4 + koorX*4 + 3] = 0;
            }
        }
    }
    else{
        int koorX;
        k = ((long double)(X1 - X0))/((long double)(Y1 - Y0));
        b = (long double)(X0) - (((long double)(X1 - X0))/((long double)(Y1 - Y0))) * (long double)(Y0);
        if(Y0 <= Y1)
            for(int koorY = Y0; koorY <= Y1; koorY++){
                koorX = k*(long double)(koorY) + b;
                if(koorX >= 0 && koorX < Ix-1 && koorY >= 0 && koorY < Iy-1){
                    image[koorY*Ix*4 + koorX*4    ] = colorB;
                    image[koorY*Ix*4 + koorX*4 + 1] = colorG;
                    image[koorY*Ix*4 + koorX*4 + 2] = colorR;
                    image[koorY*Ix*4 + koorX*4 + 3] = 0;
                }
            }
        else
            for(int koorY = Y0; koorY >= Y1; koorY--){
                koorX = k*(long double)(koorY) + b;
                if(koorX >= 0 && koorX < Ix-1 && koorY >= 0 && koorY < Iy-1){
                    image[koorY*Ix*4 + koorX*4    ] = colorB;
                    image[koorY*Ix*4 + koorX*4 + 1] = colorG;
                    image[koorY*Ix*4 + koorX*4 + 2] = colorR;
                    image[koorY*Ix*4 + koorX*4 + 3] = 0;
                }
            }
    }
}

