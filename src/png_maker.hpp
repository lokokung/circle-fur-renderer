#include <png.h>

#include "utilities.hpp"

class PNGMaker {
    public:
        int xres;
        int yres;
        
        PNGMaker(int xres, int yres);
        ~PNGMaker();

        void setPixel(int x, int y, float r, float g, float b);
        int saveImage(const char* output_file);

    private:
        Vec3f *screen;
};
