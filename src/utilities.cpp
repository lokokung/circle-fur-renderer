#include "utilities.hpp"

/* Constructs a triple of integers. */
Vec3i::Vec3i(int i, int j, int k) {
    this->i = i;
    this->j = j;
    this->k = k;
}

/* Constructs a triple of integers. */
Vec3i::Vec3i(int *t) {
    *this = *((Vec3i *) t);
}

/* Constructs a point in 3D space. */
Vec3f::Vec3f(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

/* Constructs a point in 3D space. */
Vec3f::Vec3f(float *v) {
    *this = *((Vec3f *) v);
}

/* Constructs a point in homogeneous space. */
Vec4f::Vec4f(float x, float y, float z, float w) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}

/* Constructs a point in homogeneous space. */
Vec4f::Vec4f(float *v) {
    *this = *((Vec4f *) v);
}
