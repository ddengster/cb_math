
/*
VectorMath.h - v1.0

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
*/

#pragma once

#include <iostream>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <cstring>

#define PI 3.14159265358979323846f
#define ZERO_ERROR_RANGE 1e-07f
#define PIOVER180 0.0174532925199432f
inline float DegToRad(float deg) { return (PI / 180.0f) * deg; }
inline float RadToDeg(float rad) { return (rad / PI) * 180.f; }

/******************************************************************************/
/*!
  \brief
  2d vector
*/
/******************************************************************************/
class Vector2
{
public:
  float x, y;

public:
  //Constructors
  inline Vector2() { }
  inline Vector2(float fX, float fY) : x(fX), y(fY) { }
  inline Vector2(const Vector2& rhs) : x(rhs.x), y(rhs.y) { }
  inline explicit Vector2(float scalar) : x(scalar), y(scalar) { }
  inline explicit Vector2(float coordinate[2]) : x(coordinate[0]), y(coordinate[1]) { }
  //explicit Vector2(const Vector3& rhs); use function instead { start.x, start.y};

  //access op
  inline float  operator[](unsigned int i) const { assert(i < 2); return *(&x + i); }
  inline float& operator[](unsigned int i)       { assert(i < 2); return *(&x + i); }
  inline const float* ptr() const                { return &x; }
  inline float* ptr()                            { return &x; }
  inline void Set(float fx, float fy)            { x = fx; y = fy; }

  //assignment op
  inline Vector2& operator=(Vector2 rhs)  { x = rhs.x; y = rhs.y; return *this; }
  inline Vector2& operator=(float scalar) { x = y = scalar; return *this; }

  //comparison op
  inline bool operator==(Vector2 rhs) const { return ((x == rhs.x) && (y == rhs.y)); }
  inline bool operator!=(Vector2 rhs) const { return ((x != rhs.x) || (y != rhs.y)); }
  inline bool operator<(Vector2 rhs) const  { return ((x < rhs.x) && (y < rhs.y)); }
  inline bool operator>(Vector2 rhs) const  { return ((x > rhs.x) && (y > rhs.y)); }
  inline bool IsEqual(Vector2 rhs, float error = 0) const { return ((fabs(x - rhs.x) <= error) && (fabs(y - rhs.y) <= error)); }

  //basic member arithmetic op
  inline Vector2 operator+(Vector2 rhs) const  { return Vector2(x + rhs.x, y + rhs.y); }
  inline Vector2 operator-(Vector2 rhs) const  { return Vector2(x - rhs.x, y - rhs.y); }
  inline Vector2 operator*(float scalar) const { return Vector2(x * scalar, y * scalar); }
  inline Vector2 operator/(float scalar) const { assert(scalar != 0.0); float fInv = 1.0f / scalar; return Vector2(x * fInv, y * fInv); }
  inline Vector2 operator*(Vector2 rhs) const  { return Vector2(x * rhs.x, y * rhs.y); }
  inline Vector2 operator/(Vector2 rhs) const  { return Vector2(x / rhs.x, y / rhs.y); }
  inline Vector2 operator-() const { return Vector2(-x, -y); }

  inline Vector2& operator+=(Vector2 rhs)  { x += rhs.x; y += rhs.y; return *this; }
  inline Vector2& operator-=(Vector2 rhs)  { x -= rhs.x; y -= rhs.y; return *this; }
  inline Vector2& operator*=(Vector2 rhs)  { x *= rhs.x; y *= rhs.y; return *this; }
  inline Vector2& operator/=(Vector2 rhs)  { x /= rhs.x; y /= rhs.y; return *this; }
  inline Vector2& operator+=(float scalar) { x += scalar; y += scalar; return *this; }
  inline Vector2& operator-=(float scalar) { x -= scalar; y -= scalar; return *this; }
  inline Vector2& operator*=(float scalar) { x *= scalar; y *= scalar; return *this; }
  inline Vector2& operator/=(float scalar) { assert(scalar != 0.0); x /= scalar; y /= scalar; return *this; }

  //advanced member op
  inline float Length() const                     { return sqrt(x * x + y * y); }
  inline float SquaredLength() const              { return x * x + y * y; }
  inline float Distance(Vector2 vec) const        { return (*this - vec).Length(); }
  inline float SquaredDistance(Vector2 vec) const { return (*this - vec).SquaredLength(); }
  inline float DotProduct(Vector2 rhs) const      { return x * rhs.x + y * rhs.y; }
  inline Vector2 MidPoint(Vector2 vec) const      { return Vector2((x + vec.x) * 0.5f, (y + vec.y) * 0.5f); }

  /** Normalises the vector to a length of 1. **/
  inline void Normalise() {
    float fLength = sqrt(x * x + y * y);

    // Will also work for zero-sized vectors
    if (fLength > ZERO_ERROR_RANGE)
    {
      float fInvLength = 1.0f / fLength;
      x *= fInvLength;
      y *= fInvLength;
    }
  }
  inline Vector2 NormalisedCopy(void) const { Vector2 ret = *this; ret.Normalise(); return ret;  }

  inline void MakeFloor(Vector2 cmp)           { if (cmp.x < x) x = cmp.x; if(cmp.y < y) y = cmp.y; }
  inline void MakeCeil(Vector2 cmp)            { if (cmp.x > x) x = cmp.x; if (cmp.y > y) y = cmp.y; }
  inline Vector2 Perpendicular(void) const     { return Vector2(-y, x); }  // (-y, x) 90 degrees anticlockwise
  inline Vector2 AntiPerpendicular(void) const { return Vector2(y, -x); }  // (y, -x) 90 degrees clockwise
  inline float CrossProduct(Vector2 rhs) const { return x * rhs.y - y * rhs.x; } //similar to computing determinant ad - bc
  inline bool IsZeroLength(void) const {
    float sqlen = (x * x) + (y * y);
    return (sqlen < (ZERO_ERROR_RANGE * ZERO_ERROR_RANGE));
  }
  inline Vector2 Reflect(Vector2 normal) const { return Vector2(*this - (2 * this->DotProduct(normal) * normal)); }
  inline Vector2 RotateDeg(float deg) const
  {
    float rad = DegToRad(deg);
    return RotateRad(rad);
  }
  inline Vector2 RotateRad(float rad) const
  {
    return Vector2(x * cos(rad) - y * sin(rad), x * sin(rad) + y * sin(rad));
  }

  static const Vector2 ZERO;
  static const Vector2 UNIT_X;
  static const Vector2 UNIT_Y;
  static const Vector2 NEGATIVE_UNIT_X;
  static const Vector2 NEGATIVE_UNIT_Y;
  static const Vector2 UNIT_SCALE;

  //static arithmetic func
  static inline void Add(Vector2& results, Vector2 lhs, Vector2 rhs)        { results.x = lhs.x + rhs.x; results.y = lhs.y + rhs.y; }
  static inline void Subtract(Vector2& results, Vector2 lhs, Vector2 rhs)   { results.x = lhs.x - rhs.x; results.y = lhs.y - rhs.y; }
  static inline void Multiply(Vector2& results, Vector2 lhs, float scalar)  { results.x = lhs.x * scalar; results.y = lhs.y * scalar; }
  static inline void Multiply(Vector2& results, float scalar, Vector2 rhs)  { results.x = scalar * rhs.x; results.y = scalar * rhs.y; }
  static inline void Multiply(Vector2& results, Vector2 lhs, Vector2 rhs)   { results.x = lhs.x * rhs.x; results.y = lhs.y * rhs.y; }
  static inline void Divide(Vector2& results, float scalar, Vector2 rhs)    { results.x = scalar / rhs.x; results.y = scalar / rhs.y; }
  static inline void Divide(Vector2& results, Vector2 lhs, float scalar)    { results.x = lhs.x / scalar; results.y = lhs.y / scalar; }
  static inline void Divide(Vector2& results, Vector2 lhs, Vector2 rhs)     { results.x = lhs.x / rhs.x; results.y = lhs.y / rhs.y; }

  //global friend func
  inline friend Vector2 operator+(Vector2 lhs, float rhs) { return Vector2(lhs.x + rhs, lhs.y + rhs); }
  inline friend Vector2 operator+(float lhs, Vector2 rhs) { return Vector2(lhs + rhs.x, lhs + rhs.y); }
  inline friend Vector2 operator-(Vector2 lhs, float rhs) { return Vector2(lhs.x - rhs, lhs.y - rhs); }
  inline friend Vector2 operator-(float lhs, Vector2 rhs) { return Vector2(lhs - rhs.x, lhs - rhs.y); }
  inline friend Vector2 operator*(float lhs, Vector2 rhs) { return Vector2(lhs * rhs.x, lhs * rhs.y); }
  inline friend Vector2 operator/(float lhs, Vector2 rhs) { return Vector2(lhs / rhs.x, lhs / rhs.y); }

  /** Function for writing to a stream.  */
  inline /*shared*/ friend std::ostream& operator << (std::ostream& o, Vector2 v)
  {
    o << "Vector2(" << v.x << ", " << v.y << ")";
    return o;
  }

  inline friend std::istream& operator>>(std::istream& inStream, Vector2& vector)
  {
    inStream >> vector.x;
    inStream.ignore();

    inStream >> vector.y;

    return inStream;
  }
};

/******************************************************************************/
//2D Vector with integer coords
class Vector2I
{
public:
  int x, y;

public:
  inline Vector2I() { }
  inline Vector2I(int fX, int fY) : x(fX), y(fY) { }

  inline Vector2 ToVec2F() { return Vector2((float)x, (float)y); }
  inline Vector2I& operator=(Vector2I rhs) { x = rhs.x; y = rhs.y; return *this; }
  inline Vector2I& operator=(Vector2 rhs)  { x = (int)rhs.x; y = (int)rhs.y; return *this; }
  inline void RoundFrom(Vector2 rhs) { 
    x = (int)(rhs.x > 0.f ? rhs.x + 0.5f : rhs.x - 0.5f);  y = (int)(rhs.y > 0.f ? rhs.y + 0.5f : rhs.y - 0.5f);
  }

  //comparison op
  inline bool operator==(Vector2I rhs) const { return ((x == rhs.x) && (y == rhs.y)); }
  inline bool operator!=(Vector2I rhs) const { return ((x != rhs.x) || (y != rhs.y)); }

  //basic member arithmetic op
  inline Vector2I operator+(Vector2I rhs) const { return Vector2I(x + rhs.x, y + rhs.y); }
  inline Vector2I operator-(Vector2I rhs) const { return Vector2I(x - rhs.x, y - rhs.y); }
  inline Vector2I operator*(Vector2I rhs) const { return Vector2I(x * rhs.x, y * rhs.y); }
  inline Vector2I operator/(Vector2I rhs) const { return Vector2I(x / rhs.x, y / rhs.y); }
  inline Vector2I operator-() const { return Vector2I(-x, -y); }

  inline Vector2I operator+(int rhs) const { return Vector2I(x + rhs, y + rhs); }
  inline Vector2I operator-(int rhs) const { return Vector2I(x - rhs, y - rhs); }
  inline Vector2I operator*(int rhs) const { return Vector2I(x * rhs, y * rhs); }
  inline Vector2I operator/(int rhs) const { return Vector2I(x / rhs, y / rhs); }

  inline Vector2I& operator+=(Vector2I rhs) { x += rhs.x; y += rhs.y;  return *this; }
  inline Vector2I& operator-=(Vector2I rhs) { x -= rhs.x; y -= rhs.y;  return *this; }
  inline Vector2I& operator*=(Vector2I rhs) { x *= rhs.x; y *= rhs.y;  return *this; }
  inline Vector2I& operator/=(Vector2I rhs) { x /= rhs.x; y /= rhs.y;  return *this; }

  //global friend func
  inline friend Vector2I operator+(Vector2I lhs, int rhs) { return Vector2I(lhs.x + rhs, lhs.y + rhs); }
  inline friend Vector2I operator+(int lhs, Vector2I rhs) { return Vector2I(lhs + rhs.x, lhs + rhs.y); }
  inline friend Vector2I operator-(Vector2I lhs, int rhs) { return Vector2I(lhs.x - rhs, lhs.y - rhs); }
  inline friend Vector2I operator-(int lhs, Vector2I rhs) { return Vector2I(lhs - rhs.x, lhs - rhs.y); }
  inline friend Vector2I operator*(int lhs, Vector2I rhs) { return Vector2I(lhs * rhs.x, lhs * rhs.y); }
  inline friend Vector2I operator/(int lhs, Vector2I rhs) { return Vector2I(lhs / rhs.x, lhs / rhs.y); }

  //advanced member op
  inline Vector2I Perpendicular(void) const     { return Vector2I(-y, x); }  // (-y, x) 90 degrees anticlockwise
  inline Vector2I AntiPerpendicular(void) const { return Vector2I(y, -x); }  // (y, -x) 90 degrees clockwise
  inline int DotProduct(Vector2I rhs) const     { return x * rhs.x + y * rhs.y; }
  inline int CrossProduct(Vector2I rhs) const   { return (x * rhs.y) - (y * rhs.x); }

  /** Function for writing to a stream.  */
  inline /*shared*/ friend std::ostream& operator << (std::ostream& o, Vector2I v)
  {
    o << "Vector2I(" << v.x << ", " << v.y << ")";
    return o;
  }

  inline friend std::istream& operator>>(std::istream& inStream, Vector2I& vector)
  {
    inStream >> vector.x;
    inStream.ignore();

    inStream >> vector.y;

    return inStream;
  }
};

/******************************************************************************/
/*!
  \brief
  3d vector
*/
/******************************************************************************/
class Vector3
{
public:
  float x, y, z;

public:
  //Constructors
  inline Vector3() { }
  inline Vector3(float fX, float fY, float fZ) : x(fX), y(fY), z(fZ) { }
  inline Vector3(const Vector3& rhs) : x(rhs.x), y(rhs.y), z(rhs.z) { }
  inline explicit Vector3(float scalar) : x(scalar), y(scalar), z(scalar) { }
  inline explicit Vector3(float coordinate[3]) : x(coordinate[0]), y(coordinate[1]), z(coordinate[2]) { }
  inline explicit Vector3(Vector2 rhs) : x(rhs.x), y(rhs.y), z(0.0f) { } 

  //access op
  inline float  operator[](unsigned int i) const { assert(i < 3); return *(&x + i); }
  inline float& operator[](unsigned int i) { assert(i < 3); return *(&x + i); }
  inline const float* ptr() const { return &x; }
  inline float* ptr() { return &x; }
  inline void Set(float fx, float fy, float fz) { x = fx; y = fy; z = fz; }

  //assignment op
  inline Vector3& operator=(Vector3 rhs) { x = rhs.x; y = rhs.y; z = rhs.z; return *this; }

  //comparison op
  inline bool operator<(Vector3 rhs) const { return ((x < rhs.x) && (y < rhs.y) && (z < rhs.z)); }
  inline bool operator>(Vector3 rhs) const { return ((x > rhs.x) && (y > rhs.y) && (z > rhs.z)); }
  inline bool IsEqual(Vector3 rhs, float error = 0.0f) const { return ((fabs(x - rhs.x) <= error) && (fabs(y - rhs.y) <= error) && (fabs(z - rhs.z) <= error)); }

  //basic member arithmetic op
  inline Vector3 operator+(Vector3 rhs) const  { return Vector3(x + rhs.x, y + rhs.y, z + rhs.z); }
  inline Vector3 operator-(Vector3 rhs) const  { return Vector3(x - rhs.x, y - rhs.y, z - rhs.z); }
  inline Vector3 operator*(float scalar) const { return Vector3(x * scalar, y * scalar, z * scalar); }
  inline Vector3 operator/(float scalar) const {
    assert(scalar != 0.0);
    float fInv = 1.0f / scalar;
    return Vector3(x * fInv, y * fInv, z * fInv);
  }
  inline Vector3 operator*(Vector3 rhs) const { return Vector3(x * rhs.x, y * rhs.y, z * rhs.z); }
  inline Vector3 operator/(Vector3 rhs) const { return Vector3(x / rhs.x, y / rhs.y, z * rhs.z); }
  inline Vector3 operator-() const { return Vector3(-x, -y, -z); }

  inline Vector3& operator+=(Vector3 rhs)  { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
  inline Vector3& operator-=(Vector3 rhs)  { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
  inline Vector3& operator*=(Vector3 rhs)  { x *= rhs.x; y *= rhs.y; z *= rhs.z; return *this; }
  inline Vector3& operator/=(Vector3 rhs)  { x /= rhs.x; y /= rhs.y; z /= rhs.z; return *this; }
  inline Vector3& operator+=(float scalar) { x += scalar; y += scalar; z += scalar; return *this; }
  inline Vector3& operator-=(float scalar) { x -= scalar; y -= scalar; z -= scalar; return *this; }
  inline Vector3& operator*=(float scalar) { x *= scalar; y *= scalar; z *= scalar; return *this; }
  inline Vector3& operator/=(float scalar) { assert(scalar != 0.0); x /= scalar; y /= scalar; z /= scalar; return *this; }

  //advanced member arithmetic op
  inline float Length() const                     { return sqrt(x * x + y * y + z * z); }
  inline float SquaredLength() const              { return x * x + y * y + z * z; }
  inline float Distance(Vector3 vec) const        { return (*this - vec).Length(); }
  inline float SquaredDistance(Vector3 vec) const { return (*this - vec).SquaredLength(); }
  inline float DotProduct(Vector3 rhs) const      { return x * rhs.x + y * rhs.y + z * rhs.z; }
  inline Vector3 MidPoint(Vector3 vec) const      { return Vector3((x + vec.x) * 0.5f, (y + vec.y) * 0.5f, (z + vec.z) * 0.5f); }

  /** Normalises the vector to a length of 1. **/
  inline void Normalise()  {
    float fLength = sqrt(x * x + y * y + z * z);

    // Will also work for zero-sized vectors
    if (fLength > ZERO_ERROR_RANGE)
    {
      float fInvLength = 1.0f / fLength;
      x *= fInvLength;
      y *= fInvLength;
      z *= fInvLength;
    }
  }
  inline Vector3 NormalisedCopy(void) const { Vector3 ret = *this; ret.Normalise(); return ret; }

  inline void MakeFloor(Vector3 cmp) { if (cmp.x < x) x = cmp.x; if (cmp.y < y) y = cmp.y; if (cmp.z < z) z = cmp.z; }
  inline void MakeCeil(Vector3 cmp)  { if (cmp.x > x) x = cmp.x; if (cmp.y > y) y = cmp.y; if (cmp.z > z) z = cmp.z; }
  inline Vector3 CrossProduct(Vector3 rhs) const {
    return Vector3(y * rhs.z - z * rhs.y,
      z * rhs.x - x * rhs.z,
      x * rhs.y - y * rhs.x);
  }
  inline float CrossProductXZPlane(Vector3 rhs) const { return (x * rhs.z) - (z *  rhs.x); }

  inline bool IsZeroLength(void) const {
    float sqlen = (x * x) + (y * y) + (z * z);
    return (sqlen < (ZERO_ERROR_RANGE * ZERO_ERROR_RANGE));
  }
  inline Vector3 Reflect(Vector3 normal) const { return Vector3(*this - (2 * this->DotProduct(normal) * normal)); }

  inline bool DirectionEquals(Vector3 rhs, float toleranceRad = 1e-03f) const //involves 2 normalise, inefficient
  {
    float dot = NormalisedCopy().DotProduct(rhs.NormalisedCopy());
    float angle = acos(dot);

    return (abs(angle) <= toleranceRad);
  }
  inline void Clamp(Vector3 min, Vector3 max) {
    x = x < min.x ? min.x : (x > max.x) ? max.x : x;
    y = y < min.y ? min.y : (y > max.y) ? max.y : y;
    z = z < min.z ? min.z : (z > max.z) ? max.z : z;
  }
  inline void Clamp(float min, float max) {
    x = x < min ? min : (x > max) ? max : x;
    y = y < min ? min : (y > max) ? max : y;
    z = z < min ? min : (z > max) ? max : z;
  }
  void GetRotationToVector3(float& deg, Vector3& rotationaxis, const Vector3& othervector) const;

  static const Vector3 ZERO;
  static const Vector3 UNIT_X;
  static const Vector3 UNIT_Y;
  static const Vector3 UNIT_Z;
  static const Vector3 NEGATIVE_UNIT_X;
  static const Vector3 NEGATIVE_UNIT_Y;
  static const Vector3 NEGATIVE_UNIT_Z;
  static const Vector3 UNIT_SCALE;

  //static arithmetic func
  static inline void Add(Vector3& results, Vector3 lhs, Vector3 rhs)       { results.x = lhs.x + rhs.x; results.y = lhs.y + rhs.y; results.z = lhs.z + rhs.z; }
  static inline void Subtract(Vector3& results, Vector3 lhs, Vector3 rhs)  { results.x = lhs.x + rhs.x; results.y = lhs.y + rhs.y; results.z = lhs.z + rhs.z; }
  static inline void Multiply(Vector3& results, Vector3 lhs, float scalar) { results.x = lhs.x * scalar; results.y = lhs.y * scalar; results.z = lhs.z * scalar; }
  static inline void Multiply(Vector3& results, float scalar, Vector3 rhs) { results.x = rhs.x * scalar; results.y = rhs.y * scalar; results.z = rhs.z * scalar; }
  static inline void Multiply(Vector3& results, Vector3 lhs, Vector3 rhs)  { results.x = lhs.x * rhs.x; results.y = lhs.y * rhs.y; results.z = lhs.z * rhs.z; }
  static inline void Divide(Vector3& results, float scalar, Vector3 rhs)   { results.x = scalar / rhs.x; results.y = scalar / rhs.y; results.z = scalar / rhs.z; }
  static inline void Divide(Vector3& results, Vector3 lhs, float scalar)   { results.x = lhs.x / scalar; results.y = lhs.y / scalar; results.z = lhs.z / scalar; }
  static inline void Divide(Vector3& results, Vector3 lhs, Vector3 rhs)    { results.x = lhs.x / rhs.x; results.y = lhs.y / rhs.y; results.z = lhs.z / rhs.z; }


  //global friend func
  inline friend Vector3 operator+(Vector3 lhs, float rhs) { return Vector3(lhs.x + rhs, lhs.y + rhs, lhs.z + rhs); }
  inline friend Vector3 operator+(float lhs, Vector3 rhs) { return Vector3(lhs + rhs.x, lhs + rhs.y, lhs + rhs.z); }
  inline friend Vector3 operator-(Vector3 lhs, float rhs) { return Vector3(lhs.x - rhs, lhs.y - rhs, lhs.z - rhs); }
  inline friend Vector3 operator-(float lhs, Vector3 rhs) { return Vector3(lhs - rhs.x, lhs - rhs.y, lhs - rhs.z); }
  inline friend Vector3 operator*(float lhs, Vector3 rhs) { return Vector3(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z); }
  inline friend Vector3 operator/(float lhs, Vector3 rhs) { return Vector3(lhs / rhs.x, lhs / rhs.y, lhs / rhs.z); }

  /** Function for writing to a stream.  */
  inline /*shared*/ friend std::ostream& operator << (std::ostream& o, Vector3 v)
  {
    o << "Vector3(" << v.x << ", " << v.y << ", " << v.z << ")";
    return o;
  }

  inline friend std::istream& operator>>(std::istream& inStream, Vector3& vector)
  {
    inStream >> vector.x;
    inStream.ignore();

    inStream >> vector.y;
    inStream.ignore();
    inStream >> vector.z;

    return inStream;
  }
};

inline Vector2 Vec3To2(Vector3 r) { return Vector2(r.x, r.y); }
/******************************************************************************/
/*!
  \brief
  4d vector
*/
/******************************************************************************/
class Vector4
{
public:
  //union for color compatibility
  float x, y, z, w;

public:
  //Constructors
  inline Vector4() { }
  inline Vector4(float fX, float fY, float fZ, float fW) : x(fX), y(fY), z(fZ), w(fW) { }
  inline Vector4(const Vector4& rhs) : x(rhs.x), y(rhs.y), z(rhs.z), w(rhs.w) { }
  inline explicit Vector4(float scalar) : x(scalar), y(scalar), z(scalar), w(scalar) { }
  inline explicit Vector4(float coordinate[4]) : x(coordinate[0]), y(coordinate[1]), z(coordinate[2]), w(coordinate[3]) { }
  //inline explicit Vector4(const Vector3& rhs) : x(rhs.x), y(rhs.y), z(rhs.z), w(0.0f) { }

  //access op
  inline float  operator[](unsigned int i) const { assert(i < 4); return *(&x + i); }
  inline float& operator[](unsigned int i) { assert(i < 4); return *(&x + i); }
  inline const float* ptr() const { return &x; }
  inline float* ptr() { return &x; }
  inline void Set(float fx, float fy, float fz, float fw) { x = fx; y = fy; z = fz; w = fw; }

  //assignment op
  inline Vector4& operator=(Vector4 rhs) { x = rhs.x; y = rhs.y; z = rhs.z; w = rhs.w; return *this; }
  inline Vector4& operator=(Vector3 rhs) { x = rhs.x; y = rhs.y; z = rhs.z; w = 0;     return *this; }
  inline Vector4& operator=(float scalar) { x = y = z = w = scalar; return *this; }

  //comparison op
  inline bool operator==(Vector4 rhs) const { return ((x == rhs.x) && (y == rhs.y) && (z == rhs.z) && (w == rhs.w)); }
  inline bool operator!=(Vector4 rhs) const { return ((x != rhs.x) || (y != rhs.y) || (z != rhs.z) || (w != rhs.w)); }
  inline bool operator<(Vector4 rhs) const { return ((x < rhs.x) && (y < rhs.y) && (z < rhs.z) && (w < rhs.w)); }
  inline bool operator>(Vector4 rhs) const { return ((x > rhs.x) && (y > rhs.y) && (z > rhs.z) && (w > rhs.w)); }
  inline bool IsEqual(Vector4 rhs, float error = 0) const {
    return ((fabs(x - rhs.x) <= error) && (fabs(y - rhs.y) <= error) && (fabs(z - rhs.z) <= error) && (fabs(w - rhs.w) <= error));
  }

  //basic member arithmetic op
  inline Vector4 operator+(Vector4 rhs) const { return Vector4(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w); }
  inline Vector4 operator-(Vector4 rhs) const { return Vector4(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w); }
  inline Vector4 operator*(float scalar) const { return Vector4(x * scalar, y * scalar, z * scalar, w * scalar); }
  inline Vector4 operator/(float scalar) const {
    assert(scalar != 0.0);
    float fInv = 1.0f / scalar;
    return Vector4(x * fInv, y * fInv, z * fInv, w * fInv);
  }
  inline Vector4 operator*(Vector4 rhs) const { return Vector4(x * rhs.x, y * rhs.y, z * rhs.z, w * rhs.w); }
  inline Vector4 operator/(Vector4 rhs) const { return Vector4(x / rhs.x, y / rhs.y, z * rhs.z, w * rhs.w); }
  inline Vector4 operator-() const { return Vector4(-x, -y, -z, -w); }

  inline Vector4& operator+=(Vector4 rhs) { x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w; return *this; }
  inline Vector4& operator-=(Vector4 rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; w -= rhs.w; return *this; }
  inline Vector4& operator*=(Vector4 rhs) { x *= rhs.x; y *= rhs.y; z *= rhs.z; w *= rhs.w; return *this; }
  inline Vector4& operator/=(Vector4 rhs) { x /= rhs.x; y /= rhs.y; z /= rhs.z; w /= rhs.w; return *this; }
  inline Vector4& operator+=(float scalar) { x += scalar; y += scalar; z += scalar; w += scalar; return *this; }
  inline Vector4& operator-=(float scalar) { x -= scalar; y -= scalar; z -= scalar; w -= scalar; return *this; }
  inline Vector4& operator*=(float scalar) { x *= scalar; y *= scalar; z *= scalar; w *= scalar; return *this; }
  inline Vector4& operator/=(float scalar) { x /= scalar; y /= scalar; z /= scalar; w /= scalar; return *this; }

  //advanced member arithmetic op
  inline float Length() const                     { return sqrt(x * x + y * y + z * z + w * w); }
  inline float SquaredLength() const              { return x * x + y * y + z * z + w * w; }
  inline float Distance(Vector4 vec) const        { return (*this - vec).Length(); }
  inline float SquaredDistance(Vector4 vec) const { return (*this - vec).SquaredLength(); }
  inline float DotProduct(Vector4 rhs) const      { return x * rhs.x + y * rhs.y + z * rhs.z + w * rhs.w; }
  inline Vector4 MidPoint(Vector4 vec) const      { return Vector4((x + vec.x) * 0.5f, (y + vec.y) * 0.5f, (z + vec.z) * 0.5f, (w + vec.w) * 0.5f); }
  
  /** Normalises the vector to a length of 1. **/
  inline void Normalise() {
    float fLength = sqrt(x * x + y * y + z * z + w * w);

    // Will also work for zero-sized vectors
    if (fLength > ZERO_ERROR_RANGE)
    {
      float fInvLength = 1.0f / fLength;
      x *= fInvLength;
      y *= fInvLength;
      z *= fInvLength;
      w *= fInvLength;
    }
  }
  inline Vector4 NormalisedCopy(void) const { Vector4 ret = *this; ret.Normalise(); return ret; }

  inline void MakeFloor(Vector4 cmp) {
    if (cmp.x < x) x = cmp.x;
    if (cmp.y < y) y = cmp.y;
    if (cmp.z < z) z = cmp.z;
    if (cmp.w < w) z = cmp.w;
  }
  inline void MakeCeil(Vector4 cmp) {
    if (cmp.x > x) x = cmp.x;
    if (cmp.y > y) y = cmp.y;
    if (cmp.z > z) z = cmp.z;
    if (cmp.w > w) z = cmp.w;
  }
  inline Vector4 CrossProduct(Vector4 rhs) const {
    //Vector3 cross product, with w set to 0
    return Vector4(y * rhs.z - z * rhs.y,
      z * rhs.x - x * rhs.z,
      x * rhs.y - y * rhs.x,
      0);
  }
  inline bool IsZeroLength(void) const {
    float sqlen = (x * x) + (y * y) + (z * z) + (w * w);
    return (sqlen < (ZERO_ERROR_RANGE * ZERO_ERROR_RANGE));
  }
  inline Vector4 Reflect(Vector4 normal) const {
    return Vector4(*this - (2 * this->DotProduct(normal) * normal));
  }

  inline bool DirectionEquals(Vector4 rhs, float toleranceRad = 1e-03f) const {
    float dot = NormalisedCopy().DotProduct(rhs.NormalisedCopy());
    float angle = acos(dot);

    return (abs(angle) <= toleranceRad);
  }
  inline void Clamp(const Vector4 &min, const Vector4 &max) {
    x = x < min.x ? min.x : (x > max.x) ? max.x : x;
    y = y < min.y ? min.y : (y > max.y) ? max.y : y;
    z = z < min.z ? min.z : (z > max.z) ? max.z : z;
    w = w < min.w ? min.w : (w > max.w) ? max.w : w;
  }
  inline void Clamp(float min, float max) {
    x = x < min ? min : (x > max) ? max : x;
    y = y < min ? min : (y > max) ? max : y;
    z = z < min ? min : (z > max) ? max : z;
    w = w < min ? min : (w > max) ? max : w;
  }

  static const Vector4 ZERO;
  static const Vector4 UNIT_X;
  static const Vector4 UNIT_Y;
  static const Vector4 UNIT_Z;
  static const Vector4 UNIT_W;
  static const Vector4 NEGATIVE_UNIT_X;
  static const Vector4 NEGATIVE_UNIT_Y;
  static const Vector4 NEGATIVE_UNIT_Z;
  static const Vector4 NEGATIVE_UNIT_W;
  static const Vector4 UNIT_SCALE;

  //static arithmetic func
  static inline void Add(Vector4& results, Vector4 lhs, Vector4 rhs)       { results.x = lhs.x + rhs.x; results.y = lhs.y + rhs.y; results.z = lhs.z + rhs.z; results.w = lhs.w + rhs.w; }
  static inline void Subtract(Vector4& results, Vector4 lhs, Vector4 rhs)  { results.x = lhs.x - rhs.x; results.y = lhs.y - rhs.y; results.z = lhs.z - rhs.z; results.w = lhs.w - rhs.w; }
  static inline void Multiply(Vector4& results, Vector4 lhs, float scalar) { results.x = lhs.x * scalar; results.y = lhs.y * scalar; results.z = lhs.z * scalar; results.w = lhs.w * scalar; }
  static inline void Multiply(Vector4& results, float scalar, Vector4 rhs) { results.x = rhs.x * scalar; results.y = rhs.y * scalar; results.z = rhs.z * scalar; results.w = rhs.w * scalar; }
  static inline void Multiply(Vector4& results, Vector4 lhs, Vector4 rhs)  { results.x = lhs.x * rhs.x; results.y = lhs.y * rhs.y; results.z = lhs.z * rhs.z; results.w = lhs.w * rhs.w; }
  static inline void Divide(Vector4& results, float scalar, Vector4 rhs)   { results.x = scalar / rhs.x; results.y = scalar / rhs.y; results.z = scalar / rhs.z; results.w = scalar / rhs.w; }
  static inline void Divide(Vector4& results, Vector4 lhs, float scalar)   { results.x = lhs.x / scalar; results.y = lhs.y / scalar; results.z = lhs.z / scalar; results.w = lhs.w / scalar; }
  static inline void Divide(Vector4& results, Vector4 lhs, Vector4 rhs)    { results.x = lhs.x / rhs.x; results.y = lhs.y / rhs.y; results.z = lhs.z / rhs.z; results.w = lhs.w / rhs.w; }

  //global friend func
  inline friend Vector4 operator+(Vector4 lhs, float rhs) { return Vector4(lhs.x + rhs, lhs.y + rhs, lhs.z + rhs, lhs.w + rhs); }
  inline friend Vector4 operator+(float lhs, Vector4 rhs) { return Vector4(lhs + rhs.x, lhs + rhs.y, lhs + rhs.z, lhs + rhs.w); }
  inline friend Vector4 operator-(Vector4 lhs, float rhs) { return Vector4(lhs.x - rhs, lhs.y - rhs, lhs.z - rhs, lhs.w - rhs); }
  inline friend Vector4 operator-(float lhs, Vector4 rhs) { return Vector4(lhs - rhs.x, lhs - rhs.y, lhs - rhs.z, lhs - rhs.w); }
  inline friend Vector4 operator*(float lhs, Vector4 rhs) { return Vector4(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z, lhs * rhs.w); }
  inline friend Vector4 operator/(float lhs, Vector4 rhs) { return Vector4(lhs / rhs.x, lhs / rhs.y, lhs / rhs.z, lhs / rhs.w); }

  /** Function for writing to a stream.
  */
  inline /*shared*/ friend std::ostream& operator << (std::ostream& o, Vector4 v)
  {
    o << "Vector4(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
    return o;
  }

  inline friend std::istream& operator>>(std::istream& inStream, Vector4& vector)
  {
    inStream >> vector.x;
    inStream.ignore();
    inStream >> vector.y;
    inStream.ignore();
    inStream >> vector.z;
    inStream.ignore();
    inStream >> vector.w;

    return inStream;
  }
};

/****************************************************************************/

