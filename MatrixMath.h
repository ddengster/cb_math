
#pragma once

#include "VectorMath.h"

/******************************************************************************/
/*!
  \brief
    2x2 Matrix
*/
/******************************************************************************/
class Matrix2
{
public:
  union 
  {
    float m[4];
    float m_2d[2][2]; //[row][column]
  };

  static const Matrix2 ZERO;
  static const Matrix2 IDENTITY;

  // Constructors
  inline Matrix2() { }
  inline Matrix2(const Matrix2& rhs) { memcpy(m, rhs.m, 4 * sizeof(float)); }
  inline Matrix2(const float arr[2][2]) { memcpy(m, arr, 4 * sizeof(float)); }
  inline Matrix2(float _00, float _01, float _10, float _11)
  {
    m_2d[0][0] = _00; m_2d[0][1] = _01;
    m_2d[1][0] = _10; m_2d[1][1] = _11;
  }

  inline Matrix2(Vector2 row1, Vector2 row2)
  {
    m_2d[0][0] = row1.x; m_2d[0][1] = row1.y;
    m_2d[1][0] = row2.x; m_2d[1][1] = row2.y;
  }

  // Assignment operators
  inline Matrix2& operator *= (Matrix2 rhs) {
    Matrix2 lhs(*this);

    // actual matrix multiplication 
    // can use a while loop but im going for speed
    m_2d[0][0] = lhs.m_2d[0][0] * rhs.m_2d[0][0] + lhs.m_2d[0][1] * rhs.m_2d[1][0];
    m_2d[0][1] = lhs.m_2d[0][0] * rhs.m_2d[0][1] + lhs.m_2d[0][1] * rhs.m_2d[1][1];

    m_2d[1][0] = lhs.m_2d[1][0] * rhs.m_2d[0][0] + lhs.m_2d[1][1] * rhs.m_2d[1][0];
    m_2d[1][1] = lhs.m_2d[1][0] * rhs.m_2d[0][1] + lhs.m_2d[1][1] * rhs.m_2d[1][1];
    return *this;
  }
  // Binary operators
  inline Matrix2 operator* (Matrix2 rhs) { return Matrix2(*this) *= rhs; }
  inline Vector2 operator* (Vector2 rhs) {
    return Vector2(m_2d[0][0] * rhs.x + m_2d[0][1] * rhs.y,
      m_2d[1][0] * rhs.x + m_2d[1][1] * rhs.y);
  }

  inline void Identity()              { m_2d[0][0] = m_2d[1][1] = 1.0f;  m_2d[0][1] = m_2d[1][0] = 0.0f; }
  inline void Scale(float x, float y) { m_2d[0][0] = x; m_2d[1][1] = y;  m_2d[0][1] = m_2d[1][0] = 0.0f; }
  inline void RotRad(float angle)     { m_2d[0][0] = m_2d[1][1] = cos(angle); m_2d[0][1] = -(m_2d[1][0] = sin(angle)); }
  inline void RotDeg(float angle)     { angle *= PIOVER180;  m_2d[0][0] = m_2d[1][1] = cos(angle); m_2d[0][1] = -(m_2d[1][0] = sin(angle)); }
  inline void Transpose()             { float temp = m_2d[0][1]; m_2d[0][1] = m_2d[1][0]; m_2d[1][0] = temp; }
  void Inverse();
};


/******************************************************************************/
/*!
  \brief
  3x3 Matrix. Column major 3x3 matrix. ie. XYZ Translation to the rightmost column
  Top left element has row index 0, column index 0 and bottom right has row index 2, column index 2
*/
/******************************************************************************/
class Matrix3
{
public:
  union
  {
    float m[9];
    float m_2d[3][3]; // [row] [column]
  };

  static const Matrix3 ZERO;
  static const Matrix3 IDENTITY;

public:
  //Constructors
  inline Matrix3() { }
  inline Matrix3(const Matrix3& rhs)             { memcpy(m, rhs.m, 9 * sizeof(float)); }
  inline explicit Matrix3(const float arr[3][3]) { memcpy(m, arr, 9 * sizeof(float)); }
  Matrix3(float entry00, float entry01, float entry02,
    float entry10, float entry11, float entry12,
    float entry20, float entry21, float entry22) // no forced inline
  {
    m_2d[0][0] = entry00; m_2d[0][1] = entry01; m_2d[0][2] = entry02;
    m_2d[1][0] = entry10; m_2d[1][1] = entry11; m_2d[1][2] = entry12;
    m_2d[2][0] = entry20; m_2d[2][1] = entry21; m_2d[2][2] = entry22;
  }
  inline Matrix3(Vector3 row0, Vector3 row1, Vector3 row2)
  {
    m_2d[0][0] = row0.x; m_2d[0][1] = row0.y; m_2d[0][2] = row0.z;
    m_2d[1][0] = row1.x; m_2d[1][1] = row1.y; m_2d[1][2] = row1.z;
    m_2d[2][0] = row2.x; m_2d[2][1] = row2.y; m_2d[2][2] = row2.z;
  }

  //access op
  inline float* operator[](unsigned int row) const { assert(row < 3); return (float*)m_2d[row]; }
  inline const float* ptr() const                  { return m; }
  inline float* ptr()                              { return m; }

  //assignment op
  inline Matrix3& operator=(const Matrix3& rhs) { memcpy(m, rhs.m, 9 * sizeof(float)); return *this; }

  //mutators and accessors
  inline Vector3 GetColumn(unsigned int column) const { assert(column < 3); return Vector3(m_2d[0][column], m_2d[1][column], m_2d[2][column]); }
  inline void SetColumn(Vector3 vec, unsigned int column) {
    assert(column < 3);
    m_2d[0][column] = vec.x; m_2d[1][column] = vec.y; m_2d[2][column] = vec.z;
  }
  inline Vector3 GetRow(unsigned int row) const {
    assert(row < 3);
    return Vector3(m_2d[row][0], m_2d[row][1], m_2d[row][2]);
  }
  inline void SetRow(Vector3 vec, unsigned int row) {
    assert(row < 3);
    m_2d[row][0] = vec.x; m_2d[row][1] = vec.y; m_2d[row][2] = vec.z;
  }
  inline void ResetBottomRow() { m_2d[2][0] = m_2d[2][1] = 0.f; m_2d[2][2] = 1.f; }

  //comparison op
  inline bool operator==(const Matrix3& rhs) const {
    for (unsigned int i = 0; i<9; ++i)
      if (m[i] != rhs.m[i])
        return false;
    return true;
  }
  inline bool operator!=(const Matrix3& rhs) const { return !(*this == rhs); }

  //basic member arithmetic op
  inline Matrix3 operator*(const Matrix3& rhs) const //concatenation
  {
    Matrix3 res;
    res.m[0] = m_2d[0][0] * rhs.m_2d[0][0] + m_2d[0][1] * rhs.m_2d[1][0] + m_2d[0][2] * rhs.m_2d[2][0];
    res.m[1] = m_2d[0][0] * rhs.m_2d[0][1] + m_2d[0][1] * rhs.m_2d[1][1] + m_2d[0][2] * rhs.m_2d[2][1];
    res.m[2] = m_2d[0][0] * rhs.m_2d[0][2] + m_2d[0][1] * rhs.m_2d[1][2] + m_2d[0][2] * rhs.m_2d[2][2];

    res.m[3] = m_2d[1][0] * rhs.m_2d[0][0] + m_2d[1][1] * rhs.m_2d[1][0] + m_2d[1][2] * rhs.m_2d[2][0];
    res.m[4] = m_2d[1][0] * rhs.m_2d[0][1] + m_2d[1][1] * rhs.m_2d[1][1] + m_2d[1][2] * rhs.m_2d[2][1];
    res.m[5] = m_2d[1][0] * rhs.m_2d[0][2] + m_2d[1][1] * rhs.m_2d[1][2] + m_2d[1][2] * rhs.m_2d[2][2];

    res.m[6] = m_2d[2][0] * rhs.m_2d[0][0] + m_2d[2][1] * rhs.m_2d[1][0] + m_2d[2][2] * rhs.m_2d[2][0];
    res.m[7] = m_2d[2][0] * rhs.m_2d[0][1] + m_2d[2][1] * rhs.m_2d[1][1] + m_2d[2][2] * rhs.m_2d[2][1];
    res.m[8] = m_2d[2][0] * rhs.m_2d[0][2] + m_2d[2][1] * rhs.m_2d[1][2] + m_2d[2][2] * rhs.m_2d[2][2];
    return res;
  }
  inline Matrix3 operator-() const  {
    Matrix3 res(*this);
    res.m[0] = -res.m[0]; res.m[1] = -res.m[1]; res.m[2] = -res.m[2];
    res.m[3] = -res.m[3]; res.m[4] = -res.m[4]; res.m[5] = -res.m[5];
    res.m[6] = -res.m[6]; res.m[7] = -res.m[7]; res.m[8] = -res.m[8];
    return res;
  }
  inline Vector3 operator*(Vector3 rhs) const {
    return Vector3(
      m_2d[0][0] * rhs[0] + m_2d[0][1] * rhs[1] + m_2d[0][2] * rhs[2],
      m_2d[1][0] * rhs[0] + m_2d[1][1] * rhs[1] + m_2d[1][2] * rhs[2],
      m_2d[2][0] * rhs[0] + m_2d[2][1] * rhs[1] + m_2d[2][2] * rhs[2]);
  }
  inline Matrix3& operator*=(const Matrix3& rhs) { *this = *this * rhs; return *this; }

  //static arithmetic op
  static inline void Multiply(Matrix3& result, const Matrix3& lhs, const Matrix3& rhs) { result = lhs * rhs; }
  static inline void Multiply(Vector3& result, const Matrix3& lhs, Vector3 vec)        { result = lhs * vec; }

  inline void Identity() { memset(m, 0, 9 * sizeof(float)); m_2d[0][0] = 1.0f; m_2d[1][1] = 1.0f; m_2d[2][2] = 1.0f; }
  //advanced member arithmetic op
  float Determinant() const;
  Matrix3& Inverse();
  Matrix3& Transpose();
  Matrix3& Normalise(); //per-vector normalize
  Matrix3 InverseCopy() const   { Matrix3 mat(*this); mat.Inverse(); return mat; }
  Matrix3 TransposeCopy() const { Matrix3 mat(*this); mat.Transpose(); return mat; }
  Matrix3 NormaliseCopy() const { Matrix3 mat(*this); mat.Normalise(); return mat; }

  //geometric properties
  inline void SetTranslate(Vector2 vec) { m_2d[0][2] = vec.x; m_2d[1][2] = vec.y; }
  inline Vector2 GetTranslate() const   { return Vector2(m_2d[0][2], m_2d[1][2]); }
  //warning: will not work with rotation
  inline void SetScale(Vector2 vec) { m_2d[0][0] = vec.x; m_2d[1][1] = vec.y; }
  inline Vector2 GetScale() const   { return Vector2(m_2d[0][0], m_2d[1][1]); }

  //Reset->Scale->Rotate->Translate
  inline void MakeScaleMatrix(Vector2 vec) { memset(m, 0, 9 * sizeof(float)); m_2d[0][0] = vec.x; m_2d[1][1] = vec.y; m_2d[2][2] = 1.0f; }
  inline void MakeTransform(Vector2 translation = Vector2(0, 0), Vector2 scale = Vector2(1.0f, 1.0f), float rotationrad = 0.f)
  {
    float cosineval = cos(rotationrad);
    //cosineval = (abs(cosineval) < ZERO_ERROR_RANGE) ? 0.0f : cosineval;
    float sineval = sin(rotationrad);
    //sineval = (abs(sineval) < ZERO_ERROR_RANGE) ? 0.0f : sineval;
    m_2d[0][0] = scale.x * cosineval;
    m_2d[0][1] = scale.x * -sineval;
    m_2d[1][0] = scale.y * sineval;
    m_2d[1][1] = scale.y * cosineval;

    m_2d[0][2] = translation.x;
    m_2d[1][2] = translation.y;
    m_2d[2][2] = 1.0f;

    m_2d[1][0] = m_2d[2][0] = m_2d[2][1] = 0.0f;
  }

  inline friend std::ostream& operator<<(std::ostream& outstream, const Matrix3& matrix)
  {
    outstream << "(" << matrix.m[0] << ", " << matrix.m[1] << ", " << matrix.m[2] << ") \n";
    outstream << "(" << matrix.m[3] << ", " << matrix.m[4] << ", " << matrix.m[5] << ") \n";
    outstream << "(" << matrix.m[6] << ", " << matrix.m[7] << ", " << matrix.m[8] << ") \n";

    return outstream;
  }

  static Matrix3 GetRotationMatrixToVector(const Vector3& vectortorotate, const Vector3& destvector);
};


/******************************************************************************/
/*!
  \brief
  Quaternion object that helps in rotations
*/
/******************************************************************************/
class Quaternion
{
public:
  float w, x, y, z;

  static const Quaternion ZERO;
  static const Quaternion IDENTITY;
public:
  //constructors
  inline Quaternion(float fw = 1.0f, float fx = 0.0f, float fy = 0.0f, float fz = 0.0f)
    :w(fw), x(fx), y(fy), z(fz) { }
  inline Quaternion(const Quaternion& rhs)
    : w(rhs.w), x(rhs.x), y(rhs.y), z(rhs.z) { }
  
  inline Quaternion(float anglerad, Vector3 axis, bool isradian = true) 
  { SetRotationToAxisRad(isradian ? anglerad : DegToRad(anglerad), axis); }
  inline Quaternion(float* valptr) { memcpy(&w, valptr, sizeof(float) * 4); }
  inline Quaternion(const Matrix3& mat) { GenerateFromRotationMatrix3(mat); }

  //quaternion that rotates directly to the player (with using slerp at 1.0f). See http://gamedev.stackexchange.com/questions/15070/orienting-a-model-to-face-a-target.
  //Also known in unity as quaternion.LookAt
  //up is only used in the case where a and b are exactly in opposite directions (recommended to pick a cross b for a direction)
  inline Quaternion(const Vector3& a_normalised, const Vector3& b_normalised, const Vector3& up) 
  {
    //@todo: fix this, the 2 if clauses are the same
    //http://lolengine.net/blog/2013/09/18/beautiful-maths-quaternion-from-vectors#comment-16
    float dot = a_normalised.DotProduct(b_normalised);
    if (abs(dot + 1.0f) < 0.00001f) //dot product is close to -1.0f
    {
      //a and b point in completely opposite directions (NOT simply facing away)
      SetRotationToAxisRad(PI, up);
    }
    else if (abs(dot + 1.0f) < 0.00001f)
    {
      *this = Quaternion::IDENTITY;
    }
    else
    {
      float angle = acos(dot);
      Vector3 rotaxis = a_normalised.CrossProduct(b_normalised);
      rotaxis.Normalise();
      SetRotationToAxisRad(angle, rotaxis);
    }
  }

  //access op
  inline float operator[](unsigned int i) const { assert(i < 4); return *(&w + i); }
  inline float& operator[](unsigned int i)      { assert(i < 4); return *(&w + i); }
  inline float* ptr()             { return &w; }
  inline const float* ptr() const { return &w; }

  //mutators/accessors
  Matrix3 GetRotationMatrix3() const;
  void GenerateFromRotationMatrix3(const Matrix3& rotationmat);
  inline void SetRotationToAxisRad(float anglerad, Vector3 axisvector) {
    float fHalfAngle = 0.5f * anglerad;
    float v = sin(fHalfAngle);
    w = cos(fHalfAngle);
    x = v * axisvector.x;
    y = v * axisvector.y;
    z = v * axisvector.z;
  }
  inline void SetRotationToAxisDeg(float angledeg, Vector3 axisvector) { SetRotationToAxisRad(DegToRad(angledeg), axisvector); }
  void GetRotationToAxisRad(float& anglerad, Vector3& axisvector);
  void GetRotationToAxisDeg(float& angledeg, Vector3& axisvector);

  inline void Normalise() { 
    float magnitude = sqrt(w * w + x * x + y * y + z * z); 
    w = w / magnitude; x = x / magnitude; y = y / magnitude; z = z / magnitude; 
  }

  //assignment op
  inline Quaternion& operator=(Quaternion rhs) { w = rhs.w; x = rhs.x; y = rhs.y; z = rhs.z; return *this; }

  //comparision op
  inline bool operator==(Quaternion rhs) const { return (rhs.x == x) && (rhs.y == y) && (rhs.z == z) && (rhs.w == w); }
  inline bool operator!=(Quaternion rhs) const { return !(*this == rhs); }

  //basic arithmetic op
  inline Quaternion operator+(Quaternion rhs) const { Quaternion quat; quat.w = w + rhs.w; quat.x = x + rhs.x; quat.y = y + rhs.y; quat.z = z + rhs.z; return quat; }
  inline Quaternion operator-(Quaternion rhs) const { Quaternion quat; quat.w = w - rhs.w; quat.x = x - rhs.x; quat.y = y - rhs.y; quat.z = z - rhs.z; return quat; }
  inline Quaternion operator*(Quaternion rhs) const {
    //what you want to use: combines rotations of 2 Quaternions together. Order of mult matters
    Quaternion quat;
    quat.w = w * rhs.w - x * rhs.x - y * rhs.y - z * rhs.z;
    quat.x = w * rhs.x + x * rhs.w + y * rhs.z - z * rhs.y;
    quat.y = w * rhs.y + y * rhs.w + z * rhs.x - x * rhs.z;
    quat.z = w * rhs.z + z * rhs.w + x * rhs.y - y * rhs.x;
    return quat;
  }
  inline Quaternion operator*(float scalar) const { Quaternion quat; quat.w = w * scalar; quat.x = x * scalar; quat.y = y * scalar; quat.z = z * scalar; return quat; }
  inline Vector3 operator*(Vector3 rhs) const {
    Vector3 uv, uuv;
    Vector3 qvec(x, y, z);
    uv = qvec.CrossProduct(rhs);
    uuv = qvec.CrossProduct(uv);
    uv *= (2.0f * w);
    uuv *= 2.0f;

    Vector3 res = rhs + uv + uuv;
    /*
    res.x = (abs(res.x) < ZERO_ERROR_RANGE) ? 0.0f : res.x;
    res.y = (abs(res.y) < ZERO_ERROR_RANGE) ? 0.0f : res.y;
    res.z = (abs(res.z) < ZERO_ERROR_RANGE) ? 0.0f : res.z;
    */
    return res;
  }
  inline Quaternion operator-() const { Quaternion quat; quat.w = -w; quat.x = -x; quat.y = -y; quat.z = -z; return quat; }

  inline Quaternion& operator+=(Quaternion rhs) { *this = *this + rhs; return *this; }
  inline Quaternion& operator-=(Quaternion rhs) { *this = *this - rhs; return *this; }
  inline Quaternion& operator*=(Quaternion rhs) { *this = *this * rhs; return *this; }
  inline Quaternion& operator*=(float scalar)   { *this = *this * scalar; return *this; }

  //static op
  static inline void Multiply(Quaternion& result, const Quaternion& lhs, const Quaternion& rhs) { result = lhs * rhs; }
  static inline void Multiply(Quaternion& result, const Quaternion& lhs, float scalar)          { result = lhs * scalar; }

  inline float Dot(Quaternion rhs) const { return w * rhs.w + x * rhs.x + y* rhs.y + z* rhs.z; }
  //global friend func
  friend inline Quaternion operator* (float scalar, const Quaternion& rkQ) { Quaternion quat; quat = rkQ * scalar; return quat; }

  //advanced op(INCOMPLETE)
  //Slerp, inverse normalize?
  static Quaternion Slerp(float fT, const Quaternion& rkP, const Quaternion& rkQ, bool shortestpath = false);
  static Quaternion SlerpExtraSpins(float fT, const Quaternion& rkP, const Quaternion& rkQ, int extraspins);

  // spherical quadratic interpolation
  static Quaternion Squad(float fT, const Quaternion& rkP, const Quaternion& rkA, const Quaternion& rkB, const Quaternion& rkQ, bool shortestpath = false);

  // normalised linear interpolation - faster but less accurate (non-constant rotation velocity)
  static Quaternion Nlerp(float fT, const Quaternion& rkP, const Quaternion& rkQ, bool shortestpath = false);

  inline void Conjugate()  { x = -x; y = -y; z = -z; }

  /*WARNING: different coordinate systems for roll, pitch and yaw:
  - x axis points forward(to you)
  - y axis points to the right
  - z axis points downwards
  */
  inline void GetRoll(float& degval)  //! rotation about custom x axis (x pointing towards you)
  { degval = RadToDeg(atan2(2 * (x*y + w * z), w*w + x * x - y * y - z * z)); }
  inline void GetPitch(float& degval) //! rotation about custom y axis (y pointing to right)
  { degval = RadToDeg(atan2(2*(y*z + w*x), w*w - x*x - y*y + z*z)); }
  inline void GetYaw(float& degval)   //! rotation about custom z axis (z pointing downwards)
  { degval = RadToDeg(asin(-2 * (x*z - w*y))); }

  inline friend std::ostream& operator <<(std::ostream& o, Quaternion q)
  {
    o << "Quaternion(" << q.w << ", " << q.x << ", " << q.y << ", " << q.z << ")";
    return o;
  }
};

/******************************************************************************/
/*!
  \brief
  4x4 Matrix. Column major 4x4 matrix. ie. XYZW Translation to the rightmost column
  Top left element has row index 0, column index 0 and bottom right has row index 3, column index 3
*/
/******************************************************************************/
class Matrix4
{
public:
  union
  {
    float m[16];
    float m_2d[4][4]; // [row] [column]
  };

  static const Matrix4 ZERO;
  static const Matrix4 IDENTITY;

public:
  //Constructors
  inline Matrix4() { }
  inline Matrix4(const Matrix4& rhs) { memcpy(m, rhs.m, 16 * sizeof(float)); }
  inline explicit Matrix4(const float arr[4][4]) { memcpy(m, arr, 16 * sizeof(float)); }
  Matrix4(float entry00, float entry01, float entry02, float entry03,
    float entry10, float entry11, float entry12, float entry13,
    float entry20, float entry21, float entry22, float entry23,
    float entry30, float entry31, float entry32, float entry33) // no forced inline
  {
    m_2d[0][0] = entry00; m_2d[0][1] = entry01; m_2d[0][2] = entry02; m_2d[0][3] = entry03;
    m_2d[1][0] = entry10; m_2d[1][1] = entry11; m_2d[1][2] = entry12; m_2d[1][3] = entry13;
    m_2d[2][0] = entry20; m_2d[2][1] = entry21; m_2d[2][2] = entry22; m_2d[2][3] = entry23;
    m_2d[3][0] = entry30; m_2d[3][1] = entry31; m_2d[3][2] = entry32; m_2d[3][3] = entry33;
  }
  inline Matrix4(const Vector4& row0, const Vector4& row1, const Vector4& row2, const Vector4& row3)
  {
    m_2d[0][0] = row0.x; m_2d[0][1] = row0.y; m_2d[0][2] = row0.z; m_2d[0][3] = row0.w;
    m_2d[1][0] = row1.x; m_2d[1][1] = row1.y; m_2d[1][2] = row1.z; m_2d[1][3] = row1.w;
    m_2d[2][0] = row2.x; m_2d[2][1] = row2.y; m_2d[2][2] = row2.z; m_2d[2][3] = row2.w;
    m_2d[3][0] = row3.x; m_2d[3][1] = row3.y; m_2d[3][2] = row3.z; m_2d[3][3] = row3.w;
  }
  inline Matrix4(const Matrix3& rhs)
  {
    m[0] = rhs.m[0]; m[1] = rhs.m[1]; m[2] = rhs.m[2];
    m[4] = rhs.m[3]; m[5] = rhs.m[4]; m[6] = rhs.m[5];
    m[8] = rhs.m[6]; m[9] = rhs.m[7]; m[10] = rhs.m[8];
    m[3] = m[7] = m[11] = m[12] = m[13] = m[14] = 0;
    m[15] = 1;
  }
  inline Matrix4(Quaternion rot)
  {
    Matrix3 m3x3(rot.GetRotationMatrix3());
    operator=(m3x3);
  }
  inline Matrix4(const Quaternion& rot, Vector3 origin)
  {
    Matrix3 m3x3(rot.GetRotationMatrix3());
    operator=(m3x3);
    SetTranslate(origin);
  }

  //access op
  inline float* operator[](unsigned int row) const { assert(row < 4); return (float*)m_2d[row]; }
  inline const float* ptr() const { return m; }
  inline float* ptr()             { return m; }

  //assignment op
  inline Matrix4& operator=(const Matrix4& rhs) { memcpy(m, rhs.m, 16 * sizeof(float)); return *this; }
  inline Matrix4& operator=(const Matrix3& rhs) {
    m[0] = rhs.m[0]; m[1] = rhs.m[1]; m[2] = rhs.m[2];
    m[4] = rhs.m[3]; m[5] = rhs.m[4]; m[6] = rhs.m[5];
    m[8] = rhs.m[6]; m[9] = rhs.m[7]; m[10] = rhs.m[8];
    m[3] = m[7] = m[11] = m[12] = m[13] = m[14] = 0;
    m[15] = 1;
    return *this;
  }

  //mutators and accessors
  inline Vector3 GetColumn(unsigned int column) const {
    assert(column < 4);
    return Vector3(m_2d[0][column], m_2d[1][column], m_2d[2][column]);
  }
  inline Vector4 GetColumn4(unsigned int column) const {
    assert(column < 4);
    return Vector4(m_2d[0][column], m_2d[1][column], m_2d[2][column], m_2d[3][column]);
  }
  inline void SetColumn(const Vector3& vec, unsigned int column) {
    assert(column < 4);
    m_2d[0][column] = vec.x; m_2d[1][column] = vec.y; m_2d[2][column] = vec.z;
  }
  inline Vector3 GetRow(unsigned int row) const {
    assert(row < 4);
    return Vector3(m_2d[row][0], m_2d[row][1], m_2d[row][2]);
  }
  inline Vector4 GetRow4(unsigned int row) const {
    assert(row < 4);
    return Vector4(m_2d[row][0], m_2d[row][1], m_2d[row][2], m_2d[row][3]);
  }
  inline void SetRow(const Vector3& vec, unsigned int row) {
    assert(row < 4);
    m_2d[row][0] = vec.x; m_2d[row][1] = vec.y; m_2d[row][2] = vec.z;
  }
  inline void ResetBottomRow() { m_2d[3][0] = m_2d[3][1] = m_2d[3][2] = 0.f; m_2d[3][3] = 1.f; }
  inline Matrix3 Extract3x3Matrix() const {
    return Matrix3(m_2d[0][0], m_2d[0][1], m_2d[0][2],
      m_2d[1][0], m_2d[1][1], m_2d[1][2],
      m_2d[2][0], m_2d[2][1], m_2d[2][2]);
  }
  inline Quaternion ExtractQuaternion() const { Matrix3 mat(Extract3x3Matrix()); return Quaternion(mat); }

  //comparison op
  inline bool operator==(const Matrix4& rhs) const {
    for (unsigned int i = 0; i<16; ++i)
      if (m[i] != rhs.m[i])
        return false;
    return true;
    //use memcmp?
    //memcmp(
  }
  inline bool operator!=(const Matrix4& rhs) const { return !(*this == rhs); }

  //basic member arithmetic op
  inline Matrix4 operator*(const Matrix4& rhs) const {
    Matrix4 res;
    res.m[0] = m_2d[0][0] * rhs.m_2d[0][0] + m_2d[0][1] * rhs.m_2d[1][0] + m_2d[0][2] * rhs.m_2d[2][0] + m_2d[0][3] * rhs.m_2d[3][0];
    res.m[1] = m_2d[0][0] * rhs.m_2d[0][1] + m_2d[0][1] * rhs.m_2d[1][1] + m_2d[0][2] * rhs.m_2d[2][1] + m_2d[0][3] * rhs.m_2d[3][1];
    res.m[2] = m_2d[0][0] * rhs.m_2d[0][2] + m_2d[0][1] * rhs.m_2d[1][2] + m_2d[0][2] * rhs.m_2d[2][2] + m_2d[0][3] * rhs.m_2d[3][2];
    res.m[3] = m_2d[0][0] * rhs.m_2d[0][3] + m_2d[0][1] * rhs.m_2d[1][3] + m_2d[0][2] * rhs.m_2d[2][3] + m_2d[0][3] * rhs.m_2d[3][3];

    res.m[4] = m_2d[1][0] * rhs.m_2d[0][0] + m_2d[1][1] * rhs.m_2d[1][0] + m_2d[1][2] * rhs.m_2d[2][0] + m_2d[1][3] * rhs.m_2d[3][0];
    res.m[5] = m_2d[1][0] * rhs.m_2d[0][1] + m_2d[1][1] * rhs.m_2d[1][1] + m_2d[1][2] * rhs.m_2d[2][1] + m_2d[1][3] * rhs.m_2d[3][1];
    res.m[6] = m_2d[1][0] * rhs.m_2d[0][2] + m_2d[1][1] * rhs.m_2d[1][2] + m_2d[1][2] * rhs.m_2d[2][2] + m_2d[1][3] * rhs.m_2d[3][2];
    res.m[7] = m_2d[1][0] * rhs.m_2d[0][3] + m_2d[1][1] * rhs.m_2d[1][3] + m_2d[1][2] * rhs.m_2d[2][3] + m_2d[1][3] * rhs.m_2d[3][3];

    res.m[8] = m_2d[2][0] * rhs.m_2d[0][0] + m_2d[2][1] * rhs.m_2d[1][0] + m_2d[2][2] * rhs.m_2d[2][0] + m_2d[2][3] * rhs.m_2d[3][0];
    res.m[9] = m_2d[2][0] * rhs.m_2d[0][1] + m_2d[2][1] * rhs.m_2d[1][1] + m_2d[2][2] * rhs.m_2d[2][1] + m_2d[2][3] * rhs.m_2d[3][1];
    res.m[10] = m_2d[2][0] * rhs.m_2d[0][2] + m_2d[2][1] * rhs.m_2d[1][2] + m_2d[2][2] * rhs.m_2d[2][2] + m_2d[2][3] * rhs.m_2d[3][2];
    res.m[11] = m_2d[2][0] * rhs.m_2d[0][3] + m_2d[2][1] * rhs.m_2d[1][3] + m_2d[2][2] * rhs.m_2d[2][3] + m_2d[2][3] * rhs.m_2d[3][3];

    res.m[12] = m_2d[3][0] * rhs.m_2d[0][0] + m_2d[3][1] * rhs.m_2d[1][0] + m_2d[3][2] * rhs.m_2d[2][0] + m_2d[3][3] * rhs.m_2d[3][0];
    res.m[13] = m_2d[3][0] * rhs.m_2d[0][1] + m_2d[3][1] * rhs.m_2d[1][1] + m_2d[3][2] * rhs.m_2d[2][1] + m_2d[3][3] * rhs.m_2d[3][1];
    res.m[14] = m_2d[3][0] * rhs.m_2d[0][2] + m_2d[3][1] * rhs.m_2d[1][2] + m_2d[3][2] * rhs.m_2d[2][2] + m_2d[3][3] * rhs.m_2d[3][2];
    res.m[15] = m_2d[3][0] * rhs.m_2d[0][3] + m_2d[3][1] * rhs.m_2d[1][3] + m_2d[3][2] * rhs.m_2d[2][3] + m_2d[3][3] * rhs.m_2d[3][3];
    return res;
  }
  inline Matrix4 operator-() const {
    Matrix4 res(*this);
    res.m[0] = -res.m[0]; res.m[1] = -res.m[1]; res.m[2] = -res.m[2]; res.m[3] = -res.m[3];
    res.m[4] = -res.m[4]; res.m[5] = -res.m[5]; res.m[6] = -res.m[6]; res.m[7] = -res.m[7];
    res.m[8] = -res.m[8]; res.m[9] = -res.m[9]; res.m[10] = -res.m[10]; res.m[11] = -res.m[11];
    res.m[12] = -res.m[12]; res.m[13] = -res.m[13]; res.m[14] = -res.m[14]; res.m[15] = -res.m[15];
    return res;
  }
  inline Vector4 operator*(Vector4 rhs) const {
    return Vector4(
      m_2d[0][0] * rhs.x + m_2d[0][1] * rhs.y + m_2d[0][2] * rhs.z + m_2d[0][3] * rhs.w,
      m_2d[1][0] * rhs.x + m_2d[1][1] * rhs.y + m_2d[1][2] * rhs.z + m_2d[1][3] * rhs.w,
      m_2d[2][0] * rhs.x + m_2d[2][1] * rhs.y + m_2d[2][2] * rhs.z + m_2d[2][3] * rhs.w,
      m_2d[3][0] * rhs.x + m_2d[3][1] * rhs.y + m_2d[3][2] * rhs.z + m_2d[3][3] * rhs.w);
  }

  inline Matrix4& operator*=(const Matrix4& rhs) { *this = *this * rhs; return *this; }
  
  //static arithmetic op
  static inline void Multiply(Matrix4& result, const Matrix4& lhs, const Matrix4& rhs) { result = lhs * rhs; }
  static inline void Multiply(Vector4& result, const Matrix4& lhs, Vector4 vec)        { result = lhs * vec; }

  //advanced member arithmetic op
  float Determinant() const;
  Matrix4& Inverse();
  Matrix4& Transpose();
  Matrix4 InverseCopy() const { Matrix4 mat(*this); mat.Inverse(); return mat; }
  Matrix4 TransposeCopy() const      { Matrix4 mat(*this); mat.Transpose(); return mat; }

  //geometric properties
  inline void Identity()                       { memset(m, 0, 16 * sizeof(float)); m_2d[0][0] = m_2d[1][1] = m_2d[2][2] = m_2d[3][3] = 1.f; }
  inline void SetTranslate(Vector3 vec)        { m_2d[0][3] = vec.x; m_2d[1][3] = vec.y; m_2d[2][3] = vec.z; }
  inline Vector3 GetTranslate() const          { return Vector3(m_2d[0][3], m_2d[1][3], m_2d[2][3]); }
  inline void MakeTranslateMatrix(Vector3 vec) { Identity(); m_2d[0][3] = vec.x; m_2d[1][3] = vec.y; m_2d[2][3] = vec.z; m_2d[3][3] = 1.f; }
  
  inline void SetScale(Vector3 vec)        { m_2d[0][0] = vec.x; m_2d[1][1] = vec.y; m_2d[2][2] = vec.z; } //warning: will not work with rotation
  inline Vector3 GetScale() const          { return Vector3(m_2d[0][0], m_2d[1][1], m_2d[2][2]); }
  inline void MakeScaleMatrix(Vector3 vec) { memset(m, 0, 16 * sizeof(float)); m_2d[0][0] = vec.x; m_2d[1][1] = vec.y; m_2d[2][2] = vec.z; m_2d[3][3] = 1.f; }

  //Reset->Scale->Rotate->Translate
  inline void MakeTransform(const Vector3& translation = Vector3(0, 0, 0), const Vector3& scale = Vector3(1.0f, 1.0f, 1.0f), const Quaternion& orientation = Quaternion())
  {
    Matrix4 rot4(orientation.GetRotationMatrix3());;
    Matrix4 scale4;
    scale4.MakeScaleMatrix(scale);

    Matrix4 t;
    t.MakeTranslateMatrix(translation);
    *this = t * rot4 * scale4;
  }

  //Camera Matrices
  void MakeWorldToViewMatrixRH(const Vector3& position, const Vector3& target, const Vector3& up_vec = Vector3(0, 1, 0)); //aka the lookat matrix
  void MakeViewToWorldMatrixRH(const Vector3& position, const Vector3& target, const Vector3& up_vec = Vector3(0, 1, 0));
  inline void MakePerpectiveMatrixRH0To1Depth(float fovdeg, float aspectratio, float neardist, float fardist)
  {
    memset(m_2d, 0, sizeof(float) * 16);

    m_2d[1][1] = 1.0f / tan(DegToRad(fovdeg) * 0.5f);
    m_2d[0][0] = m_2d[1][1] / aspectratio;

    //map to 0 to 1 range
    m_2d[2][2] = (fardist) / (neardist - fardist);
    m_2d[2][3] = (neardist * fardist) / (neardist - fardist);

    m_2d[3][2] = -1.0f;
  }
  inline void MakePerpectiveMatrixRHNeg1To1Depth(float fovdeg, float aspectratio, float neardist, float fardist) 
  {
    //needs testing
    memset(m_2d, 0, sizeof(float) * 16);

    m_2d[1][1] = 1.0f / tan(DegToRad(fovdeg) * 0.5f);
    m_2d[0][0] = m_2d[1][1] / aspectratio;

    //map to -1 to 1 range; based on gluPerspectiveMatrix
    m_2d[2][2] = (fardist + neardist) / (neardist - fardist);
    m_2d[2][3] = (2.0f * fardist * neardist) / (neardist - fardist);

    m_2d[3][2] = -1.0f;
  }

  //(left, bottom) lower left of nearplane, (right, top) upper right of nearplane
  inline void MakeOrthogonalMatrixRH0To1Depth(float left, float right, float bottom, float top, float neardist, float fardist)
  {
    memset(m_2d, 0, sizeof(float) * 16);

    float rlSum = static_cast<float>(right + left);
    float rlDif = static_cast<float>(right - left);
    float tbSum = static_cast<float>(top + bottom);
    float tbDif = static_cast<float>(top - bottom);
    //float fnSum = static_cast<float>(fardist + neardist);
    float fnDif = static_cast<float>(fardist - neardist);

    //from http://msdn.microsoft.com/en-us/library/windows/desktop/bb205349%28v=vs.85%29.aspx, centered
    //from https://msdn.microsoft.com/en-us/library/windows/desktop/bb205348%28v=vs.85%29.aspx, do not use D3DXMatrixOrthoRH since it does not off center

    m_2d[0][0] = 2.0f / rlDif;
    m_2d[1][1] = 2.0f / tbDif;
    m_2d[2][2] = 1.0f / -fnDif;
    m_2d[3][3] = 1.0f;

    m_2d[0][3] = -rlSum / rlDif;
    m_2d[1][3] = -tbSum / tbDif;
    m_2d[2][3] = neardist / -fnDif;
  }

  inline void MakeOrthogonalMatrixRHNeg1To1Depth(float left, float right, float bottom, float top, float neardist, float fardist) {
    memset(m_2d, 0, sizeof(float) * 16);

    float rlSum = static_cast<float>(right + left);
    float rlDif = static_cast<float>(right - left);
    float tbSum = static_cast<float>(top + bottom);
    float tbDif = static_cast<float>(top - bottom);
    float fnSum = static_cast<float>(fardist + neardist);
    float fnDif = static_cast<float>(fardist - neardist);

    //from http://www.songho.ca/opengl/gl_projectionmatrix.html
    //needs testing, need -ve sign for [2][2] and [2][3]?
    m_2d[0][0] = 2.0f / rlDif;
    m_2d[1][1] = 2.0f / tbDif;
    m_2d[2][2] = 2.0f / fnDif;
    m_2d[3][3] = 1.0f;

    m_2d[0][3] = -rlSum / rlDif;
    m_2d[1][3] = -tbSum / tbDif;
    m_2d[2][3] = fnSum / fnDif;
  }


  inline friend std::ostream& operator<<(std::ostream& outstream, const Matrix4& matrix)
  {
    outstream << "(" << matrix.m[0] << ", " << matrix.m[1] << ", " << matrix.m[2] << ", " << matrix.m[3] << ") \n";
    outstream << "(" << matrix.m[4] << ", " << matrix.m[5] << ", " << matrix.m[6] << ", " << matrix.m[7] << ") \n";
    outstream << "(" << matrix.m[8] << ", " << matrix.m[9] << ", " << matrix.m[10] << ", " << matrix.m[11] << ") \n";
    outstream << "(" << matrix.m[12] << ", " << matrix.m[13] << ", " << matrix.m[14] << ", " << matrix.m[15] << ") \n";

    return outstream;
  }
};
