
/*
MatrixMath.cpp - v1.0

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

#include "MatrixMath.h"

/******************************************************************************/

const Matrix2 Matrix2::ZERO(0.0f, 0.0f, 0.0f, 0.0f);
const Matrix2 Matrix2::IDENTITY(1.0f, 0.0f, 0.0f, 1.0f);

void Matrix2::Inverse()
{
  // store a temp matrix
  Matrix2 mtx(*this);

  // calculate the determinant of the matrix
  float det = mtx.m_2d[0][0] *mtx.m_2d[1][1] - mtx.m_2d[0][1] *mtx.m_2d[1][0];
  if (abs(det) < 0.001f)
  {
    Identity();
    return;
  }
  det = 1 / det;

  // part that remain the same
  m_2d[0][0] =  mtx.m_2d[1][1] * det;
  m_2d[1][1] =  mtx.m_2d[0][0] * det;
  m_2d[0][1] = -mtx.m_2d[1][0] * det;
  m_2d[1][0] = -mtx.m_2d[0][1] * det;
}

/******************************************************************************/

const Matrix3 Matrix3::ZERO(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
const Matrix3 Matrix3::IDENTITY(1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f);

float Matrix3::Determinant() const
{
  float factor1 = m_2d[1][1] * m_2d[2][2] - m_2d[1][2] * m_2d[2][1];
  float factor2 = m_2d[1][2] * m_2d[2][0] - m_2d[1][0] * m_2d[2][2];
  float factor3 = m_2d[1][0] * m_2d[2][1] - m_2d[1][1] * m_2d[2][0];

  float det = m[0] * factor1 + m[1] * factor2 + m[2] * factor3;

  return det;
}

Matrix3& Matrix3::Inverse()
{
  float det = Determinant();
  assert(abs(det) > ZERO_ERROR_RANGE);

  //transpose already included
  m_2d[0][0] = m_2d[1][1] * m_2d[2][2] - m_2d[1][2] * m_2d[2][1];
  m_2d[0][1] = m_2d[0][2] * m_2d[2][1] - m_2d[0][1] * m_2d[2][2];
  m_2d[0][2] = m_2d[0][1] * m_2d[1][2] - m_2d[0][2] * m_2d[1][1];

  m_2d[1][0] = m_2d[1][2] * m_2d[2][0] - m_2d[1][0] * m_2d[2][2];
  m_2d[1][1] = m_2d[0][0] * m_2d[2][2] - m_2d[0][2] * m_2d[2][0];
  m_2d[1][2] = m_2d[0][2] * m_2d[1][0] - m_2d[0][0] * m_2d[1][2];

  m_2d[2][0] = m_2d[1][0] * m_2d[2][1] - m_2d[1][1] * m_2d[2][0];
  m_2d[2][1] = m_2d[0][1] * m_2d[2][0] - m_2d[0][0] * m_2d[2][1];
  m_2d[2][2] = m_2d[0][0] * m_2d[1][1] - m_2d[0][1] * m_2d[1][0];

  float inv = 1.0f / det;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
      m_2d[i][j] *= inv;
  }
  return *this;
}

Matrix3& Matrix3::Normalise()
{
  Vector3 r0 = GetColumn(0);
  Vector3 r1 = GetColumn(1);
  Vector3 r2 = GetColumn(2);
  r0.Normalise();
  r1.Normalise();
  r2.Normalise();
  *this = Matrix3(r0, r1, r2);
  return *this;
}

Matrix3 Matrix3::GetRotationMatrixToVector(const Vector3& vectortorotate, const Vector3& destvector)
{
  Vector3 p = vectortorotate.NormalisedCopy();
  Vector3 q = destvector.NormalisedCopy();
  Vector3 u = p.CrossProduct(q);
  u.Normalise();
  float pdotq = p.DotProduct(q);
  float pdotq_sq = pdotq * pdotq;

  Matrix3 mat;
  mat.m_2d[0][0] = (u.x * u.x / (1 - pdotq_sq)) + (pdotq * (1 - (u.x * u.x / (1 - pdotq_sq))));
  mat.m_2d[0][1] = (u.x * u.y * (1 - pdotq) / (1 - pdotq_sq)) - u.z;
  mat.m_2d[0][2] = (u.x * u.z * (1 - pdotq) / (1 - pdotq_sq)) + u.y;

  mat.m_2d[1][0] = (u.x * u.y * (1 - pdotq) / (1 - pdotq_sq)) + u.z;
  mat.m_2d[1][1] = (u.y * u.y / (1 - pdotq_sq)) + (pdotq * (1 - (u.y * u.y / (1 - pdotq_sq))));
  mat.m_2d[1][2] = (u.y * u.z * (1 - pdotq) / (1 - pdotq_sq)) - u.x;

  mat.m_2d[2][0] = (u.x * u.z * (1 - pdotq) / (1 - pdotq_sq)) - u.y;
  mat.m_2d[2][1] = (u.y * u.z * (1 - pdotq) / (1 - pdotq_sq)) + u.x;
  mat.m_2d[2][2] = (u.z * u.z / (1 - pdotq_sq)) + (pdotq * (1 - (u.z * u.z / (1 - pdotq_sq))));

  return mat;
}


/******************************************************************************/

const Quaternion Quaternion::ZERO(0.0f, 0.0f, 0.0f, 0.0);
const Quaternion Quaternion::IDENTITY(1.0f, 0.0f, 0.0f, 0.0f);

Matrix3 Quaternion::GetRotationMatrix3() const
{
  float fTx = 2.0f * x;
  float fTy = 2.0f * y;
  float fTz = 2.0f * z;
  float fTwx = fTx  * w;
  float fTwy = fTy  * w;
  float fTwz = fTz  * w;
  float fTxx = fTx  * x;
  float fTxy = fTy  * x;
  float fTxz = fTz  * x;
  float fTyy = fTy  * y;
  float fTyz = fTz  * y;
  float fTzz = fTz  * z;

  Matrix3 rotationmat;
  rotationmat[0][0] = 1.0f - (fTyy + fTzz);
  rotationmat[0][1] = fTxy - fTwz;
  rotationmat[0][2] = fTxz + fTwy;
  rotationmat[1][0] = fTxy + fTwz;
  rotationmat[1][1] = 1.0f - (fTxx + fTzz);
  rotationmat[1][2] = fTyz - fTwx;
  rotationmat[2][0] = fTxz - fTwy;
  rotationmat[2][1] = fTyz + fTwx;
  rotationmat[2][2] = 1.0f - (fTxx + fTyy);

  return rotationmat;
}

void Quaternion::GenerateFromRotationMatrix3(const Matrix3& rotationmat)
{
  // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
  // article "Quaternion Calculus and Fast Animation".

  float trace = rotationmat[0][0] + rotationmat[1][1] + rotationmat[2][2];
  float root;

  if (trace > 0.0)
  {
    // |w| > 1/2, may as well choose w > 1/2
    root = sqrt(trace + 1.0f);  // 2w
    w = 0.5f * root;
    root = 0.5f / root;  // 1/(4w)
    x = (rotationmat[2][1] - rotationmat[1][2]) * root;
    y = (rotationmat[0][2] - rotationmat[2][0]) * root;
    z = (rotationmat[1][0] - rotationmat[0][1]) * root;
  }
  else
  {
    // |w| <= 1/2
    static const unsigned int s_iNext[3] = { 1, 2, 0 };
    unsigned int i = 0;
    if (rotationmat[1][1] > rotationmat[0][0])
      i = 1;
    if (rotationmat[2][2] > rotationmat[i][i])
      i = 2;
    unsigned int j = s_iNext[i];
    unsigned int k = s_iNext[j];

    root = sqrt(rotationmat[i][i] - rotationmat[j][j] - rotationmat[k][k] + 1.0f);
    float* apkQuat[3] = { &x, &y, &z };
    *apkQuat[i] = 0.5f * root;
    root = 0.5f / root;
    w = (rotationmat[k][j] - rotationmat[j][k])*root;
    *apkQuat[j] = (rotationmat[j][i] + rotationmat[i][j]) * root;
    *apkQuat[k] = (rotationmat[k][i] + rotationmat[i][k]) * root;
  }
}

void Quaternion::GetRotationToAxisRad(float& anglerad, Vector3& axisvector)
{
  float sqlength = x*x + y*y + z*z;
  if (sqlength > 0.0f)
  {
    anglerad = 2.0f * acos(w);
    float invlength = 1.f / sqrt(sqlength);
    axisvector.x = x * invlength;
    axisvector.y = y * invlength;
    axisvector.z = z * invlength;
  }
  else
  {
    // angle is 0 (mod 2*pi), so any axis will do
    anglerad = 0.0f;
    axisvector.x = 1.0f;
    axisvector.y = 0.0f;
    axisvector.z = 0.0f;
  }
}

void Quaternion::GetRotationToAxisDeg(float& angledeg, Vector3& axisvector)
{
  float sqlength = x*x + y*y + z*z;
  if (sqlength > 0.0f)
  {
    float rad = 2.0f * acos(w);
    angledeg = RadToDeg(rad);
    float invlength = 1.f / sqrt(sqlength);
    axisvector.x = x * invlength;
    axisvector.y = y * invlength;
    axisvector.z = z * invlength;
  }
  else
  {
    // angle is 0 (mod 2*pi), so any axis will do
    angledeg = 0.0f;
    axisvector.x = 1.0f;
    axisvector.y = 0.0f;
    axisvector.z = 0.0f;
  }
}

Quaternion Quaternion::Slerp(float fT, const Quaternion& rkP, const Quaternion& rkQ, bool shortestpath)
{
#if 1
  float dotp = rkP.Dot(rkQ);
  Quaternion rkT;

  // Do we need to invert rotation?
  if (dotp < 0.0f && shortestpath)
  {
    dotp = -dotp;
    rkT = -rkQ;
  }
  else
  {
    rkT = rkQ;
  }

  const float epsilon = 0.0001f;
  if (abs(dotp) < (1.0f - epsilon))
  {
    // Standard case (slerp)
#if 1
    float sn = sqrt(1.f - dotp * dotp);
    float angle = atan2(sn, dotp);
    float invsin = 1.0f / sn;
    float coeff0 = sin((1.0f - fT) * angle) * invsin;
    float coeff1 = sin(fT * angle) * invsin;
    return coeff0 * rkP + coeff1 * rkT;
#else
    float omega, sinom;
    omega = acos(dotp); // extract theta from dot product's cos theta
    sinom = sin(omega);
    float sclp = sin((1.0f - fT) * omega) / sinom;
    float sclq = sin(fT * omega) / sinom;

    Quaternion q = sclp * rkP + sclq * rkT;
    return q;
#endif
  }
  else
  {
    // There are two situations:
    // 1. "rkP" and "rkQ" are very close (dotp ~= +1), so we can do a linear
    //    interpolation safely.
    // 2. "rkP" and "rkQ" are almost inverse of each other (dotp ~= -1), there
    //    are an infinite number of possibilities interpolation. but we haven't
    //    have method to fix this case, so just use linear interpolation here.
    Quaternion t = (1.0f - fT) * rkP + fT * rkT;
    // taking the complement requires renormalisation
    //t.Normalise();
    return t;
  }
#endif
}

Quaternion Quaternion::SlerpExtraSpins(float fT, const Quaternion& rkP, const Quaternion& rkQ, int extraspins)
{
  float dotp = rkP.Dot(rkQ);
  float angle = acos(dotp);

  const float epsilon = 0.0001f;
  if (abs(angle) < epsilon)
    return rkP;

  float sn = sin(angle);
  float phase = PI * extraspins * fT;
  float invsin = 1.0f / sn;
  float coeff0 = sin((1.0f - fT) * angle - phase) * invsin;
  float coeff1 = sin(fT * angle + phase) * invsin;
  return coeff0 * rkP + coeff1 * rkQ;
}

Quaternion Quaternion::Squad(float fT, const Quaternion& rkP, const Quaternion& rkA, const Quaternion& rkB, const Quaternion& rkQ, bool shortestpath)
{
  float slerpT = 2.0f * fT * (1.0f - fT);
  Quaternion kSlerpP = Slerp(fT, rkP, rkQ, shortestpath);
  Quaternion kSlerpQ = Slerp(fT, rkA, rkB);
  return Slerp(slerpT, kSlerpP, kSlerpQ);
}

Quaternion Quaternion::Nlerp(float fT, const Quaternion& rkP, const Quaternion& rkQ, bool shortestpath)
{
  Quaternion result;
  float dotp = rkP.Dot(rkQ);
  if (dotp < 0.0f && shortestpath)
  {
    result = rkP + fT * ((-rkQ) - rkP);
  }
  else
  {
    result = rkP + fT * (rkQ - rkP);
  }
  result.Normalise();
  return result;
}

/******************************************************************************/

const Matrix4 Matrix4::ZERO(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
const Matrix4 Matrix4::IDENTITY(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);

inline static float MINOR(const Matrix4& m, unsigned int r0, unsigned int r1, unsigned int r2,
  unsigned int c0, unsigned int c1, unsigned int c2)
{
  return	m.m_2d[r0][c0] * (m.m_2d[r1][c1] * m.m_2d[r2][c2] - m.m_2d[r2][c1] * m.m_2d[r1][c2]) -
    m.m_2d[r0][c1] * (m.m_2d[r1][c0] * m.m_2d[r2][c2] - m.m_2d[r2][c0] * m.m_2d[r1][c2]) +
    m.m_2d[r0][c2] * (m.m_2d[r1][c0] * m.m_2d[r2][c1] - m.m_2d[r2][c0] * m.m_2d[r1][c1]);
}

float Matrix4::Determinant() const
{
  return	m_2d[0][0] * MINOR(*this, 1, 2, 3, 1, 2, 3) -
    m_2d[0][1] * MINOR(*this, 1, 2, 3, 0, 2, 3) +
    m_2d[0][2] * MINOR(*this, 1, 2, 3, 0, 1, 3) -
    m_2d[0][3] * MINOR(*this, 1, 2, 3, 0, 1, 2);
}

Matrix4& Matrix4::Inverse()
{
  float m00 = m_2d[0][0], m01 = m_2d[0][1], m02 = m_2d[0][2], m03 = m_2d[0][3];
  float m10 = m_2d[1][0], m11 = m_2d[1][1], m12 = m_2d[1][2], m13 = m_2d[1][3];
  float m20 = m_2d[2][0], m21 = m_2d[2][1], m22 = m_2d[2][2], m23 = m_2d[2][3];
  float m30 = m_2d[3][0], m31 = m_2d[3][1], m32 = m_2d[3][2], m33 = m_2d[3][3];

  float v0 = m20 * m31 - m21 * m30;
  float v1 = m20 * m32 - m22 * m30;
  float v2 = m20 * m33 - m23 * m30;
  float v3 = m21 * m32 - m22 * m31;
  float v4 = m21 * m33 - m23 * m31;
  float v5 = m22 * m33 - m23 * m32;

  float t00 = +(v5 * m11 - v4 * m12 + v3 * m13);
  float t10 = -(v5 * m10 - v2 * m12 + v1 * m13);
  float t20 = +(v4 * m10 - v2 * m11 + v0 * m13);
  float t30 = -(v3 * m10 - v1 * m11 + v0 * m12);

  float invDet = 1 / (t00 * m00 + t10 * m01 + t20 * m02 + t30 * m03);

  float d00 = t00 * invDet;
  float d10 = t10 * invDet;
  float d20 = t20 * invDet;
  float d30 = t30 * invDet;

  float d01 = -(v5 * m01 - v4 * m02 + v3 * m03) * invDet;
  float d11 = +(v5 * m00 - v2 * m02 + v1 * m03) * invDet;
  float d21 = -(v4 * m00 - v2 * m01 + v0 * m03) * invDet;
  float d31 = +(v3 * m00 - v1 * m01 + v0 * m02) * invDet;

  v0 = m10 * m31 - m11 * m30;
  v1 = m10 * m32 - m12 * m30;
  v2 = m10 * m33 - m13 * m30;
  v3 = m11 * m32 - m12 * m31;
  v4 = m11 * m33 - m13 * m31;
  v5 = m12 * m33 - m13 * m32;

  float d02 = +(v5 * m01 - v4 * m02 + v3 * m03) * invDet;
  float d12 = -(v5 * m00 - v2 * m02 + v1 * m03) * invDet;
  float d22 = +(v4 * m00 - v2 * m01 + v0 * m03) * invDet;
  float d32 = -(v3 * m00 - v1 * m01 + v0 * m02) * invDet;

  v0 = m21 * m10 - m20 * m11;
  v1 = m22 * m10 - m20 * m12;
  v2 = m23 * m10 - m20 * m13;
  v3 = m22 * m11 - m21 * m12;
  v4 = m23 * m11 - m21 * m13;
  v5 = m23 * m12 - m22 * m13;

  float d03 = -(v5 * m01 - v4 * m02 + v3 * m03) * invDet;
  float d13 = +(v5 * m00 - v2 * m02 + v1 * m03) * invDet;
  float d23 = -(v4 * m00 - v2 * m01 + v0 * m03) * invDet;
  float d33 = +(v3 * m00 - v1 * m01 + v0 * m02) * invDet;

  m[0] = d00; m[1] = d01; m[2] = d02; m[3] = d03;
  m[4] = d10; m[5] = d11; m[6] = d12; m[7] = d13;
  m[8] = d20; m[9] = d21; m[10] = d22; m[11] = d23;
  m[12] = d30; m[13] = d31; m[14] = d32; m[15] = d33;

  return *this;
}

Matrix4& Matrix4::Transpose()
{
  float temp;
  temp = m_2d[0][1]; m_2d[0][1] = m_2d[1][0]; m_2d[1][0] = temp;
  temp = m_2d[0][2]; m_2d[0][2] = m_2d[2][0]; m_2d[2][0] = temp;
  temp = m_2d[1][2]; m_2d[1][2] = m_2d[2][1]; m_2d[2][1] = temp;
  temp = m_2d[0][3]; m_2d[0][3] = m_2d[3][0]; m_2d[3][0] = temp;
  temp = m_2d[1][3]; m_2d[1][3] = m_2d[3][1]; m_2d[3][1] = temp;
  temp = m_2d[2][3]; m_2d[2][3] = m_2d[3][2]; m_2d[3][2] = temp;

  return *this;
}

void Matrix4::MakeWorldToViewMatrixRH(const Vector3& position, const Vector3& target, const Vector3& up_vec)
{
  Vector3 lookat = target - position;
  lookat.Normalise();
  Vector3 up = up_vec.NormalisedCopy();

  Vector3 side = (lookat.CrossProduct(up)).NormalisedCopy();
  up = side.CrossProduct(lookat); //up must be perpendicular to other 2 vectors
  up.Normalise();

  //right handed, lookat is -ve!
  m_2d[0][0] = side.x;     m_2d[0][1] = side.y;      m_2d[0][2] = side.z;      m_2d[0][3] = -side.DotProduct(position);
  m_2d[1][0] = up.x;       m_2d[1][1] = up.y;        m_2d[1][2] = up.z;        m_2d[1][3] = -up.DotProduct(position);
  m_2d[2][0] = -lookat.x;  m_2d[2][1] = -lookat.y;   m_2d[2][2] = -lookat.z;   m_2d[2][3] = lookat.DotProduct(position);
  m_2d[3][0] = 0.0f;       m_2d[3][1] = 0.0f;        m_2d[3][2] = 0.0f;        m_2d[3][3] = 1.0f;
  /*
  //internet version, RH, produces the same results. note the cross product order
  Vector3 lookat = mCamPos - mCamTarget;
  lookat.Normalise();
  Vector3 up = mCamUp.NormalisedCopy();

  Vector3 side = (up.CrossProduct(lookat)).NormalisedCopy();
  up = lookat.CrossProduct(side); //up must be perpendicular to other 2 vectors
  up.Normalise();

  //right handed, lookat is -ve!
  ret.m_2d[0][0] = side.x;     ret.m_2d[0][1] = side.y;      ret.m_2d[0][2] = side.z;      ret.m_2d[0][3] = -side.DotProduct(mCamPos);
  ret.m_2d[1][0] = up.x;       ret.m_2d[1][1] = up.y;        ret.m_2d[1][2] = up.z;        ret.m_2d[1][3] = -up.DotProduct(mCamPos);
  ret.m_2d[2][0] = lookat.x;  ret.m_2d[2][1] = lookat.y;   ret.m_2d[2][2] = lookat.z;   ret.m_2d[2][3] = -lookat.DotProduct(mCamPos);
  ret.m_2d[3][0] = 0.0f;       ret.m_2d[3][1] = 0.0f;        ret.m_2d[3][2] = 0.0f;        ret.m_2d[3][3] = 1.0f;
  */
}

void Matrix4::MakeViewToWorldMatrixRH(const Vector3& position, const Vector3& target, const Vector3& up_vec)
{
  Vector3 lookat = target - position;
  lookat.Normalise();
  Vector3 up = up_vec.NormalisedCopy();

  Vector3 side = (lookat.CrossProduct(up)).NormalisedCopy();
  up = side.CrossProduct(lookat); //up must be perpendicular to other 2 vectors
  up.Normalise();

  //right handed, lookat is -ve!
  m_2d[0][0] = side.x;   m_2d[0][1] = up.x;   m_2d[0][2] = -lookat.x;   m_2d[0][3] = position.x;
  m_2d[1][0] = side.y;   m_2d[1][1] = up.y;   m_2d[1][2] = -lookat.y;   m_2d[1][3] = position.y;
  m_2d[2][0] = side.z;   m_2d[2][1] = up.z;   m_2d[2][2] = -lookat.z;   m_2d[2][3] = position.z;
  m_2d[3][0] = 0.0f;     m_2d[3][1] = 0.0f;   m_2d[3][2] = 0.0f;        m_2d[3][3] = 1.0f;
}
