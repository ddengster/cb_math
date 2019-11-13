
/*
VectorMath.cpp - v1.0

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

#include "VectorMath.h"

/**************************************************************/

const Vector2 Vector2::ZERO(0.0f, 0.0f);
const Vector2 Vector2::UNIT_X(1.0f, 0.0f);
const Vector2 Vector2::UNIT_Y(0.0f, 1.0f);
const Vector2 Vector2::NEGATIVE_UNIT_X(-1.0f, 0.0f);
const Vector2 Vector2::NEGATIVE_UNIT_Y(0.0f, -1.0f);
const Vector2 Vector2::UNIT_SCALE(1, 1);

/**************************************************************/

const Vector3 Vector3::ZERO(0.0f, 0.0f, 0.0f);
const Vector3 Vector3::UNIT_X(1.0f, 0.0f, 0.0f);
const Vector3 Vector3::UNIT_Y(0.0f, 1.0f, 0.0f);
const Vector3 Vector3::UNIT_Z(0.0f, 0.0f, 1.0f);
const Vector3 Vector3::NEGATIVE_UNIT_X(-1.0f, 0.0f, 0.0f);
const Vector3 Vector3::NEGATIVE_UNIT_Y(0.0f, -1.0f, 0.0f);
const Vector3 Vector3::NEGATIVE_UNIT_Z(0.0f, 0.0f, -1.0f);
const Vector3 Vector3::UNIT_SCALE(1, 1, 1);

void Vector3::GetRotationToVector3(float& deg, Vector3& rotationaxis, const Vector3& othervector) const
{
  Vector3 vec = this->NormalisedCopy();
  Vector3 othervec = othervector.NormalisedCopy();

  //no rotation to 0 vector
  if (vec.IsEqual(Vector3::ZERO, 0.00005f) || othervec.IsEqual(Vector3::ZERO, 0.00005f))
  {
    deg = 0.f;
    rotationaxis = Vector3::ZERO;
    return;
  }

  rotationaxis = vec.CrossProduct(othervec);

  float radians = acos(vec.DotProduct(othervec));
  deg = RadToDeg(radians);
}

/**************************************************************/

const Vector4 Vector4::ZERO(0.0f, 0.0f, 0.0f, 0.0f);
const Vector4 Vector4::UNIT_X(1.0f, 0.0f, 0.0f, 0.0f);
const Vector4 Vector4::UNIT_Y(0.0f, 1.0f, 0.0f, 0.0f);
const Vector4 Vector4::UNIT_Z(0.0f, 0.0f, 1.0f, 0.0f);
const Vector4 Vector4::UNIT_W(0.0f, 0.0f, 0.0f, 1.0f);
const Vector4 Vector4::NEGATIVE_UNIT_X(-1.0f, 0.0f, 0.0f, 0.0f);
const Vector4 Vector4::NEGATIVE_UNIT_Y(0.0f, -1.0f, 0.0f, 0.0f);
const Vector4 Vector4::NEGATIVE_UNIT_Z(0.0f, 0.0f, -1.0f, 0.0f);
const Vector4 Vector4::NEGATIVE_UNIT_W(0.0f, 0.0f, 0.0f, -1.0f);
const Vector4 Vector4::UNIT_SCALE(1, 1, 1, 1);

/**************************************************************/

