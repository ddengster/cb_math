
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

