
/*
main.cpp - v1.0
Contains a bunch of tests.

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
#include "MatrixMath.h"
#include <iostream>

int main()
{
  using namespace std;

  //vec2
  {
    Vector2 v(1, 0);
    cout << "v: " << v << endl;

    Vector2 perp = v.Perpendicular();
    cout << "v2perp: " << perp << endl;
    assert(perp.x == 0.f);
    assert(perp.y == 1.f);
  }

  //vec3
  {
    Vector3 v(1, 0, 0);
    cout << "v: " << v << endl;

    Vector3 v2 = v + Vector3(0, 1, 1);

    Vector3 cp = v.CrossProduct(v2);
    cout << "cp: " << cp << endl;
    assert(cp.x == 0.f);
    assert(cp.y == -1.f);
    assert(cp.z == 1.f);
  }

  //mtx2
  {
    Matrix2 m(Matrix2::IDENTITY);
    m.RotDeg(90.f);
    Vector2 v(1, 0);
    Vector2 v2 = m * v;
    cout << "v2: " << v2 << endl;
    assert(v2.IsEqual(Vector2(0, 1.f), 1e-07));
  }

  //mtx3
  {
    Vector2 s(2, 2), t(1, 2);
    Matrix3 scale(Matrix3::IDENTITY), translate(Matrix3::IDENTITY);
    scale.MakeScaleMatrix(s);
    translate.SetTranslate(t);
    Vector3 pt(3, 4, 1);

    Vector3 newpt = translate * scale * pt;
    cout << "newpt: " << newpt << endl;

    Matrix3 transform;
    transform.MakeTransform(t, s);
    Vector3 newptv2 = transform * pt;
    cout << "newptv2: " << newptv2 << endl;
    assert(newpt.IsEqual(newptv2));
  }

  //mtx4
  {
    Vector3 s(2, 2, 2), t(1, 2, 3);
    Quaternion q(90.f, Vector3(1, 0, 0), false);
    
    Matrix4 scale(Matrix4::IDENTITY), translate(Matrix4::IDENTITY), rot(q);
    scale.MakeScaleMatrix(s);
    translate.SetTranslate(t);
    Vector4 pt(3, 4, 1, 1);

    Vector4 newpt = translate * rot * scale * pt;
    cout << "newpt: " << newpt << endl;

    Matrix4 transform;
    transform.MakeTransform(t, s, q);
    Vector4 newptv2 = transform * pt;
    cout << "newptv2: " << newptv2 << endl;
    assert(newpt.IsEqual(newptv2));
  }
  return 0;
}