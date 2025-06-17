#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <iostream>
using namespace std;

// TODO: implement this class and add more fields as necessary,
class Triangle: public Object3D {

public:
	Triangle() = delete;

    // a b c are three vertex positions of the triangle
	Triangle( const Vector3f& a, const Vector3f& b, const Vector3f& c, Material* m) : Object3D(m) {
		vertices[0] = a;
        vertices[1] = b;
        vertices[2] = c;
        // 计算法向量（未标准化）
        normal = Vector3f::cross(b - a, c - a);
	}

	bool intersect( const Ray& ray,  Hit& hit , float tmin) override {
		if (Vector3f::dot(normal, ray.getDirection()) == 0){
            return false;
        }

        //向量准备
        Vector3f rd=ray.getDirection();
        Vector3f e1=vertices[0]-vertices[1];
        Vector3f e2=vertices[0]-vertices[2];
        Vector3f s=vertices[0]-ray.getOrigin();

        // 利用Crammer法则计算t、beta、gamma
        float t=Matrix3f(s,e1,e2).determinant()/Matrix3f(rd,e1,e2).determinant();
        float beta=Matrix3f(rd,s,e2).determinant()/Matrix3f(rd,e1,e2).determinant();
        float gamma=Matrix3f(rd,e1,s).determinant()/Matrix3f(rd,e1,e2).determinant();

        if (t<tmin||t>=hit.getT()||beta<0||gamma<0||beta>1||gamma>1||(beta+gamma)
<0||(beta+gamma)>1){
            return false;
        }

        // 保证法向量朝外
        Vector3f dnormal=normal;
        if (Vector3f::dot(normal, ray.getDirection())>0){
            dnormal=-normal;
        }
        hit.set(t, material, dnormal);
        return true;
	}
	Vector3f normal;
	Vector3f vertices[3];
protected:

};

#endif //TRIANGLE_H
