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

        Vector3f AB = vertices[1] - vertices[0];
		Vector3f AC = vertices[2] - vertices[0];
		Vector3f crossProduct = Vector3f::cross(AB,AC);
		area = crossProduct.length() * 0.5f;
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


    PointLight* generateRandLight() override {


        // 随机生成一个点
        float u = RandNum(); // [0, 1)
        float v = RandNum() * (1.0f - u); // [0, 1-u)
        float w = 1.0f - u - v; // 确保 u + v + w = 1

        // 获取三角形顶点
        Vector3f A = vertices[0];
        Vector3f B = vertices[1];
        Vector3f C = vertices[2];

        // 插值生成随机点
        Vector3f randomPoint = A * w + B * u + C * v;

        Vector3f AB = vertices[1] - vertices[0];
        Vector3f AC = vertices[2] - vertices[0];

        float double_area = Vector3f::cross(AB,AC).length();

        // 创建并返回点光源
        return new PointLight(randomPoint, material->emissionColor, normal.normalized(), area);
    }

    // BoundingVolume generateBoxFromObject() override {
    //     // 计算三角形顶点的最小和最大坐标
    //     Vector3f minPoint = vertices[0];
    //     Vector3f maxPoint = vertices[0];

    //     for (int i = 1; i < 3; ++i) {
    //         minPoint = Vector3f(
    //             std::min(minPoint.x(), vertices[i].x()),
    //             std::min(minPoint.y(), vertices[i].y()),
    //             std::min(minPoint.z(), vertices[i].z())
    //         );
    //         maxPoint = Vector3f(
    //             std::max(maxPoint.x(), vertices[i].x()),
    //             std::max(maxPoint.y(), vertices[i].y()),
    //             std::max(maxPoint.z(), vertices[i].z())
    //         );
    //     }

    //     return BoundingVolume(minPoint, maxPoint);
    // }

	// 防止包围盒没有厚度，保证各维度有差值
	BoundingVolume generateBoxFromObject() override {
		float minx = std::min({vertices[0].x(), vertices[1].x(), vertices[2].x()});
		float miny = std::min({vertices[0].y(), vertices[1].y(), vertices[2].y()});
		float minz = std::min({vertices[0].z(), vertices[1].z(), vertices[2].z()});
		float maxx = std::max({vertices[0].x(), vertices[1].x(), vertices[2].x()});
		float maxy = std::max({vertices[0].y(), vertices[1].y(), vertices[2].y()});
		float maxz = std::max({vertices[0].z(), vertices[1].z(), vertices[2].z()});
		if (minx == maxx) {
			minx -= 0.001f;
			maxx += 0.001f;
		}
		if (miny == maxy) {
			miny -= 0.001f;
			maxy += 0.001f;
		}
		if (minz == maxz) {
			minz -= 0.001f;
			maxz += 0.001f;
		}
		Vector3f min(minx, miny, minz);
		Vector3f max(maxx, maxy, maxz);

		return BoundingVolume(min, max);
	}

	Vector3f normal;
	Vector3f vertices[3];
    float area;
protected:

};

#endif //TRIANGLE_H
