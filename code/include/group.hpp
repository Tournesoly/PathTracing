#ifndef GROUP_H
#define GROUP_H


#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include "object_tree.hpp"
#include <iostream>
#include <vector>


// TODO: Implement Group - add data structure to store a list of Object*
class Group : public Object3D {

public:
    std::vector<Object3D*> group_objects;

    Group() {

    }

    explicit Group (int num_objects) {

    }

    ~Group() override {

    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        bool intersected = false;
        for (auto obj : group_objects) {
            if (obj->intersect(r, h, tmin)) {
                intersected = true;
            }
        }
        return intersected;
    }

    bool intersectWithBVH(const Ray &r, Hit &h, float tmin) {
        return tree->intersect(r, h, tmin);
    }

    PointLight* generateRandLight() override {
        printf("Generate a random light in group\n");
        return new PointLight(Vector3f(0, 0, 0), material->emissionColor);
    }

    BoundingVolume generateBoxFromObject() override {
        //nothing
        //
        return BoundingVolume(Vector3f(-INFINITY, -INFINITY, -INFINITY), Vector3f(INFINITY, INFINITY, INFINITY));
    }


    void addObject(int index, Object3D *obj) {
        group_objects.push_back(obj);
    }

    int getGroupSize() {
        return group_objects.size();
    }

    // 获取发光物体数量
    int getNumEmissions() {
        return group_emissions.size();
    }

    // 从某一发光物体采样出一个点光源
    PointLight* getEmission(int i) {
        // 若仍未加载发光物体则加载
        PointLight * light = group_emissions[i]->generateRandLight();
        return light;
    }


    void initGroup(){
        // emission
        for (auto obj : group_objects) {
            if (obj->material->emissionColor.x() != 0 || obj->material->emissionColor.y() != 0 || obj->material->emissionColor.z() != 0){
                group_emissions.push_back(obj);
            }
        }

        // OBTree
        tree = new OBTree(group_objects, 0, group_objects.size());
    }


private:
    // 所有发光物体
    std::vector<Object3D*> group_emissions;
        // 发光物体已被加载
    OBTree *tree;
};

#endif
	
