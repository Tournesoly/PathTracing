#ifndef GROUP_H
#define GROUP_H


#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
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

    PointLight* generateRandLight() override {
        printf("Generate a random light in group\n");
        return new PointLight(Vector3f(0, 0, 0), material->emissionColor);
    }

    void addObject(int index, Object3D *obj) {
        group_objects.push_back(obj);
    }

    int getGroupSize() {
        return group_objects.size();
    }

        // 获取发光物体数量
    int getNumEmissions() {
        // 若仍未加载发光物体则加载
        if (loaded == 0) {
            loaded = 1;
            int n=0;
            for (auto obj : group_objects) {
                if (obj->material->emissionColor.x() != 0 || obj->material->emissionColor.y() != 0 || obj->material->emissionColor.z() != 0){
                    group_emissions.push_back(obj);
                }
            }
        }
        return group_emissions.size();
    }

    // 从某一发光物体采样出一个点光源
    PointLight* getEmission(int i) {
        // 若仍未加载发光物体则加载
        
        if (loaded == 0) {
            loaded = 1;
            for (auto obj : group_objects) {
                if (obj->material->emissionColor.x() != 0 || obj->material->emissionColor.y() != 0 || obj->material->emissionColor.z() != 0){
                    group_emissions.push_back(obj);
                }
            }
        }
        PointLight * light = group_emissions[i]->generateRandLight();
        return light;
    }


private:
    // 所有发光物体
    std::vector<Object3D*> group_emissions;
        // 发光物体已被加载
    bool loaded=0;
};

#endif
	
