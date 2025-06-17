#ifndef TRACER_HPP
#define TRACER_HPP

#include "scene_parser.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include "material.hpp"
#include "image.hpp"
#include "camera.hpp"

#include "group.hpp"
#include "light.hpp"

using namespace std;


Vector3f reflectDir(const Vector3f& I, const Vector3f& N) {
    // 标准化入射角度
    Vector3f newI = I.normalized();
    // 返回反射方向
    return (newI - 2 * Vector3f::dot(newI, N) * N).normalized();
}

#include <cmath>

// 计算折射方向
Vector3f refractDir_whitted(const Vector3f& I, const Vector3f& N, float ior) {
    Vector3f normalized_I = I.normalized();
    Vector3f normalized_N = N.normalized();

    // 判断是从外部进入还是从内部出射
    float cosThetaI = -Vector3f::dot(normalized_I, normalized_N);

    if (cosThetaI < 0) {
        // 光线从物体内部射出 -> 调整法线方向和折射率
        return refractDir_whitted(normalized_I, -normalized_N, 1.0f / ior);
    }

    // 折射率比值
    float eta = 1.0f / ior;

    // sin^2(theta_t) = eta^2 * (1 - cos^2(theta_i))
    float sin2ThetaT = eta * eta * (1.0f - cosThetaI * cosThetaI);

    // 全反射判断
    if (sin2ThetaT > 1.0f) {
        // 发生全反射，返回反射方向
        return reflectDir(I, N);
    }

    // 计算折射方向
    float cosThetaT = sqrt(1.0f - sin2ThetaT);
    Vector3f T = eta * normalized_I + (eta * cosThetaI - cosThetaT) * normalized_N;
    return T.normalized();
}



class RayCastTracer {
public:
    RayCastTracer(const SceneParser* sceneParser, const string& outputFile) : _sceneParser(sceneParser), _outputFile(outputFile){};

    void render() {
        Camera *camera = _sceneParser->getCamera();
        Image img(camera->getWidth(), camera->getHeight());
        // 循环屏幕空间的像素
        for ( int x = 0; x < camera->getWidth() ; ++x) {
            for ( int y = 0; y < camera->getHeight() ; ++y) {
                // 计算当前像素(x,y)处相机出射光线camRay
                Ray camRay = camera->generateRay(Vector2f(x, y)) ;
                Group* baseGroup = _sceneParser->getGroup() ;
                Hit hit ;
                // 判断camRay是否和场景有交点，并返回最近交点的数据，存储在hit 中
                bool isIntersect = baseGroup->intersect (camRay, hit , 0) ;
                if ( isIntersect ) {
                    Vector3f finalColor = Vector3f ::ZERO;
                    // 找到交点之后，累加来自所有光源的光强影响
                    for ( int li = 0; li < _sceneParser->getNumLights() ; ++li ) {
                        Light* light = _sceneParser->getLight( li );
                        Vector3f L, lightColor ;
                        float fuck_dist ;
                        // 获得光照强度
                        light->getIllumination(camRay. pointAtParameter( hit .getT() ) , L, lightColor, fuck_dist) ;
                        // 计算局部光强
                        finalColor += hit.getMaterial()->Shade(camRay, hit , L, lightColor) ;
                    }
                    img.SetPixel(x, y, finalColor);
                } else {
                // 不存在交点，返回背景色
                    img.SetPixel (x, y, _sceneParser->getBackgroundColor() ) ;
                }
            }
        }
        img.SaveBMP(_outputFile.c_str());

    }

private:
    const SceneParser* _sceneParser;
    const string &_outputFile;

};


// ========================
// Whitted 风格光线追踪器
// ========================
class WhittedStyleTracer {
public:
    WhittedStyleTracer(const SceneParser* sceneParser, int maxDepth, const string& outputFile): _sceneParser(sceneParser), _maxDepth(maxDepth), _outputFile(outputFile){};

    Vector3f trace(const Ray& ray, int depth){
        if (depth <= 0) {
            return Vector3f::ZERO; // 最大递归深度限制
        }

        Hit hit;
        if (!_sceneParser->getGroup()->intersect(ray, hit, 0.001f)) {
            return _sceneParser->getBackgroundColor(); // 没有交点，返回背景色
        }

        Material* material = hit.getMaterial();
        Vector3f hitPoint = ray.pointAtParameter(hit.getT());
        Vector3f normal = hit.getNormal().normalized();
        Vector3f finalColor = Vector3f::ZERO;

        if (material->type_d_rl_rr.x()){
            //漫反射，累加来自所有光源的光强影响
            for ( int li = 0; li < _sceneParser->getNumLights() ; ++li ) {
                // 点光源
                Light* light = _sceneParser->getLight( li );
                Vector3f lightDir, lightColor ;
                float lightDist;
                // 获得光照信息
                light->getIllumination(hitPoint, lightDir, lightColor, lightDist);
                Ray shadowRay(hitPoint + normal * 0.001, lightDir);
                Hit shadowHit;
                if ( _sceneParser->getGroup()->intersect(shadowRay, shadowHit, 0.001f) && shadowHit.getT() < lightDist ){
                    printf("shadowHit and lightDist: %f %f\n", shadowHit.getT(), lightDist);
                    // 存在阴影
                }else {
                    // 计算局部光强 phong 模型
                    finalColor += hit.getMaterial()->Shade(ray, hit , lightDir, lightColor) ;
                }

            }
        } else if (material->type_d_rl_rr.y()){
            //反射
            Vector3f refl_ray = reflectDir(ray.getDirection(),normal);
            Ray reflRay(hitPoint + normal * 0.001, refl_ray);
            finalColor += trace(reflRay, depth - 1);
        } else if (material->type_d_rl_rr.z()){
            //折射
            Vector3f refr_ray = refractDir_whitted(ray.getDirection(),normal, material->refr_rate);
            Ray refrRay(hitPoint - normal * 0.001, refr_ray);
            // 有可能全反射
            if (Vector3f::dot(refr_ray, normal) > 0){
                refrRay.setOrigin(hitPoint + normal * 0.001);
            }
            finalColor += trace(refrRay, depth - 1);
        }

        return finalColor;


    }
    void render(){
        Camera* camera = _sceneParser->getCamera();
        int width = camera->getWidth();
        int height = camera->getHeight();
        Image image = Image(width, height);

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                Vector2f point(i, j);
                Ray ray = camera->generateRay(point);
                Vector3f color = trace(ray, _maxDepth);
                image.SetPixel(i, j, color);
            }
        }

        image.SaveBMP(_outputFile.c_str());
        printf("Finished rendering %s\n", _outputFile.c_str());
    }

private:
    const SceneParser* _sceneParser;
    int _maxDepth;
    string _outputFile;
};

// ========================
// 简单路径追踪器（基础 Path Tracer）
// ========================
class SimplePathTracer {
public:
    SimplePathTracer(const SceneParser* scene, int maxDepth, int samplesPerPixel, const std::string& outputFile);

    Vector3f trace(Ray ray){

    }
    void render(){

    }

private:
    const SceneParser* scene;
    int maxDepth;
    int samplesPerPixel;
    std::string outputFile;

};

// ========================
// 带下一事件估计（NEE）的路径追踪器
// ========================
class NEEedPathTracer {
public:
    NEEedPathTracer(const SceneParser* scene, int maxDepth, int samplesPerPixel, const std::string& outputFile);

    Vector3f traced(Ray ray, int depth){

    }
    void render(){
        
    }

private:
    const SceneParser* scene;
    int maxDepth;
    int samplesPerPixel;
    std::string outputFile;


};

#endif // TRACER_HPP