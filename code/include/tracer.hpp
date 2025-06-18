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
#include <cmath>

using namespace std;

#include <cstdlib>
#include <ctime>

// 初始化随机种子（只需在程序启动时调用一次）
static bool seedInitialized = false;
inline void initRand() {
    if (!seedInitialized) {
        srand(static_cast<unsigned int>(time(0))); // 使用当前时间作为种子
        seedInitialized = true;
    }
}

// 随机数生成函数
inline double RandNum() {
    initRand(); // 确保种子初始化
    return static_cast<double>(rand()) / RAND_MAX; // 返回 [0, 1) 的随机数
}


Vector3f reflectDir(const Vector3f& I, const Vector3f& N) {
    // 标准化入射角度
    Vector3f newI = I.normalized();
    // 返回反射方向
    return (newI - 2 * Vector3f::dot(newI, N) * N).normalized();
}


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


// 计算折射方向
Vector3f refractDir(const Vector3f& I, const Vector3f& N, double ior) {
    // 标准化入射方向和法线
    Vector3f incident = I.normalized();
    Vector3f normal = N.normalized();

    // 判断光线方向：从外部进入 (cosThetaI > 0) 或内部射出 (cosThetaI < 0)
    double cosThetaI = Vector3f::dot(incident, normal);
    Vector3f adjustedNormal = normal;
    double eta = ior; // 折射率比

    if (cosThetaI < 0) {
        // 光线从内部射出，调整法线方向和折射率
        adjustedNormal = -normal;
        eta = 1.0 / ior; // 反向传播时取倒数
        cosThetaI = -cosThetaI; // 修正 cos 值
    }

    // 计算 sin^2(theta_t) 使用 Snell 定律
    double sin2ThetaI = 1.0 - cosThetaI * cosThetaI;
    double sin2ThetaT = eta * eta * sin2ThetaI;

    // 全反射检查
    if (sin2ThetaT >= 1.0) {
        // 发生全反射，返回反射方向
        return (incident + 2.0 * cosThetaI * adjustedNormal).normalized();
    }

    // 计算 cos(theta_t)
    double cosThetaT = std::sqrt(1.0 - sin2ThetaT);

    // 计算折射方向
    Vector3f refractDir = eta * incident + (eta * cosThetaI - cosThetaT) * adjustedNormal;

    // 菲涅尔反射概率 (Schlick 近似)
    double R0 = ((1.0 - ior) * (1.0 - ior)) / ((1.0 + ior) * (1.0 + ior));
    double reflectance = R0 + (1.0 - R0) * std::pow(1.0 - cosThetaI, 5.0);

    // 随机决定反射或折射
    if (RandNum() < reflectance) {
        // 选择反射
        return (incident + 2.0 * cosThetaI * adjustedNormal).normalized();
    }

    // 选择折射
    return refractDir.normalized();
}

// 漫反射方向生成函数
Vector3f diffuseDir(const Vector3f &normal) {
    // 生成两个随机数 [0, 1)
    double r1 = 2.0 * M_PI * RandNum(); // 角度 [0, 2π)
    double r2 = RandNum();              // 随机值 [0, 1)

    // 使用平方根变换生成余弦加权分布
    double r2s = std::sqrt(r2);

    // 计算球面坐标
    double x = std::cos(r1) * r2s; // x 坐标
    double y = std::sin(r1) * r2s; // y 坐标
    double z = std::sqrt(1.0 - r2); // z 坐标，余弦加权

    // 创建临时向量
    Vector3f dir(x, y, z);

    // 构造局部坐标系
    Vector3f w = normal; // 法线作为 z 轴
    Vector3f u = Vector3f::cross(std::abs(w.x()) > 0.1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0), w).normalized(); // u 轴
    Vector3f v = Vector3f::cross(w, u); // v 轴

    // 将方向转换到世界坐标系
    Vector3f worldDir = u * dir.x() + v * dir.y() + w * dir.z();

    // 确保方向在法线半球内并归一化
    if (Vector3f::dot(worldDir,normal) < 0) {
        worldDir = -worldDir; // 如果方向在下半球，反转
    }
    return worldDir.normalized();
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
        printf("WhittedTracing: Finished rendering %s\n", _outputFile.c_str());
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
    SimplePathTracer(const SceneParser* sceneParser, const string& outputFile): _sceneParser(sceneParser), _outputFile(outputFile){};

    Vector3f trace(Ray ray) {
        Vector3f color = Vector3f::ZERO;
        Vector3f throughput = Vector3f(1, 1, 1);

        while (true) {
            Hit hit;
            if (!_sceneParser->getGroup()->intersect(ray, hit, 0.001f)) {
                color += _sceneParser->getBackgroundColor() * throughput;
                break;
            }

            Material* material = hit.getMaterial();
            Vector3f hitPoint = ray.pointAtParameter(hit.getT());
            Vector3f normal = hit.getNormal().normalized();
            Vector3f input = ray.direction.normalized();
            Vector3f emitted = material->emissionColor;
            Vector3f diffuse = material->diffuseColor;

            color += diffuse * emitted * throughput;
            throughput = throughput * diffuse;

            ray.origin = hitPoint;
            double rate = RandNum();
            Vector3f type = material->type_d_rl_rr;
            Vector3f newdir;
            if (rate <= type.x()) { // 漫反射
                ray.direction = diffuseDir(normal);
            } else if (rate <= type.x() + type.y()) { // 反射
                ray.direction = reflectDir(ray.direction, normal);
            } else { // 折射
                ray.direction = refractDir(ray.direction, normal, material->refr_rate);
            }

            if (type.x() != 0) {
                throughput *= std::abs(Vector3f::dot(ray.direction, normal));
            }

            // 俄罗斯轮盘赌
            float q = 0.1f; // 固定终止概率，可动态调整
            double rouletteRand = RandNum(); // 新生成的随机数
            if (rouletteRand < q) { // 以概率 q 终止
                break;
            }
            throughput = throughput / (1 - q); // 补偿继续概率
        }

        return color;
    }

    void render(){
        Camera* camera = _sceneParser->getCamera();
        int width = camera->getWidth();
        int height = camera->getHeight();
        Image image = Image(width, height);

        int samplesPerPixel = 1000;

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                Vector3f finalColor(0, 0, 0);
                for (int s = 0; s < samplesPerPixel; s++) {
                    float u = (x + RandNum());
                    float v = (y + RandNum());
                    Ray ray = camera->generateRay(Vector2f(u, v));
                    finalColor += trace(ray);
                }
                finalColor[0] /= samplesPerPixel;
                finalColor[1] /= samplesPerPixel;
                finalColor[2] /= samplesPerPixel;
                image.SetPixel(x, y, finalColor);
            }
        }

        image.SaveBMP(_outputFile.c_str());
        printf("SimplePathTracing: Finished rendering %s\n", _outputFile.c_str());
    }

private:
    const SceneParser* _sceneParser;
    std::string _outputFile;

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