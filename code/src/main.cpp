#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>

#include "scene_parser.hpp"
#include "image.hpp"
#include "camera.hpp"
#include "group.hpp"
#include "light.hpp"

#include <string>

#include "tracer.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    for (int argNum = 1; argNum < argc; ++argNum) {
        std::cout << "Argument " << argNum << " is: " << argv[argNum] << std::endl;
    }

    if (argc != 5) {
        cout << "Usage: ./bin/PA1 <input scene file> <output bmp file> type spp" << endl;
        return 1;
    }
    string inputFile = argv[1];
    string outputFile = argv[2];  // only bmp is allowed.
    string tracerType = argv[3];
    int samplesPerPixel = atoi(argv[4]);

    // TODO: Main RayCasting Logic
    // First, parse the scene using SceneParser.
    // Then loop over each pixel in the image, shooting a ray
    // through that pixel and finding its intersection with
    // the scene.  Write the color at the intersection to that
    // pixel in your output image.
    SceneParser sceneParser(inputFile.c_str());

    int maxDepth = 10;            // 默认最大递归深度

    if (tracerType == "raycast") {
        RayCastTracer tracer(&sceneParser, outputFile);
        tracer.render();
    }else if (tracerType == "whitted") {
        WhittedStyleTracer tracer(&sceneParser, maxDepth, outputFile);
        tracer.render();
    } else if (tracerType == "simple") {
        SimplePathTracer tracer(&sceneParser, samplesPerPixel, outputFile);
        tracer.render();
    } else if (tracerType == "nee") {
        NEEedPathTracer tracer(&sceneParser, samplesPerPixel, outputFile);
        tracer.render();
    }else if (tracerType == "glossy") {
        GlossyPathTracer tracer(&sceneParser, samplesPerPixel, outputFile);
        tracer.render();
    }else if (tracerType == "bvh"){
        BVHTracer tracer(&sceneParser, samplesPerPixel, outputFile);
        tracer.render();
    } else {
        cerr << "Unknown tracer type: " << tracerType << endl;
        cout << "Available tracer types: whitted, simple, nee" << endl;
        return 1;
    }

    cout << "Hello! Computer Graphics!" << endl;
    return 0;
}



