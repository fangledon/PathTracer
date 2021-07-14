#include <iostream>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <chrono>

#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <embree3/rtcore.h>
#include <lodepng/lodepng.h>

#include "Scene.h"
#include "SceneLoader.h"
#include "Integrator.h"
#include "RenderPool.h"

#include "Engine.h"

using Clock = chrono::high_resolution_clock;
using TimePoint = chrono::time_point<Clock>;
using Duration = chrono::duration<float>;

const unsigned int WINDOW_DIM = 32;

static unsigned char convertColorChannel(float channel)
{
    return static_cast<unsigned char>(fmin(255.0f, fmax(0.0f, 255.0f * channel)));
}

static void saveImage(
    const vector<vec3>& imageData,
    uvec2 imageSize,
    const string& fileName)
{
    vector<unsigned char> imageByteData(imageSize.y * imageSize.x * 3);
    for (size_t y = 0; y < imageSize.y; y++) {
        for (size_t x = 0; x < imageSize.x; x++) {

            vec3 color = imageData[y * imageSize.x + x];

            size_t outPixelBase = 3 * (y * imageSize.x + x);
            imageByteData[outPixelBase + 0] = convertColorChannel(color.r);
            imageByteData[outPixelBase + 1] = convertColorChannel(color.g);
            imageByteData[outPixelBase + 2] = convertColorChannel(color.b);
        }
    }

    unsigned int error = lodepng::encode(
        fileName,
        imageByteData,
        imageSize.x,
        imageSize.y,
        LCT_RGB);
    if (error) {
        throw runtime_error(
            "LodePNG error (" + to_string(error) + "): "
            + lodepng_error_text(error));
    }
}

static void embreeErrorFunction(void* userPtr, RTCError error, const char* str)
{
    (void) userPtr;
    cerr << "Embree error (" << error << "): " << str << endl;
}

static void printLoadingBar(float completion, int numBars = 60)
{
    int barsComplete = static_cast<int>(floor(numBars * completion));
    int percentComplete = static_cast<int>(floor(100 * completion));

    ostringstream oss;

    oss << "\r[";
    int j = 1;
    for (; j <= barsComplete; j++) {
        oss << '#';
    }
    for (; j <= numBars; j++) {
        oss << ' ';
    }
    oss << "] " << percentComplete << "%\r";

    cout << oss.str() << flush;
}

void render(const string& sceneFilePath)
{

    cout << "Loading scene..." << endl;

    RTCDevice embreeDevice = rtcNewDevice(nullptr);
    if (!embreeDevice) throw runtime_error("Could not initialize Embree device.");

    rtcSetDeviceErrorFunction(embreeDevice, embreeErrorFunction, nullptr);

    Scene* scene;
    loadScene(sceneFilePath, embreeDevice, &scene);
    
    RayTracerIntegrator integrator;
    integrator.setScene(scene);

    cout << "Preparing render jobs..." << endl;

    int numThreads = thread::hardware_concurrency();

    vector<RenderJob*> jobs;
    for (unsigned int y = 0; y < scene->imageSize.y; y += WINDOW_DIM) {
        for (unsigned int x = 0; x < scene->imageSize.x; x += WINDOW_DIM) {
            uvec2 startPixel = uvec2(x, y);
            uvec2 windowSize = uvec2(
                fmin(x + WINDOW_DIM, scene->imageSize.x) - x,
                fmin(y + WINDOW_DIM, scene->imageSize.y) - y);
            jobs.push_back(new RenderJob(startPixel, windowSize));
        }
    }

    vector<vec3> imageData(scene->imageSize.y * scene->imageSize.x);

    cout
        << "Rendering... ("
        << jobs.size() << " jobs, "
        << numThreads << " threads)"
        << endl;

    TimePoint startTime = Clock::now();
    {
        RenderPool pool(scene, &integrator, numThreads, jobs);

        size_t numCompletedJobs = 0;
        while (numCompletedJobs < jobs.size()) {

            vector<RenderJob*> completedJobs;
            pool.getCompletedJobs(completedJobs);

            for (RenderJob* job : completedJobs) {
                vector<vec3> result = job->getResult();
                for (unsigned int wy = 0; wy < job->windowSize.y; wy++) {
                    unsigned int y = job->startPixel.y + wy;
                    for (unsigned int wx = 0; wx < job->windowSize.x; wx++) {
                        unsigned int x = job->startPixel.x + wx;
                        imageData[y * scene->imageSize.x + x] = result[wy * job->windowSize.x + wx];
                    }
                }
            }
            numCompletedJobs += completedJobs.size();

            printLoadingBar(static_cast<float>(numCompletedJobs) / jobs.size());
        }
    }
    TimePoint endTime = Clock::now();
    Duration renderTime = endTime - startTime;

    cout << endl;
    cout << "Render time: " << renderTime.count() << "s" << endl;

    rtcReleaseScene(scene->embreeScene);
    rtcReleaseDevice(embreeDevice);

    saveImage(imageData, scene->imageSize, scene->outputFileName);
    cout << "Image saved as '" << scene->outputFileName << "'" << endl;
}
