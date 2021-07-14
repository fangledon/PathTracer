#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <glm/glm.hpp>

#include "Scene.h"
#include "Integrator.h"

using namespace glm;
using namespace std;

class RenderJob {

public:

    const uvec2 startPixel;
    const uvec2 windowSize;

private:

    vector<vec3> _result;

public:

    RenderJob(uvec2 startPixel, uvec2 windowSize);
    void render(Scene* scene, Integrator* integrator);
    vector<vec3> getResult();

};

class RenderPool {

private:

    Scene* _scene;
    Integrator* _integrator;
    vector<thread> _threads;
    mutex _mutex;
    condition_variable _condition;
    size_t _nextJob;
    vector<RenderJob*>& _jobQueue;
    vector<RenderJob*> _completedJobs;

public:

    RenderPool(Scene* scene, Integrator* integrator, int numThreads, vector<RenderJob*>& jobs);
    ~RenderPool();
    void getCompletedJobs(vector<RenderJob*>& completedJobs);

private:

    static void threadMain(RenderPool* pool);

};
