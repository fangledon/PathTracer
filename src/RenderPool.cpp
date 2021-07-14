#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <glm/glm.hpp>

#include "Scene.h"
#include "Integrator.h"

#include "RenderPool.h"

RenderJob::RenderJob(uvec2 startPixel, uvec2 windowSize)
    : startPixel(startPixel),
      windowSize(windowSize),
      _result(windowSize.x * windowSize.y)
{
}

void RenderJob::render(Scene* scene, Integrator* integrator)
{
    for (size_t wy = 0; wy < windowSize.y; wy++) {
        size_t y = startPixel.y + wy;
        for (size_t wx = 0; wx < windowSize.x; wx++) {
            size_t x = startPixel.x + wx;

            vec3 target =
                scene->camera.imagePlaneTopLeft
                + (x + 0.5f) * scene->camera.pixelRight
                + (y + 0.5f) * scene->camera.pixelDown;
            vec3 direction = normalize(target - scene->camera.origin);

            _result[wy * windowSize.x + wx] = integrator->traceRay(scene->camera.origin,
                                                                   direction, target);
        }
    }
}

vector<vec3> RenderJob::getResult()
{
    return move(_result);
}

RenderPool::RenderPool(Scene* scene, Integrator* integrator, int numThreads, vector<RenderJob*>& jobs)
    : _scene(scene), _integrator(integrator), _nextJob(0), _jobQueue(jobs)
{
    for (int i = 0; i < numThreads; i++) {
        _threads.push_back(thread(threadMain, this));
    }
}

RenderPool::~RenderPool()
{
    for (thread& thread : _threads) {
        thread.join();
    }
}

void RenderPool::getCompletedJobs(vector<RenderJob*>& completedJobs)
{
    {
        unique_lock<mutex> lock(_mutex);

        _condition.wait(lock, [this]{ return _completedJobs.size() > 0; });
        completedJobs = move(_completedJobs);
    }
}

void RenderPool::threadMain(RenderPool* pool)
{
    while (true) {

        size_t jobIndex;
        {
            unique_lock<mutex> lock(pool->_mutex);

            if (pool->_nextJob >= pool->_jobQueue.size()) break;

            jobIndex = pool->_nextJob;
            pool->_nextJob++;
        }

        pool->_jobQueue[jobIndex]->render(pool->_scene, pool->_integrator);

        {
            unique_lock<mutex> lock(pool->_mutex);

            pool->_completedJobs.push_back(pool->_jobQueue[jobIndex]);
            pool->_condition.notify_all();
        }
    }
}
