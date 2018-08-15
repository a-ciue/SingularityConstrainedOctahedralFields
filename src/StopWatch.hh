/*
 * StopWatch.hh
 */

#ifndef STOPWATCH_HH
#define STOPWATCH_HH

#include <chrono>

template <typename UnitT = std::chrono::milliseconds>
class StopWatch
{
public:
    // Constructor
    StopWatch() {}
    // Destructor
    ~StopWatch() {}

    // Start time measurement
    void start()
    {
        starttime_ = std::chrono::steady_clock::now();
    }

    // Restart, return time elapsed until now.
    double restart()
    {
        double t = elapsed();
        start();
        return t;
    }

    // Stop time measurement, return time.
    double stop()
    {
        endtime_ = std::chrono::steady_clock::now();
        return elapsed();
    }

    // Get the total time in UnitT (watch has to be stopped).
    double elapsed() const
    {
        auto duration = std::chrono::duration_cast<UnitT>(endtime_ - starttime_);
        return (double) duration.count();
    }

private:
    std::chrono::steady_clock::time_point starttime_, endtime_;
};


#endif /* STOPWATCH_HH*/
