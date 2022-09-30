#ifndef INCLUDE_JET_TIMER_H_
#define INCLUDE_JET_TIMER_H_

#include <chrono>

namespace jet {

    //! Simple timer class.
    class Timer {
    public:
        //! Constructs the timer and start ticking.
        Timer();

        //! Returns the time duration since the creation or reset in seconds.
        double durationInSeconds() const;

        //! Resets the timer.
        void reset();

    private:
        std::chrono::steady_clock _clock;
        std::chrono::steady_clock::time_point _startingPoint;
    };

    Timer::Timer() {
        _startingPoint = _clock.now();
    }

    double Timer::durationInSeconds() const {
        auto end = std::chrono::steady_clock::now();
        auto count = std::chrono::duration_cast<std::chrono::microseconds>(
            end - _startingPoint).count();
        return count / 1000000.0;
    }

    void Timer::reset() {
        _startingPoint = _clock.now();
    }

}  // namespace jet

#endif  // INCLUDE_JET_TIMER_H_