#ifndef INCLUDE_JET_PHYSICS_ANIMATION_H_
#define INCLUDE_JET_PHYSICS_ANIMATION_H_

#include "animation.h"

#include "pch+.h"

#include "constants.h"
#include "timer.h"

#include <limits>

namespace jet {

    //
    //      Abstract base class for physics-based animation.
    //
    // This class represents physics-based animation by adding time-integration
    // specific functions to Animation class.
    //
    class PhysicsAnimation : public Animation {
    public:
        // Default constructor.
        PhysicsAnimation();

        // Destructor.
        virtual ~PhysicsAnimation();

        //
        //      Returns true if fixed sub-timestepping is used.
        //
        // When performing a time-integration, it is often required to take
        // sub-timestepping for better results. The sub-stepping can be either
        // fixed rate or adaptive, and this function returns which feature is
        // currently selected.
        //
        // \return     True if using fixed sub time steps, false otherwise.
        //
        bool isUsingFixedSubTimeSteps() const;

        //
        //      Sets true if fixed sub-timestepping is used.
        //
        // When performing a time-integration, it is often required to take
        // sub-timestepping for better results. The sub-stepping can be either
        // fixed rate or adaptive, and this function sets which feature should be
        // selected.
        //
        //   isUsing True to enable fixed sub-stepping.
        //
        void setIsUsingFixedSubTimeSteps(bool isUsing);

        //
        //      Returns the number of fixed sub-timesteps.
        //
        // When performing a time-integration, it is often required to take
        // sub-timestepping for better results. The sub-stepping can be either
        // fixed rate or adaptive, and this function returns the number of fixed
        // sub-steps.
        //
        // \return     The number of fixed sub-timesteps.
        //
        unsigned int numberOfFixedSubTimeSteps() const;

        //
        //      Sets the number of fixed sub-timesteps.
        //
        // When performing a time-integration, it is often required to take
        // sub-timestepping for better results. The sub-stepping can be either
        // fixed rate or adaptive, and this function sets the number of fixed
        // sub-steps.
        //
        //  numberOfSteps The number of fixed sub-timesteps.
        //
        void setNumberOfFixedSubTimeSteps(unsigned int numberOfSteps);

        // Advances a single frame.
        void advanceSingleFrame();

        //
        //      Returns current frame.
        //
        Frame currentFrame() const;

        //
        //      Sets current frame cursor (but do not invoke update()).
        //
        void setCurrentFrame(const Frame& frame);

        //
        //      Returns current time in seconds.
        //
        // This function returns the current time which is calculated by adding
        // current frame + sub-timesteps it passed.
        //
        double currentTimeInSeconds() const;

    protected:
        //
        //      Called when a single time-step should be advanced.
        //
        // When Animation::update function is called, this class will internally
        // subdivide a frame into sub-steps if needed. Each sub-step, or time-step,
        // is then taken to move forward in time. This function is called for each
        // time-step, and a subclass that inherits PhysicsAnimation class should
        // implement this function for its own physics model.
        //
        //  timeIntervalInSeconds The time interval in seconds
        //
        virtual void onAdvanceTimeStep(double timeIntervalInSeconds) = 0;

        //
        //      Returns the required number of sub-timesteps for given time
        //             interval.
        //
        // The required number of sub-timestep can be different depending on the
        // physics model behind the implementation. Override this function to
        // implement own logic for model specific sub-timestepping for given
        // time interval.
        //
        //  timeIntervalInSeconds The time interval in seconds.
        //
        // \return     The required number of sub-timesteps.
        //
        virtual unsigned int numberOfSubTimeSteps(
            double timeIntervalInSeconds) const;

        //
        //      Called at frame 0 to initialize the physics state.
        //
        // Inheriting classes can override this function to setup initial condition
        // for the simulation.
        //
        virtual void onInitialize();

    private:
        Frame _currentFrame;
        bool _isUsingFixedSubTimeSteps = true;
        unsigned int _numberOfFixedSubTimeSteps = 1;
        double _currentTime = 0.0;

        void onUpdate(const Frame& frame) final;

        void advanceTimeStep(double timeIntervalInSeconds);

        void initialize();
    };

    typedef std::shared_ptr<PhysicsAnimation> PhysicsAnimationPtr;

    PhysicsAnimation::PhysicsAnimation() { _currentFrame.index = -1; }

    PhysicsAnimation::~PhysicsAnimation() {}

    bool PhysicsAnimation::isUsingFixedSubTimeSteps() const {
        return _isUsingFixedSubTimeSteps;
    }

    void PhysicsAnimation::setIsUsingFixedSubTimeSteps(bool isUsing) {
        _isUsingFixedSubTimeSteps = isUsing;
    }

    unsigned int PhysicsAnimation::numberOfFixedSubTimeSteps() const {
        return _numberOfFixedSubTimeSteps;
    }

    void PhysicsAnimation::setNumberOfFixedSubTimeSteps(
        unsigned int numberOfSteps) {
        _numberOfFixedSubTimeSteps = numberOfSteps;
    }

    void PhysicsAnimation::advanceSingleFrame() {
        Frame f = _currentFrame;
        update(++f);
    }

    Frame PhysicsAnimation::currentFrame() const { return _currentFrame; }

    void PhysicsAnimation::setCurrentFrame(const Frame& frame) {
        _currentFrame = frame;
    }

    double PhysicsAnimation::currentTimeInSeconds() const { return _currentTime; }

    unsigned int PhysicsAnimation::numberOfSubTimeSteps(
        double timeIntervalInSeconds) const {
        UNUSED_VARIABLE(timeIntervalInSeconds);

        // Returns number of fixed sub-timesteps by default
        return _numberOfFixedSubTimeSteps;
    }

    void PhysicsAnimation::onUpdate(const Frame& frame) {
        if (frame.index > _currentFrame.index) {
            if (_currentFrame.index < 0) {
                initialize();
            }

            int32_t numberOfFrames = frame.index - _currentFrame.index;

            for (int32_t i = 0; i < numberOfFrames; ++i) {
                advanceTimeStep(frame.timeIntervalInSeconds);
            }

            _currentFrame = frame;
        }
    }

    void PhysicsAnimation::advanceTimeStep(double timeIntervalInSeconds) {
        _currentTime = _currentFrame.timeInSeconds();

        if (_isUsingFixedSubTimeSteps) {
            //JET_INFO << "Using fixed sub-timesteps: " << _numberOfFixedSubTimeSteps;

            // Perform fixed time-stepping
            const double actualTimeInterval =
                timeIntervalInSeconds /
                static_cast<double>(_numberOfFixedSubTimeSteps);

            for (unsigned int i = 0; i < _numberOfFixedSubTimeSteps; ++i) {
                //JET_INFO << "Begin onAdvanceTimeStep: " << actualTimeInterval << " (1/" << 1.0 / actualTimeInterval << ") seconds";

                Timer timer;
                onAdvanceTimeStep(actualTimeInterval);

                //JET_INFO << "End onAdvanceTimeStep (took "<< timer.durationInSeconds() << " seconds)";

                _currentTime += actualTimeInterval;
            }
        }
        else {
            //JET_INFO << "Using adaptive sub-timesteps";

            // Perform adaptive time-stepping
            double remainingTime = timeIntervalInSeconds;
            while (remainingTime > kEpsilonD) {
                unsigned int numSteps = numberOfSubTimeSteps(remainingTime);
                double actualTimeInterval =
                    remainingTime / static_cast<double>(numSteps);

                //JET_INFO << "Number of remaining sub-timesteps: " << numSteps;

                //JET_INFO << "Begin onAdvanceTimeStep: " << actualTimeInterval << " (1/" << 1.0 / actualTimeInterval << ") seconds";

                Timer timer;
                onAdvanceTimeStep(actualTimeInterval);

                //JET_INFO << "End onAdvanceTimeStep (took "<< timer.durationInSeconds() << " seconds)";

                remainingTime -= actualTimeInterval;
                _currentTime += actualTimeInterval;
            }
        }
    }

    void PhysicsAnimation::initialize() { onInitialize(); }

    void PhysicsAnimation::onInitialize() {
        // Do nothing
    }

}  // namespace jet

#endif  // INCLUDE_JET_PHYSICS_ANIMATION_H_