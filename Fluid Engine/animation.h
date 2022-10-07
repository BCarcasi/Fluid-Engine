#ifndef INCLUDE_JET_ANIMATION_H_
#define INCLUDE_JET_ANIMATION_H_

#include "macros.h"
#include <memory>

#include "timer.h"

#include "pch+.h"

#include "private_helpers.h"

namespace jet {

    //!
    //! \brief Representation of an animation frame.
    //!
    //! This struct holds current animation frame index and frame interval in
    //! seconds.
    //!
    struct Frame final {
        //! Frame index.
        int index = 0;

        //! Time interval in seconds between two adjacent frames.
        double timeIntervalInSeconds = 1.0 / 60.0;

        //! Constructs Frame instance with 1/60 seconds time interval.
        Frame();

        //! Constructs Frame instance with given time interval.
        Frame(int newIndex, double newTimeIntervalInSeconds);

        //! Returns the elapsed time in seconds.
        double timeInSeconds() const;

        //! Advances single frame.
        void advance();

        //! Advances multiple frames.
        //! \param delta Number of frames to advance.
        void advance(int delta);

        //! Advances single frame (prefix).
        Frame& operator++();

        //! Advances single frame (postfix).
        Frame operator++(int);
    };

    Frame::Frame() {
    }

    Frame::Frame(int newIndex, double newTimeIntervalInSeconds)
        : index(newIndex)
        , timeIntervalInSeconds(newTimeIntervalInSeconds) {
    }

    double Frame::timeInSeconds() const {
        return index * timeIntervalInSeconds;
    }

    void Frame::advance() {
        ++index;
    }

    void Frame::advance(int delta) {
        index += delta;
    }

    Frame& Frame::operator++() {
        advance();
        return *this;
    }

    Frame Frame::operator++(int i) {
        UNUSED_VARIABLE(i);

        Frame result = *this;
        advance();
        return result;
    }

    //!
    //! \brief Abstract base class for animation-related class.
    //!
    //! This class represents the animation logic in very abstract level.
    //! Generally animation is a function of time and/or its previous state.
    //! This base class provides a virtual function update() which can be
    //! overriden by its sub-classes to implement their own state update logic.
    //!
    class Animation {
    public:
        Animation();

        virtual ~Animation();

        //!
        //! \brief Updates animation state for given \p frame.
        //!
        //! This function updates animation state by calling Animation::onUpdate
        //! function.
        //!
        void update(const Frame& frame);

    protected:
        //!
        //! \brief The implementation of this function should update the animation
        //!     state for given Frame instance \p frame.
        //!
        //! This function is called from Animation::update when state of this class
        //! instance needs to be updated. Thus, the inherited class should overrride
        //! this function and implement its logic for updating the animation state.
        //!
        virtual void onUpdate(const Frame& frame) = 0;
    };

    //! Shared pointer for the Animation type.
    typedef std::shared_ptr<Animation> AnimationPtr;

    Animation::Animation() {
    }

    Animation::~Animation() {
    }

    void Animation::update(const Frame& frame) {
        Timer timer;

        //JET_INFO << "Begin updating frame: " << frame.index << " timeIntervalInSeconds: " << frame.timeIntervalInSeconds << " (1/" << 1.0 / frame.timeIntervalInSeconds << ") seconds";

        onUpdate(frame);

        //JET_INFO << "End updating frame (took " << timer.durationInSeconds() << " seconds)";
    }

}  // namespace jet

#endif  // INCLUDE_JET_ANIMATION_H_