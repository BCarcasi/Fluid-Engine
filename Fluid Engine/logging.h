#ifndef INCLUDE_JET_LOGGING_H_
#define INCLUDE_JET_LOGGING_H_


#include "logging.h"
#include "macros.h"

#include <chrono>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>

#include <iostream>
#include <sstream>
#include <string>

namespace jet {

    //! Level of the logging.
    //! All < Debug < Info < Warn < Error < Off.
    enum class LoggingLevel : uint8_t {
        All = 0,
        Debug = 1,
        Info = 2,
        Warn = 3,
        Error = 4,
        Off = 5
    };

    //!
    //! \brief Super simple logger implementation.
    //!
    //! This is a super simple logger implementation that has minimal logging
    //! capability. Currently, the class doesn't support multi-thread logging.
    //!
    class Logger final {
    public:
        //! Constructs a logger with logging level.
        explicit Logger(LoggingLevel level);

        //! Destructor.
        ~Logger();

        //! Writes a value to the buffer stream.
        template <typename T>
        const Logger& operator<<(const T& x) const {
            _buffer << x;
            return *this;
        }

    private:
        LoggingLevel _level;
        mutable std::stringstream _buffer;
    };

    //! Helper class for logging.
    class Logging {
    public:
        //! Sets the output stream for the info level logs.
        static void setInfoStream(std::ostream* strm);

        //! Sets the output stream for the warning level logs.
        static void setWarnStream(std::ostream* strm);

        //! Sets the output stream for the error level logs.
        static void setErrorStream(std::ostream* strm);

        //! Sets the output stream for the debug level logs.
        static void setDebugStream(std::ostream* strm);

        //! Sets the output stream for all the log levelss.
        static void setAllStream(std::ostream* strm);

        //! Returns the header string.
        static std::string getHeader(LoggingLevel level);

        //! Sets the logging level.
        static void setLevel(LoggingLevel level);

        //! Mutes the logger.
        static void mute();

        //! Un-mutes the logger.
        static void unmute();
    };

    //! Info-level logger.
    extern Logger infoLogger;

    //! Warn-level logger.
    extern Logger warnLogger;

    //! Error-level logger.
    extern Logger errorLogger;

    //! Debug-level logger.
    extern Logger debugLogger;

#define JET_INFO \
    (Logger(LoggingLevel::Info) << Logging::getHeader(LoggingLevel::Info) \
     << "[" << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")] ")
#define JET_WARN \
    (Logger(LoggingLevel::Warn) << Logging::getHeader(LoggingLevel::Warn) \
     << "[" << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")] ")
#define JET_ERROR \
    (Logger(LoggingLevel::Error) << Logging::getHeader(LoggingLevel::Error) \
     << "[" << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")] ")
#define JET_DEBUG \
    (Logger(LoggingLevel::Debug) << Logging::getHeader(LoggingLevel::Debug) \
     << "[" << __FILE__ << ":" << __LINE__ << " (" << __func__ << ")] ")


    static std::mutex critical;

    static std::ostream* infoOutStream = &std::cout;
    static std::ostream* warnOutStream = &std::cout;
    static std::ostream* errorOutStream = &std::cerr;
    static std::ostream* debugOutStream = &std::cout;
    static LoggingLevel sLoggingLevel = LoggingLevel::All;

    inline std::ostream* levelToStream(LoggingLevel level) {
        switch (level) {
        case LoggingLevel::Info:
            return infoOutStream;
        case LoggingLevel::Warn:
            return warnOutStream;
        case LoggingLevel::Error:
            return errorOutStream;
        case LoggingLevel::Debug:
            return debugOutStream;
        default:
            return infoOutStream;
        }
    }

    inline std::string levelToString(LoggingLevel level) {
        switch (level) {
        case LoggingLevel::Info:
            return "INFO";
        case LoggingLevel::Warn:
            return "WARN";
        case LoggingLevel::Error:
            return "ERROR";
        case LoggingLevel::Debug:
            return "DEBUG";
        default:
            return "";
        }
    }

    inline bool isLeq(LoggingLevel a, LoggingLevel b) {
        return (uint8_t)a <= (uint8_t)b;
    }

    Logger::Logger(LoggingLevel level) : _level(level) {}

    Logger::~Logger() {
        std::lock_guard<std::mutex> lock(critical);
        if (isLeq(sLoggingLevel, _level)) {
            auto strm = levelToStream(_level);
            (*strm) << _buffer.str() << std::endl;
            strm->flush();
        }
    }

    void Logging::setInfoStream(std::ostream* strm) {
        std::lock_guard<std::mutex> lock(critical);
        infoOutStream = strm;
    }

    void Logging::setWarnStream(std::ostream* strm) {
        std::lock_guard<std::mutex> lock(critical);
        warnOutStream = strm;
    }

    void Logging::setErrorStream(std::ostream* strm) {
        std::lock_guard<std::mutex> lock(critical);
        errorOutStream = strm;
    }

    void Logging::setDebugStream(std::ostream* strm) {
        std::lock_guard<std::mutex> lock(critical);
        debugOutStream = strm;
    }

    void Logging::setAllStream(std::ostream* strm) {
        setInfoStream(strm);
        setWarnStream(strm);
        setErrorStream(strm);
        setDebugStream(strm);
    }

    std::string Logging::getHeader(LoggingLevel level) {
        auto now =
            std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        char timeStr[20];
#ifdef JET_WINDOWS
        tm time;
        localtime_s(&time, &now);
#ifdef _MSC_VER
        strftime(timeStr, sizeof(timeStr), "%F %T", &time);
#else
        // Such as MinGW - https://sourceforge.net/p/mingw-w64/bugs/793/
        strftime(timeStr, sizeof(timeStr), "%Y-%m-%d %H:%M:%S", &time);
#endif
#else
        strftime(timeStr, sizeof(timeStr), "%F %T", std::localtime(&now));
#endif
        char header[256];
        snprintf(header, sizeof(header), "[%s] %s ", levelToString(level).c_str(),
            timeStr);
        return header;
    }

    void Logging::setLevel(LoggingLevel level) {
        std::lock_guard<std::mutex> lock(critical);
        sLoggingLevel = level;
    }

    void Logging::mute() { setLevel(LoggingLevel::Off); }

    void Logging::unmute() { setLevel(LoggingLevel::All); }
}  // namespace jet

#endif  // INCLUDE_JET_LOGGING_H_