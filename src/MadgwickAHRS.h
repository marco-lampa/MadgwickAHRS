//=============================================================================================
// MadgwickAHRS.h
//=============================================================================================
//
// Implementation of Madgwick's IMU and AHRS algorithms.
// See: http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/
//
// From the x-io website "Open-source resources available on this website are
// provided under the GNU General Public Licence unless an alternative licence
// is provided in source."
//
// Date			Author          Notes
// 29/09/2011	SOH Madgwick    Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
//
//=============================================================================================
#ifndef MadgwickAHRS_h
#define MadgwickAHRS_h
#include <math.h>

//--------------------------------------------------------------------------------------------
// Variable declaration
class Madgwick{
public:
// Helper types for orientation return
    struct OrientationQuat {
        float x{0.f};  // Often called q1
        float y{0.f};  // Often called q2
        float z{0.f};  // Often called q3
        float w{1.f};  // Often called q0
        void normalize()
        {
            const float recipNorm = invSqrt(w * w + x * x + y * y + z * z);
	        w *= recipNorm;
	        x *= recipNorm;
	        y *= recipNorm;
	        z *= recipNorm;
        }
    };

    struct OrientationRPY {
        const float roll;
        const float pitch;
        const float yaw;

        OrientationRPY() = delete;
        OrientationRPY(const OrientationQuat& q) :
            roll{atan2f(q.w*q.x + q.y*q.z, 0.5f - q.x*q.x - q.y*q.y)},
	        pitch{asinf(-2.0f * (q.x*q.z - q.w*q.y))},
	        yaw{atan2f(q.x*q.y + q.w*q.z, 0.5f - q.y*q.y - q.z*q.z)}
        {}
        OrientationRPY(const float r, const float p, const float y):
            roll{r}, pitch{p}, yaw{y}
        {}
    };
private:
    static constexpr float radToDeg = 57.29578f;
    float beta;				// algorithm gain
    float invSampleFreq;
    OrientationQuat quat;

    static inline float invSqrt(float x) {
        // Fast inverse square-root
        // See: http://en.wikipedia.org/wiki/Fast_inverse_square_root
        float halfx = 0.5f * x;
        float y = x;
        long i = *(long*)&y;
        i = 0x5f3759df - (i>>1);
        y = *(float*)&i;
        y = y * (1.5f - (halfx * y * y));
        y = y * (1.5f - (halfx * y * y));
        return y;
    }

//-------------------------------------------------------------------------------------------
// Function declarations
public:
    Madgwick(void);
    void begin(float sampleFrequency) { invSampleFreq = 1.0f / sampleFrequency; }
    void update(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz);
    void updateIMU(float gx, float gy, float gz, float ax, float ay, float az);
    OrientationRPY getOrientationRPY(const bool deg = false)
    {
        static const float radToDeg = 57.29578f;
        const OrientationRPY rpyRad = OrientationRPY(quat);
        if (deg) {
            return OrientationRPY(rpyRad.roll * radToDeg, rpyRad.pitch * radToDeg, rpyRad.yaw * radToDeg);
        }
        return rpyRad;
    }
    OrientationQuat getOrientationQuat()
    {
        return quat;
    }
};
#endif
