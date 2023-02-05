//=============================================================================================
// MadgwickAHRS.c
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
// 19/02/2012	SOH Madgwick	Magnetometer measurement is normalised
//
//=============================================================================================

//-------------------------------------------------------------------------------------------
// Header files

#include "MadgwickAHRS.h"
#include <math.h>

//-------------------------------------------------------------------------------------------
// Definitions

#define sampleFreqDef   512.0f          // sample frequency in Hz
#define betaDef         0.1f            // 2 * proportional gain


//============================================================================================
// Functions

float invSqrt(float x)
{
	return 1.f / sqrtf(x);
}

//-------------------------------------------------------------------------------------------
// AHRS algorithm update

Madgwick::Madgwick() {
	beta = betaDef;
	quat = OrientationQuat();
	invSampleFreq = 1.0f / sampleFreqDef;
}

void Madgwick::update(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz) {
	float recipNorm;
	float s0, s1, s2, s3;
	float qDot1, qDot2, qDot3, qDot4;
	float hx, hy;
	float _2q0mx, _2q0my, _2q0mz, _2q1mx, _2bx, _2bz, _4bx, _4bz, _2q0, _2q1, _2q2, _2q3, _2q0q2, _2q2q3, q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;

	// Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
	if((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f)) {
		updateIMU(gx, gy, gz, ax, ay, az);
		return;
	}

	// Convert gyroscope degrees/sec to radians/sec
	gx *= 0.0174533f;
	gy *= 0.0174533f;
	gz *= 0.0174533f;

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-quat.x * gx - quat.y * gy - quat.z * gz);
	qDot2 = 0.5f * (quat.w * gx + quat.y * gz - quat.z * gy);
	qDot3 = 0.5f * (quat.w * gy - quat.x * gz + quat.z * gx);
	qDot4 = 0.5f * (quat.w * gz + quat.x * gy - quat.y * gx);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;

		// Normalise magnetometer measurement
		recipNorm = invSqrt(mx * mx + my * my + mz * mz);
		mx *= recipNorm;
		my *= recipNorm;
		mz *= recipNorm;

		// Auxiliary variables to avoid repeated arithmetic
		_2q0mx = 2.0f * quat.w * mx;
		_2q0my = 2.0f * quat.w * my;
		_2q0mz = 2.0f * quat.w * mz;
		_2q1mx = 2.0f * quat.x * mx;
		_2q0 = 2.0f * quat.w;
		_2q1 = 2.0f * quat.x;
		_2q2 = 2.0f * quat.y;
		_2q3 = 2.0f * quat.z;
		_2q0q2 = 2.0f * quat.w * quat.y;
		_2q2q3 = 2.0f * quat.y * quat.z;
		q0q0 = quat.w * quat.w;
		q0q1 = quat.w * quat.x;
		q0q2 = quat.w * quat.y;
		q0q3 = quat.w * quat.z;
		q1q1 = quat.x * quat.x;
		q1q2 = quat.x * quat.y;
		q1q3 = quat.x * quat.z;
		q2q2 = quat.y * quat.y;
		q2q3 = quat.y * quat.z;
		q3q3 = quat.z * quat.z;

		// Reference direction of Earth's magnetic field
		hx = mx * q0q0 - _2q0my * quat.z + _2q0mz * quat.y + mx * q1q1 + _2q1 * my * quat.y + _2q1 * mz * quat.z - mx * q2q2 - mx * q3q3;
		hy = _2q0mx * quat.z + my * q0q0 - _2q0mz * quat.x + _2q1mx * quat.y - my * q1q1 + my * q2q2 + _2q2 * mz * quat.z - my * q3q3;
		_2bx = sqrtf(hx * hx + hy * hy);
		_2bz = -_2q0mx * quat.y + _2q0my * quat.x + mz * q0q0 + _2q1mx * quat.z - mz * q1q1 + _2q2 * my * quat.z - mz * q2q2 + mz * q3q3;
		_4bx = 2.0f * _2bx;
		_4bz = 2.0f * _2bz;

		// Gradient decent algorithm corrective step
		s0 = -_2q2 * (2.0f * q1q3 - _2q0q2 - ax) + _2q1 * (2.0f * q0q1 + _2q2q3 - ay) - _2bz * quat.y * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * quat.z + _2bz * quat.x) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * quat.y * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s1 = _2q3 * (2.0f * q1q3 - _2q0q2 - ax) + _2q0 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * quat.x * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + _2bz * quat.z * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * quat.y + _2bz * quat.w) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * quat.z - _4bz * quat.x) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s2 = -_2q0 * (2.0f * q1q3 - _2q0q2 - ax) + _2q3 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * quat.y * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + (-_4bx * quat.y - _2bz * quat.w) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * quat.x + _2bz * quat.z) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * quat.w - _4bz * quat.y) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s3 = _2q1 * (2.0f * q1q3 - _2q0q2 - ax) + _2q2 * (2.0f * q0q1 + _2q2q3 - ay) + (-_4bx * quat.z + _2bz * quat.x) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * quat.w + _2bz * quat.y) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * quat.x * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= beta * s0;
		qDot2 -= beta * s1;
		qDot3 -= beta * s2;
		qDot4 -= beta * s3;
	}

	// Integrate rate of change of quaternion to yield quaternion
	quat.w += qDot1 * invSampleFreq;
	quat.x += qDot2 * invSampleFreq;
	quat.y += qDot3 * invSampleFreq;
	quat.z += qDot4 * invSampleFreq;

	// Normalise quaternion
	quat.normalize();
}

//-------------------------------------------------------------------------------------------
// IMU algorithm update

void Madgwick::updateIMU(float gx, float gy, float gz, float ax, float ay, float az) {
	float recipNorm;
	float s0, s1, s2, s3;
	float qDot1, qDot2, qDot3, qDot4;
	float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3;

	// Convert gyroscope degrees/sec to radians/sec
	gx *= 0.0174533f;
	gy *= 0.0174533f;
	gz *= 0.0174533f;

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-quat.x * gx - quat.y * gy - quat.z * gz);
	qDot2 = 0.5f * (quat.w * gx + quat.y * gz - quat.z * gy);
	qDot3 = 0.5f * (quat.w * gy - quat.x * gz + quat.z * gx);
	qDot4 = 0.5f * (quat.w * gz + quat.x * gy - quat.y * gx);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;

		// Auxiliary variables to avoid repeated arithmetic
		_2q0 = 2.0f * quat.w;
		_2q1 = 2.0f * quat.x;
		_2q2 = 2.0f * quat.y;
		_2q3 = 2.0f * quat.z;
		_4q0 = 4.0f * quat.w;
		_4q1 = 4.0f * quat.x;
		_4q2 = 4.0f * quat.y;
		_8q1 = 8.0f * quat.x;
		_8q2 = 8.0f * quat.y;
		q0q0 = quat.w * quat.w;
		q1q1 = quat.x * quat.x;
		q2q2 = quat.y * quat.y;
		q3q3 = quat.z * quat.z;

		// Gradient decent algorithm corrective step
		s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
		s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * quat.x - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
		s2 = 4.0f * q0q0 * quat.y + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
		s3 = 4.0f * q1q1 * quat.z - _2q1 * ax + 4.0f * q2q2 * quat.z - _2q2 * ay;
		recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= beta * s0;
		qDot2 -= beta * s1;
		qDot3 -= beta * s2;
		qDot4 -= beta * s3;
	}

	// Integrate rate of change of quaternion to yield quaternion
	quat.w += qDot1 * invSampleFreq;
	quat.x += qDot2 * invSampleFreq;
	quat.y += qDot3 * invSampleFreq;
	quat.z += qDot4 * invSampleFreq;

	// Normalise quaternion
	recipNorm = invSqrt(quat.w * quat.w + quat.x * quat.x + quat.y * quat.y + quat.z * quat.z);
	quat.w *= recipNorm;
	quat.x *= recipNorm;
	quat.y *= recipNorm;
	quat.z *= recipNorm;
}
