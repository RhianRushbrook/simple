#pragma once
#include <cmath>

namespace simple {

	namespace math {
		// Math ops
		constexpr float Pi = 3.14159265358979323846f;
		constexpr float Tau = 6.28318530717958647692f;
		constexpr float Epsilon = 0.0001f;
		constexpr float Root2 = 1.41421356237309504880f;

		inline float Sin(float rad) { return std::sinf(rad); }
		inline float Cos(float rad) { return std::cosf(rad); }
		inline float DegToRad(float deg) { return (deg * Tau / 360); }
		inline float RadToDeg(float rad) { return (rad * 360 / Tau); }
		inline float Min(float a, float b) { return ((a) < (b) ? (a) : (b)); }
		inline float Max(float a, float b) { return ((a) > (b) ? (a) : (b)); }
		inline float Abs(float val) { return ((val) < 0 ? -(val) : (val)); }
		inline float Clamp(float val, float min, float max) { return Max(min, Min(val, max)); }
		inline float Sign(float val) { return (val < 0 ? -1.0f : 1.0f); }
		inline float Sqr(float val) { return (val * val); }
		inline float Sqrt(float val) { return std::sqrt(val); }

		// Structures
		struct vec2 {

			float x = 0.0f;
			float y = 0.0f;

			vec2() : x(0.0f), y(0.0f) {	}
			vec2(float x_, float y_) : x(x_), y(y_) {	}
			vec2(const vec2& v) : x(v.x), y(v.y) {	}

			void Set(float x_, float y_) {
				x = x_;
				y = y_;
			}

			float Len() const {
				return Sqrt((x * x) + (y * y));
			}

			vec2 operator -  () const { return vec2(-x, -y); }
			vec2 operator +  (const vec2& v) const { return vec2(x + v.x, y + v.y); }
			vec2 operator +  (const float& f) const { return vec2(x + f, y + f); }
			vec2 operator -  (const vec2& v) const { return vec2(x - v.x, y - v.y); }
			vec2 operator -  (const float& f) const { return vec2(x - f, y - f); }
			vec2 operator *  (const float& f) const { return vec2(x * f, y * f); }
			vec2 operator *  (const vec2& v)  const { return vec2(x * v.x, y * v.y); }
			vec2 operator /  (const float& f) const { return vec2(x / f, y / f); }
			//vec2 operator ^  (const float& i) const { return vec2(std::powf(x, i), std::pow(y,  i)); }
			void operator += (const vec2& v) { x += v.x; y += v.y; }
			void operator -= (const vec2& v) { x -= v.x; y -= v.y; }
			void operator *= (const float& f) { x *= f; y *= f; }
			void operator *= (const vec2& v) { x *= v.x; y *= v.y; }
			void operator /= (const float& f) { x /= f; y /= f; }
		};


		// Vector ops
		inline float VecDot(const vec2& a, const vec2& b) { return (a.x * b.x) + (a.y * b.y); }
		inline float VecSqr(const vec2& v) { return (v.x * v.x) + (v.y * v.y); }
		inline vec2 VecCross(const vec2 v) { return vec2(v.y, -v.x); }
		inline vec2 VecAbs(const vec2& v) { return vec2(Abs(v.x), Abs(v.y)); }
		inline vec2 VecMin(const vec2& a, const vec2& b) { return vec2(Min(a.x, b.x), Min(a.y, b.y)); }
		inline vec2 VecMax(const vec2& a, const vec2& b) { return vec2(Max(a.x, b.x), Max(a.y, b.y)); }
		inline vec2 VecClamp(const vec2& v, const vec2& min, const vec2& max) { return VecMin(max, VecMax(min, v)); }
		inline vec2 VecNorm(const vec2& v) {
			if (v.Len() < Epsilon) { return vec2(); }

			float r = 1 / v.Len();
			return v * r;
		}


		struct mat3 {

			float val[3][3];

			mat3() {
				val[0][0] = 1; val[1][0] = 0; val[2][0] = 0;
				val[0][1] = 0; val[1][1] = 1; val[2][1] = 0;
				val[0][2] = 0; val[1][2] = 0; val[2][2] = 1;
			}
			mat3(const mat3& m) {
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						val[i][j] = m.val[i][j];
			}

			mat3 operator + (const mat3& m) const {
				mat3 out;
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						out.val[i][j] = val[i][j] + m.val[i][j];
				return out;
			}
			mat3 operator - (const mat3& m) const {
				mat3 out;
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						out.val[i][j] = val[i][j] - m.val[i][j];
				return out;
			}
			mat3 operator * (const float& f) const {
				mat3 out;
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						out.val[i][j] = val[i][j] * f;
				return out;
			}
			mat3 operator * (const mat3& m) const {
				mat3 out;
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++) {
						float v = val[0][j] * m.val[i][0] +
							val[1][j] * m.val[i][1] +
							val[2][j] * m.val[i][2];
						out.val[i][j] = v;
					}
				return out;
			}

		};

	} // Math end
}
