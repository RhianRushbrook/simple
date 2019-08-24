#pragma once
#include <cmath>
namespace simple {

	namespace math {
		// Primitives
		constexpr float Pi = 3.141592741f;
		constexpr float Epsilon = 0.0001f;
		constexpr float Root2 = 1.41421356237f;

		inline float Sin(float rad) { return sinf(rad); }
		inline float Cos(float rad) { return cosf(rad); }
		inline float DegToRad(float deg) { return (deg * Pi / 180); }
		inline float RadToDeg(float rad) { return (rad * 180 / Pi); }
		inline float Min(float a, float b) {return ((a) < (b) ? (a) : (b)); }
		inline float Max(float a, float b) { return ((a) > (b) ? (a) : (b)); }
		inline float Abs(float val) { return ((val) < 0 ? -(val) : (val)); }
		inline float Clamp(float val, float min, float max) { return Max(min, Min(val, max)); }
		inline float Sign(float val) { return (val < 0 ? -1.0f : 1.0f); }
		inline float Sqr(float val) { return (val * val); }
		inline float Sqrt(float val) { return sqrt(val); }

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

			void Rotate(float radians) {
				float c = Cos(radians);
				float s = Sin(radians);

				float xp = x * c - y * s;
				float yp = x * s + y * c;

				x = xp;
				y = yp;
			}

			float LenSqr() const {
				return x * x + y * y;
			}

			float Len() const {
				return sqrt(x * x + y * y);
			}

			vec2 operator + (const vec2& v) const {
				return vec2(x + v.x, y + v.y);
			}

			vec2 operator + (float f) const {
				return vec2(x + f, y + f);
			}

			vec2 operator - (const vec2& v) const {
				return vec2(x - v.x, y - v.y);
			}

			vec2 operator * (const float& f) const {
				return vec2(x * f, y * f);
			}

			vec2 operator / (const float& f) const {
				return vec2(x / f, y / f);
			}

			vec2& operator += (const vec2& v) {
				this->x += v.x;	this->y += v.y;
				return *this;
			}

			vec2& operator -= (const vec2& v) {
				this->x -= v.x;	this->y -= v.y;
				return *this;
			}

			vec2& operator *= (const float f) {
				this->x *= f;	this->y *= f;
				return *this;
			}

			vec2& operator /= (const float f) {
				this->x /= f;	this->y /= f;
				return *this;
			}
		};

		// Vector ops
		inline vec2 VecAbs(vec2 v) { return vec2(Abs(v.x), Abs(v.y)); }
		inline vec2 VecMin(vec2 a, vec2 b) { return vec2(Min(a.x, b.x), Min(a.y, b.y)); }
		inline vec2 VecMax(vec2 a, vec2 b) { return vec2(Max(a.x, b.x), Max(a.y, b.y)); }
		inline vec2 VecClamp(vec2 v, vec2 min, vec2 max) { return VecMin(max, VecMax(min, v)); }

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


		inline float VecDot(const vec2& a, const vec2& b) {
			return (a.x * b.x) + (a.y * b.y);
		};

		inline float VecSqr(const vec2& v) {
			return (v.x * v.x) + (v.y * v.y);
		};

		inline vec2 VecNorm(const vec2& v) {
			if (v.Len() > Epsilon) {
				float r = 1 / v.Len();
				return v * r;
			}
			return vec2();
		};

	}


	namespace camera {
	// TO-DO
	// Local to World ect
	}

	namespace structure {

		struct rotation {
			float cos;
			float sin;
		};

		struct transform {
			math::vec2 pos;
			rotation rot;
		};
	}

	namespace collision {

		constexpr int MaxVertexCount = 8;

		namespace body {
			struct halfspace {
				math::vec2 direction;
				float distance;
			};

			struct circle {
				math::vec2 position;
				float radius;
			};

			struct box {
				math::vec2 min;
				math::vec2 max;
			};

			struct capsule {
				math::vec2 start;
				math::vec2 end;
				float radius;
			};

			struct poly {
				int nVertex;
				math::vec2 vertex[MaxVertexCount];
				math::vec2 norms[MaxVertexCount];
			};
		};

		struct manifold {
			int count;
			float depths[2];

			math::vec2 contact_points[2];
		};

		
		// Box ops

		body::poly BoxToPoly(body::box& box) {
			body::poly out;
			out.nVertex = 4;
			out.vertex[0] = box.min;
			out.vertex[1] = math::vec2(box.max.x, box.min.y);
			out.vertex[2] = box.max;
			out.vertex[3] = math::vec2(box.min.x, box.max.y);
			return out;
		};


		// Poly math

		inline void poly_add(body::poly& poly, float x, float y) {
			for (int i = 0; i < poly.nVertex; i++) {
				poly.vertex[i].x += x;
				poly.vertex[i].y += y;
			}
		}


		int PointVsCircle(math::vec2& point, body::circle& circle);
		int PointVsBox(math::vec2& point, body::box& box);
		int PointVsCapsule(math::vec2& point, body::capsule& capsule);
		int PointVsPoly(math::vec2& point, body::poly& poly);

		int CircleVsCircle(body::circle& circle1, body::circle& circle2);
		int CircleVsBox(body::circle& circle, body::box& box);
		int CircleVsCapsule(body::circle& circle, body::capsule& capsule);
		int CircleVsPoly(body::circle& circle, body::poly& poly);

		int BoxVsBox(body::box& a, body::box& b);
		int BoxVsCapsule(body::box& box, body::capsule& capsule);
		int BoxVsPoly(body::box& box, body::poly& poly);


#ifndef SIMPLE_COLLISION_FUNCTIONS

		int PointVsCircle(math::vec2& point, body::circle& circle) {
			math::vec2 adjusted_point = point - circle.position;
			float point_distance = math::VecSqr(adjusted_point);
			float circle_radius = circle.radius * circle.radius;

			return (point_distance <= circle_radius);
		}

		int PointVsBox(math::vec2& point, body::box& box) {
			return (point.x >= box.min.x && point.x <= box.max.x
				 && point.y >= box.min.y && point.y <= box.max.y);
		}

		int PointVsCapsule(math::vec2& point, body::capsule& capsule) {

			float length_sqr = math::VecSqr(capsule.end - capsule.start);
			if (length_sqr == 0.0f) 
				return math::VecSqr( point - capsule.start) <= math::Sqr(capsule.radius);
			float t = math::Clamp(math::VecDot(point - capsule.start, capsule.end - capsule.start) / length_sqr, 0.0f, 1.0f);
			math::vec2 projection = capsule.start + ((capsule.end - capsule.start) * t);
			float dist_sqr = math::VecSqr(projection - point);

			return dist_sqr <= math::Sqr(capsule.radius);
		}

		int PointVsPoly(math::vec2& point, body::poly& poly) {
			for (int i = 0; i < poly.nVertex; i++) {
				math::vec2 start = poly.vertex[i];
				math::vec2 end = poly.vertex[(i + 1) % poly.nVertex];

				math::vec2 face_direction = math::VecNorm(end - start);
				math::vec2 face_normal = { face_direction.y * -1, face_direction.x };

				float point_along_normal = math::VecDot(point, face_normal);
				float face_along_normal = math::VecDot(start, face_normal);

				if (point_along_normal <= face_along_normal)
					return 0;
			}
			return 1;
		}


		// Circle collisions

		int CircleVsCircle(body::circle& circle1, body::circle& circle2) {
			math::vec2 local_space = circle2.position - circle1.position;
			float circle_distance = math::VecSqr(local_space);
			float combined_radius = math::Sqr(circle1.radius + circle2.radius);

			return (circle_distance <= combined_radius);
		}


		int CircleVsBox(body::circle& circle, body::box& box) {
			math::vec2 a = circle.position;
			math::vec2 b = math::VecClamp(a, box.min, box.max);
			math::vec2 ab = b - a;

			return math::VecSqr(ab) <= math::Sqr(circle.radius);
		}

		int CircleVsCapsule(body::circle& circle, body::capsule& capsule) {
			float length_sqr = math::VecSqr(capsule.end - capsule.start);
			if (length_sqr == 0.0f)
				return math::VecSqr(circle.position - capsule.start) <= math::Sqr(capsule.radius + circle.radius);
			float t = math::Clamp(math::VecDot(circle.position - capsule.start, capsule.end - capsule.start) / length_sqr, 0.0f, 1.0f);
			math::vec2 projection = capsule.start + ((capsule.end - capsule.start) * t);
			float dist_sqr = math::VecSqr(projection - circle.position);

			return dist_sqr <= math::Sqr(capsule.radius + circle.radius);
		}

		int CircleVsPoly(body::circle& circle, body::poly& poly) {
			for (int i = 0; i < poly.nVertex; i++) {
				math::vec2 start = poly.vertex[i];
				math::vec2 end = poly.vertex[(i + 1) % poly.nVertex];

				math::vec2 face_direction = math::VecNorm(end - start);
				math::vec2 face_normal = { face_direction.y * -1, face_direction.x };

				float point_along_normal = math::VecDot(circle.position, face_normal);
				float face_along_normal = math::VecDot(start, face_normal);

				if (point_along_normal + circle.radius <= face_along_normal)
					return 0;
			}
			return 1;
		}


		// box collisions

		int BoxVsBox(body::box& a, body::box& b) {
			float mLeft = a.min.x - b.max.x;
			float mRight = a.max.x - b.min.x;
			float mTop = a.min.y - b.max.y;
			float mBot = a.max.y - b.min.y;

			return (mLeft <= 0 && mRight >= 0 && mTop <= 0 && mBot >= 0);
		}

		// util

		void build_regular_poly(int nSides, float side_length, body::poly& out) {

			float ang = (2 * math::Pi) / nSides;

			out.nVertex = nSides;

			for (int i = 0; i < nSides; i++) {
				out.vertex[i] = { side_length * math::Cos(ang * i), side_length * math::Sin(ang * i) };
			}
		}

#define SIMPLE_COLLISION_FUNCTIONS
#endif
	}

}
