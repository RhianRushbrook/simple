#pragma once
#include <cmath>
namespace simple {

	namespace math {
		// Math ops
		constexpr float Pi = 3.14159265358979323846f;
		constexpr float Epsilon = 0.0001f;
		constexpr float Root2 = 1.41421356237309504880f;

		inline float Sin(float rad) { return std::sinf(rad); }
		inline float Cos(float rad) { return std::cosf(rad); }
		inline float DegToRad(float deg) { return (deg * Pi / 180); }
		inline float RadToDeg(float rad) { return (rad * 180 / Pi); }
		inline float Min(float a, float b) {return ((a) < (b) ? (a) : (b)); }
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

			vec2  operator +  (const vec2&  v) const { return vec2(x + v.x, y + v.y); }
			vec2  operator +  (const float& f) const { return vec2(x + f, y + f); }
			vec2  operator -  (const vec2&  v) const { return vec2(x - v.x, y - v.y); }
			vec2  operator -  (const float& f) const { return vec2(x - f, y - f); }
			vec2  operator *  (const float& f) const { return vec2(x * f, y * f); }
			vec2  operator /  (const float& f) const { return vec2(x / f, y / f); }
			vec2& operator += (const vec2& v) { x += v.x;	y += v.y; return *this; }
			vec2& operator -= (const vec2& v) { x -= v.x;	y -= v.y; return *this; }
			vec2& operator *= (const float f) { x *= f;	y *= f;	return *this; }
			vec2& operator /= (const float f) { x /= f;	y /= f;	return *this; }
		};

		// Vector ops
		inline float VecDot(const vec2& a, const vec2& b) { return (a.x * b.x) + (a.y * b.y); }
		inline float VecSqr(const vec2& v) { return (v.x * v.x) + (v.y * v.y); }

		inline vec2 VecAbs(vec2 v) { return vec2(Abs(v.x), Abs(v.y)); }
		inline vec2 VecMin(vec2 a, vec2 b) { return vec2(Min(a.x, b.x), Min(a.y, b.y)); }
		inline vec2 VecMax(vec2 a, vec2 b) { return vec2(Max(a.x, b.x), Max(a.y, b.y)); }
		inline vec2 VecClamp(vec2 v, vec2 min, vec2 max) { return VecMin(max, VecMax(min, v)); }
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


	namespace time {
		float DeltaTime();
	}


	namespace camera {

		struct view_port {
			math::vec2 offset;
			float width, height;
			float scale = 1;

			math::vec2 Min() const {
				return offset;
			}
			math::vec2 Max() const {
				return (math::vec2(width, height) / scale) + offset;
			}

			void Set(const math::vec2& vec);
			void Set(float x, float y);

			void Pan(const math::vec2& vec);
			void Pan(float x, float y);

			void Zoom(float val, const math::vec2& vec = { 0.0f, 0.0f });
			void Zoom(float val, int x = 0, int y = 0);
		};

		math::vec2 WorldToPort(const view_port& port, const math::vec2& vec);
		math::vec2 WorldToPort(const view_port& port, const float x, const float y);

		math::vec2 PortToWorld(const view_port& port, const math::vec2& vec);
		math::vec2 PortToWorld(const view_port& port, const float x, const float y);

	} // Camera end


	namespace structure {
		// Move to math?
		struct rotation {
			float cos;
			float sin;
		};

		struct transform {
			math::vec2 pos;
			rotation rot;
		};
	} // Structure end

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
				float radius = 0.0f;
			};

			struct poly {
				int nVertex = 0;
				math::vec2 vertex[MaxVertexCount];
				math::vec2 norms[MaxVertexCount];
			};
		}

		struct manifold {
			int count;
			float depths[2];

			math::vec2 contact_points[2];
		};


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
		
		int CapsuleVsCapsule(body::capsule& a, body::capsule& b);
		int CapsuleVsPoly(body::capsule& capsule, body::poly& poly);

		int PolyVsPoly(body::poly& a, body::poly& b);
	} // Collision end

} // simple end


namespace simple {

	namespace camera {

		math::vec2 WorldToPort(const view_port& port, const math::vec2& vec) {
			return (vec - port.offset) * port.scale;
		}

		math::vec2 WorldToPort(const view_port& port, const float x, const float y) {
			return (math::vec2(x, y) - port.offset) * port.scale;
		}

		math::vec2 PortToWorld(const view_port& port, const math::vec2& vec) {
			return (vec / port.scale) + port.offset;
		}

		math::vec2 PortToWorld(const view_port& port, const float x, const float y) {
			return (math::vec2(x, y) / port.scale) + port.offset;
		}

		void view_port::Set(const math::vec2& vec) {
			offset = vec;
		}
		
		void view_port::Set(float x, float y) {
			offset.x = x;
			offset.y = y;
		}

		void view_port::Pan(const math::vec2& vec) {
			offset += vec / scale;
		}
		
		void view_port::Pan(float x, float y) {
			offset.x += x / scale;
			offset.y += y / scale;
		}

		void view_port::Zoom(float val, const math::vec2& vec) {
			math::vec2 initial = PortToWorld(*this, vec);
			scale *= val;
			math::vec2 finish = PortToWorld(*this, vec);
			Set(offset + (initial - finish));
		}

		void view_port::Zoom(float val, int x, int y) {
			math::vec2 initial = PortToWorld(*this, x, y);
			scale *= val;
			math::vec2 finish = PortToWorld(*this, x, y);
			Set(offset + (initial - finish));
		}
		
	} // Camera end


	namespace collision {

		int PointVsCircle(math::vec2& point, body::circle& circle) {
			float sqr_distance = math::VecSqr(point - circle.position);

			return (sqr_distance <= math::Sqr(circle.radius));
		}

		int PointVsBox(math::vec2& point, body::box& box) {
			return (point.x >= box.min.x && point.x <= box.max.x
				&& point.y >= box.min.y && point.y <= box.max.y);
		}

		int PointVsCapsule(math::vec2& point, body::capsule& capsule) {
			float sqr_length = math::VecSqr(capsule.end - capsule.start);
			if (sqr_length == 0.0f)
				return (math::VecSqr(point - capsule.start) <= math::Sqr(capsule.radius));

			float t = math::Clamp(math::VecDot(point - capsule.start, capsule.end - capsule.start) / sqr_length, 0.0f, 1.0f);
			math::vec2 projection = capsule.start + ((capsule.end - capsule.start) * t);
			float sqr_distance = math::VecSqr(projection - point);

			return (sqr_distance <= math::Sqr(capsule.radius));
		}

		int PointVsPoly(math::vec2& point, body::poly& poly) {
			for (int i = 0; i < poly.nVertex; i++) {
				float point_along_normal = math::VecDot(point, poly.norms[i]);
				float face_along_normal = math::VecDot(poly.vertex[i], poly.norms[i]);

				if (point_along_normal <= face_along_normal)
					return 0;
			}
			return 1;
		}


		// Circle collisions

		int CircleVsCircle(body::circle& circle1, body::circle& circle2) {
			float sqr_distance = math::VecSqr(circle2.position - circle1.position);
			float sqr_radius = math::Sqr(circle1.radius + circle2.radius);

			return (sqr_distance <= sqr_radius);
		}

		int CircleVsBox(body::circle& circle, body::box& box) {
			math::vec2 projected_circle = math::VecClamp(circle.position, box.min, box.max);
			float sqr_distance = math::VecSqr(projected_circle - circle.position);

			return (sqr_distance <= math::Sqr(circle.radius));
		}

		int CircleVsCapsule(body::circle& circle, body::capsule& capsule) {
			float sqr_length = math::VecSqr(capsule.end - capsule.start);
			float sqr_radius = math::Sqr(capsule.radius + circle.radius);
			if (sqr_length == 0.0f)
				return (math::VecSqr(circle.position - capsule.start) <= sqr_radius);

			float t = math::Clamp(math::VecDot(circle.position - capsule.start, capsule.end - capsule.start) / sqr_length, 0.0f, 1.0f);
			math::vec2 projection = capsule.start + ((capsule.end - capsule.start) * t);
			float sqr_distance = math::VecSqr(projection - circle.position);

			return sqr_distance <= sqr_radius;
		}


		// box collisions

		int BoxVsBox(body::box& a, body::box& b) {
			float mLeft = a.min.x - b.max.x;
			float mRight = a.max.x - b.min.x;
			float mTop = a.min.y - b.max.y;
			float mBot = a.max.y - b.min.y;

			return (mLeft <= 0 && mRight >= 0 && mTop <= 0 && mBot >= 0);
		}


		// Temp util

		void build_regular_poly(int nSides, float side_length, body::poly& out) {

			float ang = (2 * math::Pi) / nSides;

			out.nVertex = nSides;

			for (int i = 0; i < nSides; i++) {
				out.vertex[i] = { side_length * math::Cos(ang * i), side_length * math::Sin(ang * i) };
			}
			for (int i = 0; i < nSides; i++) {
				math::vec2 plane = out.vertex[(i + 1) % nSides] - out.vertex[i];
				out.norms[i] = { -plane.y, plane.x };
		}	}

		body::poly BoxToPoly(simple::collision::body::box& box) {
			simple::collision::body::poly out;
			out.nVertex = 4;
			out.vertex[0] = box.min;
			out.vertex[1] = simple::math::vec2(box.max.x, box.min.y);
			out.vertex[2] = box.max;
			out.vertex[3] = simple::math::vec2(box.min.x, box.max.y);
			return out;
		}

		inline void poly_add(body::poly& poly, float x, float y) {
			for (int i = 0; i < poly.nVertex; i++) {
				poly.vertex[i].x += x;
				poly.vertex[i].y += y;
		}	}
		
	} // Collision end

} // Simple end
