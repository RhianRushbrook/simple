#pragma once
#include "simple_math.h"

namespace simple {

	namespace collision {

		constexpr int MaxVertexCount = 8;

		namespace body {

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


		struct rotation {
			float cos;
			float sin;
		};

		struct transform {
			math::vec2 pos;
			rotation rot;
		};

		struct ray {
			math::vec2 position;
			math::vec2 direction;
			float max_distance;
		};

		struct manifold {
			int count;
			math::vec2 contact_normal;
			math::vec2 contact_point;
			float distance;
		};

		struct simplex_vert {
			math::vec2 support1;
			math::vec2 support2;
			math::vec2 vector = support2 - support1;
			float u;
			int support_index1;
			int support_index2;
		};

		struct simplex {
			int simplex_count;

			simplex_vert vert_1, vert_2, vert_3;

			math::vec2 GetSearchDirection() const;
		};


		// Structure ops

		rotation DegToRotation(float ang);
		math::vec2 ApplyTransform(const math::vec2& v, const transform& t);


		// Bool collision check

		int RayVsCircle(const ray& ray, const body::circle& circle, manifold& out);
		int RayVsBox(const ray& ray, const body::box& box, manifold& out);
		int RayVsCapsule(const ray& ray, const body::capsule& capsule, manifold& out);
		int RayVsPoly(const ray& ray, const body::poly& poly, manifold& out);

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

		int PolyVsPoly(body::poly& p1, body::poly& p2);


		// Temp poly utility

		void build_regular_poly(int nSides, float side_length, body::poly& out);
		body::poly BoxToPoly(simple::collision::body::box& box);
		void poly_add(body::poly& poly, float x, float y);

	} // Collision end
}