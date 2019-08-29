#include "simple_collision.h"

namespace simple {

	namespace collision {

		//math::vec2 simplex::GetSearchDirection() const {
		//	switch (simplex_count) {
		//	case 1:
		//		return -vert_1.vector;
		//	case 2:
		//		return math::vec2();
		//
		//	}
		//}
		rotation DegToRotation(float ang) {
			rotation out = { math::Cos(ang), math::Sin(ang) };
			return out;
		}

		math::vec2 ApplyTransform(const math::vec2& v, const transform& t) {
			float xx = (v.x * t.rot.cos) + (v.y * t.rot.sin);
			float yy = (v.x * t.rot.sin) - (v.y * t.rot.cos);
			return math::vec2(xx + t.pos.x, yy + t.pos.y);
		}


		// Circle collision

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


		// Poly collisions

		int PolyVsPoly(body::poly& p1, body::poly& p2) {

			body::poly* a = &p1;
			body::poly* b = &p2;

			for (int shape = 0; shape < 2; shape++) {

				if (shape == 1) {
					a = &p2;
					b = &p1;
				}

				for (unsigned i = 0; i < a->nVertex; i++) {

					math::vec2 EdgeNormal = a->norms[i];

					float aMin = INFINITY, aMax = -INFINITY;
					for (unsigned j = 0; j < a->nVertex; j++) {
						float Dot = math::VecDot(a->vertex[j], EdgeNormal);
						aMax = math::Max(aMax, Dot);
						aMin = math::Min(aMin, Dot);
					}

					float bMin = INFINITY, bMax = -INFINITY;
					for (unsigned k = 0; k < b->nVertex; k++) {
						float Dot = math::VecDot(b->vertex[k], EdgeNormal);
						bMax = math::Max(bMax, Dot);
						bMin = math::Min(bMin, Dot);
					}

					if (!(bMax >= aMin && aMax >= bMin))
						return 0;
				}
			}
			return 1;
		}


		// Ray collision 

		int RayVsCircle(const ray& ray, const body::circle& circle, manifold& out) {
			math::vec2 local_space = ray.position - circle.position;
			float sqr_distance = math::VecSqr(local_space) - math::Sqr(circle.radius);
			float b = math::VecDot(local_space, ray.direction);
			float disc = math::Sqr(b) - sqr_distance;

			if (disc < 0)
				return 0;

			float collision_distance = -b - math::Sqrt(disc);

			if (collision_distance >= 0 && collision_distance <= ray.max_distance) {
				out.count = 1;
				out.contact_point = ray.position + (ray.direction * collision_distance);
				out.contact_normal = math::VecNorm(out.contact_point - circle.position);
				return 1;
			}
			return 0;
		}

		int RayVsBox(const ray& ray, const body::box& box, manifold& out) {
			math::vec2 inv = { 1 / ray.direction.x, 1 / ray.direction.y };
			math::vec2 d0 = (box.min - ray.position) * inv;
			math::vec2 d1 = (box.max - ray.position) * inv;
			math::vec2 v0 = math::VecMin(d0, d1);
			math::vec2 v1 = math::VecMax(d0, d1);
			float lo = math::Max(v0.x, v0.y);
			float hi = math::Min(v1.x, v1.y);

			if (hi >= 0 && hi >= lo && lo <= ray.max_distance) {

				out.count = 1;
				out.contact_point = ray.position + (ray.direction * lo);

				math::vec2 box_center = (box.min + box.max) / 2;
				box_center = out.contact_point - box_center;
				math::vec2 abs_c = math::VecAbs(box_center);
				if (abs_c.x > abs_c.y)
					out.contact_normal = { math::Sign(box_center.x), 0 };
				else
					out.contact_normal = { 0, math::Sign(box_center.y) };

				return 1;
			}
			return 0;
		}


		// Utility

		void build_regular_poly(int nSides, float side_length, body::poly& out) {

			float ang = (2 * math::Pi) / nSides;

			out.nVertex = nSides;

			for (int i = 0; i < nSides; i++) {
				out.vertex[i] = { side_length * math::Cos(ang * i), side_length * math::Sin(ang * i) };
			}
			for (int i = 0; i < nSides; i++) {
				math::vec2 plane = out.vertex[(i + 1) % nSides] - out.vertex[i];
				out.norms[i] = { -plane.y, plane.x };
			}
		}

		body::poly BoxToPoly(simple::collision::body::box& box) {
			simple::collision::body::poly out;
			out.nVertex = 4;
			out.vertex[0] = box.min;
			out.vertex[1] = simple::math::vec2(box.max.x, box.min.y);
			out.vertex[2] = box.max;
			out.vertex[3] = simple::math::vec2(box.min.x, box.max.y);
			return out;
		}

		void poly_add(body::poly& poly, float x, float y) {
			for (int i = 0; i < poly.nVertex; i++) {
				poly.vertex[i].x += x;
				poly.vertex[i].y += y;
			}
		}

	} // Collision end
}