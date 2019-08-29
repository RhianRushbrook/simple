#include "simple_camera.h"

namespace simple {

	namespace camera {

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
			math::vec2 initial = PortToWorld(vec);
			scale *= val;
			math::vec2 finish = PortToWorld(vec);
			Set(offset + (initial - finish));
		}

		void view_port::Zoom(float val, int x, int y) {
			math::vec2 initial = PortToWorld(x, y);
			scale *= val;
			math::vec2 finish = PortToWorld(x, y);
			Set(offset + (initial - finish));
		}


		math::vec2 view_port::WorldToPort(const math::vec2& vec) {
			return (vec - offset) * scale;
		}

		math::vec2 view_port::WorldToPort(const float x, const float y) {
			return (math::vec2(x, y) - offset) * scale;
		}

		math::vec2 view_port::PortToWorld(const math::vec2& vec) {
			return (vec / scale) + offset;
		}

		math::vec2 view_port::PortToWorld(const float x, const float y) {
			return (math::vec2(x, y) / scale) + offset;
		}

	} // Camera end
}