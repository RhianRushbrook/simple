#pragma once
#include "simple_math.h"

namespace simple {

	namespace camera {

		class view_port {

		private:
			math::vec2 offset;

			float width, height;
			float scale = 1;

		public:
			view_port() = delete;
			view_port(float w, float h) : width(w), height(h) { }

		public:
			math::vec2 Min() const {
				return offset;
			}
			math::vec2 Max() const {
				return (math::vec2(width, height) / scale) + offset;
			}

			float Scale() const {
				return scale;
			}

		public:
			void Set(const math::vec2& vec);
			void Set(float x, float y);

			void Pan(const math::vec2& vec);
			void Pan(float x, float y);

			void Zoom(float val, const math::vec2& vec = { 0.0f, 0.0f });
			void Zoom(float val, int x = 0, int y = 0);

		public:
			math::vec2 WorldToPort(const math::vec2& vec);
			math::vec2 WorldToPort(const float x, const float y);

			math::vec2 PortToWorld(const math::vec2& vec);
			math::vec2 PortToWorld(const float x, const float y);
		};

	} // Camera end
}